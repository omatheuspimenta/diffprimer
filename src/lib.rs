use pyo3::prelude::*;
use pyo3::types::PyModule;
use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use bio::alignment::distance::simd::bounded_levenshtein;
use bio::pattern_matching::myers::long::Myers;
use bio::alignment::pairwise::Aligner;
use bio::alignment::AlignmentOperation;

// 2-bit encoding for DNA bases
const A: u64 = 0b00;
const C: u64 = 0b01;
const G: u64 = 0b10;
const T: u64 = 0b11;

#[derive(Clone, Copy, Debug)]
struct Kmer {
    value: u64,
    k: usize,
}

impl Kmer {
    fn new(k: usize) -> Self {
        assert!(k <= 31, "k must be <= 31 for u64 representation");
        Kmer { value: 0, k }
    }

    fn push(&mut self, base: u8) -> Option<()> {
        let encoded = match base.to_ascii_uppercase() {
            b'A' => A,
            b'C' => C,
            b'G' => G,
            b'T' => T,
            _ => return None, // Skip N or invalid bases
        };

        let mask = (1u64 << (2 * self.k)) - 1;
        self.value = ((self.value << 2) | encoded) & mask;
        Some(())
    }

    fn as_u64(&self) -> u64 {
        self.value
    }
}

// Enhanced FASTA reader for multi-FASTA files
// Enhanced FASTA reader for multi-FASTA files
struct MultiFastaReader {
    reader: BufReader<File>,
    line_buffer: String,
}

impl MultiFastaReader {
    fn new(path: &Path) -> std::io::Result<Self> {
        let file = File::open(path)?;
        Ok(MultiFastaReader {
            reader: BufReader::with_capacity(1024 * 1024, file),
            line_buffer: String::new(),
        })
    }
}

impl Iterator for MultiFastaReader {
    type Item = std::io::Result<(String, String)>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut current_header = String::new();
        let mut current_sequence = String::new();

        // If buffer is empty, try to read first line
        if self.line_buffer.is_empty() {
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                     // Check if it's a header
                     let trimmed = self.line_buffer.trim();
                     if !trimmed.starts_with('>') {
                         // First line must be header, or skip empty lines until header
                         // For simplicity, if it's not a header at start properly, we might error or skip
                         // But let's assume valid FASTA or we handle leading whitespace
                     }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        // We expect line_buffer to contain the Header (starting with >) from previous iteration
        // OR the very first line read above.
        let trimmed_header = self.line_buffer.trim();
        if !trimmed_header.starts_with('>') {
             // If we have content but it's not a header, it's unexpected for a new record start
             // But maybe end of file?
             if trimmed_header.is_empty() {
                 return None;
             }
             // Should probably consume until next header?
        } else {
             current_header = trimmed_header[1..].to_string();
        }

        self.line_buffer.clear();

        // Read sequence lines until next header or EOF
        loop {
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    let trimmed = self.line_buffer.trim();
                    if trimmed.starts_with('>') {
                        // Found next header, stop reading sequence
                        // line_buffer now holds the NEXT header for the NEXT iteration
                        break;
                    } else if !trimmed.is_empty() {
                         current_sequence.push_str(trimmed);
                    }
                    self.line_buffer.clear();
                }
                Err(e) => return Some(Err(e)),
            }
        }

        if current_header.is_empty() && current_sequence.is_empty() {
            return None;
        }

        Some(Ok((current_header, current_sequence)))
    }
}

// Updated k-mer extraction function that handles multifasta files
// Updated k-mer extraction function that handles multifasta files (using iterator)
fn extract_kmers_from_sequence(path: &Path, k: usize) -> std::io::Result<FxHashSet<u64>> {
    let mut kmers = FxHashSet::default();
    let reader = MultiFastaReader::new(path)?;

    // Process each sequence in the file via streaming
    for result in reader {
        let (_header, sequence) = result?;
        let seq_bytes = sequence.as_bytes();
        let mut current_kmer = Kmer::new(k);
        let mut valid_count = 0;

        for &base in seq_bytes.iter() {
            if current_kmer.push(base).is_some() {
                valid_count += 1;
                if valid_count >= k {
                    kmers.insert(current_kmer.as_u64());
                }
            } else {
                // Reset on invalid base (N, etc.)
                valid_count = 0;
                current_kmer = Kmer::new(k);
            }
        }
    }

    Ok(kmers)
}

fn count_kmers_from_reference(path: &Path, k: usize) -> std::io::Result<FxHashMap<u64, u32>> {
    let mut kmer_counts = FxHashMap::default();
    let reader = MultiFastaReader::new(path)?;

    for result in reader {
        let (_header, sequence) = result?;
        let seq_bytes = sequence.as_bytes();
        let mut current_kmer = Kmer::new(k);
        let mut valid_count = 0;

        for &base in seq_bytes.iter() {
            if current_kmer.push(base).is_some() {
                valid_count += 1;
                if valid_count >= k {
                    *kmer_counts.entry(current_kmer.as_u64()).or_default() += 1;
                }
            } else {
                valid_count = 0;
                current_kmer = Kmer::new(k);
            }
        }
    }
    Ok(kmer_counts)
}

pub fn get_exclusive_kmers(
    reference_path: &Path,
    sequences_dir: &Path,
    k: usize,
    num_threads: usize,
    max_abundance: u32,
) -> std::io::Result<FxHashSet<u64>> {
    // Configure Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

    // Extract k-mers from reference sequence with counting
    let spinner = indicatif::ProgressBar::new_spinner();
    spinner.set_style(
        indicatif::ProgressStyle::default_spinner()
            .template("{spinner} {msg}")
            .unwrap()
            .tick_strings(&[
                "🔍 A ",
                "🔍 C ",
                "🔍 G ",
                "🔍 T ",
            ]),
    );
    spinner.set_message("Extracting and counting k-mers from reference sequence...");
    spinner.enable_steady_tick(std::time::Duration::from_millis(100));

    let reference_counts = count_kmers_from_reference(reference_path, k)?;
    spinner.finish_with_message(format!("Extracted {} unique k-mers from reference", reference_counts.len()));
    
    // Filter k-mers that exceed max abundance in reference
    let reference_kmers: FxHashSet<u64> = reference_counts
        .into_iter()
        .filter(|&(_, count)| count <= max_abundance)
        .map(|(kmer, _)| kmer)
        .collect();

    println!(
        "Found {} suitable k-mers in reference (abundance <= {})", 
        reference_kmers.len(), 
        max_abundance
    );

    // Get all FASTA files
    let fasta_files: Vec<PathBuf> = std::fs::read_dir(sequences_dir)?
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            let ext = path.extension()?.to_str()?;
            if ext == "fasta" || ext == "fna" || ext == "fa" || ext == "fas" {
                Some(path)
            } else {
                None
            }
        })
        .collect();

    println!(
        "Processing {} sequence files with {} threads...",
        fasta_files.len(),
        num_threads
    );

    let progress_bar = indicatif::ProgressBar::new(fasta_files.len() as u64);
    progress_bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta}) {msg}")
            .unwrap()
            .progress_chars("#>-"),
    );

    // Shared HashSet for all k-mers found in other sequences
    let all_other_kmers = Arc::new(Mutex::new(FxHashSet::default()));

    // Process all files in parallel, each thread adds to shared HashSet
    fasta_files.par_iter().for_each(|seq_path| {
        match extract_kmers_from_sequence(seq_path, k) {
            Ok(file_kmers) => {
                // Lock and insert k-mers from this file
                {
                    let mut shared_kmers = all_other_kmers.lock().unwrap();
                    shared_kmers.extend(file_kmers);
                }
                progress_bar.inc(1);
            }
            Err(e) => {
                progress_bar.println(format!("Error processing {:?}: {}", seq_path, e));
            }
        }
    });
    progress_bar.finish_with_message("Done processing files");

    // Extract final HashSet and compute difference
    let other_kmers = Arc::try_unwrap(all_other_kmers)
        .unwrap()
        .into_inner()
        .unwrap();
    println!(
        "Found {} total k-mers in other sequences",
        other_kmers.len()
    );

    // Compute exclusive k-mers (reference - others)
    let exclusive_kmers: FxHashSet<u64> =
        reference_kmers.difference(&other_kmers).copied().collect();

    println!("Found {} exclusive k-mers", exclusive_kmers.len());
    Ok(exclusive_kmers)
}

// Find positions of k-mers in sequence
#[pyclass]
#[derive(Debug, Clone)]
pub struct Region {
    #[pyo3(get, set)]
    pub start: usize,
    #[pyo3(get, set)]
    pub end: usize,
}

#[pymethods]
impl Region {
    #[new]
    fn new(start: usize, end: usize) -> Self {
        Region { start, end }
    }

    fn __repr__(&self) -> String {
        format!("Region(start={}, end={})", self.start, self.end)
    }
}

#[pyclass]
#[derive(Debug, Clone)]
pub struct ContigRegions {
    #[pyo3(get, set)]
    pub header: String,
    #[pyo3(get, set)]
    pub sequence: String,
    #[pyo3(get, set)]
    pub regions: Vec<Region>,
}

#[pymethods]
impl ContigRegions {
    #[new]
    fn new(header: String, sequence: String, regions: Vec<Region>) -> Self {
        ContigRegions {
            header,
            sequence,
            regions,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ContigRegions(header='{}', sequence_length={}, num_regions={})",
            self.header,
            self.sequence.len(),
            self.regions.len()
        )
    }
}

// Find k-mer positions in a single sequence string
fn find_kmer_positions_in_sequence(
    sequence: &str,
    exclusive_kmers: &FxHashSet<u64>,
    k: usize,
) -> Vec<usize> {
    let mut positions = Vec::new();
    let seq_bytes = sequence.as_bytes();
    let mut current_kmer = Kmer::new(k);
    let mut valid_count = 0;

    for (pos, &base) in seq_bytes.iter().enumerate() {
        if current_kmer.push(base).is_some() {
            valid_count += 1;
            if valid_count >= k {
                if exclusive_kmers.contains(&current_kmer.as_u64()) {
                    // Position is the start of the k-mer
                    let kmer_start = pos + 1 - k;
                    positions.push(kmer_start);
                }
            }
        } else {
            // Reset on invalid base (N, etc.)
            valid_count = 0;
            current_kmer = Kmer::new(k);
        }
    }

    positions
}

// Update the merge_positions_to_regions function to accept min_length:
fn merge_positions_to_regions(positions: Vec<usize>, k: usize, min_length: usize) -> Vec<Region> {
    if positions.is_empty() {
        return Vec::new();
    }

    let mut sorted_positions = positions;
    sorted_positions.sort_unstable();

    let mut regions = Vec::new();
    let mut current_start = sorted_positions[0];
    let mut current_end = sorted_positions[0] + k;

    for &pos in &sorted_positions[1..] {
        if pos <= current_end {
            // Extend current region
            current_end = current_end.max(pos + k);
        } else {
            // Save current region if it meets minimum length
            let region_length = current_end - current_start;
            if region_length >= min_length {
                regions.push(Region {
                    start: current_start,
                    end: current_end,
                });
            }
            current_start = pos;
            current_end = pos + k;
        }
    }

    // Add the last region if it meets minimum length
    let region_length = current_end - current_start;
    if region_length >= min_length {
        regions.push(Region {
            start: current_start,
            end: current_end,
        });
    }

    regions
}

// Main function to find regions for all contigs
pub fn find_kmer_occurrence_regions_multi_fasta(
    reference_path: &Path,
    exclusive_kmers: &FxHashSet<u64>,
    k: usize,
    min_length: usize,
) -> std::io::Result<Vec<ContigRegions>> {
    // NOTE: This now reads all sequences into memory to parallelize.
    // Ideally, we would stream + parallelize, but par_iter requires the collection to exist.
    // For now, let's keep the read_all logic here locally or re-use MultiFastaReader 
    // to collect into a Vec, then par_iter.
    
    // We cannot use the streaming iterator directly with Rayon's par_iter easily without collecting.
    // So we collect (which is same memory usage as before, but validation logic is better).
    // Optimization: In future we could use Rayon's par_bridge() on the iterator!
    
    let mut sequences = Vec::new();
    let reader = MultiFastaReader::new(reference_path)?;
    for result in reader {
        sequences.push(result?);
    }

    println!(
        "Processing {} contigs for exclusive regions (min length: {})...",
        sequences.len(),
        min_length
    );
    
    let progress_bar = indicatif::ProgressBar::new(sequences.len() as u64);
    progress_bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta}) {msg}")
            .unwrap()
            .progress_chars("#>-"),
    );

    let all_contig_regions: Vec<ContigRegions> = sequences
        .par_iter()
        .map(|(header, sequence)| {
            // Use thread-local increment? indicatif handles this.
            progress_bar.inc(1);
            
            // Find k-mer positions in this sequence
            let positions = find_kmer_positions_in_sequence(sequence, exclusive_kmers, k);

            // Merge positions into regions with minimum length filter
            let regions = merge_positions_to_regions(positions, k, min_length);

            ContigRegions {
                header: header.clone(),
                sequence: sequence.clone(),
                regions,
            }
        })
        .collect();

    progress_bar.finish_with_message("Done processing contigs");

    Ok(all_contig_regions)
}

// Write regions summary to file
pub fn write_regions_file(
    contig_regions: &[ContigRegions],
    output_path: &Path,
) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;

    for contig in contig_regions {
        writeln!(file, "HEADER: {}", contig.header)?;
        writeln!(
            file,
            "REGIONS: {} exclusive regions found",
            contig.regions.len()
        )?;

        for (idx, region) in contig.regions.iter().enumerate() {
            writeln!(
                file,
                "\tRegion {}: {}-{} (length: {})",
                idx + 1,
                region.start,
                region.end,
                region.end - region.start
            )?;
        }
        writeln!(file)?; // Empty line between contigs
    }

    println!("Regions summary written to: {:?}", output_path);
    Ok(())
}

// Write exclusive regions as FASTA sequences
pub fn write_exclusive_regions_fasta(
    contig_regions: &[ContigRegions],
    output_path: &Path,
) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;
    let mut total_regions = 0;

    for contig in contig_regions {
        for (region_idx, region) in contig.regions.iter().enumerate() {
            // Extract sequence for this region
            let region_sequence = &contig.sequence[region.start..region.end];

            // Create header for this region
            let region_header = format!("{}_region_{}_{}:{}", contig.header, region_idx + 1, region.start, region.end);

            // Write FASTA entry
            writeln!(file, ">{}", region_header)?;

            // Write sequence in lines of 80 characters (standard FASTA format)
            for chunk in region_sequence.as_bytes().chunks(80) {
                writeln!(file, "{}", String::from_utf8_lossy(chunk))?;
            }

            total_regions += 1;
        }
    }

    println!(
        "Extracted {} exclusive regions to FASTA: {:?}",
        total_regions, output_path
    );
    Ok(())
}

#[pyfunction]
fn process_seqs<'py>(
    py: Python<'py>,
    reference_path: &str,
    sequences_dir: &str,
    k: usize,
    num_threads: usize,
    min_region_length: usize,
    regions_output_path: &str,
    fasta_output_path: &str,
    max_abundance: u32,
) -> PyResult<Bound<'py, pyo3::types::PyList>> {
    let reference_path = PathBuf::from(reference_path);
    let sequences_dir = PathBuf::from(sequences_dir);
    let regions_output_path = PathBuf::from(regions_output_path);
    let fasta_output_path = PathBuf::from(fasta_output_path);

    // Release GIL for expensive computations
    #[allow(deprecated)]
    let contig_regions = py.allow_threads(move || {
        // Step 1: Get exclusive k-mers
        let exclusive_kmers = get_exclusive_kmers(
            &reference_path,
            &sequences_dir,
            k,
            num_threads,
            max_abundance,
        ).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        // Step 2: Find k-mer occurrence regions in reference
        let contig_regions = find_kmer_occurrence_regions_multi_fasta(
            &reference_path,
            &exclusive_kmers,
            k,
            min_region_length,
        ).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        // Step 3: Write regions summary to file
        write_regions_file(&contig_regions, &regions_output_path)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        // Step 4: Write exclusive regions as FASTA
        write_exclusive_regions_fasta(&contig_regions, &fasta_output_path)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        Ok::<_, std::io::Error>(contig_regions)
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;



    // Convert to Python list
    Ok(pyo3::types::PyList::new(py, contig_regions)?)
}

// ============================================================================
// PRIMER SPECIFICITY ANALYSIS
// ============================================================================

#[pyclass]
#[derive(Debug, Clone)]
pub struct PrimerCandidate {
    #[pyo3(get, set)]
    pub header: String,
    #[pyo3(get, set)]
    pub region_sequence: String,
    #[pyo3(get, set)]
    pub start: usize,
    #[pyo3(get, set)]
    pub end: usize,
    #[pyo3(get, set)]
    pub left_primer_seq: String,
    #[pyo3(get, set)]
    pub right_primer_seq: String,
    #[pyo3(get, set)]
    pub left_primer_offset: usize,
    #[pyo3(get, set)]
    pub right_primer_offset: usize,
}

#[pymethods]
impl PrimerCandidate {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        header: String, 
        region_sequence: String, 
        start: usize, 
        end: usize,
        left_primer_seq: String, 
        right_primer_seq: String,
        left_primer_offset: usize,
        right_primer_offset: usize,
    ) -> Self {
        PrimerCandidate {
            header, region_sequence, start, end, 
            left_primer_seq, right_primer_seq,
            left_primer_offset, right_primer_offset
        }
    }
}

/// Tag indicating primer specificity classification
#[pyclass]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_camel_case_types)]
pub enum PrimerSpecificityTag {
    /// Global similarity below threshold - primer is likely specific (safe)
    Unique_LowSim,
    /// High global similarity AND primer region is conserved - risk of false positive
    NonSpecific_HighSim,
    /// High global similarity BUT primer region has significant gaps/mismatches (safe)
    Specific_In_SimRegion,
}

#[pymethods]
impl PrimerSpecificityTag {
    fn __repr__(&self) -> String {
        match self {
            PrimerSpecificityTag::Unique_LowSim => "Unique_LowSim".to_string(),
            PrimerSpecificityTag::NonSpecific_HighSim => "NonSpecific_HighSim".to_string(),
            PrimerSpecificityTag::Specific_In_SimRegion => "Specific_In_SimRegion".to_string(),
        }
    }
    
    fn __str__(&self) -> String {
        self.__repr__()
    }
}

/// Result of primer specificity analysis for a single region
#[pyclass]
#[derive(Debug, Clone)]
pub struct PrimerSpecificityResult {
    #[pyo3(get, set)]
    pub region_header: String,
    #[pyo3(get, set)]
    pub region_start: usize,
    #[pyo3(get, set)]
    pub region_end: usize,
    #[pyo3(get)]
    pub tag: PrimerSpecificityTag,
    #[pyo3(get, set)]
    pub max_similarity: f64,
    #[pyo3(get, set)]
    pub most_similar_target: String,
    #[pyo3(get, set)]
    pub local_distance: u32,
}

#[pymethods]
impl PrimerSpecificityResult {
    #[new]
    fn new(
        region_header: String,
        region_start: usize,
        region_end: usize,
        tag: PrimerSpecificityTag,
        max_similarity: f64,
        most_similar_target: String,
        local_distance: u32,
    ) -> Self {
        PrimerSpecificityResult {
            region_header,
            region_start,
            region_end,
            tag,
            max_similarity,
            most_similar_target,
            local_distance,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PrimerSpecificityResult(region='{}', start={}, end={}, tag={:?}, max_sim={:.2}%, target='{}', local_dist={})",
            self.region_header,
            self.region_start,
            self.region_end,
            self.tag,
            self.max_similarity,
            self.most_similar_target,
            self.local_distance,
        )
    }
}

/// Calculate percentage similarity from Levenshtein distance
#[inline]
fn calculate_similarity(distance: u32, len_a: usize, len_b: usize) -> f64 {
    let max_len = len_a.max(len_b);
    if max_len == 0 {
        return 100.0;
    }
    (1.0 - (distance as f64 / max_len as f64)) * 100.0
}

/// Check primer specificity for candidates against a database of sequences.
pub fn check_primer_specificity_candidates(
    candidates: &[PrimerCandidate],
    target_sequences: &[(String, Vec<u8>)], // Changed to Vec<u8>
    similarity_threshold: f64,
    local_mismatch_threshold: u32,
) -> Vec<PrimerSpecificityResult> {

    
    println!(
        "Checking specificity for {} primer candidates against {} targets...",
        candidates.len(),
        target_sequences.len()
    );

    let progress_bar = indicatif::ProgressBar::new(candidates.len() as u64);
    progress_bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.green/white} {pos}/{len} ({eta}) {msg}")
            .unwrap()
            .progress_chars("█▓░"),
    );

    let results: Vec<PrimerSpecificityResult> = candidates
        .par_iter()
        .map(|candidate| {
            progress_bar.inc(1);
            let region_seq = candidate.region_sequence.as_bytes();
            let region_len = region_seq.len();
            
            // 1. Calculate max allowed distance (Global Similarity Threshold)
            let max_allowed_distance = ((region_len as f64) * (1.0 - similarity_threshold / 100.0)) as usize;
            
            // 2. Myers Bit-Vector Algorithm for Global Search
            // Using long::Myers for arbitrary length patterns
            let myers = Myers::<u64>::new(region_seq);
            
            let mut min_distance: usize = usize::MAX;
            let mut best_target_match: Option<(&str, usize)> = None; // (target_header, distance)
            let mut found_nonspecific = false;

            // Iterate ALL targets (do not break on error, to find best statistic)
            for (target_header, target_bytes) in target_sequences {
                
                // Fast check: find matches within max_allowed_distance
                let mut best_in_target = usize::MAX;
                let mut best_end_in_target = 0;
                
                // Note: we track the BEST match in this target
                for (end_pos, dist) in myers.find_all_end(target_bytes, max_allowed_distance) {
                    if dist < best_in_target {
                        best_in_target = dist;
                        best_end_in_target = end_pos;
                    }
                }
                
                // If we found a match better than current global best, update stats
                if best_in_target < min_distance {
                    min_distance = best_in_target;
                    best_target_match = Some((target_header, best_in_target));
                }

                // If this target is a "Hit" (similarity > threshold), we MUST check for specificity
                if best_in_target <= max_allowed_distance {
                    // Case A: Global Similarity > Threshold (Match Found)
                    // Now check Local Specificity using robust Alignment
                    
                    // We found a match ending at `best_end_in_target`
                    // Extract a generous window around it to allow for indel drift
                    let target_len = target_bytes.len();
                    
                    // Window calculation: 
                    // Region len + 50bp buffer on each side.
                    // Start approx: end - region_len - buffer.
                    let buffer = 50; 
                    let window_end = (best_end_in_target + buffer).min(target_len);
                    let window_start = best_end_in_target.saturating_sub(region_len + buffer);
                    
                    let target_window = &target_bytes[window_start..window_end];
                    
                    // Perform Semiglobal Alignment
                    // Region vs TargetWindow
                    // Scores: Match=1, Mismatch=-1, GapOpen=-1, GapExtend=-1 (Unit cost approx)
                    // Or Levenshtein-like: Match=0, Mismatch=-1, Gap=-1? 
                    // Aligner maximizes score. Levenshtein minimizes distance.
                    // To approximate Levenshtein: Match=0, Subst=-1, Gap=-1.
                    // But Aligner usually uses positive match. Match=1, Subst=-1, Gap=-1.
                    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
                    let mut aligner = Aligner::with_capacity(region_len, target_window.len(), -1, -1, score);
                    let alignment = aligner.semiglobal(region_seq, target_window);
                    
                    // Map primer coordinates using alignment path
                    let check_primer = |p_offset: usize, p_seq: &str| -> bool {
                        let p_len = p_seq.len();
                        let p_end_offset = p_offset + p_len;
                        
                        let mut t_start_idx = None;
                        let mut t_end_idx = None;
                        
                        let mut x = alignment.xstart; // region index
                        let mut y = alignment.ystart; // window index
                        
                        // If alignment starts after primer (clipping), we consider it mismatch?
                        if x > p_offset { return false; }

                        let mut covered = false;

                        for op in &alignment.operations {
                             if x == p_offset { t_start_idx = Some(y); }
                             if x == p_end_offset { t_end_idx = Some(y); covered = true; break; }
                             
                             match op {
                                 AlignmentOperation::Match | AlignmentOperation::Subst => {
                                     x += 1; y += 1;
                                 },
                                 AlignmentOperation::Ins => { y += 1; },
                                 AlignmentOperation::Del => { x += 1; },
                                 _ => {}
                             }
                        }
                        if x == p_end_offset { t_end_idx = Some(y); covered = true; }
                        
                        if let (Some(ts), Some(te)) = (t_start_idx, t_end_idx) {
                             if te > target_window.len() { return false; }
                             let t_slice = &target_window[ts..te];
                             // Calculate mismatch distance
                             // Use bounded_levenshtein for exact distance check
                             if let Some(dist) = bounded_levenshtein(p_seq.as_bytes(), t_slice, local_mismatch_threshold as u32 + 1) {
                                 return dist < local_mismatch_threshold as u32; // Specific if dist >= threshold
                                 // Return TRUE if NonSpecific (match found)
                             }
                        }
                        false // Default to Specific if mapping fails
                    };
                    
                    let l_nonspecific = check_primer(candidate.left_primer_offset, &candidate.left_primer_seq);
                    let r_nonspecific = check_primer(candidate.right_primer_offset, &candidate.right_primer_seq);
                    
                    if l_nonspecific && r_nonspecific {
                         found_nonspecific = true;
                         // We continue loop to improve min_distance statistics, 
                         // but we already know tag is NonSpecific.
                    }
                }
            }
            
            // Determine result tag
            let max_similarity = calculate_similarity(min_distance as u32, region_len, region_len); // Approx len
            
            let tag = if max_similarity < similarity_threshold {
                PrimerSpecificityTag::Unique_LowSim
            } else if found_nonspecific {
                PrimerSpecificityTag::NonSpecific_HighSim
            } else {
                PrimerSpecificityTag::Specific_In_SimRegion
            };
            
            let most_similar_target = best_target_match.map(|(h, _)| h.to_string()).unwrap_or_default();
            let local_distance = if min_distance == usize::MAX { u32::MAX } else { min_distance as u32 };

            PrimerSpecificityResult {
                region_header: candidate.header.clone(),
                region_start: candidate.start,
                region_end: candidate.end,
                tag,
                max_similarity,
                most_similar_target,
                local_distance,
            }
        })
        .collect();
        
    progress_bar.finish_with_message("Specificity check complete");
    
    // Summary output
    let unique = results.iter().filter(|r| r.tag == PrimerSpecificityTag::Unique_LowSim).count();
    let specific_in_sim = results.iter().filter(|r| r.tag == PrimerSpecificityTag::Specific_In_SimRegion).count();
    let nonspecific = results.iter().filter(|r| r.tag == PrimerSpecificityTag::NonSpecific_HighSim).count();
    println!("Results: {} Unique, {} Specific in Sim, {} NonSpecific", unique, specific_in_sim, nonspecific);
    
    results
}

/// Load sequences from a directory of FASTA files
pub fn load_sequences_from_directory(dir_path: &Path) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    let fasta_files: Vec<PathBuf> = std::fs::read_dir(dir_path)?
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            let ext = path.extension()?.to_str()?;
            if ext == "fasta" || ext == "fna" || ext == "fa" || ext == "fas" {
                Some(path)
            } else {
                None
            }
        })
        .collect();

    let mut all_sequences = Vec::new();
    for fasta_path in fasta_files {
        let reader = MultiFastaReader::new(&fasta_path)?;
        for result in reader {
            let (h, s) = result?;
            all_sequences.push((h, s.into_bytes()));
        }
    }
    
    Ok(all_sequences)
}

#[pyfunction]
fn check_specificity<'py>(
    py: Python<'py>,
    candidates: Vec<PrimerCandidate>,
    sequences_dir: &str,
    similarity_threshold: f64,
    local_mismatch_threshold: u32,
    num_threads: usize,
) -> PyResult<Bound<'py, pyo3::types::PyList>> {
    let sequences_dir = PathBuf::from(sequences_dir);

    // Release GIL for expensive computations
    #[allow(deprecated)]
    let results = py.allow_threads(move || {
        // Configure Rayon thread pool
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global();

        // Load target sequences
        let target_sequences = load_sequences_from_directory(&sequences_dir)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        println!("Loaded {} target sequences for specificity check", target_sequences.len());

        // Run specificity check
        let results = check_primer_specificity_candidates(
            &candidates,
            &target_sequences,
            similarity_threshold,
            local_mismatch_threshold,
        );

        Ok::<_, std::io::Error>(results)
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;

    // Convert to Python list
    Ok(pyo3::types::PyList::new(py, results)?)
}

#[pymodule]
fn kmer_extractor(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Region>()?;
    m.add_class::<ContigRegions>()?;
    m.add_class::<PrimerCandidate>()?;
    m.add_class::<PrimerSpecificityTag>()?;
    m.add_class::<PrimerSpecificityResult>()?;
    m.add_function(wrap_pyfunction!(process_seqs, m)?)?;
    m.add_function(wrap_pyfunction!(check_specificity, m)?)?;
    Ok(())
}
