use pyo3::prelude::*;
use pyo3::types::PyModule;
use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::cell::RefCell;
use bio::alignment::distance::simd::bounded_levenshtein;
use bio::pattern_matching::myers::long::Myers;
use bio::alignment::pairwise::Aligner;
use bio::alignment::AlignmentOperation;
use console::style;

/// Scoring function for semiglobal DNA alignment.
/// Defined as a named fn so it can be used as a `fn` pointer in `thread_local!`.
fn dna_score(a: u8, b: u8) -> i32 {
    if a == b { 1 } else { -1 }
}

thread_local! {
    /// Per-thread Aligner, reused across every alignment on a given Rayon worker.
    /// Avoids repeated heap allocation of the DP scoring matrix (typ. ~120 KB each).
    static THREAD_ALIGNER: RefCell<Aligner<fn(u8, u8) -> i32>> =
        RefCell::new(Aligner::with_capacity(200, 300, -1, -1, dna_score));
}

// 2-bit encoding for DNA bases
const A: u64 = 0b00;
const C: u64 = 0b01;
const G: u64 = 0b10;
const T: u64 = 0b11;

/// A struct representing a k-mer (substring of length k) encoded as a 2-bit integer.
#[derive(Clone, Copy, Debug)]
struct Kmer {
    /// The 2-bit encoded value of the k-mer.
    value: u64,
    /// The length of the k-mer.
    k: usize,
}

impl Kmer {
    /// Creates a new empty Kmer.
    ///
    /// # Arguments
    ///
    /// * `k` - The size of the k-mer. Must be <= 31.
    fn new(k: usize) -> Self {
        assert!(k <= 31, "k must be <= 31 for u64 representation");
        Kmer { value: 0, k }
    }

    /// Pushes a new base into the k-mer, shifting the previous value.
    ///
    /// Returns `Some(())` if the base is valid (A, C, G, T), or `None` otherwise.
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

    /// Returns the k-mer value as a u64 integer.
    fn as_u64(&self) -> u64 {
        self.value
    }
}

// Enhanced FASTA reader for multi-FASTA files
// Enhanced FASTA reader for multi-FASTA files
/// A buffered reader for efficiently processing multi-FASTA files.
struct MultiFastaReader {
    reader: BufReader<File>,
    line_buffer: String,
}

impl MultiFastaReader {
    /// Creates a new MultiFastaReader for the given file path.
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
/// Extracts unique k-mers from a FASTA file.
///
/// Reads the file efficiently using `MultiFastaReader` and returns a set of unique k-mer hash values.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file.
/// * `k` - The k-mer size.
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

/// Counts k-mer occurrences in a reference FASTA file.
///
/// # Arguments
///
/// * `path` - Path to the reference FASTA file.
/// * `k` - The k-mer size.
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

/// Identifies k-mers that are exclusive to the reference genome.
///
/// This function computes the set of k-mers present in the reference genome (with abundance <= `max_abundance`)
/// and removes any k-mers found in other sequence files in `sequences_dir`.
///
/// # Arguments
///
/// * `reference_path` - Path to the reference FASTA file.
/// * `sequences_dir` - Directory containing other FASTA files to compare against.
/// * `k` - The k-mer size.
/// * `num_threads` - Number of threads for parallel processing.
/// * `max_abundance` - Maximum allowed abundance in the reference for a k-mer to be considered.
///
/// # Returns
///
/// A `Result` containing a `HashSet` of exclusive k-mers.
pub fn get_exclusive_kmers(
    reference_path: &Path,
    sequences_dir: &Path,
    k: usize,
    num_threads: usize,
    max_abundance: u32,
) -> std::io::Result<FxHashSet<u64>> {
    // Configure Rayon thread pool
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global();

    // -- Step 1: scan reference and count k-mers ----------------------------
    let spinner = indicatif::ProgressBar::new_spinner();
    spinner.set_style(
        indicatif::ProgressStyle::default_spinner()
            .template("{spinner:.cyan} {msg}")
            .unwrap()
            .tick_strings(&[
                "-", "\\", "|", "/",
            ]),
    );
    spinner.set_message(format!(
        "{} Scanning reference for {}-mers...",
        style("[1/3]").bold().cyan(),
        style(k).bold()
    ));
    spinner.enable_steady_tick(std::time::Duration::from_millis(80));

    let reference_counts = count_kmers_from_reference(reference_path, k)?;
    spinner.finish_with_message(format!(
        "{} {}  --  {} unique {}-mers found in reference",
        style("[1/3]").bold().cyan(),
        style("Done").bold().green(),
        style(reference_counts.len()).bold().white(),
        k
    ));
    
    // Filter k-mers that exceed max abundance in reference
    let mut exclusive_kmers: FxHashSet<u64> = reference_counts
        .into_iter()
        .filter(|&(_, count)| count <= max_abundance)
        .map(|(kmer, _)| kmer)
        .collect();
    exclusive_kmers.shrink_to_fit();

    println!(
        "      After abundance filter (<= {})  --  {} candidate k-mers retained",
        style(max_abundance).bold(),
        style(exclusive_kmers.len()).bold().white()
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
        "\n{} Subtracting k-mers  --  {} comparison file(s)  |  {} thread(s)",
        style("[2/3]").bold().cyan(),
        style(fasta_files.len()).bold(),
        style(num_threads).bold()
    );

    let progress_bar = indicatif::ProgressBar::new(fasta_files.len() as u64);
    progress_bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("    {prefix:.bold}  [{elapsed_precise}] {bar:35.cyan/blue} {pos:>3}/{len} files  {eta} remaining")
            .unwrap()
            .progress_chars("=>-"),
    );
    progress_bar.set_prefix("Subtracting");

    // LOCK-FREE PIPELINE: Rayon extracts k-mers from all files in parallel on a
    // dedicated thread. The calling thread drains the channel and subtracts in-place,
    // pipelined with extraction. The bounded channel caps memory to
    // O(num_threads × max_file_kmers) — no Mutex contention.
    let (tx, rx) = mpsc::sync_channel::<FxHashSet<u64>>(num_threads.max(4));
    let bar_clone = progress_bar.clone();

    let extraction_handle = std::thread::spawn(move || {
        fasta_files.par_iter().for_each_with(tx, |tx, seq_path| {
            match extract_kmers_from_sequence(seq_path, k) {
                Ok(file_kmers) => {
                    bar_clone.inc(1);
                    // Blocks when channel is full — natural backpressure limits memory
                    let _ = tx.send(file_kmers);
                }
                Err(e) => {
                    bar_clone.println(format!(
                        "    {}: Error processing {:?}: {}",
                        style("WARNING").bold().yellow(),
                        seq_path,
                        e
                    ));
                }
            }
        });
        // tx is dropped here — channel closes, the rx loop below terminates
    });

    // Sequential subtraction, pipelined with parallel extraction above
    for file_kmers in rx {
        for kmer in &file_kmers {
            exclusive_kmers.remove(kmer);
        }
        // file_kmers dropped here — no accumulation
    }
    extraction_handle.join().expect("k-mer extraction thread panicked");
    progress_bar.finish_with_message("");
    progress_bar.finish_and_clear();

    println!(
        "      {} Subtraction done  --  {} exclusive k-mers remaining",
        style("Done.").bold().green(),
        style(exclusive_kmers.len()).bold().white()
    );
    Ok(exclusive_kmers)
}

// Find positions of k-mers in reference contigs, merge into regions, and extract subsequences for output
/// A struct representing a genomic region.
#[pyclass]
#[derive(Debug, Clone)]
pub struct Region {
    #[pyo3(get, set)]
    pub start: usize,
    #[pyo3(get, set)]
    pub end: usize,
    /// The extracted subsequence for this region (avoids storing the full contig).
    #[pyo3(get, set)]
    pub subsequence: String,
}

#[pymethods]
impl Region {
    /// Creates a new Region.
    #[new]
    fn new(start: usize, end: usize, subsequence: String) -> Self {
        Region { start, end, subsequence }
    }

    fn __repr__(&self) -> String {
        format!("Region(start={}, end={}, len={})", self.start, self.end, self.subsequence.len())
    }
}

/// A struct collecting all exclusive regions found in a contig.
#[pyclass]
#[derive(Debug, Clone)]
pub struct ContigRegions {
    #[pyo3(get, set)]
    pub header: String,
    #[pyo3(get, set)]
    pub regions: Vec<Region>,
}

#[pymethods]
impl ContigRegions {
    /// Creates a new ContigRegions object.
    #[new]
    fn new(header: String, regions: Vec<Region>) -> Self {
        ContigRegions {
            header,
            regions,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ContigRegions(header='{}', num_regions={})",
            self.header,
            self.regions.len()
        )
    }
}

// Find k-mer positions in a single sequence string
/// Finds all occurrences of exclusive k-mers in a sequence.
///
/// Scans the sequence using a sliding window and checks if each k-mer is present in the `exclusive_kmers` set.
///
/// # Arguments
///
/// * `sequence` - The DNA sequence to scan.
/// * `exclusive_kmers` - The set of exclusive k-mers to look for.
/// * `k` - The k-mer size.
///
/// # Returns
///
/// A vector of starting positions (0-based) for each found exclusive k-mer.
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

// Merge adjacent/overlapping k-mer positions into (start, end) coordinate pairs.
// Returns lightweight tuples — full Region structs with subsequences are built later.
/// Merges adjacent or overlapping k-mer positions into continuous regions.
///
/// # Arguments
///
/// * `positions` - A vector of k-mer start positions.
/// * `k` - The k-mer size (used to determine the end of each k-mer).
/// * `min_length` - Minimum length for a merged region to be kept.
///
/// # Returns
///
/// A vector of `(start, end)` tuples representing the merged regions.
fn merge_positions_to_regions(positions: Vec<usize>, k: usize, min_length: usize) -> Vec<(usize, usize)> {
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
                regions.push((current_start, current_end));
            }
            current_start = pos;
            current_end = pos + k;
        }
    }

    // Add the last region if it meets minimum length
    let region_length = current_end - current_start;
    if region_length >= min_length {
        regions.push((current_start, current_end));
    }

    regions
}

// MEMORY-OPTIMIZED: Process reference contigs one-at-a-time, extract region
// subsequences immediately, write output files inline, and drop full sequences.
// This avoids collecting all (header, sequence) pairs into memory.
/// Processes reference contigs to identify and extract exclusive regions.
///
/// This function streams the reference file, finds exclusive kmers in each contig, 
/// merges them into regions, and writes the results to output files immediately.
///
/// # Arguments
///
/// * `reference_path` - Path to the reference FASTA file.
/// * `exclusive_kmers` - Set of exclusive k-mers.
/// * `k` - K-mer size.
/// * `min_length` - Minimum region length.
/// * `regions_output_path` - Path to write the regions summary text file.
/// * `fasta_output_path` - Path to write the regions FASTA file.
///
/// # Returns
///
/// A `Result` containing a vector of `ContigRegions` (with subsequences).
pub fn process_reference_contigs(
    reference_path: &Path,
    exclusive_kmers: &FxHashSet<u64>,
    k: usize,
    min_length: usize,
    regions_output_path: &Path,
    fasta_output_path: &Path,
) -> std::io::Result<Vec<ContigRegions>> {
    // Count contigs for progress bar (lightweight pass)
    let contig_count = {
        let reader = MultiFastaReader::new(reference_path)?;
        reader.count() as u64
    };

    println!(
        "\n{} Scanning {} contig(s) for exclusive regions  (k={}, min-len={} bp)...",
        style("[2/3]").bold().cyan(),
        style(contig_count).bold(),
        style(k).bold(),
        style(min_length).bold()
    );

    let progress_bar = indicatif::ProgressBar::new(contig_count);
    progress_bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("    {prefix:.bold}  [{elapsed_precise}] {bar:35.cyan/blue} {pos:>5}/{len} contigs  {eta} remaining")
            .unwrap()
            .progress_chars("=>-"),
    );
    progress_bar.set_prefix("Scanning  ");

    // Open output files once
    let mut regions_file = File::create(regions_output_path)?;
    let mut fasta_file = File::create(fasta_output_path)?;
    let mut total_fasta_regions = 0usize;

    // Stream contigs one at a time — only one full sequence in memory at a time
    let reader = MultiFastaReader::new(reference_path)?;
    let mut all_contig_regions = Vec::new();

    for result in reader {
        let (header, sequence) = result?;

        // Find k-mer positions in this sequence
        let positions = find_kmer_positions_in_sequence(&sequence, exclusive_kmers, k);

        // Merge positions into regions with minimum length filter
        let raw_regions = merge_positions_to_regions(positions, k, min_length);

        // Extract subsequences and build Region structs with embedded subsequences
        let regions: Vec<Region> = raw_regions
            .iter()
            .map(|&(start, end)| {
                let subseq = sequence[start..end].to_string();
                Region {
                    start,
                    end,
                    subsequence: subseq,
                }
            })

            .collect();

        // Write regions summary inline
        writeln!(regions_file, "HEADER: {}", header)?;
        writeln!(regions_file, "REGIONS: {} exclusive regions found", regions.len())?;
        for (idx, region) in regions.iter().enumerate() {
            writeln!(
                regions_file,
                "\tRegion {}: {}-{} (length: {})",
                idx + 1,
                region.start,
                region.end,
                region.end - region.start
            )?;
        }
        writeln!(regions_file)?;

        // Write FASTA entries inline
        for (region_idx, region) in regions.iter().enumerate() {
            let region_header = format!(
                "{}_region_{}_{}:{}",
                header, region_idx + 1, region.start, region.end
            );
            writeln!(fasta_file, ">{}", region_header)?;
            for chunk in region.subsequence.as_bytes().chunks(80) {
                writeln!(fasta_file, "{}", String::from_utf8_lossy(chunk))?;
            }
            total_fasta_regions += 1;
        }

        // Store lightweight ContigRegions (no full sequence)
        all_contig_regions.push(ContigRegions {
            header,
            regions,
        });

        // `sequence` is dropped here — memory freed immediately
        progress_bar.inc(1);
    }

    progress_bar.finish_and_clear();
    println!(
        "      {} Contig scan done  --  {} exclusive region(s) extracted",
        style("Done.").bold().green(),
        style(total_fasta_regions).bold().white()
    );
    println!("      Regions report  -->  {}", style(regions_output_path.display()).cyan());
    println!("      Regions FASTA   -->  {}", style(fasta_output_path.display()).cyan());

    Ok(all_contig_regions)
}

/// Python-exposed function to find regions with exclusive k-mers.
///
/// This is the main entry point for the k-mer analysis part of the pipeline.
/// It wraps the Rust logic and releases the Python GIL during heavy computations.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `reference_path` - Path to reference FASTA.
/// * `sequences_dir` - Directory with other sequence files.
/// * `k` - K-mer size.
/// * `num_threads` - Number of threads.
/// * `min_region_length` - Minimum region length.
/// * `regions_output_path` - Output path for text report.
/// * `fasta_output_path` - Output path for FASTA file.
/// * `max_abundance` - Max abundance in reference.
///
/// # Returns
///
/// A Python list of `ContigRegions` objects.
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
        // Step 1: Get exclusive k-mers (subtract-in-place, memory efficient)
        let exclusive_kmers = get_exclusive_kmers(
            &reference_path,
            &sequences_dir,
            k,
            num_threads,
            max_abundance,
        ).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        // Step 2+3+4: Stream reference contigs, find regions, write output files
        // inline — only one contig sequence in memory at a time
        let contig_regions = process_reference_contigs(
            &reference_path,
            &exclusive_kmers,
            k,
            min_region_length,
            &regions_output_path,
            &fasta_output_path,
        ).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        Ok::<_, std::io::Error>(contig_regions)
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{}", e)))?;

    // Convert to Python list
    Ok(pyo3::types::PyList::new(py, contig_regions)?)
}

// ============================================================================
// PRIMER SPECIFICITY ANALYSIS
// ============================================================================

/// A struct representing a candidate region for primer design.
///
/// This struct holds information about a specific genomic region and the primers
/// designed for it, which will be evaluated for specificity.
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
    /// Creates a new PrimerCandidate.
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

/// Result of primer specificity analysis for a single region.
#[pyclass]
#[derive(Debug, Clone)]
pub struct PrimerSpecificityResult {
    /// The header of the region being tested.
    #[pyo3(get, set)]
    pub region_header: String,
    /// start position of the region.
    #[pyo3(get, set)]
    pub region_start: usize,
    /// End position of the region.
    #[pyo3(get, set)]
    pub region_end: usize,
    /// The specificity classification tag.
    #[pyo3(get)]
    pub tag: PrimerSpecificityTag,
    /// The maximum similarity percentage found in the database.
    #[pyo3(get, set)]
    pub max_similarity: f64,
    /// The header of the target sequence with the highest similarity.
    #[pyo3(get, set)]
    pub most_similar_target: String,
    /// The local mismatch distance (number of edits) for the best match.
    #[pyo3(get, set)]
    pub local_distance: u32,
}

#[pymethods]
impl PrimerSpecificityResult {
    /// Creates a new PrimerSpecificityResult.
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

/// Calculates percentage similarity from Levenshtein distance.
///
/// # Arguments
///
/// * `distance` - The Levenshtein edit distance.
/// * `len_a` - Length of the first string.
/// * `len_b` - Length of the second string.
///
/// # Returns
///
/// The similarity percentage (0.0 to 100.0).
#[inline]
fn calculate_similarity(distance: u32, len_a: usize, len_b: usize) -> f64 {
    let max_len = len_a.max(len_b);
    if max_len == 0 {
        return 100.0;
    }
    (1.0 - (distance as f64 / max_len as f64)) * 100.0
}

/// Checks primer specificity for candidates against a database of sequences.
///
/// This function iterates over all candidates and all target sequences, finding the best match
/// for each candidate region in the target sequences. It computes global similarity and,
/// if a match is found, checks for local primer specificity.
///
/// # Arguments
///
/// * `candidates` - A slice of `PrimerCandidate`s to check.
/// * `target_sequences` - A slice of target sequences (header, sequence bytes).
/// * `similarity_threshold` - Minimum similarity percentage to consider a match "non-specific".
/// * `local_mismatch_threshold` - Minimum mismatches required in the primer binding site to be considered specific.
///
/// # Returns
///
/// A vector of `PrimerSpecificityResult`s.
pub fn check_primer_specificity_candidates(
    candidates: &[PrimerCandidate],
    target_sequences: &[(String, Vec<u8>)], // Placeholder for tuple type
    similarity_threshold: f64,
    local_mismatch_threshold: u32,
) -> Vec<PrimerSpecificityResult> {

    println!(
        "\n{} Specificity check  --  {} candidate(s)  x  {} target sequence(s)",
        style("[3/3]").bold().cyan(),
        style(candidates.len()).bold(),
        style(target_sequences.len()).bold()
    );

    let progress_bar = indicatif::ProgressBar::new(candidates.len() as u64);
    progress_bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("    {prefix:.bold}  [{elapsed_precise}] {bar:35.green/white} {pos:>5}/{len} primers  {eta} remaining  ({per_sec})")
            .unwrap()
            .progress_chars("#>-"),
    );
    progress_bar.set_prefix("Checking  ");

    let results: Vec<PrimerSpecificityResult> = candidates
        .par_iter()
        .map(|candidate| {
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
                    // Perform Semiglobal Alignment using the thread-local Aligner.
                    // Avoids allocating a fresh DP scoring matrix for every candidate-target pair.
                    // check_primer is defined inside the with() callback so that `alignment`
                    // (a reference into the Aligner's internal state) stays valid for both calls.
                    let (l_nonspecific, r_nonspecific) = THREAD_ALIGNER.with(|cell| {
                        let mut aligner = cell.borrow_mut();
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

                        let l = check_primer(candidate.left_primer_offset, &candidate.left_primer_seq);
                        let r = check_primer(candidate.right_primer_offset, &candidate.right_primer_seq);
                        (l, r)
                    });

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

            let result = PrimerSpecificityResult {
                region_header: candidate.header.clone(),
                region_start: candidate.start,
                region_end: candidate.end,
                tag,
                max_similarity,
                most_similar_target,
                local_distance,
            };

            progress_bar.inc(1);
            result
        })
        .collect();
        
    progress_bar.finish_and_clear();

    // -- Summary -------------------------------------------------------------
    let unique          = results.iter().filter(|r| r.tag == PrimerSpecificityTag::Unique_LowSim).count();
    let specific_in_sim = results.iter().filter(|r| r.tag == PrimerSpecificityTag::Specific_In_SimRegion).count();
    let nonspecific     = results.iter().filter(|r| r.tag == PrimerSpecificityTag::NonSpecific_HighSim).count();
    println!("      {} Specificity check done:", style("Done.").bold().green());
    println!("        Unique (low similarity)      {}", style(unique).bold().white());
    println!("        Specific in similar region   {}", style(specific_in_sim).bold().white());
    println!("        Non-specific (high sim hit)  {}", style(nonspecific).bold().white());
    
    results
}

/// Loads sequences from a directory of FASTA files.
///
/// Reads all FASTA files in the given directory and returns a collected vector of (header, sequence) tuples.
///
/// # Arguments
///
/// * `dir_path` - Path to the directory containing FASTA files.
#[allow(clippy::type_complexity)]
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

/// Python-exposed function to check primer specificity.
///
/// This function wraps the Rust `check_primer_specificity_candidates` logic,
/// handling loading of sequences and parallel execution.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `candidates` - List of `PrimerCandidate` objects.
/// * `sequences_dir` - Directory containing target sequences to check against.
/// * `similarity_threshold` - Global similarity threshold.
/// * `local_mismatch_threshold` - Local mismatch threshold for primer binding.
/// * `num_threads` - Number of threads to use.
///
/// # Returns
///
/// A Python list of `PrimerSpecificityResult` objects.
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

        println!(
            "      Loaded {} target sequence(s) from disk",
            style(target_sequences.len()).bold().white()
        );

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
