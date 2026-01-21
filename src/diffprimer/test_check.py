from diffprimer.kmer_extractor import process_seqs, check_specificity, PrimerCandidate

print("Step 1: Processing sequences to find unique regions...")
contig_regions = process_seqs(
    reference_path="/home/matheus/Documents/ESALQ/pedro_vila/data/manisopliae/ref/ref5.fasta",
    sequences_dir="/home/matheus/Documents/ESALQ/pedro_vila/data/manisopliae/selected_fasta/",
    k=21,
    num_threads=4,
    min_region_length=200,
    regions_output_path="/home/matheus/Documents/ESALQ/pedro_vila/data/manisopliae/results/exclusive_regions.txt",
    fasta_output_path="/home/matheus/Documents/ESALQ/pedro_vila/data/manisopliae/results/exclusive_sequences.fasta",
    max_abundance=1
)

print(f"Found {len(contig_regions)} contigs with unique regions.")

# Create mock candidates
candidates = []
for contig in contig_regions:
    # Use rust Region object properties
    # Note: contig.regions is a list of Region objects
    for region in contig.regions:
        # Extract sequence
        # Rust ContigRegions has full sequence
        # Region has start/end
        seq = contig.sequence[region.start : region.end]
        
        if len(seq) < 50:
            continue
            
        # Mock primers: first 20bp and last 20bp
        left_seq = seq[:20]
        right_seq = seq[-20:]
        
        cand = PrimerCandidate(
            header=f"{contig.header}_mock_{region.start}",
            region_sequence=seq,
            start=region.start,
            end=region.end,
            left_primer_seq=left_seq,
            right_primer_seq=right_seq,
            left_primer_offset=0,
            right_primer_offset=len(seq)-20
        )
        candidates.append(cand)

print(f"Created {len(candidates)} mock primer candidates. Checking specificity...")

# Check specificity
results = check_specificity(
    candidates=candidates,
    sequences_dir="/home/matheus/Documents/ESALQ/pedro_vila/data/manisopliae/selected_fasta/",
    similarity_threshold=80.0,
    local_mismatch_threshold=3,
    num_threads=4
)

print("Results:")
for r in results:
    print(f"{r.region_header}: {r.tag} (MaxSim: {r.max_similarity:.1f}%, Target: {r.most_similar_target})")