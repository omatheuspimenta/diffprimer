import pytest
import os
import sys
from pathlib import Path

# Ensure we can import the package
# This assumes the package is installed in the environment or available via source
# For the rust extension, it must be compiled and installed.

def test_exclusive_regions_identification(tmp_path):
    """
    Verify that unique regions are correctly identified.
    Reference: Has Contig_1 and Contig_2.
    Target: Has only Contig_1.
    Expected: Exclusive regions found in Contig_2. Contig_1 should have 0 exclusive regions.
    """
    # Create synthetic data setup
    base_dir = tmp_path / "data"
    base_dir.mkdir()
    reference_path = base_dir / "reference.fasta"
    sequences_dir = base_dir / "targets"
    sequences_dir.mkdir()
    target_path = sequences_dir / "target.fasta"

    # Define sequences
    # Contig 1: Shared exactly
    # 60bp of random-ish sequence
    contig1_seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    
    # Contig 2: Has a unique core
    # Flanking: 20bp shared with contig1 (to test boundaries)
    # Core: 60bp UNIQUE non-repetitive sequence
    # Flanking: 20bp shared
    flank = "ATGCATGCATGCATGCATGC"
    
    unique_core = (
        "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
    ) # 70bp
    
    contig2_seq = flank + unique_core + flank
    
    # Write Reference
    with open(reference_path, "w") as f:
        f.write(f">ref_contig_1\n{contig1_seq}\n")
        f.write(f">ref_contig_2\n{contig2_seq}\n")
        
    # Write Target (Only Contig 1 equivalent)
    with open(target_path, "w") as f:
        f.write(f">target_contig_1\n{contig1_seq}\n")
    
    # Output paths
    regions_out = tmp_path / "regions.txt"
    fasta_out = tmp_path / "exclusive.fasta"
    
    k = 21
    min_len = 50 
    
    try:
        from diffprimer.kmer_extractor import process_seqs
    except ImportError:
        pytest.fail("Could not import diffprimer.kmer_extractor. Ensure the extension is built/installed.")

    # Run extraction
    # Signature: 
    # process_seqs(py, reference_path, sequences_dir, k, num_threads, min_region_length, regions_output_path, fasta_output_path)
    # Note: 'py' argument is handled by pyo3 automatically when called from Python
    
    results = process_seqs(
        str(reference_path),
        str(sequences_dir),
        k,
        1, # num_threads
        min_len,
        str(regions_out),
        str(fasta_out),
        1 # max_abundance
    )
    
    # Validate Results
    # Results is a list of ContigRegions
    # ContigRegions(header, sequence, regions)
    
    headers = {r.header: r for r in results}
    
    # Contig 1 should have NO exclusive regions because it is identical in target
    if "ref_contig_1" in headers:
         c1 = headers["ref_contig_1"]
         assert len(c1.regions) == 0, f"Contig 1 should not have exclusive regions, found: {c1.regions}"
         
    # Contig 2 should have exclusive regions (the core)
    assert "ref_contig_2" in headers
    c2 = headers["ref_contig_2"]
    assert len(c2.regions) > 0, "Contig 2 should have exclusive regions"
    
    # Verify the exclusive region contains part of the unique core
    found_unique = False
    for reg in c2.regions:
        # Region now carries its own subsequence — no need to slice from full sequence
        if "TCTGATAGCAGC" in reg.subsequence:
            found_unique = True
            break
    
    assert found_unique, "Did not find the expected unique core sequence in extracted regions"

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
