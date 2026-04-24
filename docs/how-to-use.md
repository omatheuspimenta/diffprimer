# How to Use DiffPrimer

DiffPrimer coordinates the entire computational pipeline: discovering taxonomic exclusive genomic regions, structurally designing optimal primers, and executing rigorous *in silico* specificity verifications.

## Basic Execution

To launch a standard analysis, provide the target reference sequence and the directory containing the background comparison sequences:

```bash
diffprimer \
    --reference-file reference.fasta \
    --sequences-path genomes_directory/
```

## Advanced Execution

For a comprehensive analysis encompassing gene annotation mapping, customized thermodynamic constraints, and complete specificity verification, use the following configuration:

```bash
diffprimer \
    --reference-file data/target_species.fasta \
    --sequences-path data/background_species/ \
    --annotation-path data/target_annotations.gff3 \
    --config-file primer3_config.ini \
    --min-region-length 200 \
    --cpus 8 \
    --check-specificity
```

[Download Example Files](#){ .md-button .md-button--primary }

> [!NOTE]
> The example files download link is currently a placeholder. Comprehensive datasets will be provided via the official GitHub repository prior to release.

## Command-Line Arguments

| Argument | Short | Description |
| :--- | :--- | :--- |
| `--reference-file` | `-r` | **Required.** The path to the target reference genome (in FASTA format) for which diagnostic markers are to be discovered. |
| `--sequences-path` | `-s` | **Required.** The path to the directory containing the background genomes (in FASTA format). The algorithm will explicitly exclude any $k$-mers found within these genomes. |
| `--annotation-path` | `-a` | *(Optional)* The path to a GFF3 annotation file corresponding to the reference genome. Enables structural annotation of the discovered regions, mapping them directly to gene locus and product information. |
| `--config-file` | `-c` | *(Optional)* A custom Primer3 configuration file to explicitly define thermodynamic parameters such as melting temperature ($T_m$), GC content limits, and amplicon size constraints. |
| `--reference-max-abundance` | | The maximum allowable intra-genomic frequency of a $k$-mer within the reference genome to be retained as a candidate marker. A value of `1` enforces strict single-copy unique markers. Elevating this threshold permits multi-copy markers. (Default: `1`) |
| `--kmer-size` | `-k` | The fundamental $k$-mer word size utilized for the alignment-free uniqueness subtraction. (Default: `21`) |
| `--min-region-length` | `-m` | The minimum contiguous length (in base pairs) required for a unique block of $k$-mers to be assembled into a valid candidate amplicon. (Default: `200`) |
| `--check-specificity` | | **Highly Recommended.** Triggers the comprehensive *in silico* specificity verification module. If omitted, primer pairs are designed strictly based on region uniqueness without secondary cross-hybridization checks against the background data. |
| `--similarity-threshold` | | The global sequence similarity percentage threshold utilized by the Myers bit-vector algorithm to flag highly homologous off-target candidate regions. (Default: `80.0`) |
| `--local-mismatch-threshold` | | The critical positional mismatch score boundary for localized specificity verification. Candidate primer alignments that accumulate a penalty score lower than this threshold are flagged as unsafe (non-specific). (Default: `7.0`) |
| `--penalty` | | A comma-separated numeric array of six discrete values dictating the severity of alignment mismatches. The primary value defines the default penalty for the 5' primer region, while the subsequent five values dictate penalties along the critical 3' terminal end. (Default: `1.0,3.0,3.0,3.0,3.0,3.0`) |

## Output Format

The execution outputs a structured CSV matrix detailing one row per successfully designed and verified primer pair. The core data columns include:

-   `Sequence_Header`: The specific contig or chromosome identifier from the reference genome.
-   `Region_Start`/`End`: The spatial coordinates demarcating the discovered unique genomic region.
-   `Forward_Primer` / `Reverse_Primer`: The generated oligonucleotide sequences, supplemented with thermodynamic statistics (`_Tm`, `_GC`).
-   **`Specificity_Tag`**: The ultimate classification assigned by the specificity verification engine (requires `--check-specificity`):
    -   `Specific_LowGlobalSim`: The amplicon diverges significantly from all background sequences. The primer pair is highly specific.
    -   `Specific_PositionalMismatches`: The amplicon displays background homology, but critical positional mismatches (primarily at the 3' terminus) biochemically inhibit off-target amplification. **Specific and safe for deployment.**
    -   `NonSpecific_Amplification`: The primer pair strongly hybridizes to a background sequence with insufficient mismatch penalties to prevent extension. **Cross-reactive; do not use.**
    -   `Not_Checked`: The specificity verification module was explicitly bypassed.
-   `Gene_Name` / `Product`: Structural annotation data extracted from the provided GFF3 file.

## Primer3 Thermodynamic Configuration

To customize the internal Primer3 thermodynamic engine, supply an initialization file via the `--config-file` argument. A standard parameterization template is formatted as follows:

```
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=27
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=57.0
PRIMER_MAX_TM=63.0
PRIMER_MIN_GC=20.0
PRIMER_MAX_GC=80.0
```
