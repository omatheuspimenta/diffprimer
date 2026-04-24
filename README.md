# DiffPrimer

![diffprimer_logo](docs/logo.png)

**DiffPrimer** is a bioinformatics tool designed to identify unique genomic regions in a reference genome compared to a set of other "background" genomes and automatically design specific PCR primers for these markers.

It combines high-performance k-mer analysis (written in Rust) with standard primer design tools (Primer3) to generate diagnostic markers that are exclusive to your target organism.

## Key Features

-   **Exclusive Region Discovery**: Identifies genomic regions present in your reference but absent in a database of other genomes.
-   **Automated Primer Design**: Integrated with **Primer3** to design optimal primer pairs for identified unique regions.
-   **Specificity Check**: 
    -   Validates designed primers against *all* input background genomes.
    -   Uses a hybrid **Myers Bit-Vector Algorithm** (Global) and **Semiglobal Alignment** (Local) to detect potential off-target amplification, even with insertions/deletions.
-   **Annotation Integration**: Cross-references unique regions with GFF3 annotation files to identify which genes (if any) the markers overlap with.
-   **Parallel Processing**: Fully parallelized core for fast execution on large datasets.

---

## Installation

### Prerequisites

-   **Operating System**: Linux (Recommended) or macOS.
-   **Python**: Version **3.13** or higher (Strict requirement).
-   **Rust**: Required to build the core engine. Install via [rustup](https://rustup.rs/).

### Quick Install (Recommended)

We recommend using [`uv`](https://github.com/astral-sh/uv) or `pip` in a clean environment.

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/yourusername/diffprimer.git
    cd diffprimer
    ```

2.  **Create a Virtual Environment**:
    ```bash
    # Using uv (Recommended)
    uv venv --python 3.13
    source .venv/bin/activate
    
    # OR using standard python
    python3.13 -m venv .venv
    source .venv/bin/activate
    ```

3.  **Install**:
    ```bash
    uv pip install .
    # OR
    pip install .
    ```

### Troubleshooting (Conda/Anaconda Users)

If you are working inside a Conda environment, the system linker may conflict with Rust, causing build errors (e.g., `undefined symbol: getauxval`). To fix this:

1.  **Deactivate Conda** completely for the build step.
2.  Use the system's standard PATH.

```bash
# In your terminal
unset CONDA_PREFIX
export PATH=/usr/local/bin:/usr/bin:/bin:$HOME/.cargo/bin
source .venv/bin/activate
uv pip install --force-reinstall .
```

---

## Usage

Running `diffprimer` executes the full pipeline: finding unique regions -> designing primers -> verifying specificity.

### Basic Command

```bash
diffprimer \
    --reference-file reference.fasta \
    --sequences-path genomes_directory/ \
```

### Complete Example

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

### Arguments Explained

| Argument | Short | Description |
| :--- | :--- | :--- |
| `--reference-file` | `-r` | **Required.** Path to the target reference genome (FASTA). |
| `--sequences-path` | `-s` | **Required.** Directory containing background genomes (FASTA) to compare against. Regions found in these genomes will be excluded. |
| `--annotation-path` | `-a` | (Optional) GFF3 file for the reference. Used to annotate output regions with gene names + product info. |
| `--config-file` | `-c` | (Optional) Primer3 configuration file defining Tm, GC%, and size constraints. |
| `--reference-max-abundance` | | Maximum allowed frequency of a k-mer in the reference genome to be considered a candidate marker. Use 1 for strictly unique markers (default). Increasing this value allows markers that are repeated a few times in the reference. (Default: 1) |
| `--kmer-size` | `-k` | K-mer size for uniqueness check (Default: 21). |
| `--min-region-length` | `-m` | Minimum length of unique regions to keep (Default: 200 bp). |
| `--check-specificity` | | **Highly Recommended.** Enables the rigorous cryptographic check of primer specificity. Without this, primers are only designed on unique regions but not physically verified against off-targets. |
| `--similarity-threshold` | | Global similarity threshold (%) for off-target flagging (Default: 80.0). |
| `--local-mismatch-threshold` | | Positional mismatch score threshold for the local specificity check. Standard mismatches add 1, 3' end mismatches add 3 (configurable via `--penalty-array`). Default is 7, allowing up to 2 mismatches in the 3' region to still be considered non-specific (capable of amplifying). |
| `--penalty-array` | | Penalty values for mismatches. Provide a comma-separated list of 6 values. The first value is the penalty outside the 3' region. The next 5 values are the penalties for the 5 nucleotides of the 3' region. (Default: 1,3,3,3,3,3) |

---

## Output Format

The output is a CSV file containing one row per designed primer pair. Key columns include:

-   `Sequence_Header`: Contig/Region name.
-   `Region_Start`/`End`: Coordinates of the unique region.
-   `Forward_Primer` / `Reverse_Primer` (and `_Tm`, `_GC`, etc.): Primer details.
-   **`Specificity_Tag`**: The result of the specificity analysis (if `--check-specificity` is used).
    -   `Specific_LowGlobalSim`: The region is globally unique; off-targets are extremely dissimilar, so primers are specific.
    -   `Specific_PositionalMismatches`: The region has high global similarity to a background genome, BUT the primers *themselves* are specific because they have biological mismatches in critical locations (like the 3' end) that prevent binding.
    -   `NonSpecific_Amplification`: The primers bind perfectly or near-perfectly to a similar region in a background genome resulting in cross-reactivity. **Do NOT use.**
    -   `Not_Checked`: Specificity check was skipped.
-   `Gene_Name` / `Product`: Annotation info (if GFF3 provided).

## Configuration

You can customize [Primer3 settings](https://primer3.org/manual.html) by providing a file with `--config-file`. Example format:
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
