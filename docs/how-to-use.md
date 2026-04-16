# How to Use DiffPrimer

The main command is `run`. This executes the full pipeline: finding unique regions -> designing primers -> verifying specificity.

## Basic Usage

```bash
diffprimer run \
    --reference-file reference.fasta \
    --sequences-path genomes_directory/ \
```

## Complete Example

```bash
diffprimer run \
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
> The example files download link is currently a placeholder. GitHub will provide example files in the future.

## Arguments Explained

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
| `--local-mismatch-threshold` | | The score threshold for local mismatch positional penalty. Below this threshold, an off-target is considered non-specific. Penalty is 3.0 at 3' end and 1.0 elsewhere. (Default: 7.0) |

## Output Format

The output is a CSV file containing one row per designed primer pair. Key columns include:

-   `Sequence_Header`: Contig/Region name.
-   `Region_Start`/`End`: Coordinates of the unique region.
-   `Forward_Primer` / `Reverse_Primer` (and `_Tm`, `_GC`, etc.): Primer details.
-   **`Specificity_Tag`**: The result of the specificity analysis (if `--check-specificity` is used).
    -   `Specific_LowGlobalSim`: The region is globally unique; primers are safe.
    -   `Specific_PositionalMismatches`: The region has some similarity to background, but the primers *themselves* are specific due to sufficient positional mismatch penalties. **Safe to use.**
    -   `NonSpecific_Amplification`: The primers bind closely to a similar region in a background genome (mismatch score < threshold). **Do NOT use.**
    -   `Not_Checked`: Specificity check was skipped.
-   `Gene_Name` / `Product`: Annotation info (if GFF3 provided).

## Primer3 Configuration

You can customize Primer3 settings by providing a file with `--config-file`. Example format:

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

