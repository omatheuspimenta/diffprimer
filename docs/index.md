# Welcome to DiffPrimer

<p align="center">
  <img src="logo.png" alt="DiffPrimer Logo" width="300"/>
</p>

**DiffPrimer** is a high-performance bioinformatics pipeline engineered for the discovery of taxonomically exclusive genomic regions and the automated design of highly specific Polymerase Chain Reaction (PCR) diagnostic primers. 

By integrating alignment-free $k$-mer set arithmetic with the industry-standard Primer3 engine, DiffPrimer systematically isolates genetic markers unique to a target organism against a comprehensive background of off-target genomes.

## Key Features

- **Exclusive Region Discovery**: Rapidly identifies genomic subsequences present strictly within your target reference while being entirely absent across a designated background genomic database.
- **Automated Primer Design**: Seamlessly integrates with **Primer3** to generate thermodynamically optimal primer pairs tailored specifically for the discovered exclusive regions.
- **Rigorous Specificity Verification**:
  - Validates all candidate primer pairs exhaustively against every provided background genome.
  - Employs a robust two-stage verification algorithm, combining the **Myers Bit-Vector Algorithm** for rapid global homology filtering with **Needleman-Wunsch Semiglobal Alignment** to detect complex off-target hybridizations involving indels.
  - Implements a configurable **Positional Mismatch Penalty Array**, mathematically modeling the biochemical necessity of 3' terminal matching to prevent false-positive amplification.
- **Annotation Integration**: Cross-references discovered unique amplicons with standard GFF3 annotation files, directly mapping diagnostic markers to their corresponding gene products.
- **High-Performance Parallelization**: Leverages a highly concurrent Rust core, ensuring scalable and rapid execution across large-scale genomic datasets.

## Pipeline Architecture

The DiffPrimer computational workflow operates in three distinct phases:

1. **Alignment-Free Region Identification**: The pipeline scans the target reference genome, dynamically subtracting any $k$-mers shared with the background set. Contiguous blocks of unshared $k$-mers are assembled into continuous unique amplicons.
2. **Thermodynamic Primer Design**: For each assembled unique region, the system interfaces with Primer3 to compute and retrieve the optimal primer pairs adhering to strict thermodynamic constraints.
3. **In Silico Specificity Verification**: *(Highly Recommended)* The pipeline exhaustively tests the designed primers against the background genomes. By penalizing local mismatches, particularly at the 3' end, the system definitively tags each primer pair, ensuring that cross-reactive false positives are strictly filtered out.

## Installation

DiffPrimer is distributed as a Python package. We strongly recommend installing it within an isolated environment using standard package managers such as `pip` or [`uv`](https://github.com/astral-sh/uv).

```bash
# Standard installation via pip
pip install diffprimer
```

For advanced installation guidelines, environment configuration, and troubleshooting, please refer to the main repository [README](https://github.com/omatheuspimenta/diffprimer).

## Citation

> [!NOTE]  
> Official citation information is pending publication.

## Contributing

We welcome contributions from the computational biology community! Should you encounter bugs, require technical support, or wish to propose new features, please submit an issue via our [GitHub repository](https://github.com/omatheuspimenta/diffprimer). 

Pull requests are encouraged. Please ensure that all proposed changes are accompanied by appropriate test coverage.
