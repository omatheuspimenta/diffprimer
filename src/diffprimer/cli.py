from typing import Annotated

from rich.console import Console
from rich.text import Text
from typer import Exit, Option, Typer

from diffprimer import __version__
from diffprimer.main import main as _main

app = Typer(rich_markup_mode="rich")
console = Console()


def show_version(flag: bool) -> None:
    if flag:
        console.print(
            Text("diffprimer version: ", style="bold cyan")
            + Text(__version__, style="bold green")
        )
        raise Exit(code=0)


@app.command(epilog=f"diffprimer version: v{__version__}")
def main(
    reference_file: Annotated[
        str,
        Option(
            "--reference-file",
            "-r",
            help=(
                "Reference genome file in FASTA format.\n\n"
                "This sequence is used as the baseline for identifying "
                "exclusive regions across the input sequences."
            ),
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel="Input/Output",
        ),
    ],
    sequences_path: Annotated[
        str,
        Option(
            "--sequences-path",
            "-s",
            help=(
                "Directory containing FASTA files of sequences to compare "
                "against the reference.\n\n"
                "Each file should represent one genome or sequence set."
            ),
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            resolve_path=True,
            rich_help_panel="Input/Output",
        ),
    ],
    config_file: Annotated[
        str | None,
        Option(
            "--config-file",
            "-c",
            help=(
                "Primer3 configuration file used for primer design.\n\n"
                "This file defines parameters such as primer length, melting "
                "temperature, GC content, and other Primer3 settings.\n\n"
                "For more details on configuring Primer3, refer to the "
                "Primer3 documentation: https://primer3.org/manual.html"
            ),
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel="Input/Output",
        ),
    ] = None,
    annotation_path: Annotated[
        str | None,
        Option(
            "--annotation-path",
            "-a",
            help=(
                "Genome annotation file associated with the reference sequence.\n\n"
                "Typically provided in GFF3 format and used to map "
                "exclusive regions to annotated features."
            ),
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            rich_help_panel="Input/Output",
        ),
    ] = None,
    k: Annotated[
        int,
        Option(
            "--kmer-size",
            "-k",
            help=(
                "K-mer size used during sequence comparison.\n\n"
                "Default: 21."
            ),
            rich_help_panel="Configuration",
        ),
    ] = 21,
    cpus: Annotated[
        int | None,
        Option(
            "--cpus",
            "-p",
            help=(
                "Number of CPU cores to use for parallel processing.\n\n"
                "If not provided, all available CPUs will be used."
            ),
            rich_help_panel="Configuration",
        ),
    ] = None,
    min_region_length: Annotated[
        int,
        Option(
            "--min-region-length",
            "-m",
            help=(
                "Minimum length (in base pairs) of exclusive regions to consider.\n\n"
                "This filter applies only to exclusive regions and helps remove "
                "short, potentially non-informative segments. Default: 200."
            ),
            rich_help_panel="Configuration",
        ),
    ] = 200,
    reference_max_abundance: Annotated[
        int,
        Option(
            "--reference-max-abundance",
            help=(
                "Maximum allowed frequency of a k-mer in the reference genome "
                "to be considered a candidate marker.\n\n"
                "Use 1 for strictly unique markers (default). Increasing this value "
                "allows markers that are repeated a few times in the reference."
            ),
            rich_help_panel="Configuration",
        ),
    ] = 1,
    check_specificity: Annotated[
        bool,
        Option(
            "--check-specificity",
            help=(
                "Check if the designed primers are specific to the target region."
            ),
            rich_help_panel="Configuration",
        ),
    ] = False,
    similarity_threshold: Annotated[
        float,
        Option(
            "--similarity-threshold",
            help=(
                "Similarity threshold (as a percentage) for specificity checking.\n\n"
                "Primers with similarity above this threshold to non-target regions "
                "will be flagged. Default: 80.0."
            ),
            rich_help_panel="Configuration",
        ),
    ] = 80.0,
    local_mismatch_threshold: Annotated[
        int,
        Option(
            "--local-mismatch-threshold",
            help=(
                "Local mismatch score threshold for primer binding.\n\n"
                "Mismatches in the 5' end add the penalty to the score, while 3' end mismatches add the penalty from the penalty array. "
                "A score below this threshold means the primer is considered to bind (non-specific hit). "
                "The default is 7, which allows up to 2 mismatches in the 3' region, using default penalty values."
            ),
            rich_help_panel="Configuration",
        ),
    ] = 7,
    penalty_array: Annotated[
        str,
        Option(
            "--penalty",
            help=(
                "Penalty values for mismatches. Provide a comma-separated list of 6 values. "
                "The first value is the penalty outside the 3' region. "
                "The next 5 values are the penalties for the 5 nucleotides of the 3' region. "
                "Default: '1,3,3,3,3,3'."
            ),
            rich_help_panel="Configuration",
        ),
    ] = "1,3,3,3,3,3",
    version: Annotated[
        bool,
        Option(
            "--version",
            "-v",
            help="Show the version and exit.",
            is_eager=True,
            callback=show_version,
            rich_help_panel="Configuration",
        ),
    ] = False,
) -> None:
    """
    Run the diffprimer analysis pipeline.

    This command compares a reference genome against a collection of sequences
    to identify exclusive genomic regions and design primers targeting those regions.

    The workflow includes:
    - K-mer–based comparison of sequences
    - Identification of exclusive regions
    - Mapping regions to genome annotations
    - Primer design using Primer3

    Example usage:

        diffprimer \
            --reference-file reference.fasta \
            --sequences-path genomes/ \
            --annotation-path reference.gff3 \
            --config-file config.ini \
            --kmer-size 21 \
            --cpus 8

    All input files must be readable, and paths are resolved to absolute paths.
    """

    _main(
        reference_file=reference_file,
        sequences_path=sequences_path,
        annotation_path=annotation_path,
        config_file=config_file,
        k=k,
        cpus=cpus,
        min_region_length=min_region_length,
        reference_max_abundance=reference_max_abundance,
        check_specificity=check_specificity,
        similarity_threshold=similarity_threshold,
        local_mismatch_threshold=local_mismatch_threshold,
        penalty_array=penalty_array,
    )
