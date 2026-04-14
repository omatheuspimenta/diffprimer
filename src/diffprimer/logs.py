"""Log class for diffprimer.
Provides logging functionalities for the diffprimer application.
"""

import logging

from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
    MofNCompleteColumn,
    TransferSpeedColumn,
)

# ---------------------------------------------------------------------------
# Logging configuration
# ---------------------------------------------------------------------------
FORMAT = "%(message)s"
logging.basicConfig(
    format=FORMAT,
    level="INFO",
    handlers=[
        RichHandler(
            show_time=False,      # Rust output has no timestamps; keep both sides uniform
            show_level=False,     # Rust output has no level badges
            show_path=False,
            markup=True,
            rich_tracebacks=True,
            tracebacks_show_locals=False,
        )
    ],
)

logtext = logging.getLogger("rich")
logtext.setLevel(logging.INFO)


def _make_progress(total: bool = False) -> Progress:
    """
    Create a configured Rich Progress instance.

    Args:
        total (bool): If True, shows M-of-N completion and a rate column
                      (used for primer design, specificity check, CSV writing).
                      If False, uses a simpler bar (scanning, loading).

    Returns:
        Progress: A configured Rich Progress object.
    """
    # Single, unobtrusive spinner that works well in scientific contexts
    SPINNER = "line"

    if total:
        return Progress(
            SpinnerColumn(SPINNER, style="cyan"),
            TextColumn("[bold cyan]{task.description}[/bold cyan]"),
            TimeElapsedColumn(),   # elapsed BEFORE bar -- matches Rust layout
            BarColumn(bar_width=35, style="cyan", complete_style="green"),
            MofNCompleteColumn(),
            TimeRemainingColumn(),
            TransferSpeedColumn(),
            transient=True,
        )

    return Progress(
        SpinnerColumn(SPINNER, style="cyan"),
        TextColumn("[bold cyan]{task.description}[/bold cyan]"),
        TimeElapsedColumn(),       # elapsed BEFORE bar -- matches Rust layout
        BarColumn(bar_width=35, style="cyan", complete_style="green"),
        TimeRemainingColumn(),
        transient=True,
    )


class diffprimerLog:
    """diffprimer logging class.
    Provides methods for logging messages at different levels.
    """

    @staticmethod
    def info(message: str) -> None:
        """Log an info message."""
        logtext.info(message)

    @staticmethod
    def warning(message: str) -> None:
        """Log a warning message."""
        logtext.warning(message)

    @staticmethod
    def error(message: str) -> None:
        """Log an error message."""
        logtext.error(message)
