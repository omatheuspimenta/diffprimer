"""Log class for diffprimer.
Provides logging functionalities for the diffprimer application.
"""

import logging
from random import choice

from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    MofNCompleteColumn,
    RenderableColumn,
)

# Logging configuration
FORMAT = "%(message)s"
# format="[%(levelname)s]: "
logging.basicConfig(
    format=FORMAT,
    level="INFO",
    handlers=[RichHandler(show_time=True, show_path=False, markup=True)],
)

logtext = logging.getLogger("rich")
logtext.setLevel(20)


def _make_progress(total: bool = False) -> Progress:
    """
    Create a configured Rich Progress instance.

    Args:
        total (bool): If True, includes a column for M of N completion (checking specificity, etc).
                      If False, uses a simpler progress bar (scanning, etc).

    Returns:
        Progress: A configured Rich Progress object.
    """
    spinners = [
        "aesthetic",
        "shark",
        "dots",
        "line",
        "bouncingBall",
        "moon",
        "earth",
        "monkey",
        "runner",
        "pong",
        "weather",
        "clock",
    ]
    if total:
        return Progress(
            SpinnerColumn(choice(spinners)),
            TaskProgressColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            RenderableColumn(),
            transient=True,
        )
    return Progress(
        SpinnerColumn(choice(spinners)),
        TaskProgressColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
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
