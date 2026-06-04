from __future__ import annotations

import logging

from rich.console import Console
from rich.logging import RichHandler

CONSOLE = Console(stderr=True)
LOG = logging.getLogger("ivsBLASTn")


def setup_logging(verbose: bool) -> None:
    """Configure rich colored logging."""

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=CONSOLE, rich_tracebacks=True, show_path=False)],
    )
