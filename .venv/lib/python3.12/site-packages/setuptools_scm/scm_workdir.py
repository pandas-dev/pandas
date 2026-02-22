from __future__ import annotations

import logging

from dataclasses import dataclass
from datetime import date
from datetime import datetime
from datetime import timezone
from pathlib import Path

from ._config import Configuration
from .version import ScmVersion

log = logging.getLogger(__name__)


def get_latest_file_mtime(changed_files: list[str], base_path: Path) -> date | None:
    """Get the latest modification time of the given files.

    Args:
        changed_files: List of relative file paths
        base_path: Base directory path to resolve relative paths

    Returns:
        The date of the most recently modified file, or None if no valid files found
    """
    if not changed_files or changed_files == [""]:
        return None

    latest_mtime = 0.0
    for filepath in changed_files:
        full_path = base_path / filepath
        try:
            file_stat = full_path.stat()
            latest_mtime = max(latest_mtime, file_stat.st_mtime)
        except OSError:
            # File might not exist or be accessible, skip it
            log.debug("Failed to get mtime for %s", full_path)
            continue

    if latest_mtime > 0:
        # Convert to UTC date
        dt = datetime.fromtimestamp(latest_mtime, timezone.utc)
        return dt.date()

    return None


@dataclass()
class Workdir:
    path: Path

    def run_describe(self, config: Configuration) -> ScmVersion:
        raise NotImplementedError(self.run_describe)
