from __future__ import annotations

import os
from stat import S_IXGRP, S_IXOTH, S_IXUSR
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path


def make_exe(filename: Path) -> None:
    original_mode = filename.stat().st_mode
    levels = [S_IXUSR, S_IXGRP, S_IXOTH]
    for at in range(len(levels), 0, -1):
        try:
            mode = original_mode
            for level in levels[:at]:
                mode |= level
            filename.chmod(mode)
            break
        except OSError:
            continue


def set_tree(folder: Path, stat: int) -> None:
    for root, _, files in os.walk(str(folder)):
        for filename in files:
            os.chmod(os.path.join(root, filename), stat)


__all__ = (
    "make_exe",
    "set_tree",
)
