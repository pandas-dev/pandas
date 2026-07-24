from __future__ import annotations

import logging
import os
import shutil
import sys
from stat import S_IWUSR
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

LOGGER = logging.getLogger(__name__)


def ensure_dir(path: Path) -> None:
    if not path.exists():
        LOGGER.debug("create folder %s", path)
        os.makedirs(str(path))


def ensure_safe_to_do(src: Path, dest: Path) -> None:
    if src == dest:
        msg = f"source and destination is the same {src}"
        raise ValueError(msg)
    if not dest.exists():
        return
    if dest.is_dir() and not dest.is_symlink():
        LOGGER.debug("remove directory %s", dest)
        safe_delete(dest)
    else:
        LOGGER.debug("remove file %s", dest)
        dest.unlink()


def symlink(src: Path, dest: Path) -> None:
    ensure_safe_to_do(src, dest)
    LOGGER.debug("symlink %s", _Debug(src, dest))
    dest.symlink_to(src, target_is_directory=src.is_dir())


def copy(src: Path, dest: Path) -> None:
    ensure_safe_to_do(src, dest)
    is_dir = src.is_dir()
    method = copytree if is_dir else shutil.copy
    LOGGER.debug("copy %s", _Debug(src, dest))
    method(str(src), str(dest))


def copytree(src: str, dest: str) -> None:
    for root, _, files in os.walk(src):
        dest_dir = os.path.join(dest, os.path.relpath(root, src))
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        for name in files:
            src_f = os.path.join(root, name)
            dest_f = os.path.join(dest_dir, name)
            shutil.copy(src_f, dest_f)


def safe_delete(dest: Path) -> None:
    def onerror(func: object, path: str, exc_info: object) -> None:  # ruff:ignore[unused-function-argument]
        if not os.access(path, os.W_OK):
            os.chmod(path, S_IWUSR)
            func(path)  # ty: ignore[call-non-callable]
        else:
            raise  # ruff:ignore[misplaced-bare-raise]

    if sys.version_info >= (3, 12):
        shutil.rmtree(str(dest), ignore_errors=True, onexc=onerror)
    else:
        shutil.rmtree(str(dest), ignore_errors=True, onerror=onerror)


class _Debug:
    def __init__(self, src: Path, dest: Path) -> None:
        self.src = src
        self.dest = dest

    def __str__(self) -> str:
        return f"{'directory ' if self.src.is_dir() else ''}{self.src!s} to {self.dest!s}"


__all__ = [
    "copy",
    "copytree",
    "ensure_dir",
    "safe_delete",
    "symlink",
]
