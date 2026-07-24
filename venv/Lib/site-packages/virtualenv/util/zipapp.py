from __future__ import annotations

import logging
import zipfile
from pathlib import Path

from virtualenv.info import ROOT

LOGGER = logging.getLogger(__name__)


def read(full_path: str | Path) -> str:
    sub_file = _get_path_within_zip(full_path)
    with zipfile.ZipFile(ROOT, "r") as zip_file, zip_file.open(sub_file) as file_handler:
        return file_handler.read().decode("utf-8")


def extract(full_path: str | Path, dest: Path) -> None:
    LOGGER.debug("extract %s to %s", full_path, dest)
    sub_file = _get_path_within_zip(full_path)
    with zipfile.ZipFile(ROOT, "r") as zip_file:
        info = zip_file.getinfo(sub_file)
        info.filename = dest.name
        zip_file.extract(info, str(dest.parent))


def _get_path_within_zip(full_path: str | Path) -> str:
    # Use Path.relative_to so symlinks and ``..`` segments cannot slip through a string ``startswith`` check. The zipapp
    # root is a real file we own so ``resolve`` is safe; anything that does not resolve under ROOT is a bug or an
    # attempt to escape the archive and we refuse it.
    resolved = Path(full_path).resolve()
    root = Path(ROOT).resolve()
    try:
        relative = resolved.relative_to(root)
    except ValueError as exc:
        msg = f"full_path={resolved} should be within ROOT={root}"
        raise RuntimeError(msg) from exc
    # Zip entries always use forward slashes regardless of platform.
    return relative.as_posix()


__all__ = [
    "extract",
    "read",
]
