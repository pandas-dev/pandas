from __future__ import annotations

import os.path
from typing import TYPE_CHECKING

from virtualenv.util.lock import NoOpFileLock

from .via_disk_folder import AppDataDiskFolder, PyInfoStoreDisk

if TYPE_CHECKING:
    from pathlib import Path
    from typing import NoReturn


class ReadOnlyAppData(AppDataDiskFolder):
    can_update = False

    def __init__(self, folder: str) -> None:
        if not os.path.isdir(folder):
            msg = f"read-only app data directory {folder} does not exist"
            raise RuntimeError(msg)
        super().__init__(folder)
        self.lock = NoOpFileLock(folder)

    def reset(self) -> None:
        msg = "read-only app data does not support reset"
        raise RuntimeError(msg)

    def py_info_clear(self) -> None:
        raise NotImplementedError

    def py_info(self, path: Path) -> _PyInfoStoreDiskReadOnly:
        return _PyInfoStoreDiskReadOnly(self.py_info_at, path)

    def embed_update_log(self, distribution: str, for_py_version: str) -> NoReturn:
        raise NotImplementedError


class _PyInfoStoreDiskReadOnly(PyInfoStoreDisk):
    def write(self, content: str) -> NoReturn:  # ruff:ignore[unused-method-argument]
        msg = "read-only app data python info cannot be updated"
        raise RuntimeError(msg)


__all__ = [
    "ReadOnlyAppData",
]
