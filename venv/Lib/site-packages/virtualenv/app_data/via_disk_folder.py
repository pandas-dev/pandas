r"""A rough layout of the current storage goes as:

::

    virtualenv-app-data
    ├── py - <version> <cache information about python interpreters>
    │  └── *.json/lock
    ├── wheel <cache wheels used for seeding>
    │   ├── house
    │   │   └── *.whl <wheels downloaded go here>
    │   └── <python major.minor> -> 3.9
    │       ├── img-<version>
    │       │   └── image
    │       │           └── <install class> -> CopyPipInstall / SymlinkPipInstall
    │       │               └── <wheel name> -> pip-20.1.1-py2.py3-none-any
    │       └── embed
    │           └── 3 -> json format versioning
    │               └── *.json -> for every distribution contains data about newer embed versions and releases
    └─── unzip <in zip app we cannot refer to some internal files, so first extract them>
         └── <virtualenv version>
             ├── py_info.py
             ├── debug.py
             └── _virtualenv.py

"""  # ruff:ignore[missing-terminal-punctuation]

from __future__ import annotations

import json
import logging
from abc import ABC
from contextlib import contextmanager, suppress
from hashlib import sha256
from typing import TYPE_CHECKING, Any

from virtualenv.util.lock import ReentrantFileLock
from virtualenv.util.path import safe_delete
from virtualenv.util.zipapp import extract
from virtualenv.version import __version__

from .base import AppData, ContentStore

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

LOGGER = logging.getLogger(__name__)


class AppDataDiskFolder(AppData):
    """Store the application data on the disk within a folder layout."""

    transient = False
    can_update = True

    def __init__(self, folder: str) -> None:
        self.lock = ReentrantFileLock(folder)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.lock.path})"

    def __str__(self) -> str:
        return str(self.lock.path)

    def reset(self) -> None:
        LOGGER.debug("reset app data folder %s", self.lock.path)
        safe_delete(self.lock.path)

    def close(self) -> None:
        """Do nothing."""

    @contextmanager
    def locked(self, path: Path) -> Generator[None]:
        path_lock = self.lock / path  # ty: ignore[unsupported-operator]
        with path_lock:
            yield path_lock.path

    @contextmanager
    def extract(self, path: Path, to_folder: Path | None) -> Generator[Path]:
        root = ReentrantFileLock(to_folder()) if to_folder is not None else self.lock / "unzip" / __version__  # ty: ignore[call-non-callable]
        with root.lock_for_key(path.name):
            dest = root.path / path.name
            if not dest.exists():
                extract(path, dest)
            yield dest

    @property
    def py_info_at(self) -> ReentrantFileLock:
        return self.lock / "py_info" / "5"  # ty: ignore[invalid-return-type]

    def py_info(self, path: Path) -> PyInfoStoreDisk:
        return PyInfoStoreDisk(self.py_info_at, path)

    def py_info_clear(self) -> None:
        """clear py info."""
        py_info_folder = self.py_info_at
        with py_info_folder:
            for filename in py_info_folder.path.iterdir():
                if filename.suffix == ".json":
                    with py_info_folder.lock_for_key(filename.stem):
                        if filename.exists():
                            filename.unlink()

    def embed_update_log(self, distribution: str, for_py_version: str) -> EmbedDistributionUpdateStoreDisk:
        return EmbedDistributionUpdateStoreDisk(self.lock / "wheel" / for_py_version / "embed" / "3", distribution)  # ty: ignore[invalid-argument-type]

    @property
    def house(self) -> Path:
        path = self.lock.path / "wheel" / "house"
        path.mkdir(parents=True, exist_ok=True)
        return path

    def wheel_image(self, for_py_version: str, name: str) -> Path:
        return self.lock.path / "wheel" / for_py_version / "image" / "1" / name


class JSONStoreDisk(ContentStore, ABC):
    def __init__(self, in_folder: ReentrantFileLock, key: str, msg_args: tuple[str, ...]) -> None:
        self.in_folder = in_folder
        self.key = key
        self.msg_args = (*msg_args, self.file)

    @property
    def file(self) -> Path:
        return self.in_folder.path / f"{self.key}.json"

    def exists(self) -> bool:
        return self.file.exists()

    def read(self) -> Any:  # ruff:ignore[any-type]
        data, bad_format = None, False
        try:
            data = json.loads(self.file.read_text(encoding="utf-8"))
        except ValueError:
            bad_format = True
        except Exception:  # ruff:ignore[blind-except, try-except-pass]
            pass
        else:
            LOGGER.debug("got %s %s from %s", *self.msg_args)
            return data
        if bad_format:
            with suppress(OSError):  # reading and writing on the same file may cause race on multiple processes
                self.remove()
        return None

    def remove(self) -> None:
        self.file.unlink()
        LOGGER.debug("removed %s %s at %s", *self.msg_args)

    @contextmanager
    def locked(self) -> Generator[None]:
        with self.in_folder.lock_for_key(self.key):
            yield

    def write(self, content: Any) -> None:  # ruff:ignore[any-type]
        folder = self.file.parent
        folder.mkdir(parents=True, exist_ok=True)
        self.file.write_text(json.dumps(content, sort_keys=True, indent=2), encoding="utf-8")
        LOGGER.debug("wrote %s %s at %s", *self.msg_args)


class PyInfoStoreDisk(JSONStoreDisk):
    def __init__(self, in_folder: ReentrantFileLock, path: Path) -> None:
        key = sha256(str(path).encode("utf-8")).hexdigest()
        super().__init__(in_folder, key, ("python info of", path))  # ty: ignore[invalid-argument-type]


class EmbedDistributionUpdateStoreDisk(JSONStoreDisk):
    def __init__(self, in_folder: ReentrantFileLock, distribution: str) -> None:
        super().__init__(
            in_folder,
            distribution,
            ("embed update of distribution", distribution),
        )


__all__ = [
    "AppDataDiskFolder",
    "JSONStoreDisk",
    "PyInfoStoreDisk",
]
