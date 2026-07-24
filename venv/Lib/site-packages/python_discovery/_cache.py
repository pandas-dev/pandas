"""Cache Protocol and built-in implementations for Python interpreter discovery."""

from __future__ import annotations

import json
import logging
from contextlib import contextmanager, suppress
from hashlib import sha256
from typing import TYPE_CHECKING, Final, Protocol, runtime_checkable

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

_LOGGER: Final[logging.Logger] = logging.getLogger(__name__)


@runtime_checkable
class ContentStore(Protocol):
    """A store for reading and writing cached content."""

    def exists(self) -> bool:
        """Return whether the cached content exists."""
        ...

    def read(self) -> dict | None:
        """Read the cached content, or ``None`` if unavailable or corrupt."""
        ...

    def write(self, content: dict) -> None:
        """
        Persist *content* to the store.

        :param content: interpreter metadata to cache.
        """
        ...

    def remove(self) -> None:
        """Delete the cached content."""
        ...

    @contextmanager
    def locked(self) -> Generator[None]:
        """Context manager that acquires an exclusive lock on this store."""
        ...


@runtime_checkable
class PyInfoCache(Protocol):
    """Cache interface for Python interpreter information."""

    def py_info(self, path: Path) -> ContentStore:
        """
        Return the content store for the interpreter at *path*.

        :param path: absolute path to a Python executable.
        """
        ...

    def py_info_clear(self) -> None:
        """Remove all cached interpreter information."""
        ...


class DiskContentStore:
    """JSON file-based content store with file locking."""

    def __init__(self, folder: Path, key: str) -> None:
        self._folder = folder
        self._key = key

    @property
    def _file(self) -> Path:
        return self._folder / f"{self._key}.json"

    def exists(self) -> bool:
        return self._file.exists()

    def read(self) -> dict | None:
        data, bad_format = None, False
        try:
            data = json.loads(self._file.read_text(encoding="utf-8"))
        except ValueError:
            bad_format = True
        except OSError:
            _LOGGER.debug("failed to read %s", self._file, exc_info=True)
        else:
            _LOGGER.debug("got python info from %s", self._file)
            return data
        if bad_format:
            with suppress(OSError):
                self.remove()
        return None

    def write(self, content: dict) -> None:
        self._folder.mkdir(parents=True, exist_ok=True)
        self._file.write_text(json.dumps(content, sort_keys=True, indent=2), encoding="utf-8")
        _LOGGER.debug("wrote python info at %s", self._file)

    def remove(self) -> None:
        with suppress(OSError):
            self._file.unlink()
        _LOGGER.debug("removed python info at %s", self._file)

    @contextmanager
    def locked(self) -> Generator[None]:
        from filelock import FileLock  # ruff:ignore[import-outside-top-level]

        lock_path = self._folder / f"{self._key}.lock"
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        with FileLock(str(lock_path)):
            yield


class DiskCache:
    """
    File-system based Python interpreter info cache (``<root>/py_info/4/<sha256>.json``).

    :param root: root directory for the on-disk cache.
    """

    def __init__(self, root: Path) -> None:
        self._root = root

    @property
    def _py_info_dir(self) -> Path:
        return self._root / "py_info" / "4"

    def py_info(self, path: Path) -> DiskContentStore:
        """
        Return the content store for the interpreter at *path*.

        :param path: absolute path to a Python executable.
        """
        key = sha256(str(path).encode("utf-8")).hexdigest()
        return DiskContentStore(self._py_info_dir, key)

    def py_info_clear(self) -> None:
        """Remove all cached interpreter information."""
        folder = self._py_info_dir
        if folder.exists():
            for entry in folder.iterdir():
                if entry.suffix == ".json":
                    with suppress(OSError):
                        entry.unlink()


class NoOpContentStore(ContentStore):
    """Content store that does nothing -- implements ContentStore protocol."""

    def exists(self) -> bool:  # ruff:ignore[no-self-use]
        return False

    def read(self) -> dict | None:  # ruff:ignore[no-self-use]
        return None

    def write(self, content: dict) -> None:
        pass

    def remove(self) -> None:
        pass

    @contextmanager
    def locked(self) -> Generator[None]:  # ruff:ignore[no-self-use]
        yield


class NoOpCache(PyInfoCache):
    """Cache that does nothing -- implements PyInfoCache protocol."""

    def py_info(self, path: Path) -> NoOpContentStore:  # ruff:ignore[unused-method-argument, no-self-use]
        return NoOpContentStore()

    def py_info_clear(self) -> None:
        pass


__all__ = [
    "ContentStore",
    "DiskCache",
    "DiskContentStore",
    "NoOpCache",
    "NoOpContentStore",
    "PyInfoCache",
]
