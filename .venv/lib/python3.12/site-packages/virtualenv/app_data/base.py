"""Application data stored by virtualenv."""

from __future__ import annotations

from abc import ABC, abstractmethod
from contextlib import contextmanager
from typing import TYPE_CHECKING

from virtualenv.info import IS_ZIPAPP

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path


class AppData(ABC):
    """Abstract storage interface for the virtualenv application."""

    @abstractmethod
    def close(self) -> None:
        """Called before virtualenv exits."""

    @abstractmethod
    def reset(self) -> None:
        """Called when the user passes in the reset app data."""

    @abstractmethod
    def py_info(self, path: Path) -> ContentStore:
        """Return a content store for cached interpreter information at the given path.

        :param path: the interpreter executable path

        :returns: a content store for the cached data

        """
        raise NotImplementedError

    @abstractmethod
    def py_info_clear(self) -> None:
        """Clear all cached interpreter information."""
        raise NotImplementedError

    @property
    def can_update(self) -> bool:
        """``True`` if this app data store supports updating cached content."""
        raise NotImplementedError

    @abstractmethod
    def embed_update_log(self, distribution: str, for_py_version: str) -> ContentStore:
        """Return a content store for the embed update log of a distribution.

        :param distribution: the package name (e.g. ``pip``)
        :param for_py_version: the target Python version string

        :returns: a content store for the update log

        """
        raise NotImplementedError

    @property
    def house(self) -> Path:
        """The root directory of the application data store."""
        raise NotImplementedError

    @property
    def transient(self) -> bool:
        """``True`` if this app data store is transient and does not persist across runs."""
        raise NotImplementedError

    @abstractmethod
    def wheel_image(self, for_py_version: str, name: str) -> Path:
        """Return the path to a cached wheel image.

        :param for_py_version: the target Python version string
        :param name: the package name

        :returns: the path to the cached wheel

        """
        raise NotImplementedError

    @contextmanager
    def ensure_extracted(self, path: Path, to_folder: Path | None = None) -> Generator[Path]:
        """Ensure a path is available on disk, extracting from zipapp if needed.

        :param path: the path to ensure is available
        :param to_folder: optional target directory for extraction

        :returns: yields the usable path on disk

        """
        if IS_ZIPAPP:
            with self.extract(path, to_folder) as result:
                yield result
        else:
            yield path

    @abstractmethod
    @contextmanager
    def extract(self, path: Path, to_folder: Path | None) -> Generator[Path]:
        """Extract a path from the zipapp to a location on disk.

        :param path: the path to extract
        :param to_folder: optional target directory

        :returns: yields the extracted path

        """
        raise NotImplementedError

    @abstractmethod
    @contextmanager
    def locked(self, path: Path) -> Generator[None]:
        """Acquire an exclusive lock on the given path.

        :param path: the path to lock

        """
        raise NotImplementedError


class ContentStore(ABC):
    """A store for reading and writing cached content."""

    @abstractmethod
    def exists(self) -> bool:
        """Check if the stored content exists.

        :returns: ``True`` if content exists

        """
        raise NotImplementedError

    @abstractmethod
    def read(self) -> str:
        """Read the stored content.

        :returns: the stored content as a string

        """
        raise NotImplementedError

    @abstractmethod
    def write(self, content: str) -> None:
        """Write content to the store.

        :param content: the content to write

        """
        raise NotImplementedError

    @abstractmethod
    def remove(self) -> None:
        """Remove the stored content."""
        raise NotImplementedError

    @abstractmethod
    @contextmanager
    def locked(self) -> Generator[None]:
        """Acquire an exclusive lock on this content store."""


__all__ = [
    "AppData",
    "ContentStore",
]
