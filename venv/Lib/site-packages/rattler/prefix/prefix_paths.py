from __future__ import annotations
import os
from typing import List, TYPE_CHECKING, Literal, Optional

from rattler.package.paths_json import FileMode
from rattler.rattler import PyPrefixPaths, PyPrefixPathsEntry, PyPrefixPathType

# `os.PathLike` started to be generic from Python 3.9
if TYPE_CHECKING:
    BasePathLike = os.PathLike[str]
else:
    BasePathLike = os.PathLike


class PrefixPathType:
    _inner: PyPrefixPathType

    def __init__(
        self,
        path_type: Literal[
            "hardlink",
            "softlink",
            "directory",
            "pyc_file",
            "windows_python_entry_point_script",
            "windows_python_entry_point_exe",
            "unix_python_entry_point",
        ],
    ) -> None:
        """
        Create a new PrefixPathType instance.

        Parameters
        ----------
        path_type : str
            The type of path. Must be one of: "hardlink", "softlink", "directory"

        Examples
        --------
        ```python
        >>> path_type = PrefixPathType("hardlink")
        >>> path_type.hardlink
        True
        >>>
        ```
        """
        self._inner = PyPrefixPathType(path_type)

    @classmethod
    def from_py_path_type(cls, py_path_type: PyPrefixPathType) -> PrefixPathType:
        """Construct Rattler PathType from FFI PyPathType object."""
        path_type = cls.__new__(cls)
        path_type._inner = py_path_type
        return path_type

    @property
    def hardlink(self) -> bool:
        """
        Whether the path should be hardlinked (the default) (once installed)
        """
        return self._inner.hardlink

    @property
    def softlink(self) -> bool:
        """
        Whether the path should be softlinked (once installed)
        """
        return self._inner.softlink

    @property
    def directory(self) -> bool:
        """
        This is a directory
        """
        return self._inner.directory

    @property
    def pyc_file(self) -> bool:
        """
        This is a file compiled from Python code when a noarch package was installed
        """
        return self._inner.pyc_file

    @property
    def windows_python_entry_point_script(self) -> bool:
        """
        A Windows entry point python script (a <entrypoint>-script.py Python script file)
        """
        return self._inner.windows_python_entry_point_script

    @property
    def windows_python_entry_point_exe(self) -> bool:
        """
        A Windows entry point python script (a <entrypoint>.exe executable)
        """
        return self._inner.windows_python_entry_point_exe

    @property
    def unix_python_entry_point(self) -> bool:
        """
        A Unix entry point python script (a <entrypoint> Python script file)
        """
        return self._inner.unix_python_entry_point


class PrefixPathsEntry(BasePathLike):
    _inner: PyPrefixPathsEntry

    def __init__(
        self,
        relative_path: os.PathLike[str],
        path_type: PrefixPathType,
        prefix_placeholder: Optional[str] = None,
        file_mode: Optional[FileMode] = None,
        sha256: Optional[bytes] = None,
        sha256_in_prefix: Optional[bytes] = None,
        size_in_bytes: Optional[int] = None,
        original_path: Optional[os.PathLike[str]] = None,
    ) -> None:
        """
        Create a new PrefixPathsEntry instance.

        Parameters
        ----------
        relative_path : os.PathLike[str]
            The relative path from the root of the package
        path_type : PrefixPathType
            Determines how to include the file when installing the package
        prefix_placeholder : Optional[str], optional
            The placeholder prefix used in the file, by default None
        file_mode : Optional[FileMode], optional
            The file mode of the path, by default None
        sha256 : Optional[bytes], optional
            The sha256 of the path, by default None
        sha256_in_prefix : Optional[bytes], optional
            The sha256 of the path in the prefix, by default None
        size_in_bytes : Optional[int], optional
            The size of the path in bytes, by default None
        original_path : Optional[os.PathLike[str]], optional
            The original path of the file, by default None
        """
        self._inner = PyPrefixPathsEntry(
            relative_path,
            path_type._inner,
            prefix_placeholder,
            file_mode._inner if file_mode else None,
            sha256,
            sha256_in_prefix,
            size_in_bytes,
            original_path,
        )

    def __fspath__(self) -> str:
        return str(self._inner.path)

    @classmethod
    def _from_py_paths_entry(cls, py_paths_entry: PyPrefixPathsEntry) -> PrefixPathsEntry:
        """Construct Rattler PathsEntry from FFI PyPathsEntry object."""
        entry = cls.__new__(cls)
        entry._inner = py_paths_entry
        return entry

    @property
    def relative_path(self) -> os.PathLike[str]:
        """
        The relative path of the entry.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/tk-8.6.12-h8ffe710_0.json"
        ... )
        >>> paths = r.paths_data
        >>> relative_path = paths.paths[0].relative_path
        >>>
        ...
        ```
        """
        return self._inner.relative_path

    @relative_path.setter
    def relative_path(self, path: os.PathLike[str]) -> None:
        self._inner.set_relative_path(path)

    @property
    def no_link(self) -> bool:
        """
        Whether this file should be linked

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> no_link = paths.paths[0].no_link
        >>>
        ```
        """

    @no_link.setter
    def no_link(self, no_link: bool) -> None:
        self._inner.set_no_link(no_link)

    @property
    def path_type(self) -> PrefixPathType:
        """
        The type of the path.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> path_type = paths.paths[0].path_type
        >>>
        ```
        """
        return PrefixPathType.from_py_path_type(self._inner.path_type)

    @path_type.setter
    def path_type(self, path_type: PrefixPathType) -> None:
        self._inner.set_path_type(path_type._inner)

    @property
    def prefix_placeholder(self) -> str | None:
        """
        The prefix placeholder for the path.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> prefix_placeholder = paths.paths[0].prefix_placeholder
        >>>
        ```
        """
        return self._inner.prefix_placeholder

    @prefix_placeholder.setter
    def prefix_placeholder(self, placeholder: Optional[str]) -> None:
        self._inner.set_prefix_placeholder(placeholder)

    @property
    def file_mode(self) -> FileMode:
        """
        The file mode of the path.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> file_mode = paths.paths[0].file_mode
        >>>
        ```
        """
        return FileMode._from_py_file_mode(self._inner.file_mode)

    @file_mode.setter
    def file_mode(self, file_mode: Optional[FileMode]) -> None:
        self._inner.set_file_mode(file_mode._inner if file_mode else None)

    @property
    def sha256(self) -> bytes:
        """
        The sha256 of the path.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> sha256 = paths.paths[0].sha256
        >>>
        ```
        """
        return self._inner.sha256

    @sha256.setter
    def sha256(self, sha256: Optional[bytes]) -> None:
        self._inner.set_sha256(sha256)

    @property
    def sha256_in_prefix(self) -> bytes:
        """
        The sha256 of the path in the prefix.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> sha256_in_prefix = paths.paths[0].sha256_in_prefix
        >>>
        ```
        """
        return self._inner.sha256_in_prefix

    @sha256_in_prefix.setter
    def sha256_in_prefix(self, sha256: Optional[bytes]) -> None:
        self._inner.set_sha256_in_prefix(sha256)

    @property
    def size_in_bytes(self) -> int:
        """
        The size of the path in bytes.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> size_in_bytes = paths.paths[0].size_in_bytes
        >>>
        ```
        """
        return self._inner.size_in_bytes

    @size_in_bytes.setter
    def size_in_bytes(self, size: Optional[int]) -> None:
        self._inner.set_size_in_bytes(size)


class PrefixPaths:
    _paths: PyPrefixPaths

    @classmethod
    def _from_py_prefix_paths(cls, py_prefix_paths: PyPrefixPaths) -> PrefixPaths:
        """Construct Rattler PrefixRecord from FFI PyPrefixRecord object."""
        paths = cls.__new__(cls)
        paths._paths = py_prefix_paths
        return paths

    def __init__(self, paths_version: int = 1) -> None:
        """
        Create a new PrefixPaths instance.

        Parameters
        ----------
        paths_version : int, optional
            The version of the paths file format, by default 1

        Examples
        --------
        ```python
        >>> paths = PrefixPaths()
        >>> paths.paths_version
        1
        >>>
        ```
        """
        self._paths = PyPrefixPaths(paths_version)

    @property
    def paths_version(self) -> int:
        """
        The version of the file.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> paths.paths_version
        1
        >>>
        ```
        """
        return self._paths.paths_version

    @paths_version.setter
    def paths_version(self, version: int) -> None:
        self._paths.paths_version = version

    @property
    def paths(self) -> List[PrefixPathsEntry]:
        """
        All entries included in the package.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> paths = r.paths_data
        >>> paths.paths # doctest:+ELLIPSIS
        [...]
        >>>
        ```
        """
        return [PrefixPathsEntry._from_py_paths_entry(path) for path in self._paths.paths]

    @paths.setter
    def paths(self, paths: List[PrefixPathsEntry]) -> None:
        self._paths.paths = [path._inner for path in paths]

    def __repr__(self) -> str:
        """
        Returns a representation of the PrefixPaths.

        Examples
        --------
        ```python
        >>> from rattler.prefix.prefix_record import PrefixRecord
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> r.paths_data # doctest:+ELLIPSIS
        PrefixPaths(paths=[...])
        >>>
        ```
        """
        return f"PrefixPaths(paths={self.paths})"
