from __future__ import annotations
import os
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Literal
from rattler.rattler import (
    PyPathsJson,
    PyPathsEntry,
    PyPathType,
    PyPrefixPlaceholder,
    PyFileMode,
)

if TYPE_CHECKING:
    from rattler.networking.client import Client


class PathsJson:
    """
    A representation of the `paths.json` file found in package archives.
    The `paths.json` file contains information about every file included with the package.
    """

    _inner: PyPathsJson

    @staticmethod
    def from_package_archive(path: os.PathLike[str]) -> PathsJson:
        """
        Parses the package file from archive.
        Note: If you want to extract multiple `info/*` files then this will be slightly
              slower than manually iterating over the archive entries with
              custom logic as this skips over the rest of the archive


        """
        return PathsJson._from_py_paths_json(PyPathsJson.from_package_archive(path))

    @staticmethod
    def from_path(path: os.PathLike[str]) -> PathsJson:
        """
        Parses the object from a file specified by a `path`, using a format
        appropriate for the file type.

        For example, if the file is in JSON format, this function reads the data
        from the file at the specified path, parse the JSON string and return the
        resulting object. If the file is not in a parsable format or if the file
        could not read, this function returns an error.

         Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> paths_json
        PathsJson()
        >>>
        ```
        """
        return PathsJson._from_py_paths_json(PyPathsJson.from_path(Path(path)))

    @staticmethod
    def from_package_directory(path: os.PathLike[str]) -> PathsJson:
        """
        Parses the object by looking up the appropriate file from the root of the
        specified Conda archive directory, using a format appropriate for the file
        type.

        For example, if the file is in JSON format, this function reads the
        appropriate file from the archive, parse the JSON string and return the
        resulting object. If the file is not in a parsable format or if the file
        could not be read, this function returns an error.
        """
        return PathsJson._from_py_paths_json(PyPathsJson.from_package_directory(Path(path)))

    @staticmethod
    def from_str(string: str) -> PathsJson:
        """
        Parses the object from a string, using a format appropriate for the file
        type.

        For example, if the file is in JSON format, this function parses the JSON
        string and returns the resulting object. If the file is not in a parsable
        format, this function returns an error.
        """
        return PathsJson._from_py_paths_json(PyPathsJson.from_str(string))

    @classmethod
    async def from_remote_url(cls, client: Client, url: str) -> Optional[PathsJson]:
        """
        Fetches `info/paths.json` from a remote package archive URL.
        """
        py_paths_json = await PyPathsJson.from_remote_url(client._client, url)
        if py_paths_json is None:
            return None
        return cls._from_py_paths_json(py_paths_json)

    @staticmethod
    def package_path() -> Path:
        """
        Returns the path to the file within the Conda archive.

        The path is relative to the root of the archive and includes any necessary
        directories.
        """
        return PathsJson.package_path()

    @staticmethod
    def from_deprecated_package_directory(path: os.PathLike[str]) -> PathsJson:
        """
        Constructs a new instance by reading older (deprecated) files from a package directory.

        In older package archives the `paths.json` file does not exist. These packages contain the
        information normally present in the `paths.json` file spread over different files in the
        archive.

        This function reads the different files and tries to reconstruct a `paths.json` from it.
        """
        return PathsJson._from_py_paths_json(PyPathsJson.from_deprecated_package_directory(path))

    @staticmethod
    def from_package_directory_with_deprecated_fallback(
        path: os.PathLike[str],
    ) -> PathsJson:
        """
        Reads the file from a package archive directory. If the `paths.json` file could not be found
        use the `from_deprecated_package_directory` method as a fallback.
        """
        return PathsJson._from_py_paths_json(PyPathsJson.from_package_directory_with_deprecated_fallback(path))

    @property
    def paths(self) -> List[PathsEntry]:
        """
        All entries included in the package.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> paths_json.paths
        [PathsEntry(relative_path="Lib/site-packages/conda-22.9.0-py3.8.egg-info/PKG-INFO", no_link=False, path_type=PathType(hardlink=True), prefix_placeholder="None", sha256="1323efbd9b3abb527b06435392b39de11710eb3a814e87a8174230c8f5a0826a", size_in_bytes=1229), ...]
        >>> paths_json.paths = [PathsEntry(relative_path="new/path", no_link=True, path_type=PathType("softlink"), prefix_placeholder=None, sha256=None, size_in_bytes=None)]
        >>> len(paths_json.paths)
        1
        >>>
        ```
        """
        return [PathsEntry._from_py_paths_entry(path) for path in self._inner.paths]

    @paths.setter
    def paths(self, paths: List[PathsEntry]) -> None:
        self._inner.paths = [entry._inner for entry in paths]

    @property
    def paths_version(self) -> int:
        """
        The version of the file.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> paths_json.paths_version
        1
        >>> paths_json.paths_version = 2
        >>> paths_json.paths_version
        2
        >>>
        ```
        """
        return self._inner.paths_version

    @paths_version.setter
    def paths_version(self, version: int) -> None:
        self._inner.paths_version = version

    @classmethod
    def _from_py_paths_json(cls, py_paths_json: PyPathsJson) -> PathsJson:
        paths_json = cls.__new__(cls)
        paths_json._inner = py_paths_json

        return paths_json

    def __repr__(self) -> str:
        """
        Returns a representation of the PathsJson.
        """
        return "PathsJson()"


class PathsEntry:
    """
    A single entry in the `paths.json` file.
    """

    _inner: PyPathsEntry

    def __init__(
        self,
        relative_path: str,
        no_link: bool,
        path_type: PathType,
        prefix_placeholder: Optional[PrefixPlaceholder],
        sha256: Optional[bytes],
        size_in_bytes: Optional[int],
    ) -> None:
        """
        Create a new paths entry.

        Parameters
        ----------
        relative_path : str
            The relative path from the root of the package
        no_link : bool
            Whether or not this file should be linked when installing the package
        path_type : PathType
            How to include the file when installing the package (hardlink, softlink, or directory)
        prefix_placeholder : Optional[PrefixPlaceholder]
            The placeholder prefix used in the file, if any
        sha256 : Optional[bytes]
            The SHA256 hash of the file contents (only used in paths.json version 1)
        size_in_bytes : Optional[int]
            The size of the file in bytes (only used in paths.json version 1)

        Examples
        --------
        ```python
        >>> # Create a basic file entry
        >>> entry = PathsEntry(
        ...     relative_path="lib/file.txt",
        ...     no_link=False,
        ...     path_type=PathType("hardlink"),
        ...     prefix_placeholder=None,
        ...     sha256=None,
        ...     size_in_bytes=None
        ... )
        >>> str(entry.relative_path)
        'lib/file.txt'
        >>> entry.no_link
        False
        >>> entry.path_type.hardlink
        True
        >>>
        >>> # Create an entry with prefix placeholder
        >>> placeholder = PrefixPlaceholder(FileMode("text"), "/old/prefix")
        >>> sha256 = bytes.fromhex("c609c2f1a8594abf959388e559d76241e51b0216faa7b37f529255eb1fc2c5eb")
        >>> entry = PathsEntry(
        ...     relative_path="bin/script",
        ...     no_link=False,
        ...     path_type=PathType("hardlink"),
        ...     prefix_placeholder=placeholder,
        ...     sha256=sha256,
        ...     size_in_bytes=1234
        ... )
        >>> entry.prefix_placeholder.placeholder
        '/old/prefix'
        >>> entry.size_in_bytes
        1234
        >>>
        ```
        """
        if prefix_placeholder is not None:
            prefix_placeholder = prefix_placeholder._inner
        self._inner = PyPathsEntry(relative_path, no_link, path_type._inner, prefix_placeholder, sha256, size_in_bytes)

    @property
    def relative_path(self) -> Path:
        """
        The relative path from the root of the package.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> str(entry.relative_path)
        'Lib/site-packages/conda-22.9.0-py3.8.egg-info/PKG-INFO'
        >>> entry.relative_path = "new/path"
        >>> str(entry.relative_path)
        'new/path'
        >>>
        ```
        """
        return self._inner.relative_path

    @relative_path.setter
    def relative_path(self, path: str) -> None:
        self._inner.relative_path = path

    @property
    def no_link(self) -> bool:
        """
        Whether or not this file should be linked or not when installing the package.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> entry.no_link
        False
        >>> entry.no_link = True
        >>> entry.no_link
        True
        >>>
        ```
        """
        return self._inner.no_link

    @no_link.setter
    def no_link(self, no_link: bool) -> None:
        self._inner.no_link = no_link

    @property
    def path_type(self) -> PathType:
        """
        Determines how to include the file when installing the package.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> entry.path_type
        PathType(hardlink=True)
        >>> new_type = PathType("softlink")
        >>> entry.path_type = new_type
        >>> entry.path_type
        PathType(softlink=True)
        >>>
        ```
        """
        return PathType._from_py_path_type(self._inner.path_type)

    @path_type.setter
    def path_type(self, path_type: "PathType") -> None:
        self._inner.path_type = path_type._inner

    @property
    def prefix_placeholder(self) -> Optional[PrefixPlaceholder]:
        """
        Optionally the placeholder prefix used in the file. If this value is `None`
        the prefix is not present in the file.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> entry.prefix_placeholder
        >>> new_placeholder = PrefixPlaceholder(FileMode("text"), "placeholder")
        >>> entry.prefix_placeholder = new_placeholder
        >>> entry.prefix_placeholder
        PrefixPlaceholder(file_mode=FileMode("text"), placeholder="placeholder")
        >>>
        ```
        """
        if placeholder := self._inner.prefix_placeholder:
            return PrefixPlaceholder._from_py_prefix_placeholder(placeholder)

        return None

    @prefix_placeholder.setter
    def prefix_placeholder(self, placeholder: Optional[PrefixPlaceholder]) -> None:
        if placeholder is None:
            self._inner.prefix_placeholder = None
        else:
            self._inner.prefix_placeholder = placeholder._inner

    @property
    def sha256(self) -> Optional[bytes]:
        """
        A hex representation of the SHA256 hash of the contents of the file.
        This entry is only present in version 1 of the paths.json file.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> entry.sha256.hex()
        '1323efbd9b3abb527b06435392b39de11710eb3a814e87a8174230c8f5a0826a'
        >>> entry.sha256 = bytes.fromhex('058016a01bb3845320c81755882a367e03a449c1898a3de4f3ea54112fb3eba4')
        >>> entry.sha256.hex()
        '058016a01bb3845320c81755882a367e03a449c1898a3de4f3ea54112fb3eba4'
        >>>
        ```
        """
        return self._inner.sha256

    @sha256.setter
    def sha256(self, sha: Optional[bytes]) -> None:
        self._inner.sha256 = sha

    @property
    def size_in_bytes(self) -> Optional[int]:
        """
        The size of the file in bytes.
        This entry is only present in version 1 of the paths.json file.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> entry.size_in_bytes
        1229
        >>> entry.size_in_bytes = 42
        >>> entry.size_in_bytes
        42
        >>>
        ```
        """
        if size := self._inner.size_in_bytes:
            return size

        return None

    @size_in_bytes.setter
    def size_in_bytes(self, size: Optional[int]) -> None:
        self._inner.size_in_bytes = size

    @classmethod
    def _from_py_paths_entry(cls, py_paths_entry: PyPathsEntry) -> PathsEntry:
        paths_entry = cls.__new__(cls)
        paths_entry._inner = py_paths_entry

        return paths_entry

    def __repr__(self) -> str:
        """
        Returns a representation of the PathsEntry.
        """
        sha256_str = self.sha256.hex() if self.sha256 else None
        return f'PathsEntry(relative_path="{self.relative_path}", no_link={self.no_link}, path_type={self.path_type}, prefix_placeholder="{self.prefix_placeholder}", sha256="{sha256_str}", size_in_bytes={self.size_in_bytes})'


class PathType:
    """
    The path type of the path entry
    """

    _inner: PyPathType

    def __init__(self, path_type: Literal["hardlink", "softlink", "directory"]) -> None:
        self._inner = PyPathType(path_type)

    @property
    def hardlink(self) -> bool:
        """
        The path should be hard linked (the default).

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> path_type = entry.path_type
        >>> path_type.hardlink
        True
        >>>
        ```
        """
        return self._inner.hardlink

    @property
    def softlink(self) -> bool:
        """
        The path should be soft linked.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> path_type = entry.path_type
        >>> path_type.softlink
        False
        >>>
        ```
        """
        return self._inner.softlink

    @property
    def directory(self) -> bool:
        """
        This should explicitly create an empty directory.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[0]
        >>> path_type = entry.path_type
        >>> path_type.directory
        False
        >>>
        ```
        """
        return self._inner.directory

    @classmethod
    def _from_py_path_type(cls, py_paths_type: PyPathType) -> PathType:
        path_type = cls.__new__(cls)
        path_type._inner = py_paths_type

        return path_type

    def __repr__(self) -> str:
        """
        Returns a representation of the PathType.
        """
        if self._inner.hardlink:
            return "PathType(hardlink=True)"
        elif self._inner.softlink:
            return "PathType(softlink=True)"
        else:
            return "PathType(directory=True)"


class PrefixPlaceholder:
    """
    Description off a placeholder text found in a file that must be replaced
    when installing the file into the prefix.
    """

    _inner: PyPrefixPlaceholder

    def __init__(self, file_mode: FileMode, placeholder: str) -> None:
        """
        Create a new prefix placeholder.

        Parameters
        ----------
        file_mode: FileMode
            The file mode of the entry.
        placeholder: str
            The placeholder prefix used in the file.

        Examples
        --------
        ```python
        >>> placeholder = PrefixPlaceholder(FileMode("text"), "placeholder")
        >>> placeholder
        PrefixPlaceholder(file_mode=FileMode("text"), placeholder="placeholder")
        >>>
        ```
        """
        self._inner = PyPrefixPlaceholder(file_mode._inner, placeholder)

    @property
    def file_mode(self) -> FileMode:
        """
        The type of the file, either binary or text.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[-1]
        >>> entry.prefix_placeholder.file_mode
        FileMode("text")
        >>>
        ```
        """
        return FileMode._from_py_file_mode(self._inner.file_mode)

    @property
    def placeholder(self) -> str:
        """
        The placeholder prefix used in the file. This is the path of the
        prefix when the package was build.

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[-1]
        >>> entry.prefix_placeholder.placeholder
        'D:\\\\bld\\\\conda_1667595064120\\\\_h_env'
        >>>
        ```
        """
        return self._inner.placeholder

    @classmethod
    def _from_py_prefix_placeholder(cls, py_prefix_placeholder: PyPrefixPlaceholder) -> PrefixPlaceholder:
        prefix_placeholder = cls.__new__(cls)
        prefix_placeholder._inner = py_prefix_placeholder

        return prefix_placeholder

    def __repr__(self) -> str:
        """
        Returns a representation of the PrefixPlaceholder.
        """
        return f'PrefixPlaceholder(file_mode={self.file_mode}, placeholder="{self.placeholder}")'


class FileMode:
    """
    The file mode of the entry.
    """

    _inner: PyFileMode | None = None

    def __init__(self, file_mode: Literal["binary", "text"]) -> None:
        self._inner = PyFileMode(file_mode)

    @property
    def binary(self) -> bool:
        """
        The file is a binary file (needs binary prefix replacement).

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[-1]
        >>> file_mode = entry.prefix_placeholder.file_mode
        >>> file_mode.binary
        False
        >>>
        ```
        """
        return self._inner.binary if self._inner else False

    @property
    def text(self) -> bool:
        """
        The file is a text file (needs text prefix replacement).

        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[-1]
        >>> file_mode = entry.prefix_placeholder.file_mode
        >>> file_mode.text
        True
        >>>
        ```
        """
        return self._inner.text if self._inner else False

    @property
    def unknown(self) -> bool:
        """
        The file mode is unknown/unspecified
        Examples
        --------
        ```python
        >>> paths_json = PathsJson.from_path(
        ...     "../test-data/conda-22.9.0-py38haa244fe_2-paths.json"
        ... )
        >>> entry = paths_json.paths[-1]
        >>> file_mode = entry.prefix_placeholder.file_mode
        >>> file_mode.unknown
        False
        >>>
        """
        return self._inner is None

    @classmethod
    def _from_py_file_mode(cls, py_file_mode: PyFileMode) -> FileMode:
        file_mode = cls.__new__(cls)
        file_mode._inner = py_file_mode

        return file_mode

    def __repr__(self) -> str:
        """
        Returns a representation of the FileMode.
        """
        if self.binary:
            return 'FileMode("binary")'
        elif self.text:
            return 'FileMode("text")'
        else:
            return "FileMode()"
