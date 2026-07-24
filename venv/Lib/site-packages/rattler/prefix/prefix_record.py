from __future__ import annotations
import os
from enum import Enum
from typing import List, Optional

from rattler.rattler import PyRecord, PyLink
from rattler.prefix.prefix_paths import PrefixPaths
from rattler.repo_data.record import RepoDataRecord
from pathlib import Path


class LinkType(Enum):
    HARDLINK = ("hardlink",)
    COPY = ("copy",)
    SOFTLINK = ("softlink",)
    DIRECTORY = ("directory",)


class Link:
    _inner: PyLink

    def __init__(self, path: os.PathLike[str], type: Optional[LinkType]) -> None:
        self._inner = PyLink(path, type.value if type else None)


class PrefixRecord(RepoDataRecord):
    @classmethod
    def _from_py_record(cls, py_record: PyRecord) -> PrefixRecord:
        """Construct Rattler PrefixRecord from FFI PyRecord object."""

        # quick sanity check
        assert py_record.is_prefix_record
        record = cls.__new__(cls)
        record._record = py_record
        return record

    def __init__(
        self,
        repodata_record: RepoDataRecord,
        paths_data: PrefixPaths,
        link: Optional[Link] = None,
        package_tarball_full_path: Optional[os.PathLike[str]] = None,
        extracted_package_dir: Optional[os.PathLike[str]] = None,
        requested_spec: Optional[str] = None,
        requested_specs: Optional[List[str]] = None,
        files: Optional[List[os.PathLike[str]]] = None,
    ) -> None:
        record = PyRecord.create_prefix_record(
            repodata_record=repodata_record._record,
            paths_data=paths_data._paths,
            link=link._inner if link else None,
            package_tarball_full_path=package_tarball_full_path,
            extracted_package_dir=extracted_package_dir,
            files=files,
            requested_spec=requested_spec,
            requested_specs=requested_specs,
        )
        self._record = record

    @staticmethod
    def from_path(path: os.PathLike[str]) -> PrefixRecord:
        """
        Parses a PrefixRecord from a file.

        Examples
        --------
        ```python
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> assert isinstance(r, PrefixRecord)
        >>>
        ```
        """
        return PrefixRecord._from_py_record(PyRecord.from_path(path))

    def write_to_path(self, path: os.PathLike[str], pretty: bool) -> None:
        """
        Writes the contents of this instance to the file at the specified location.
        """
        self._record.write_to_path(path, pretty)

    @property
    def package_tarball_full_path(self) -> Optional[Path]:
        """
        The path to where the archive of the package was stored on disk.

        Examples
        --------
        ```python
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> str(r.package_tarball_full_path)
        'C:\\\\Users\\\\bas\\\\micromamba\\\\pkgs\\\\requests-2.28.2-pyhd8ed1ab_0.tar.bz2'
        >>>
        ```
        """
        return self._record.package_tarball_full_path

    @package_tarball_full_path.setter
    def package_tarball_full_path(self, value: Optional[os.PathLike[str]]) -> None:
        self._record.package_tarball_full_path = value

    @property
    def extracted_package_dir(self) -> Optional[Path]:
        """
        The path that contains the extracted package content.

        Examples
        --------
        ```python
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> str(r.extracted_package_dir)
        'C:\\\\Users\\\\bas\\\\micromamba\\\\pkgs\\\\requests-2.28.2-pyhd8ed1ab_0'
        >>>
        ```
        """
        return self._record.extracted_package_dir

    @extracted_package_dir.setter
    def extracted_package_dir(self, value: Optional[os.PathLike[str]]) -> None:
        self._record.extracted_package_dir = value

    @property
    def files(self) -> List[Path]:
        """
        A sorted list of all files included in this package

        Examples
        --------
        ```python
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> r.files # doctest:+ELLIPSIS
        [...]
        >>>
        ```
        """
        return self._record.files

    @files.setter
    def files(self, value: List[os.PathLike[str]]) -> None:
        self._record.files = value

    @property
    def paths_data(self) -> PrefixPaths:
        """
        Information about how files have been linked when installing the package.

        Examples
        --------
        ```python
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> r.paths_data # doctest:+ELLIPSIS
        PrefixPaths(paths=[...])
        >>>
        ```
        """
        return PrefixPaths._from_py_prefix_paths(self._record.paths_data)

    @paths_data.setter
    def paths_data(self, value: PrefixPaths) -> None:
        self._record.paths_data = value._paths

    @property
    def requested_spec(self) -> Optional[str]:
        """
        The spec that was used when this package was installed (deprecated).
        Use requested_specs instead.

        Examples
        --------
        ```python
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> r.requested_spec is None
        True
        >>>
        ```
        """
        return self._record.requested_spec

    @requested_spec.setter
    def requested_spec(self, value: Optional[str]) -> None:
        self._record.requested_spec = value

    @property
    def requested_specs(self) -> List[str]:
        """
        The specs that were used when this package was installed.
        If this package was not directly requested by the user but was instead
        installed as a dependency of another package an empty list will be
        returned.

        Examples
        --------
        ```python
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> r.requested_specs
        []
        >>>
        ```
        """
        return self._record.requested_specs

    @requested_specs.setter
    def requested_specs(self, value: List[str]) -> None:
        self._record.requested_specs = value

    def __repr__(self) -> str:
        """
        Returns a representation of the version

        Examples
        --------
        ```python
        >>> r = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/requests-2.28.2-pyhd8ed1ab_0.json"
        ... )
        >>> r
        PrefixRecord(name="requests", version="2.28.2")
        >>>
        ```
        """
        return f'PrefixRecord(name="{self.name.normalized}", version="{self.version}")'
