from __future__ import annotations
import os
from pathlib import Path
from typing import List, Optional, Type, Literal, Iterable
from types import TracebackType

from rattler.match_spec.match_spec import MatchSpec
from rattler.channel.channel import Channel
from rattler.package.package_name import PackageName
from enum import Enum

from rattler.rattler import PySparseRepoData, PyPackageFormatSelection
from rattler.repo_data.record import RepoDataRecord


class PackageFormatSelection(Enum):
    """
    Enum that describes what to do if both a `.tar.bz2` and a `.conda` package is available.
    """

    ONLY_TAR_BZ2 = PyPackageFormatSelection.OnlyTarBz2
    """
    Only use the `.tar.bz2` packages, ignore all `.conda` packages.
    """

    ONLY_CONDA = PyPackageFormatSelection.OnlyConda
    """
    Only use the `.conda` packages, ignore all `.tar.bz2` packages.
    """

    PREFER_CONDA = PyPackageFormatSelection.PreferConda
    """
    Only use the `.conda` packages if there are both a `.tar.bz2` and a `.conda` package available.
    """

    PREFER_CONDA_WITH_WHL = PyPackageFormatSelection.PreferCondaWithWhl
    """
    Only use the `.conda` packages if there are both a `.tar.bz2` and a `.conda` package available.
    Also adds `.whl` files if available.
    """

    BOTH = PyPackageFormatSelection.Both
    """
    Use both the `.tar.bz2` and the `.conda` packages.
    """


class SparseRepoData:
    """
    A class to enable loading records from a `repodata.json` file on demand.
    Since most of the time you don't need all the records from the `repodata.json`
    this can help provide some significant speedups.
    """

    def __init__(
        self,
        channel: Channel,
        subdir: str,
        path: os.PathLike[str] | str,
    ) -> None:
        if not isinstance(channel, Channel):
            raise TypeError(
                "SparseRepoData constructor received unsupported type "
                f" {type(channel).__name__!r} for the `channel` parameter"
            )
        if not isinstance(subdir, str):
            raise TypeError(
                "SparseRepoData constructor received unsupported type "
                f" {type(subdir).__name__!r} for the `subdir` parameter"
            )
        if not isinstance(path, (str, Path)):
            raise TypeError(
                "SparseRepoData constructor received unsupported type "
                f" {type(path).__name__!r} for the `path` parameter"
            )
        self._sparse = PySparseRepoData(channel._channel, subdir, str(path))

    def close(self) -> None:
        """
        Closes any mapped resources associated with this `SparseRepoData`
        instance. It is good practice to call this method when you are done
        with it. This is especially important if you want to modify or delete
        the file from which this instance was created.

        This method will release all resources associated with this instance,
        including those that are currently being used on another thread. This
        method will block until all resources are released.

        This method has no effect if the file is already closed. Once the
        instance is closed, any operation on the instance will raise a
        `ValueError`.

        As a convenience, it is allowed to call this method more than once;
        only the first call, however, will have an effect.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> sparse_data = SparseRepoData(channel, "linux-64", path)
        >>> sparse_data.close()
        >>> sparse_data.package_names() # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ValueError: I/O operation on closed file.
        >>>
        ```
        """
        self._sparse.close()

    def package_names(
        self, package_format_selection: PackageFormatSelection = PackageFormatSelection.PREFER_CONDA
    ) -> List[str]:
        """
        Returns a list over all package names in this repodata file.
        This works by iterating over all elements in the `packages` and
        `conda_packages` fields of the repodata and returning the unique
        package names.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> sparse_data = SparseRepoData(channel, "linux-64", path)
        >>> package_names = sparse_data.package_names()
        >>> package_names
        [...]
        >>> isinstance(package_names[0], str)
        True
        >>>
        ```
        """
        return self._sparse.package_names(package_format_selection.value)

    def record_count(
        self, package_format_selection: PackageFormatSelection = PackageFormatSelection.PREFER_CONDA
    ) -> int:
        """
        Returns the total number of packages in this repodata file.
        :return:
        """
        return self._sparse.record_count(package_format_selection.value)

    def load_records(
        self,
        package_name: str | PackageName,
        package_format_selection: PackageFormatSelection = PackageFormatSelection.PREFER_CONDA,
    ) -> List[RepoDataRecord]:
        """
        Returns all the records for the specified package name.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig, RepoDataRecord, PackageName
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> sparse_data = SparseRepoData(channel, "linux-64", path)
        >>> package_name = PackageName(sparse_data.package_names()[0])
        >>> records = sparse_data.load_records(package_name)
        >>> records
        [...]
        >>> isinstance(records[0], RepoDataRecord)
        True
        >>>
        ```
        """
        if not isinstance(package_name, PackageName):
            package_name = PackageName(package_name)
        return [
            RepoDataRecord._from_py_record(record)
            for record in self._sparse.load_records(package_name._name, package_format_selection.value)
        ]

    def load_all_records(
        self, package_format_selection: PackageFormatSelection = PackageFormatSelection.PREFER_CONDA
    ) -> List[RepoDataRecord]:
        """
        Returns all the records for the specified package name.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig, RepoDataRecord, PackageName
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> sparse_data = SparseRepoData(channel, "linux-64", path)
        >>> records = sparse_data.load_all_records()
        >>> records
        [...]
        >>> isinstance(records[0], RepoDataRecord)
        True
        >>>
        ```
        """
        # maybe change package_name to Union[str, PackageName]
        return [
            RepoDataRecord._from_py_record(record)
            for record in self._sparse.load_all_records(package_format_selection.value)
        ]

    def load_matching_records(
        self,
        specs: Iterable[MatchSpec],
        package_format_selection: PackageFormatSelection = PackageFormatSelection.PREFER_CONDA,
    ) -> List[RepoDataRecord]:
        """
        Returns all the records that match any of the specified MatchSpecs.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig, RepoDataRecord, PackageName
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> sparse_data = SparseRepoData(channel, "linux-64", path)
        >>> [record.file_name for record in sparse_data.load_matching_records([MatchSpec("cuda-version 12.5")])]
        ['cuda-version-12.5-hd4f0392_3.conda']
        >>>
        ```
        """
        return [
            RepoDataRecord._from_py_record(record)
            for record in self._sparse.load_matching_records(
                [spec._match_spec for spec in specs], package_format_selection.value
            )
        ]

    @property
    def subdir(self) -> str:
        """
        Returns the subdirectory from which this repodata was loaded.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> sparse_data = SparseRepoData(channel, "linux-64", path)
        >>> sparse_data.subdir
        'linux-64'
        >>>
        ```
        """
        return self._sparse.subdir

    @staticmethod
    def load_records_recursive(
        repo_data: List[SparseRepoData],
        package_names: List[PackageName],
        package_format_selection: PackageFormatSelection = PackageFormatSelection.PREFER_CONDA,
    ) -> List[List[RepoDataRecord]]:
        """
        Given a set of [`SparseRepoData`]s load all the records
        for the packages with the specified names and all the packages
        these records depend on. This will parse the records for the
        specified packages as well as all the packages these records
        depend on.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig, PackageName
        >>> channel = Channel("dummy")
        >>> subdir = "test-data/dummy/linux-64"
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> sparse_data = SparseRepoData(channel, "linux-64", path)
        >>> package_name = PackageName("python")
        >>> SparseRepoData.load_records_recursive([sparse_data], [package_name])
        [...]
        >>>
        ```
        """
        return [
            [RepoDataRecord._from_py_record(record) for record in list_of_records]
            for list_of_records in PySparseRepoData.load_records_recursive(
                [r._sparse for r in repo_data], [p._name for p in package_names], package_format_selection.value
            )
        ]

    @classmethod
    def _from_py_sparse_repo_data(cls, py_sparse_repo_data: PySparseRepoData) -> SparseRepoData:
        """
        Construct Rattler SparseRepoData from FFI PySparseRepoData object.
        """
        sparse_repo_data = cls.__new__(cls)
        sparse_repo_data._sparse = py_sparse_repo_data
        return sparse_repo_data

    def __repr__(self) -> str:
        """
        Returns a representation of the SparseRepoData.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> sparse_data = SparseRepoData(channel, "linux-64", path)
        >>> sparse_data
        SparseRepoData(subdir="linux-64")
        >>>
        ```
        """
        return f'SparseRepoData(subdir="{self.subdir}")'

    def __enter__(self) -> SparseRepoData:
        """
        Returns the `SparseRepoData` instance itself. This is used to
        enable the use of the `with` statement to automatically close
        the instance when done.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> with SparseRepoData(channel, "linux-64", path) as sparse_data:
        ...     print(sparse_data)
        ...
        SparseRepoData(subdir="linux-64")
        >>>
        ```
        """
        return self

    def __exit__(
        self,
        exctype: Optional[Type[BaseException]],
        excinst: Optional[BaseException],
        exctb: Optional[TracebackType],
    ) -> Literal[False]:
        """
        Closes the `SparseRepoData` instance when exiting the `with` statement.

        Examples
        --------
        ```python
        >>> from rattler import Channel, ChannelConfig
        >>> channel = Channel("dummy", ChannelConfig())
        >>> path = "../test-data/channels/dummy/linux-64/repodata.json"
        >>> with SparseRepoData(channel, "linux-64", path) as sparse_data:
        ...     print(sparse_data)
        ...
        SparseRepoData(subdir="linux-64")
        >>>
        ```
        """
        self.close()
        return False
