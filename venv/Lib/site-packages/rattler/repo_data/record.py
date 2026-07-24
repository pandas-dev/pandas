from __future__ import annotations

from rattler.rattler import PyRecord
from rattler.repo_data.package_record import PackageRecord


class RepoDataRecord(PackageRecord):
    def __init__(self, package_record: PackageRecord, file_name: str, url: str, channel: str) -> None:
        record = PyRecord.create_repodata_record(
            package_record._record,
            file_name,
            url,
            channel,
        )
        self._record = record

    @property
    def url(self) -> str:
        """
        The canonical URL from where to get this package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.url
        'https://conda.anaconda.org/conda-forge/win-64/libsqlite-3.40.0-hcfcfb64_0.tar.bz2'
        >>>
        ```
        """
        return self._record.url

    @url.setter
    def url(self, value: str) -> None:
        """
        Set the canonical URL for this package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.url = "https://example.com/package.tar.bz2"
        >>> record.url
        'https://example.com/package.tar.bz2'
        >>>
        ```
        """
        self._record.url = value

    @property
    def channel(self) -> str:
        """
        String representation of the channel where the
        package comes from. This could be a URL but it
        could also be a channel name.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.channel
        'https://conda.anaconda.org/conda-forge/win-64'
        >>>
        ```
        """
        return self._record.channel

    @channel.setter
    def channel(self, value: str) -> None:
        """
        Set the channel for this package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.channel = "conda-forge"
        >>> record.channel
        'conda-forge'
        >>>
        ```
        """
        self._record.channel = value

    @property
    def file_name(self) -> str:
        """
        The filename of the package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.file_name
        'libsqlite-3.40.0-hcfcfb64_0.tar.bz2'
        >>>
        ```
        """
        return self._record.file_name

    @file_name.setter
    def file_name(self, value: str) -> None:
        """
        Set the filename of the package.

        Examples
        --------
        ```python
        >>> from rattler import PrefixRecord
        >>> record = PrefixRecord.from_path(
        ...     "../test-data/conda-meta/libsqlite-3.40.0-hcfcfb64_0.json"
        ... )
        >>> record.file_name = "new-package-1.0.tar.bz2"
        >>> record.file_name
        'new-package-1.0.tar.bz2'
        >>>
        ```
        """
        self._record.file_name = value

    @classmethod
    def _from_py_record(cls, py_record: PyRecord) -> RepoDataRecord:
        """
        Construct Rattler RepoDataRecord from FFI PyRecord object.
        """

        # quick sanity check
        assert py_record.is_repodata_record
        record = cls.__new__(cls)
        record._record = py_record
        return record

    def __repr__(self) -> str:
        """
        Returns a representation of the RepoDataRecord.

        Examples
        --------
        ```python
        >>> from rattler import RepoData, Channel
        >>> repo_data = RepoData.from_path(
        ...     "../test-data/test-server/repo/noarch/repodata.json"
        ... )
        >>> repo_data.into_repo_data(Channel("test"))[0]
        RepoDataRecord(url="https://conda.anaconda.org/test/noarch/test-package-0.1-0.tar.bz2")
        >>>
        ```
        """
        return f'{type(self).__name__}(url="{self.url}")'
