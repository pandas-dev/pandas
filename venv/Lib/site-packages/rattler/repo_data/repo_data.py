from __future__ import annotations
from pathlib import Path
from typing import Optional, Union, List, TYPE_CHECKING


if TYPE_CHECKING:
    from os import PathLike
    from rattler.channel import Channel
    from rattler.repo_data import PatchInstructions, RepoDataRecord

from rattler.rattler import PyChannelInfo, PyChannelRelations, PyRepoData


class ChannelRelations:
    """
    Relationships to other channels declared inside ``repodata.json`` under
    ``info.channel_relations`` as specified by
    [CEP-42](https://github.com/conda/ceps/blob/main/cep-0042.md).

    `ChannelRelations` is a read-only view obtained via
    ``RepoData.from_path(...).info.channel_relations``; it cannot be
    constructed directly.
    """

    _inner: PyChannelRelations

    @classmethod
    def _from_inner(cls, py_channel_relations: PyChannelRelations) -> ChannelRelations:
        instance = cls.__new__(cls)
        instance._inner = py_channel_relations
        return instance

    @property
    def base(self) -> Optional[str]:
        """A reference to a channel with higher priority than the declaring channel."""
        return self._inner.base

    @property
    def overrides(self) -> Optional[str]:
        """A reference to a channel with lower priority than the declaring channel."""
        return self._inner.overrides

    def __repr__(self) -> str:
        return f"ChannelRelations(base={self.base!r}, overrides={self.overrides!r})"


class ChannelInfo:
    """
    The ``info`` section of a ``repodata.json`` file.

    `ChannelInfo` is a read-only view obtained via
    ``RepoData.from_path(...).info``; it cannot be constructed directly.
    """

    _inner: PyChannelInfo

    @classmethod
    def _from_inner(cls, py_channel_info: PyChannelInfo) -> ChannelInfo:
        instance = cls.__new__(cls)
        instance._inner = py_channel_info
        return instance

    @property
    def subdir(self) -> Optional[str]:
        """The channel's subdirectory (e.g. ``"linux-64"``)."""
        return self._inner.subdir

    @property
    def base_url(self) -> Optional[str]:
        """The base URL for all package URLs in this channel, if any."""
        return self._inner.base_url

    @property
    def channel_relations(self) -> Optional[ChannelRelations]:
        """
        Channel relations declared by this channel, see
        [CEP-42](https://github.com/conda/ceps/blob/main/cep-0042.md).

        ``None`` when the channel does not declare any relations.
        """
        relations = self._inner.channel_relations
        if relations is None:
            return None
        return ChannelRelations._from_inner(relations)

    def __repr__(self) -> str:
        return (
            f"ChannelInfo(subdir={self.subdir!r}, base_url={self.base_url!r}, "
            f"channel_relations={self.channel_relations!r})"
        )


class RepoData:
    """
    Repository metadata as deserialized from a ``repodata.json`` file.

    `RepoData` is obtained via `RepoData.from_path(path)`.
    """

    _repo_data: PyRepoData

    @classmethod
    def from_path(cls, path: Union[str, PathLike[str]]) -> RepoData:
        """
        Load a `RepoData` from a ``repodata.json`` file on disk.

        Examples
        --------
        ```python
        >>> repo_data = RepoData.from_path("../test-data/test-server/repo/noarch/repodata.json")
        >>> repo_data
        RepoData()
        >>>
        ```
        """
        if not isinstance(path, (str, Path)):
            raise TypeError(
                f"RepoData.from_path received unsupported type {type(path).__name__!r} for the `path` parameter"
            )
        return cls._from_py_repo_data(PyRepoData.from_path(path))

    @property
    def info(self) -> Optional[ChannelInfo]:
        """Returns the channel info contained in the repodata, if any."""
        info = self._repo_data.info
        if info is None:
            return None
        return ChannelInfo._from_inner(info)

    @property
    def version(self) -> Optional[int]:
        """Returns the repodata format version, if any."""
        return self._repo_data.version

    def apply_patches(self, instructions: PatchInstructions) -> None:
        """
        Apply a patch to a repodata file.
        Note that we currently do not handle revoked instructions.
        """
        self._repo_data.apply_patches(instructions._patch_instructions)

    def into_repo_data(self, channel: Channel) -> List[RepoDataRecord]:
        """
        Builds a `List[RepoDataRecord]` from the packages in a
        `RepoData` given the source of the data.

        Examples
        --------
        ```python
        >>> from rattler import Channel
        >>> repo_data = RepoData.from_path("../test-data/test-server/repo/noarch/repodata.json")
        >>> repo_data.into_repo_data(Channel("test"))
        [...]
        >>>
        ```
        """
        from rattler.repo_data import RepoDataRecord

        return [
            RepoDataRecord._from_py_record(record)
            for record in PyRepoData.repo_data_to_records(self._repo_data, channel._channel)
        ]

    @classmethod
    def _from_py_repo_data(cls, py_repo_data: PyRepoData) -> RepoData:
        """
        Construct Rattler RepoData from FFI PyRepoData object.
        """
        repo_data = cls.__new__(cls)
        repo_data._repo_data = py_repo_data
        return repo_data

    def __repr__(self) -> str:
        """
        Returns a representation of the RepoData.

        Examples
        --------
        ```python
        >>> repo_data = RepoData.from_path("../test-data/test-server/repo/noarch/repodata.json")
        >>> repo_data
        RepoData()
        >>>
        ```
        """
        return f"{type(self).__name__}()"
