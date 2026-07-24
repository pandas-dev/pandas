from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Iterable, List, Optional, Union

from rattler.channel.channel import Channel
from rattler.match_spec.match_spec import MatchSpec
from rattler.networking.client import Client
from rattler.networking.fetch_repo_data import CacheAction
from rattler.package.package_name import PackageName
from rattler.platform.platform import Platform, PlatformLiteral
from rattler.rattler import PyGateway, PyMatchSpec, PySourceConfig
from rattler.repo_data.record import RepoDataRecord

if TYPE_CHECKING:
    from rattler.repo_data.source import RepoDataSource


class _RepoDataSourceAdapter:
    """Adapter that wraps a user's RepoDataSource and converts FFI types.

    This adapter receives raw PyPlatform and PyPackageName from Rust,
    converts them to the proper Python wrapper types (Platform, PackageName),
    and calls the user's implementation.
    """

    def __init__(self, source: RepoDataSource) -> None:
        self._source = source

    async def fetch_package_records(self, py_platform: Any, py_name: Any) -> List[RepoDataRecord]:
        """Convert FFI types and delegate to the wrapped source."""
        # Wrap raw FFI types in Python wrapper classes
        platform = Platform._from_py_platform(py_platform)
        name = PackageName._from_py_package_name(py_name)

        # Call the user's implementation with proper Python types
        return await self._source.fetch_package_records(platform, name)

    def package_names(self, py_platform: Any) -> List[str]:
        """Convert FFI types and delegate to the wrapped source."""
        platform = Platform._from_py_platform(py_platform)
        return self._source.package_names(platform)


@dataclass
class SourceConfig:
    """
    Describes properties about a channel.

    This can be used to configure the Gateway to handle channels in a certain
    way.
    """

    zstd_enabled: bool = True
    """Whether the ZSTD compression is enabled or not."""

    bz2_enabled: bool = True
    """Whether the BZ2 compression is enabled or not."""

    sharded_enabled: bool = True
    """Whether sharded repodata is enabled or not."""

    jlap_enabled: Optional[bool] = None
    """Deprecated: JLAP support has been removed. This field is ignored."""

    cache_action: CacheAction = "cache-or-fetch"
    """How to interact with the cache.

    * `'cache-or-fetch'` (default): Use the cache if its up to date or fetch from the URL if there is no valid cached value.
    * `'use-cache-only'`: Only use the cache, but error out if the cache is not up to date
    * `'force-cache-only'`: Only use the cache, ignore whether or not it is up to date.
    * `'no-cache'`: Do not use the cache even if there is an up to date entry
    """

    def __post_init__(self) -> None:
        if self.jlap_enabled is not None:
            warnings.warn(
                "The 'jlap_enabled' option is deprecated and has no effect. JLAP support has been removed.",
                DeprecationWarning,
                stacklevel=2,
            )

    def _into_py(self) -> PySourceConfig:
        """
        Converts this object into a type that can be used by the Rust code.

        Examples
        --------
        ```python
        >>> SourceConfig()._into_py() # doctest: +ELLIPSIS
        <builtins.PySourceConfig object at 0x...>
        >>>
        ```
        """
        return PySourceConfig(
            zstd_enabled=self.zstd_enabled,
            bz2_enabled=self.bz2_enabled,
            sharded_enabled=self.sharded_enabled,
            cache_action=self.cache_action,
        )


class Gateway:
    """
    The gateway manages all the quircks and complex bits of efficiently acquiring
    repodata. It implements all the necessary logic to fetch the repodata from a
    remote server, cache it locally and convert it into python objects.

    The gateway can also easily be used concurrently, as it is designed to be
    thread-safe. When two threads are querying the same channel at the same time,
    their requests are coalesced into a single request. This is done to reduce the
    number of requests made to the remote server and reduce the overall memory usage.

    The gateway caches the repodata internally, so if the same channel is queried
    multiple times the records will only be fetched once. However, the conversion
    of the records to a python object is done every time the query method is called.
    Therefor, instead of requesting records directly, its more efficient to pass the
    gateway itself to methods that accepts it.
    """

    def __init__(
        self,
        cache_dir: Optional[os.PathLike[str]] = None,
        default_config: Optional[SourceConfig] = None,
        per_channel_config: Optional[dict[str, SourceConfig]] = None,
        max_concurrent_requests: int = 100,
        client: Optional[Client] = None,
        show_progress: bool = False,
    ) -> None:
        """
        Arguments:
            cache_dir: The directory where the repodata should be cached. If not specified the
                       default cache directory is used.
            default_config: The default configuration for channels.
            per_channel_config: Source configuration on a per-URL basis. This URL is used as a
                                prefix, so any channel that starts with the URL uses the configuration.
                                The configuration with the longest matching prefix is used.
            max_concurrent_requests: The maximum number of concurrent requests that can be made.
            client: An authenticated client to use for acquiring repodata. If not specified a default
                    client will be used.
            show_progress: Whether to show progress bars when fetching repodata.

        Examples
        --------
        ```python
        >>> Gateway()
        Gateway()
        >>>
        ```
        """
        default_config = default_config or SourceConfig()

        self._gateway = PyGateway(
            cache_dir=cache_dir,
            default_config=default_config._into_py(),
            per_channel_config={channel: config._into_py() for channel, config in (per_channel_config or {}).items()},
            max_concurrent_requests=max_concurrent_requests,
            client=client._client if client is not None else None,
            show_progress=show_progress,
        )

    async def query(
        self,
        sources: Iterable[Union[Channel, str, RepoDataSource]],
        platforms: Iterable[Platform | PlatformLiteral],
        specs: Iterable[MatchSpec | PackageName | str],
        recursive: bool = True,
    ) -> List[List[RepoDataRecord]]:
        """Queries the gateway for repodata from channels and custom sources.

        If `recursive` is `True` the gateway will recursively fetch the dependencies of the
        encountered records. If `recursive` is `False` only the records with the package names
        specified in `specs` are returned.

        The `specs` can either be a `MatchSpec`, `PackageName` or a string. If a string or a
        `PackageName` is provided it will be converted into a MatchSpec that matches any record
        with the given name. If a `MatchSpec` is provided all records that match the name
        specified in the spec will be returned, but only the dependencies of the records
        that match the entire spec are recursively fetched.

        The gateway caches records from channels internally, so if the same channel is queried
        multiple times the records will only be fetched once. However, the conversion of the
        records to a python object is done every time the query method is called.

        Note: Custom RepoDataSource implementations are **not cached** by the gateway. If caching
        is needed for custom sources, it must be implemented within the source itself.

        Arguments:
            sources: The sources to query. Can be channels (by name, URL, or Channel object)
                     or custom RepoDataSource implementations.
            platforms: The platforms to query.
            specs: The specs to query.
            recursive: Whether recursively fetch dependencies or not.

        Returns:
            A list of lists of `RepoDataRecord`s. The outer list contains the results for each
            source in the same order they are provided in the `sources` argument.

        Examples
        --------
        ```python
        >>> import asyncio
        >>> gateway = Gateway()
        >>> records = asyncio.run(gateway.query(["conda-forge"], ["linux-aarch64"], ["python"]))
        >>> assert len(records) == 1
        >>>
        ```
        """
        py_records = await self._gateway.query(
            sources=_convert_sources(sources),
            platforms=[
                platform._inner if isinstance(platform, Platform) else Platform(platform)._inner
                for platform in platforms
            ],
            specs=[
                spec._match_spec if isinstance(spec, MatchSpec) else PyMatchSpec(str(spec), True, True)
                for spec in specs
            ],
            recursive=recursive,
        )

        # Convert the records into python objects
        return [[RepoDataRecord._from_py_record(record) for record in records] for records in py_records]

    async def names(
        self,
        sources: Iterable[Union[Channel, str, RepoDataSource]],
        platforms: Iterable[Platform | PlatformLiteral],
    ) -> List[PackageName]:
        """Queries all the names of packages in channels or custom sources.

        Arguments:
            sources: The sources to query. Can be channels (by name, URL, or Channel object)
                     or custom RepoDataSource implementations.
            platforms: The platforms to query.

        Returns:
            A list of package names that are present in the given subdirectories.

        Examples
        --------
        ```python
        >>> import asyncio
        >>> gateway = Gateway()
        >>> records = asyncio.run(gateway.names(["conda-forge"], ["linux-64"]))
        >>> PackageName("python") in records
        True
        >>>
        ```
        """

        py_package_names = await self._gateway.names(
            sources=_convert_sources(sources),
            platforms=[
                platform._inner if isinstance(platform, Platform) else Platform(platform)._inner
                for platform in platforms
            ],
        )

        # Convert the records into python objects
        return [PackageName._from_py_package_name(package_name) for package_name in py_package_names]

    def clear_repodata_cache(
        self,
        channel: Channel | str,
        subdirs: Optional[Iterable[Platform | PlatformLiteral]] = None,
        clear_disk: bool = False,
    ) -> None:
        """
        Clears the cache for the given channel.

        Any subsequent query will re-fetch any required data from the source.

        Arguments:
            channel: The channel to clear the cache for.
            subdirs: A selection of subdirectories to clear, if `None` is specified
                     all subdirectories of the channel are cleared.
            clear_disk: If `True`, also clears the on-disk cache. By default only the
                        in-memory cache is cleared.

        Examples
        --------
        ```python
        >>> gateway = Gateway()
        >>> gateway.clear_repodata_cache("conda-forge", ["linux-64"])
        >>> gateway.clear_repodata_cache("robostack")
        >>> gateway.clear_repodata_cache("conda-forge", clear_disk=True)
        >>>
        ```
        """
        self._gateway.clear_repodata_cache(
            channel._channel if isinstance(channel, Channel) else Channel(channel)._channel,
            {subdir._inner if isinstance(subdir, Platform) else Platform(subdir)._inner for subdir in subdirs}
            if subdirs is not None
            else None,
            clear_disk,
        )

    def __repr__(self) -> str:
        """
        Returns a representation of the Gateway.

        Examples
        --------
        ```python
        >>> Gateway()
        Gateway()
        >>>
        ```
        """
        return f"{type(self).__name__}()"


def _convert_sources(sources: Iterable[Any]) -> List[Any]:
    """Convert an iterable of sources to a list suitable for the Rust gateway.

    Channels are converted to their internal PyChannel representation.
    Custom RepoDataSource implementations are wrapped in an adapter that
    converts between FFI types and Python wrapper types.

    Raises:
        TypeError: If a source doesn't implement the required interface.
    """
    from rattler.repo_data.source import RepoDataSource

    converted = []
    for source in sources:
        if isinstance(source, str):
            # String channel name/URL - convert to PyChannel
            converted.append(Channel(source)._channel)
        elif isinstance(source, Channel):
            # Channel object - extract PyChannel
            converted.append(source._channel)
        elif isinstance(source, RepoDataSource):
            # Wrap RepoDataSource in adapter for FFI type conversion
            converted.append(_RepoDataSourceAdapter(source))
        else:
            raise TypeError(
                f"Expected Channel, str, or object implementing RepoDataSource protocol, "
                f"got {type(source).__name__}. "
                f"See rattler.RepoDataSource for the required interface."
            )
    return converted
