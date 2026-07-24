from __future__ import annotations
import warnings
from dataclasses import dataclass
from typing import Callable, List, Literal, Optional, Union, TYPE_CHECKING

from rattler.networking.client import Client
from rattler.rattler import py_fetch_repo_data, PyFetchRepoDataOptions
from rattler.repo_data.sparse import SparseRepoData

if TYPE_CHECKING:
    import os
    from rattler.channel import Channel
    from rattler.platform import Platform


CacheAction = Literal["cache-or-fetch", "use-cache-only", "force-cache-only", "no-cache"]
Variant = Literal["after-patches", "from-packages", "current"]


@dataclass
class FetchRepoDataOptions:
    cache_action: CacheAction = "cache-or-fetch"
    """How to interact with the cache.

    * `'cache-or-fetch'` (default): Use the cache if its up to date or fetch from the URL if there is no valid cached value.
    * `'use-cache-only'`: Only use the cache, but error out if the cache is not up to date
    * `'force-cache-only'`: Only use the cache, ignore whether or not it is up to date.
    * `'no-cache'`: Do not use the cache even if there is an up to date entry
    """

    variant: Variant = "after-patches"
    """Which type of repodata to download

    * `'after-patches'` (default): Fetch the `repodata.json` file. This `repodata.json` has repodata patches applied.
    * `'from-packages'` Fetch the `repodata_from_packages.json` file
    * `'current'`: Fetch `current_repodata.json` file. This file contains only the latest version of each package.
    """

    zstd_enabled: bool = True
    """Whether the ZSTD compression is enabled or not."""

    bz2_enabled: bool = True
    """Whether the BZ2 compression is enabled or not."""

    jlap_enabled: Optional[bool] = None
    """Deprecated: JLAP support has been removed. This field is ignored."""

    def __post_init__(self) -> None:
        if self.jlap_enabled is not None:
            warnings.warn(
                "The 'jlap_enabled' option is deprecated and has no effect. JLAP support has been removed.",
                DeprecationWarning,
                stacklevel=2,
            )

    def _into_py(self) -> PyFetchRepoDataOptions:
        """
        Converts this object into a type that can be used by the Rust code.

        Examples
        --------
        ```python
        >>> FetchRepoDataOptions()._into_py() # doctest: +ELLIPSIS
        <builtins.PyFetchRepoDataOptions object at 0x...>
        >>>
        ```
        """
        return PyFetchRepoDataOptions(
            cache_action=self.cache_action,
            variant=self.variant,
            zstd_enabled=self.zstd_enabled,
            bz2_enabled=self.bz2_enabled,
        )


async def fetch_repo_data(
    *,
    channels: List[Channel],
    platforms: List[Platform],
    cache_path: Union[str, os.PathLike[str]],
    callback: Optional[Callable[[int, int], None]],
    client: Optional[Client] = None,
    fetch_options: Optional[FetchRepoDataOptions] = None,
) -> List[SparseRepoData]:
    """
    Returns a list of RepoData for given channels and platform.

    Arguments:
        channels: A list of `Channel`s to fetch repo data.
        platforms: A list of `Platform`s for which the repo data
                   should be fetched.
        cache_path: A `os.PathLike[str]` where the repo data should
                    be downloaded.
        callback: A `Callable[[int, int], None]` to report the download
                  progress of repo data.
        client: A `Client` to use for fetching the repo data.

    Returns:
        A list of `SparseRepoData` for requested channels and platforms.
    """
    fetch_options = fetch_options or FetchRepoDataOptions()
    repo_data_list = await py_fetch_repo_data(
        [channel._channel for channel in channels],
        [platform._inner for platform in platforms],
        cache_path,
        callback,
        client._client if client else None,
        fetch_options._into_py(),
    )

    return [SparseRepoData._from_py_sparse_repo_data(repo_data) for repo_data in repo_data_list]
