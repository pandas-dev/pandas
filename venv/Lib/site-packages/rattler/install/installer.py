from __future__ import annotations
import os
from typing import List, Optional, Protocol, runtime_checkable

from rattler.match_spec import MatchSpec
from rattler.networking.client import Client
from rattler.platform.platform import Platform
from rattler.prefix.prefix_record import PrefixRecord
from rattler.repo_data.record import RepoDataRecord

from rattler.rattler import py_install


@runtime_checkable
class InstallerReporter(Protocol):
    """
    Protocol for custom installation progress reporters.

    Implement this protocol (by subclassing or structurally) and override the
    methods we care about to receive progress callbacks during package
    installation. All methods have no-op defaults so we only need to implement
    the ones relevant to our use case.

    Methods that return an ``int`` act as opaque tokens: the value we return
    will be passed back to the corresponding ``*_complete`` callback so we can
    correlate start/finish events.

    Example
    -------
    ```python
    from rattler.install import InstallerReporter

    class MyReporter:  # no need to inherit — structural subtyping is sufficient
        def on_transaction_start(self, total_operations: int) -> None:
            print(f"Starting installation of {total_operations} operations")

        def on_download_progress(
            self,
            download_idx: int,
            progress: int,
            total: int | None,
        ) -> None:
            pct = f"{progress}/{total}" if total else str(progress)
            print(f"  downloading [{download_idx}]: {pct} bytes")

        def on_transaction_complete(self) -> None:
            print("Installation complete!")
    ```
    """

    def on_transaction_start(self, total_operations: int) -> None:
        """Called once when the installation transaction begins.

        Parameters
        ----------
        total_operations:
            The total number of install/unlink operations in this transaction.
        """

    def on_transaction_operation_start(self, operation: int) -> None:
        """Called when a single transaction operation starts.

        Parameters
        ----------
        operation:
            Index of the operation within the transaction.
        """

    def on_populate_cache_start(self, operation: int, package_name: str) -> int:
        """Called when rattler begins populating the local cache for a package.

        Parameters
        ----------
        operation:
            Index of the operation within the transaction.
        package_name:
            Normalised name of the package being cached.

        Returns
        -------
        int
            An opaque token passed back to :meth:`on_populate_cache_complete`.
        """
        return 0

    def on_validate_start(self, cache_entry: int) -> int:
        """Called when cache-entry validation begins.

        Parameters
        ----------
        cache_entry:
            Token returned by :meth:`on_populate_cache_start`.

        Returns
        -------
        int
            An opaque token passed back to :meth:`on_validate_complete`.
        """
        return 0

    def on_validate_complete(self, validate_idx: int) -> None:
        """Called when cache-entry validation finishes.

        Parameters
        ----------
        validate_idx:
            Token returned by :meth:`on_validate_start`.
        """

    def on_download_start(self, cache_entry: int) -> int:
        """Called just before a package download begins.

        Parameters
        ----------
        cache_entry:
            Token returned by :meth:`on_populate_cache_start`.

        Returns
        -------
        int
            An opaque token passed to :meth:`on_download_progress` and
            :meth:`on_download_completed`.
        """
        return 0

    def on_download_progress(self, download_idx: int, progress: int, total: Optional[int]) -> None:
        """Called periodically with download byte progress.

        Parameters
        ----------
        download_idx:
            Token returned by :meth:`on_download_start`.
        progress:
            Bytes downloaded so far.
        total:
            Total expected bytes, or ``None`` if unknown.
        """

    def on_download_completed(self, download_idx: int) -> None:
        """Called when a download finishes.

        Parameters
        ----------
        download_idx:
            Token returned by :meth:`on_download_start`.
        """

    def on_populate_cache_complete(self, cache_entry: int) -> None:
        """Called when the cache has been fully populated for a package.

        Parameters
        ----------
        cache_entry:
            Token returned by :meth:`on_populate_cache_start`.
        """

    def on_unlink_start(self, operation: int, package_name: str) -> int:
        """Called when rattler begins unlinking (removing) a package.

        Parameters
        ----------
        operation:
            Index of the operation within the transaction.
        package_name:
            Normalised name of the package being unlinked.

        Returns
        -------
        int
            An opaque token passed back to :meth:`on_unlink_complete`.
        """
        return 0

    def on_unlink_complete(self, index: int) -> None:
        """Called when a package has been fully unlinked.

        Parameters
        ----------
        index:
            Token returned by :meth:`on_unlink_start`.
        """

    def on_link_start(self, operation: int, package_name: str) -> int:
        """Called when rattler begins linking (installing) a package.

        Parameters
        ----------
        operation:
            Index of the operation within the transaction.
        package_name:
            Normalised name of the package being linked.

        Returns
        -------
        int
            An opaque token passed back to :meth:`on_link_complete`.
        """
        return 0

    def on_link_complete(self, index: int) -> None:
        """Called when a package has been fully linked into the prefix.

        Parameters
        ----------
        index:
            Token returned by :meth:`on_link_start`.
        """

    def on_transaction_operation_complete(self, operation: int) -> None:
        """Called when a single transaction operation finishes.

        Parameters
        ----------
        operation:
            Index of the operation within the transaction.
        """

    def on_transaction_complete(self) -> None:
        """Called once when the entire transaction has finished."""

    def on_post_link_start(self, package_name: str, script_path: str) -> int:
        """Called when a post-link script begins execution.

        Parameters
        ----------
        package_name:
            Name of the package whose post-link script is running.
        script_path:
            Relative path to the script within the prefix.

        Returns
        -------
        int
            An opaque token passed back to :meth:`on_post_link_complete`.
        """
        return 0

    def on_post_link_complete(self, index: int, success: bool) -> None:
        """Called when a post-link script finishes.

        Parameters
        ----------
        index:
            Token returned by :meth:`on_post_link_start`.
        success:
            ``True`` if the script exited successfully.
        """

    def on_pre_unlink_start(self, package_name: str, script_path: str) -> int:
        """Called when a pre-unlink script begins execution.

        Parameters
        ----------
        package_name:
            Name of the package whose pre-unlink script is running.
        script_path:
            Relative path to the script within the prefix.

        Returns
        -------
        int
            An opaque token passed back to :meth:`on_pre_unlink_complete`.
        """
        return 0

    def on_pre_unlink_complete(self, index: int, success: bool) -> None:
        """Called when a pre-unlink script finishes.

        Parameters
        ----------
        index:
            Token returned by :meth:`on_pre_unlink_start`.
        success:
            ``True`` if the script exited successfully.
        """


async def install(
    records: List[RepoDataRecord],
    target_prefix: str | os.PathLike[str],
    cache_dir: Optional[os.PathLike[str]] = None,
    installed_packages: Optional[List[PrefixRecord]] = None,
    reinstall_packages: Optional[set[str]] = None,
    ignored_packages: Optional[set[str]] = None,
    platform: Optional[Platform] = None,
    execute_link_scripts: bool = False,
    show_progress: bool = True,
    client: Optional[Client] = None,
    requested_specs: Optional[List[MatchSpec]] = None,
    reporter: Optional[InstallerReporter] = None,
) -> None:
    """
    Create an environment by downloading and linking the `dependencies` in
    the `target_prefix` directory.

    !!! warning

        When `execute_link_scripts` is set to `True` the post-link and pre-unlink scripts of
        packages will be executed. These scripts are not sandboxed and can be used to execute
        arbitrary code. It is therefor discouraged to enable executing link scripts.

    Example
    -------
    ```python
    >>> import asyncio
    >>> from rattler import solve, install
    >>> from tempfile import TemporaryDirectory
    >>> temp_dir = TemporaryDirectory()
    >>>
    >>> async def main():
    ...     # Solve an environment with python 3.9 for the current platform
    ...     records = await solve(sources=["conda-forge"], specs=["python 3.9.*"])
    ...
    ...     # Link the environment in a temporary directory (you can pass any kind of path here).
    ...     await install(records, target_prefix=temp_dir.name)
    ...
    ...     # That's it! The environment is now created.
    ...     # You will find Python under `f"{temp_dir.name}/bin/python"` or `f"{temp_dir.name}/python.exe"` on Windows.
    >>> asyncio.run(main())

    ```

    Arguments:
        records: A list of solved `RepoDataRecord`s.
        target_prefix: Path to the directory where the environment should
                be created.
        cache_dir: Path to directory where the dependencies will be
                downloaded and cached.
        installed_packages: A list of `PrefixRecord`s which are
                already installed in the `target_prefix`. This can be obtained by loading
                `PrefixRecord`s from `{target_prefix}/conda-meta/`.
                If `None` is specified then the `target_prefix` will be scanned for installed
                packages.
        reinstall_packages: A list of package names that should be reinstalled.
        ignored_packages: A list of package names that should be ignored (left untouched).
                These packages will not be removed, installed, or updated.
        platform: Target platform to create and link the
                environment. Defaults to current platform.
        execute_link_scripts: whether to execute the post-link and pre-unlink scripts
                that may be part of a package. Defaults to False.
        show_progress: If set to `True` a progress bar will be shown on the CLI.
                Ignored when `reporter` is provided.
        client: An authenticated client to use for downloading packages. If not specified a default
                client will be used.
        requested_specs: A list of `MatchSpec`s that were originally requested. These will be used
                to populate the `requested_specs` field in the generated `conda-meta/*.json` files.
                If `None`, the `requested_specs` field will remain empty.
        reporter: An optional :class:`InstallerReporter` instance that receives progress
                callbacks during installation. When provided, `show_progress` is ignored.
                Subclass :class:`InstallerReporter` and override the methods you need.
    """

    await py_install(
        records=records,
        target_prefix=str(target_prefix),
        cache_dir=cache_dir,
        installed_packages=installed_packages,
        reinstall_packages=reinstall_packages,
        ignored_packages=ignored_packages,
        platform=platform._inner if platform is not None else None,
        client=client._client if client is not None else None,
        execute_link_scripts=execute_link_scripts,
        show_progress=show_progress,
        requested_specs=requested_specs,
        reporter=reporter,
    )
