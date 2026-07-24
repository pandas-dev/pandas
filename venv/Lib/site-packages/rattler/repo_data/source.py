"""Protocol definition for custom repodata sources."""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Protocol, runtime_checkable

if TYPE_CHECKING:
    from rattler.platform.platform import Platform
    from rattler.package.package_name import PackageName
    from rattler.repo_data.record import RepoDataRecord


@runtime_checkable
class RepoDataSource(Protocol):
    """Protocol for custom repodata sources.

    Implement this protocol to provide repodata records from custom sources
    (databases, APIs, in-memory caches, etc.) that can be used alongside
    traditional conda channels.

    Note
    ----
    Unlike channels, custom sources are **not cached** by the gateway. The gateway's
    internal caching mechanisms only apply to channel data. If caching is needed for
    your custom source, you must implement it yourself within your source implementation.

    **Performance:** Custom sources are slower than channels because data must be
    marshalled between Python and Rust for each request. For performance-critical
    applications with large amounts of repodata, consider using channels when possible.

    Example
    -------
    ```python
    from rattler import Platform, PackageName, RepoDataRecord

    class MyCustomSource:
        async def fetch_package_records(
            self, platform: Platform, name: PackageName
        ) -> List[RepoDataRecord]:
            # Fetch records from your custom source
            return [...]

        def package_names(self, platform: Platform) -> List[str]:
            # Return all available package names for the platform
            return ["numpy", "pandas", ...]

    # Usage with Gateway
    gateway = Gateway()
    records = await gateway.query(
        sources=[channel, MyCustomSource()],  # Mix channels and custom sources
        platforms=["linux-64"],
        specs=["numpy"],
    )
    ```
    """

    async def fetch_package_records(self, platform: Platform, name: PackageName) -> List[RepoDataRecord]:
        """Fetch records for a specific package name and platform.

        This method is called by the gateway when it needs repodata records
        for a particular package. The platform parameter indicates which
        subdirectory the gateway is querying for.

        Args:
            platform: The platform to fetch records for (e.g., linux-64, noarch)
            name: The package name to fetch records for

        Returns:
            List of RepoDataRecord objects for the package
        """
        ...

    def package_names(self, platform: Platform) -> List[str]:
        """Return all available package names for the given platform.

        This is used by the gateway to know which packages are available
        in this source for a given platform/subdirectory.

        Args:
            platform: The platform to list packages for

        Returns:
            List of package name strings
        """
        ...
