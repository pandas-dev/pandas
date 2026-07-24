from __future__ import annotations

from rattler.rattler import PyRecord
from rattler.repo_data.package_record import PackageRecord


class WhlPackageRecord(PackageRecord):
    """
    A wheel package record pairing a PackageRecord with its URL or path.

    Used to build repodata from PyPI/wheel metadata without conda archives.
    """

    def __init__(self, package_record: PackageRecord, url: str) -> None:
        self._record = PyRecord.create_whl_record(package_record._record, url)

    @property
    def url(self) -> str:
        return self._record.url

    def __repr__(self) -> str:
        """Returns a representation of the WhlPackageRecord."""
        return f'WhlPackageRecord(url="{self.url}")'
