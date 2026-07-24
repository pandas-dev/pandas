from __future__ import annotations
from typing import Optional
from rattler.rattler import PyPackageHashes


class PackageHashes:
    _hashes: PyPackageHashes

    @property
    def md5(self) -> Optional[bytes]:
        """
        Returns the Sha256 hash.
        """
        return self._hashes.md5

    @property
    def sha256(self) -> Optional[bytes]:
        """
        Returns the Sha256 hash.
        """
        return self._hashes.sha256

    @classmethod
    def _from_py_package_hashes(cls, pkg_hashes: PyPackageHashes) -> PackageHashes:
        """
        Construct Rattler PackageHashes from FFI PyPackageHashes object.
        """
        hashes = cls.__new__(cls)
        hashes._hashes = pkg_hashes
        return hashes

    def __repr__(self) -> str:
        return "PackageHashes()"
