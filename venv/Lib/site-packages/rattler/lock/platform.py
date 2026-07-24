from __future__ import annotations
from typing import List, Optional

from rattler.platform.platform import Platform
from rattler.rattler import PyLockPlatform


class LockPlatform:
    """
    Represents a platform in a lock file.

    This provides access to the platform name, the underlying conda subdir,
    and any virtual packages associated with the platform.
    """

    _inner: PyLockPlatform

    def __init__(
        self,
        name: str,
        subdir: Optional[Platform] = None,
        virtual_packages: Optional[List[str]] = None,
    ) -> None:
        """
        Create a new lock platform.

        Args:
            name: The name of the platform (e.g., "linux-64", "osx-arm64").
            subdir: The underlying conda platform. If not provided, it is
                    automatically determined from the name.
            virtual_packages: The list of virtual packages for this platform.
                              Defaults to an empty list.

        Examples
        --------
        ```python
        >>> from rattler import LockPlatform
        >>> platform = LockPlatform("linux-64")
        >>> platform.name
        'linux-64'
        >>>
        ```
        """
        self._inner = PyLockPlatform(
            name,
            subdir._inner if subdir else None,
            virtual_packages,
        )

    @property
    def name(self) -> str:
        """
        The name of the platform (e.g., "linux-64", "osx-arm64").

        Examples
        --------
        ```python
        >>> from rattler import LockPlatform
        >>> platform = LockPlatform("linux-64")
        >>> platform.name
        'linux-64'
        >>>
        ```
        """
        return self._inner.name

    @property
    def subdir(self) -> Platform:
        """
        The underlying conda subdir/platform.

        Examples
        --------
        ```python
        >>> from rattler import LockPlatform
        >>> platform = LockPlatform("linux-64")
        >>> platform.subdir
        Platform(linux-64)
        >>>
        ```
        """
        return Platform._from_py_platform(self._inner.subdir)

    @property
    def virtual_packages(self) -> List[str]:
        """
        The list of virtual packages for this platform.

        Examples
        --------
        ```python
        >>> from rattler import LockPlatform
        >>> platform = LockPlatform("linux-64", virtual_packages=["__glibc=2.17"])
        >>> platform.virtual_packages
        ['__glibc=2.17']
        >>>
        ```
        """
        return self._inner.virtual_packages

    @classmethod
    def _from_py_lock_platform(cls, py_lock_platform: PyLockPlatform) -> LockPlatform:
        """
        Construct a LockPlatform from FFI PyLockPlatform object.
        """
        platform = cls.__new__(cls)
        platform._inner = py_lock_platform
        return platform

    def __repr__(self) -> str:
        """
        Returns a representation of the LockPlatform.
        """
        return f"LockPlatform(name='{self.name}')"

    def __str__(self) -> str:
        """
        Returns a string representation of the LockPlatform.
        """
        return self.name
