from __future__ import annotations
from collections.abc import Iterator
from typing import Any, Dict, Literal, Tuple, Optional

from rattler.rattler import PyPlatform
from rattler.platform.arch import Arch

PlatformLiteral = Literal[
    "noarch",
    "unknown",
    "linux-32",
    "linux-64",
    "linux-aarch64",
    "linux-armv6l",
    "linux-armv7l",
    "linux-loongarch64",
    "linux-ppc64le",
    "linux-ppc64",
    "linux-ppc",
    "linux-s390x",
    "linux-riscv32",
    "linux-riscv64",
    "freebsd-32",
    "freebsd-64",
    "freebsd-arm64",
    "osx-64",
    "osx-arm64",
    "win-32",
    "win-64",
    "win-arm64",
    "emscripten-wasm32",
    "wasi-wasm32",
    "zos-z",
]


class PlatformSingleton(type):
    _instances: Dict[str, Platform]

    def __init__(cls, *args: Tuple[Any], **kwargs: Dict[Any, Any]) -> None:
        cls._instances = {}

    def __call__(cls, platform: str, *args: Tuple[Any], **kwargs: Dict[Any, Any]) -> Platform:
        try:
            return cls._instances[platform]
        except KeyError:
            pass

        instance = super().__call__(platform, *args, **kwargs)
        cls._instances[platform] = instance
        return instance


class Platform(metaclass=PlatformSingleton):
    def __init__(self, value: PlatformLiteral | str):
        self._inner = PyPlatform(value)

    @classmethod
    def _from_py_platform(cls, py_platform: PyPlatform) -> Platform:
        """Construct Rattler version from FFI PyArch object."""
        try:
            platform = cls._instances[py_platform.name]
        except KeyError:
            platform = cls.__new__(cls)
            platform._inner = py_platform
            cls._instances[str(platform)] = platform
        return platform

    def __str__(self) -> str:
        """
        Returns a string representation of the platform.

        Examples
        --------
        ```python
        >>> str(Platform("linux-64"))
        'linux-64'
        >>>
        ```
        """
        return self._inner.name

    def __repr__(self) -> str:
        """
        Returns a representation of the platform.

        Examples
        --------
        ```python
        >>> Platform("linux-64")
        Platform(linux-64)
        >>>
        ```
        """
        return f"Platform({self._inner.name})"

    @classmethod
    def current(cls) -> Platform:
        """
        Returns the current platform.
        """
        return cls._from_py_platform(PyPlatform.current())

    @classmethod
    def all(cls) -> Iterator[Platform]:
        """
        Returns all supported platforms.

        Examples
        --------
        ```python
        >>> next(Platform.all())
        Platform(noarch)
        >>> len(list(Platform.all()))
        25
        >>>
        """
        return (cls._from_py_platform(p) for p in PyPlatform.all())

    @property
    def is_linux(self) -> bool:
        """
        Return True if the platform is linux.

        Examples
        --------
        ```python
        >>> Platform("linux-64").is_linux
        True
        >>> Platform("osx-64").is_linux
        False
        >>>
        ```
        """
        return self._inner.is_linux

    @property
    def is_osx(self) -> bool:
        """
        Return True if the platform is osx.

        Examples
        --------
        ```python
        >>> Platform("osx-64").is_osx
        True
        >>> Platform("linux-64").is_osx
        False
        >>>
        ```
        """
        return self._inner.is_osx

    @property
    def is_windows(self) -> bool:
        """
        Return True if the platform is win.

        Examples
        --------
        ```python
        >>> Platform("win-64").is_windows
        True
        >>> Platform("linux-64").is_windows
        False
        >>>
        ```
        """
        return self._inner.is_windows

    @property
    def is_unix(self) -> bool:
        """
        Return True if the platform is unix.

        Examples
        --------
        ```python
        >>> Platform("linux-64").is_unix
        True
        >>> Platform("win-64").is_unix
        False
        >>>
        ```
        """
        return self._inner.is_unix

    @property
    def arch(self) -> Optional[Arch]:
        """
        Return the architecture of the platform.

        Examples
        --------
        ```python
        >>> Platform("linux-64").arch
        Arch(x86_64)
        >>> Platform("linux-aarch64").arch
        Arch(aarch64)
        >>>
        ```
        """
        arch = self._inner.arch()
        return Arch._from_py_arch(arch) if arch is not None else None

    @property
    def only_platform(self) -> Optional[str]:
        """
        Return the platform without the architecture.

        Examples
        --------
        ```python
        >>> Platform("linux-64").only_platform
        'linux'
        >>>
        ```
        """
        return self._inner.only_platform
