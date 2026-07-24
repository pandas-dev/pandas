from __future__ import annotations

from rattler.rattler import PyArch

from typing import Literal

ArchLiteral = Literal[
    "x86",
    "x86_64",
    "aarch64",
    "armv6l",
    "armv7l",
    "loongarch64",
    "ppc64le",
    "ppc64",
    "s390x",
    "riscv32",
    "riscv64",
]


class Arch:
    def __init__(self, value: ArchLiteral) -> None:
        self._inner = PyArch(value)

    @classmethod
    def _from_py_arch(cls, py_arch: PyArch) -> Arch:
        """Construct Rattler version from FFI PyArch object."""
        arch = cls.__new__(cls)
        arch._inner = py_arch
        return arch

    def __str__(self) -> str:
        """
        Returns a string representation of the architecture.

        Examples
        --------
        ```python
        >>> str(Arch("x86_64"))
        'x86_64'
        >>>
        ```
        """
        return self._inner.as_str()

    def __repr__(self) -> str:
        """
        Returns a representation of the architecture.

        Examples
        --------
        ```python
        >>> Arch("aarch64")
        Arch(aarch64)
        >>>
        ```
        """
        return f"Arch({self._inner.as_str()})"
