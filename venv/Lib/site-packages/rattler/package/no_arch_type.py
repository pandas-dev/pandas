from __future__ import annotations
from typing import Literal, Optional

from rattler.rattler import PyNoArchType

NoArchLiteral = Optional[Literal["python", "generic", True]]


class NoArchType:
    _noarch: PyNoArchType

    def __init__(self, noarch: NoArchLiteral = None) -> None:
        if noarch is None:
            self._noarch = PyNoArchType.none()
        elif noarch == "python":
            self._noarch = PyNoArchType.python()
        elif noarch == "generic" or noarch is True:
            self._noarch = PyNoArchType.generic()
        else:
            raise ValueError(f"NoArchType constructor received unsupported value {noarch} for the `noarch` parameter")

    @classmethod
    def _from_py_no_arch_type(cls, py_no_arch_type: PyNoArchType) -> NoArchType:
        """Construct Rattler NoArchType from FFI PyNoArchType object."""
        no_arch_type = cls.__new__(cls)
        no_arch_type._noarch = py_no_arch_type
        return no_arch_type

    @property
    def generic(self) -> bool:
        """
        Return whether this NoArchType is 'generic'
        >>> NoArchType('generic').generic
        True
        >>> NoArchType('generic').python
        False
        >>>
        """
        return self._noarch.is_generic

    @property
    def none(self) -> bool:
        """
        Return whether this NoArchType is set
        >>> NoArchType(None).none
        True
        >>> NoArchType(None).python
        False
        >>>
        """
        return self._noarch.is_none

    @property
    def python(self) -> bool:
        """
        Return whether this NoArchType is 'python'
        >>> NoArchType('python').python
        True
        >>> NoArchType('python').generic
        False
        >>>
        """
        return self._noarch.is_python

    def __hash__(self) -> int:
        """
        Computes the hash of this instance.

        Examples
        --------
        ```python
        >>> hash(NoArchType("python")) == hash(NoArchType("python"))
        True
        >>> hash(NoArchType("python")) == hash(NoArchType("generic"))
        False
        >>>
        ```
        """
        return self._noarch.__hash__()

    def __eq__(self, other: object) -> bool:
        """
        Returns True if this instance represents the same NoArchType as `other`.

        Examples
        --------
        ```python
        >>> NoArchType("python") == NoArchType("generic")
        False
        >>> NoArchType("python") == NoArchType("python")
        True
        >>> NoArchType("generic") == NoArchType("generic")
        True
        >>> NoArchType("python") == "python"
        False
        >>>
        ```
        """
        if not isinstance(other, NoArchType):
            return False

        return self._noarch == other._noarch

    def __ne__(self, other: object) -> bool:
        """
        Returns True if this instance does not represents the same NoArchType as `other`.

        Examples
        --------
        ```python
        >>> NoArchType("python") != NoArchType("python")
        False
        >>> NoArchType("python") != "python"
        True
        >>>
        ```
        """
        if not isinstance(other, NoArchType):
            return True

        return self._noarch != other._noarch

    def __repr__(self) -> str:
        """
        Returns a representation of the NoArchType.

        Examples
        --------
        ```python
        >>> p = NoArchType("python")
        >>> p
        NoArchType("python")
        >>>
        ```
        """

        if self._noarch.is_python:
            return 'NoArchType("python")'
        elif self._noarch.is_generic:
            return 'NoArchType("generic")'
        else:
            return "NoArchType(None)"
