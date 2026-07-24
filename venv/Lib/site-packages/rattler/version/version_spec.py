from __future__ import annotations

from typing import TYPE_CHECKING

from rattler.rattler import PyVersionSpec

if TYPE_CHECKING:
    from rattler.version.version import Version


class VersionSpec:
    """
    A version specification represents a constraint on a version.

    Version specifications can be simple like ">=1.2.3" or complex
    with multiple constraints like ">=1.2.3,<2.0.0" or ">=1.2.3|<1.0.0".

    Examples
    --------
    ```python
    >>> from rattler import VersionSpec, Version
    >>> spec = VersionSpec(">=1.2.3,<2.0.0")
    >>> version = Version("1.5.0")
    >>> spec.matches(version)
    True
    >>> version2 = Version("2.1.0")
    >>> spec.matches(version2)
    False
    >>>
    ```
    """

    _spec: PyVersionSpec

    def __init__(self, spec: str, strict: bool = False) -> None:
        """
        Construct a new VersionSpec from a string.

        Parameters
        ----------
        spec : str
            The version specification string (e.g., ">=1.2.3,<2.0.0")
        strict : bool, optional
            If True, use strict parsing mode (default: False)

        Raises
        ------
        InvalidVersionSpecError
            If the version specification string is invalid

        Examples
        --------
        ```python
        >>> spec = VersionSpec(">=1.2.3,<2.0.0")
        >>> spec_strict = VersionSpec(">=1.2.3", strict=True)
        >>> spec_any = VersionSpec("*")
        >>>
        ```
        """
        if isinstance(spec, str):
            self._spec = PyVersionSpec(spec, strict)
        else:
            raise TypeError(
                f"VersionSpec constructor received unsupported type {type(spec).__name__!r} for the `spec` parameter"
            )

    @classmethod
    def _from_py_version_spec(cls, py_spec: PyVersionSpec) -> VersionSpec:
        """Construct Rattler VersionSpec from FFI PyVersionSpec object."""
        spec = cls.__new__(cls)
        spec._spec = py_spec
        return spec

    def matches(self, version: Version) -> bool:
        """
        Check if a version matches this version specification.

        Parameters
        ----------
        version : Version
            The version to check

        Returns
        -------
        bool
            True if the version matches the specification, False otherwise

        Examples
        --------
        ```python
        >>> from rattler import VersionSpec, Version
        >>> spec = VersionSpec(">=1.2.3,<2.0.0")
        >>> spec.matches(Version("1.5.0"))
        True
        >>> spec.matches(Version("2.1.0"))
        False
        >>> spec.matches(Version("1.2.3"))
        True
        >>> spec.matches(Version("2.0.0"))
        False
        >>>
        ```
        """
        return self._spec.matches(version._version)

    def __str__(self) -> str:
        """
        Returns the string representation of the version specification.

        Examples
        --------
        ```python
        >>> str(VersionSpec(">=1.2.3,<2.0.0"))
        '>=1.2.3,<2.0.0'
        >>>
        ```
        """
        return self._spec.as_str()

    def __repr__(self) -> str:
        """
        Returns a representation of the version specification.

        Examples
        --------
        ```python
        >>> VersionSpec(">=1.2.3,<2.0.0")
        VersionSpec(">=1.2.3,<2.0.0")
        >>>
        ```
        """
        return f'VersionSpec("{self._spec.as_str()}")'

    def __hash__(self) -> int:
        """
        Computes the hash of this instance.

        Examples
        --------
        ```python
        >>> hash(VersionSpec(">=1.2.3")) == hash(VersionSpec(">=1.2.3"))
        True
        >>> hash(VersionSpec(">=1.2.3")) == hash(VersionSpec(">=2.0.0"))
        False
        >>>
        ```
        """
        return self._spec.__hash__()

    def __eq__(self, other: object) -> bool:
        """
        Returns True if this instance represents the same version spec as `other`.

        Examples
        --------
        ```python
        >>> VersionSpec(">=1.2.3") == VersionSpec(">=1.2.3")
        True
        >>> VersionSpec(">=1.2.3") == VersionSpec(">=2.0.0")
        False
        >>>
        ```
        """
        if not isinstance(other, VersionSpec):
            return NotImplemented
        return self._spec == other._spec

    def __ne__(self, other: object) -> bool:
        """
        Returns False if this instance represents the same version spec as `other`.

        Examples
        --------
        ```python
        >>> VersionSpec(">=1.2.3") != VersionSpec(">=1.2.3")
        False
        >>> VersionSpec(">=1.2.3") != VersionSpec(">=2.0.0")
        True
        >>>
        ```
        """
        if not isinstance(other, VersionSpec):
            return NotImplemented
        return self._spec != other._spec
