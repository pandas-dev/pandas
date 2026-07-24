from __future__ import annotations

from typing import Union, Optional

from rattler.rattler import PyVersion

from rattler.version import Version


class VersionWithSource(Version):
    """
    Holds a version and the string it was created from. This is useful if
    you want to retain the original string the version was created from.
    This might be useful in cases where you have multiple strings that
    are represented by the same [`Version`] but you still want to be able to
    distinguish them.

    The string `1.1` and `1.01` represent the same version. When you print
    the parsed version though it will come out as `1.1`. You loose the
    original representation. This class stores the original source string,
    which can be accessed by `source` property.

    This class derives from `Version` which allows you to compare it with
    other `Version` objects.

    Comparison will be done based on the parsed version, not the source string.
    """

    _source: str

    def __init__(self, version: Union[str, Version]):
        if not isinstance(version, (str, Version)):
            raise TypeError(
                "VersionWithSource constructor received unsupported type "
                f" {type(version).__name__!r} for the `version` parameter"
            )

        if isinstance(version, str):
            self._source = version
            self._version = PyVersion(version)

        if isinstance(version, Version):
            self._source = str(version)
            self._version = version._version

    @classmethod
    def _from_py_version(cls, py_version: PyVersion, source: Optional[str] = None) -> VersionWithSource:
        """Construct Rattler version from FFI PyVersion object."""
        version = cls.__new__(cls)
        version._version = py_version
        version._source = source or py_version.as_str()
        return version

    def __str__(self) -> str:
        """
        Returns the string representation of the version

        Examples
        --------
        ```python
        >>> str(VersionWithSource("1.02.3"))
        '1.02.3'
        >>>
        ```
        """
        return self._source

    def __repr__(self) -> str:
        """
        Returns a representation of the version

        Examples
        --------
        ```python
        >>> VersionWithSource("1.02.3")
        VersionWithSource(version="1.2.3", source="1.02.3")
        >>>
        ```
        """
        return f'{type(self).__name__}(version="{self._version.as_str()}", source="{self._source}")'
