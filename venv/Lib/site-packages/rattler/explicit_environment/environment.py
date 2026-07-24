from pathlib import Path
from typing import List, Optional

from rattler.rattler import PyExplicitEnvironmentSpec as _PyExplicitEnvironmentSpec
from rattler.rattler import PyExplicitEnvironmentEntry as _PyExplicitEnvironmentEntry
from rattler.platform import Platform


class ExplicitEnvironmentEntry:
    """A wrapper around an explicit environment entry which represents a URL to a package"""

    def __init__(self, inner: _PyExplicitEnvironmentEntry) -> None:
        self._inner = inner

    @property
    def url(self) -> str:
        """Returns the URL of the package"""
        return self._inner.url()


class ExplicitEnvironmentSpec:
    """The explicit environment (e.g. env.txt) file that contains a list of all URLs in a environment"""

    def __init__(self, inner: _PyExplicitEnvironmentSpec) -> None:
        self._inner = inner

    @classmethod
    def from_path(cls, path: Path) -> "ExplicitEnvironmentSpec":
        """Parses the object from a file specified by a `path`, using a format appropriate for the file type.

        For example, if the file is in text format, this function reads the data from the file at
        the specified path, parses the text and returns the resulting object. If the file is
        not in a parsable format or if the file could not be read, this function raises an error.
        """
        return cls(_PyExplicitEnvironmentSpec.from_path(path))

    @classmethod
    def from_str(cls, content: str) -> "ExplicitEnvironmentSpec":
        """
        Parses the object from a string containing the explicit environment specification

        Examples:

        ```python
        >>> spec = ExplicitEnvironmentSpec.from_str('''@EXPLICIT
        ... # platform: linux-64
        ... http://repo.anaconda.com/pkgs/main/linux-64/python-3.9.0-h3.tar.bz2
        ... ''')
        >>> spec.platform
        Platform(linux-64)
        >>> spec.packages[0].url
        'http://repo.anaconda.com/pkgs/main/linux-64/python-3.9.0-h3.tar.bz2'
        >>>
        ```
        """
        return cls(_PyExplicitEnvironmentSpec.from_str(content))

    @property
    def platform(self) -> Optional[Platform]:
        """Returns the platform specified in the explicit environment specification"""
        platform = self._inner.platform()
        if platform is not None:
            return Platform._from_py_platform(platform)
        return None

    @property
    def packages(self) -> List[ExplicitEnvironmentEntry]:
        """Returns the environment entries (URLs) specified in the explicit environment specification"""
        return [ExplicitEnvironmentEntry(p) for p in self._inner.packages()]
