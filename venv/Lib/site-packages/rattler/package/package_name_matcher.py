from __future__ import annotations

from typing import Union

from rattler.package.package_name import PackageName
from rattler.rattler import PyPackageNameMatcher


class PackageNameMatcher:
    """
    A class representing a package name matcher.

    Examples
    --------
    ```python
    >>> PackageNameMatcher("rattler")
    PackageNameMatcher("rattler", exact)
    >>> PackageNameMatcher("jupyter-*")
    PackageNameMatcher("jupyter-*", glob)
    >>> PackageNameMatcher("^jupyter-.*$")
    PackageNameMatcher("^jupyter-.*$", regex)
    >>>
    ```
    """

    _package_name_matcher: PyPackageNameMatcher

    def __init__(self, package_name_matcher: str):
        self._package_name_matcher = PyPackageNameMatcher(package_name_matcher)

    def __repr__(self) -> str:
        inner = self._package_name_matcher.display_inner()
        return f"{type(self).__name__}({inner})"

    @classmethod
    def _from_py_package_name_matcher(cls, py_package_name_matcher: PyPackageNameMatcher) -> PackageNameMatcher:
        """Construct Rattler PackageNameMatcher from FFI PyPackageName object."""
        package_name_matcher = cls.__new__(cls)
        package_name_matcher._package_name_matcher = py_package_name_matcher
        return package_name_matcher

    @property
    def normalized(self) -> str:
        """
        Returns the normalized string representation of the matcher.

        For exact matches, returns the normalized package name.
        For glob and regex patterns, returns the pattern string.

        Examples
        --------
        ```python
        >>> PackageNameMatcher("rattler").normalized
        'rattler'
        >>> PackageNameMatcher("jupyter-*").normalized
        'jupyter-*'
        >>> PackageNameMatcher("^jupyter-.*$").normalized
        '^jupyter-.*$'
        >>>
        ```
        """
        return self._package_name_matcher.normalized

    def as_package_name(self) -> Union[PackageName, None]:
        """
        Converts a PackageNameMatcher to a PackageName if it is an exact matcher.

        Examples
        --------
        ```python
        >>> PackageNameMatcher("rattler").as_package_name()
        PackageName("rattler")
        >>> PackageNameMatcher("jupyter-*").as_package_name()
        >>> PackageNameMatcher("^jupyter-.*$").as_package_name()
        >>>
        ```
        """
        py_package_name = self._package_name_matcher.as_package_name()
        if py_package_name is None:
            return None
        return PackageName._from_py_package_name(py_package_name)
