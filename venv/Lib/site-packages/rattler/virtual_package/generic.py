from __future__ import annotations

from rattler.version import Version
from rattler.package import PackageName

from rattler.rattler import PyGenericVirtualPackage


class GenericVirtualPackage:
    _generic_virtual_package: PyGenericVirtualPackage

    def __init__(self, name: PackageName, version: Version, build_string: str) -> None:
        if not isinstance(name, PackageName):
            raise TypeError(
                "GenericVirtualPackage constructor received unsupported type "
                f" {type(name).__name__!r} for the `name` parameter"
            )

        if not isinstance(version, Version):
            raise TypeError(
                "GenericVirtualPackage constructor received unsupported type "
                f" {type(version).__name__!r} for the `version` parameter"
            )

        if not isinstance(build_string, str):
            raise TypeError(
                "GenericVirtualPackage constructor received unsupported type "
                f" {type(build_string).__name__!r} for the `build_string` parameter"
            )

        self._generic_virtual_package = PyGenericVirtualPackage(name._name, version._version, build_string)

    @property
    def name(self) -> PackageName:
        """
        Returns the name of the package

        Examples
        --------
        ```python
        >>> from rattler.package.package_name import PackageName
        >>> from rattler.version.version import Version
        >>> gvp = GenericVirtualPackage(PackageName("__archspec"), Version("1"), "x86_64")
        >>> gvp.name
        PackageName("__archspec")
        >>> gvp.name.source
        '__archspec'
        >>> gvp.name.normalized
        '__archspec'
        >>>
        ```
        """
        return PackageName._from_py_package_name(self._generic_virtual_package.name)

    @property
    def version(self) -> Version:
        """
        Returns the version of the package

        Examples
        --------
        ```python
        >>> from rattler.package.package_name import PackageName
        >>> from rattler.version.version import Version
        >>> gvp = GenericVirtualPackage(PackageName("__archspec"), Version("1"), "x86_64")
        >>> gvp.version
        Version("1")
        >>>
        ```
        """
        return Version._from_py_version(self._generic_virtual_package.version)

    @property
    def build_string(self) -> str:
        """
        Returns the build identifier of the package.

        Examples
        --------
        ```python
        >>> from rattler.package.package_name import PackageName
        >>> from rattler.version.version import Version
        >>> gvp = GenericVirtualPackage(PackageName("__archspec"), Version("1"), "x86_64")
        >>> gvp.build_string
        'x86_64'
        >>>
        ```
        """
        return self._generic_virtual_package.build_string

    @classmethod
    def _from_py_generic_virtual_package(
        cls, py_generic_virtual_package: PyGenericVirtualPackage
    ) -> GenericVirtualPackage:
        """
        Construct Rattler GenericVirtualPackage from FFI
        PyGenericVirtualPackage object.
        """
        generic_virtual_package = cls.__new__(cls)
        generic_virtual_package._generic_virtual_package = py_generic_virtual_package
        return generic_virtual_package

    def __str__(self) -> str:
        """
        Returns the string representation of the GenericVirtualPackage

        Examples
        --------
        ```python
        >>> from rattler.package.package_name import PackageName
        >>> from rattler.version.version import Version
        >>> gvp = GenericVirtualPackage(PackageName("__archspec"), Version("1"), "x86_64")
        >>> str(gvp)
        '__archspec=1=x86_64'
        >>>
        ```
        """

        return self._generic_virtual_package.as_str()

    def __repr__(self) -> str:
        """
        Returns a representation of the GenericVirtualPackage

        Examples
        --------
        ```python
        >>> from rattler.package.package_name import PackageName
        >>> from rattler.version.version import Version
        >>> gvp = GenericVirtualPackage(PackageName("__archspec"), Version("1"), "x86_64")
        >>> gvp
        GenericVirtualPackage("__archspec=1=x86_64")
        >>>
        ```
        """
        return f'GenericVirtualPackage("{self._generic_virtual_package.as_str()}")'
