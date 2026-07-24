from __future__ import annotations

from rattler.rattler import PyPackageName


class PackageName:
    def __init__(self, source: str) -> None:
        if not isinstance(source, str):
            raise TypeError(
                "PackageName constructor received unsupported type "
                f" {type(source).__name__!r} for the `source` parameter"
            )
        self._name = PyPackageName(source)

    @staticmethod
    def unchecked(normalized: str) -> PackageName:
        """
        Constructs a new `PackageName` from a string without checking if the string is actually a
        valid or normalized conda package name. This should only be used if you are sure that the
        input string is valid.

        Examples
        --------
        ```python
        >>> p = PackageName.unchecked("test_xyz")
        >>>
        ```
        """
        return PackageName._from_py_package_name(PyPackageName.new_unchecked(normalized))

    @staticmethod
    def from_matchspec_str(spec: str) -> PackageName:
        """
        Parses the package name part from a matchspec string without parsing the entire matchspec.

        This extracts the package name by splitting on whitespace or version constraint characters
        (`>`, `<`, `=`, `!`, `~`, `;`).

        Examples
        --------
        ```python
        >>> p = PackageName.from_matchspec_str("numpy>=1.0,<2.0")
        >>> p.source
        'numpy'
        >>> p = PackageName.from_matchspec_str("pillow >=10")
        >>> p.source
        'pillow'
        >>>
        ```
        """
        return PackageName._from_py_package_name(PyPackageName.from_matchspec_str(spec))

    @staticmethod
    def from_matchspec_str_unchecked(spec: str) -> PackageName:
        """
        Parses the package name part from a matchspec string without parsing the entire matchspec.
        This function assumes the matchspec string is a valid matchspec.

        This extracts the package name by splitting on whitespace or version constraint characters
        (`>`, `<`, `=`, `!`, `~`, `;`). The original capitalization is preserved in the source,
        while the normalized version is lowercase.

        Examples
        --------
        ```python
        >>> p = PackageName.from_matchspec_str_unchecked("Pillow >=10")
        >>> p.source
        'Pillow'
        >>> p.normalized
        'pillow'
        >>>
        ```
        """
        return PackageName._from_py_package_name(PyPackageName.from_matchspec_str_unchecked(spec))

    @classmethod
    def _from_py_package_name(cls, py_package_name: PyPackageName) -> PackageName:
        """Construct Rattler PackageName from FFI PyPackageName object."""
        package_name = cls.__new__(cls)
        package_name._name = py_package_name
        return package_name

    @property
    def source(self) -> str:
        """
        Returns the source representation of the package name.
        This is the string from which this instance was created.

        Examples
        --------
        ```python
        >>> p = PackageName("test-xyz")
        >>> p.source
        'test-xyz'
        >>>
        ```
        """
        return self._name.source

    @property
    def normalized(self) -> str:
        """
        Returns the normalized version of the package name.
        The normalized string is guaranteed to be a valid conda package name.

        Examples
        --------
        ```python
        >>> p = PackageName("test-xyz")
        >>> p.normalized
        'test-xyz'
        >>>
        ```
        """
        return self._name.normalized

    def __hash__(self) -> int:
        """
        Computes the hash of this instance.

        Examples
        --------
        ```python
        >>> hash(PackageName("test-abc")) == hash(PackageName("test-abc"))
        True
        >>> hash(PackageName("test-abc")) == hash(PackageName("test-ABC"))
        True
        >>> hash(PackageName("test-abc")) == hash(PackageName("abc-test"))
        False
        >>>
        ```
        """
        return self._name.__hash__()

    def __eq__(self, other: object) -> bool:
        """
        Returns True if this instance represents the same PackageName as `other`.

        Examples
        --------
        ```python
        >>> PackageName("test-abc") == PackageName("abc-test")
        False
        >>> PackageName("test-abc") == PackageName("test-abc")
        True
        >>> PackageName("test-abc") == PackageName("test-ABC")
        True
        >>> PackageName("test-abc") == "test-abc"
        True
        >>> PackageName("test-abc") == "not-test-abc"
        False
        >>>
        ```
        """
        if isinstance(other, str):
            return self._name == PyPackageName(other)

        if not isinstance(other, PackageName):
            return False

        return self._name == other._name

    def __ne__(self, other: object) -> bool:
        """
        Returns True if this instance does not represents the same PackageName as `other`.

        Examples
        --------
        ```python
        >>> PackageName("test-abc") != PackageName("test-abc")
        False
        >>> PackageName("test-abc") != PackageName("test-ABC")
        False
        >>> PackageName("test-abc") != PackageName("abc-test")
        True
        >>> PackageName("test-abc") != "test-abc"
        False
        >>> PackageName("test-abc") != "not-test-abc"
        True
        >>>
        ```
        """
        if isinstance(other, str):
            return self._name != PyPackageName(other)

        if not isinstance(other, PackageName):
            return True

        return self._name != other._name

    def __repr__(self) -> str:
        """
        Returns a representation of the PackageName.

        Examples
        --------
        ```python
        >>> p = PackageName("test-xyz")
        >>> p
        PackageName("test-xyz")
        >>>
        ```
        """
        return f'PackageName("{self.source}")'
