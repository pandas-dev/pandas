from __future__ import annotations

from typing import List, Optional, Tuple, Union

from rattler.rattler import PyVersion, InvalidVersionError


class Version:
    """
    This class implements an order relation between version strings.
    Version strings can contain the usual alphanumeric characters
    (A-Za-z0-9), separated into segments by dots and underscores.
    Empty segments (i.e. two consecutive dots, a leading/trailing
    underscore) are not permitted. An optional epoch number - an
    integer followed by '!' - can precede the actual version string
    (this is useful to indicate a change in the versioning scheme itself).
    Version comparison is case-insensitive.
    """

    _version: PyVersion

    def __init__(self, version: str) -> None:
        if isinstance(version, str):
            self._version = PyVersion(version)
        else:
            raise TypeError(
                f"Version constructor received unsupported type  {type(version).__name__!r} for the `version` parameter"
            )

    @classmethod
    def _from_py_version(cls, py_version: PyVersion) -> Version:
        """Construct Rattler version from FFI PyVersion object."""
        version = cls.__new__(cls)
        version._version = py_version
        return version

    @property
    def epoch(self) -> Optional[str]:
        """
        Gets the epoch of the version or `None` if the epoch was not defined.

        Examples
        --------
        ```python
        >>> v = Version('2!1.0')
        >>> v.epoch
        2
        >>>
        ```
        """
        return self._version.epoch()

    def bump_major(self) -> Version:
        """
        Returns a new version where the major segment of this version has
        been bumped.

        Examples
        --------
        ```python
        >>> v = Version('1.0')
        >>> v.bump_major()
        Version("2.0")
        >>> v = Version('9d')
        >>> v.bump_major()
        Version("10a")
        >>>
        ```
        """
        return Version._from_py_version(self._version.bump_major())

    def bump_minor(self) -> Version:
        """
        Returns a new version where the minor segment of this version has
        been bumped.

        Examples
        --------
        ```python
        >>> v = Version('1.0')
        >>> v.bump_minor()
        Version("1.1")
        >>>
        >>> Version("1").bump_minor()
        Version("1.1")
        >>>
        ```
        """
        return Version._from_py_version(self._version.bump_minor())

    def bump_patch(self) -> Version:
        """
        Returns a new version where the patch segment of this version has
        been bumped.

        Examples
        --------
        ```python
        >>> v = Version('1.0.5')
        >>> v.bump_patch()
        Version("1.0.6")
        >>> v = Version('1.1.1e')
        >>> v.bump_patch()
        Version("1.1.2a")
        >>>
        >>> Version("1.5").bump_patch()
        Version("1.5.1")
        >>>
        ```
        """
        return Version._from_py_version(self._version.bump_patch())

    def bump_last(self) -> Version:
        """
        Returns a new version where the last segment of this version has
        been bumped.

        Examples
        --------
        ```python
        >>> v = Version('1.0')
        >>> v.bump_last()
        Version("1.1")
        >>>
        ```
        """
        return Version._from_py_version(self._version.bump_last())

    def bump_segment(self, index: int) -> Version:
        """
        Returns a new version where the last segment of this version has
        been bumped.

        Examples
        --------
        ```python
        >>> v = Version('1.0')
        >>> v.bump_segment(index=1)
        Version("1.1")
        >>>
        >>> Version("1.5").bump_segment(-5) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        exceptions.VersionBumpError
        >>> Version("1.5").bump_segment(5)
        Version("1.5.0.0.0.1")
        >>>
        ```
        """
        return Version._from_py_version(self._version.bump_segment(index))

    def with_alpha(self) -> Version:
        """
        Returns a new version where the last segment of this version has
        been bumped with an alpha character. If the last segment contains a
        character, nothing is added.

        Examples
        --------
        ```python
        >>> v = Version('1.0')
        >>> v.with_alpha()
        Version("1.0.0a0")
        >>> v = Version('1.0.f')
        >>> v.with_alpha()
        Version("1.0.f")
        >>>
        ```
        """
        return Version._from_py_version(self._version.with_alpha())

    def remove_local(self) -> Version:
        """
        Returns a new version where the local segment of the version has been removed.
        Leaves the version unchanged if it does not have a local segment.

        Examples
        --------
        ```python
        >>> v = Version('1.0+3.4')
        >>> v.remove_local()
        Version("1.0")
        >>> v = Version('1.0')
        >>> v.remove_local()
        Version("1.0")
        >>>
        ```
        """
        return Version._from_py_version(self._version.remove_local())

    def extend_to_length(self, length: int) -> Version:
        """
        Returns a new version that is extended with `0s` to the specified length.

        Examples
        --------
        ```python
        >>> v = Version('1')
        >>> v.extend_to_length(3)
        Version("1.0.0")
        >>> v = Version('4!1.2+3.4')
        >>> v.extend_to_length(4)
        Version("4!1.2.0.0+3.4")
        >>>
        ```
        """
        return Version._from_py_version(self._version.extend_to_length(length))

    @property
    def has_local(self) -> bool:
        """
        Returns true if this version has a local segment defined.
        The local part of a version is the part behind the (optional) `+`.

        Examples
        --------
        ```python
        >>> v = Version('1.0+3.2-alpha0')
        >>> v.has_local
        True
        >>> v2 = Version('1.0')
        >>> v2.has_local
        False
        >>>
        ```
        """
        return self._version.has_local()

    def segments(self) -> List[List[Union[str, int]]]:
        """
        Returns a list of segments of the version. It does not contain
        the local segment of the version.

        Examples
        --------
        ```python
        >>> v = Version("1.2dev.3-alpha4.5+6.8")
        >>> v.segments()
        [[1], [2, 'dev'], [3], [0, 'alpha', 4], [5]]
        >>>
        ```
        """
        return self._version.segments()

    def local_segments(self) -> List[List[Union[str, int]]]:
        """
        Returns a list of local segments of the version. It does not
        contain the non-local segment of the version.

        Examples
        --------
        ```python
        >>> v = Version("1.2dev.3-alpha4.5+6.8")
        >>> v.local_segments()
        [[6], [8]]
        >>>
        ```
        """
        return self._version.local_segments()

    def as_major_minor(self) -> Optional[Tuple[int, int]]:
        """
        Returns the major and minor segments from the version.
        Requires a minimum of 2 segments in version to be split
        into major and minor, returns `None` otherwise.

        Examples
        --------
        ```python
        >>> v = Version('1.0')
        >>> v.as_major_minor()
        (1, 0)
        >>>
        ```
        """
        return self._version.as_major_minor()

    @property
    def is_dev(self) -> bool:
        """
        Returns true if the version contains a component name "dev",
        dev versions are sorted before non-dev version.

        Examples
        --------
        ```python
        >>> v = Version('1.0.1dev')
        >>> v.is_dev
        True
        >>> v_non_dev = Version('1.0.1')
        >>> v_non_dev >= v
        True
        >>>
        ```
        """
        return self._version.is_dev()

    def starts_with(self, other: Version) -> bool:
        """
        Checks if the version and local segment start
        same as other version.

        Examples
        --------
        ```python
        >>> v1 = Version('1.0.1')
        >>> v2 = Version('1.0')
        >>> v1.starts_with(v2)
        True
        >>>
        ```
        """
        return self._version.starts_with(other._version)

    def compatible_with(self, other: Version) -> bool:
        """
        Checks if this version is compatible with other version.
        Minor versions changes are compatible with older versions,
        major version changes are breaking and will not be compatible.

        Examples
        --------
        ```python
        >>> v1 = Version('1.0')
        >>> v2 = Version('1.2')
        >>> v_major = Version('2.0')
        >>> v1.compatible_with(v2)
        False
        >>> v2.compatible_with(v1)
        True
        >>> v_major.compatible_with(v2)
        False
        >>> v2.compatible_with(v_major)
        False
        >>>
        ```
        """
        return self._version.compatible_with(other._version)

    def pop_segments(self, n: int = 1) -> Version:
        """
        Pops `n` number of segments from the version and returns
        the new version. Raises `InvalidVersionError` if version
        becomes invalid due to the operation.

        Examples
        --------
        ```python
        >>> v = Version('2!1.0.1')
        >>> v.pop_segments() # `n` defaults to 1 if left empty
        Version("2!1.0")
        >>> v.pop_segments(2) # old version is still usable
        Version("2!1")
        >>> v.pop_segments(3) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        exceptions.InvalidVersionError: new Version must have atleast 1 valid
        segment
        >>>
        ```
        """
        new_py_version = self._version.pop_segments(n)
        if new_py_version:
            return self._from_py_version(new_py_version)
        else:
            raise InvalidVersionError("new Version must have atleast 1 valid segment")

    def with_segments(self, start: int, stop: int) -> Version:
        """
        Returns new version with with segments ranging from `start` to `stop`.
        `stop` is exclusive. Raises `InvalidVersionError` if the provided range
        is invalid.

        Examples
        --------
        ```python
        >>> v = Version('2!1.2.3')
        >>> v.with_segments(0, 2)
        Version("2!1.2")
        >>>
        ```
        """
        new_py_version = self._version.with_segments(start, stop)
        if new_py_version:
            return self._from_py_version(new_py_version)
        else:
            raise InvalidVersionError("Invalid segment range provided")

    @property
    def segment_count(self) -> int:
        """
        Returns the number of segments in the version.
        This does not include epoch or local segment of the version

        Examples
        --------
        ```python
        >>> v = Version('2!1.2.3')
        >>> v.segment_count
        3
        >>>
        ```
        """
        return self._version.segment_count()

    def strip_local(self) -> Version:
        """
        Returns a new version with local segment stripped.

        Examples
        --------
        ```python
        >>> v = Version('1.2.3+4.alpha-5')
        >>> v.strip_local()
        Version("1.2.3")
        >>>
        ```
        """
        return self._from_py_version(self._version.strip_local())

    def __str__(self) -> str:
        """
        Returns the string representation of the version

        Examples
        --------
        ```python
        >>> str(Version("1.2.3"))
        '1.2.3'
        >>>
        ```
        """
        return self._version.as_str()

    def __repr__(self) -> str:
        """
        Returns a representation of the version

        Examples
        --------
        ```python
        >>> Version("1.2.3")
        Version("1.2.3")
        >>>
        ```
        """
        return f'Version("{self._version.as_str()}")'

    def __hash__(self) -> int:
        """
        Computes the hash of this instance.

        Examples
        --------
        ```python
        >>> hash(Version("1.2.3")) == hash(Version("1.2.3"))
        True
        >>> hash(Version("1.2.3")) == hash(Version("3.2.1"))
        False
        >>> hash(Version("1")) == hash(Version("1.0.0"))
        True
        >>>
        ```
        """
        return self._version.__hash__()

    def __eq__(self, other: Version) -> bool:  # type: ignore[override]
        """
        Returns True if this instance represents the same version as `other`.

        Examples
        --------
        ```python
        >>> Version("1.2.3") == Version("1.2.3")
        True
        >>> Version("3.2.1") == Version("1.2.3")
        False
        >>> Version("1") == Version("1.0.0")
        True
        >>>
        ```
        """
        return self._version == other._version

    def __ne__(self, other: Version) -> bool:  # type: ignore[override]
        """
        Returns False if this instance represents the same version as `other`.

        Examples
        --------
        ```python
        >>> Version("1.2.3") != Version("1.2.3")
        False
        >>> Version("3.2.1") != Version("1.2.3")
        True
        >>> Version("1") != Version("1.0.0")
        False
        >>>
        ```
        """
        return self._version != other._version

    def __gt__(self, other: Version) -> bool:
        """
        Returns True if this instance should be ordered *after* `other`.

        Examples
        --------
        ```python
        >>> Version("1.2.3") > Version("1.2.3")
        False
        >>> Version("1.2.4") > Version("1.2.3")
        True
        >>> Version("1.2.3.1") > Version("1.2.3")
        True
        >>> Version("3.2.1") > Version("1.2.3")
        True
        >>> Version("1") > Version("1.0.0")
        False
        >>>
        ```
        """
        return self._version > other._version

    def __lt__(self, other: Version) -> bool:
        """
        Returns True if this instance should be ordered *before* `other`.

        Examples
        --------
        ```python
        >>> Version("1.2.3") < Version("1.2.3")
        False
        >>> Version("1.2.3") < Version("1.2.4")
        True
        >>> Version("1.2.3") < Version("1.2.3.1")
        True
        >>> Version("3.2.1") < Version("1.2.3")
        False
        >>> Version("1") < Version("1.0.0")
        False
        >>>
        ```
        """
        return self._version < other._version

    def __ge__(self, other: Version) -> bool:
        """
        Returns True if this instance should be ordered *after* or at the same location
        as `other`.

        Examples
        --------
        ```python
        >>> Version("1.2.3") >= Version("1.2.3")
        True
        >>> Version("1.2.4") >= Version("1.2.3")
        True
        >>> Version("1.2.3.1") >= Version("1.2.3")
        True
        >>> Version("3.2.1") >= Version("1.2.3")
        True
        >>> Version("1.2.3") >= Version("3.2.1")
        False
        >>> Version("1") >= Version("1.0.0")
        True
        >>>
        ```
        """
        return self._version >= other._version

    def __le__(self, other: Version) -> bool:
        """
        Returns True if this instance should be ordered *before* or at the same
        location as `other`.

        Examples
        --------
        ```python
        >>> Version("1.2.3") <= Version("1.2.3")
        True
        >>> Version("1.2.3") <= Version("1.2.4")
        True
        >>> Version("1.2.3") <= Version("1.2.3.1")
        True
        >>> Version("3.2.1") <= Version("1.2.3")
        False
        >>> Version("1") <= Version("1.0.0")
        True
        >>>
        ```
        """
        return self._version <= other._version
