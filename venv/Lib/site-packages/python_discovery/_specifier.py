"""Version specifier support using only standard library (PEP 440 compatible)."""

from __future__ import annotations

import contextlib
import operator
import re
import sys
from dataclasses import dataclass
from typing import TYPE_CHECKING, Final

_DC_KW = {"frozen": True, "kw_only": True, "slots": True} if sys.version_info >= (3, 10) else {"frozen": True}

if TYPE_CHECKING:
    from collections.abc import Iterator

_VERSION_RE: Final[re.Pattern[str]] = re.compile(
    r"""
    ^
    (\d+)               # major
    (?:\.(\d+))?        # optional minor
    (?:\.(\d+))?        # optional micro
    (?:(a|b|rc)(\d+))?  # optional pre-release suffix
    $
    """,
    re.VERBOSE,
)
_SPECIFIER_RE: Final[re.Pattern[str]] = re.compile(
    r"""
    ^
    (===|==|~=|!=|<=|>=|<|>)  # operator
    \s*
    (.+)                       # version string
    $
    """,
    re.VERBOSE,
)
_PRE_ORDER: Final[dict[str, int]] = {"a": 1, "b": 2, "rc": 3}


@dataclass(**_DC_KW)
class SimpleVersion:
    """
    Simple PEP 440-like version parser using only standard library.

    :param version_str: the original version string.
    :param major: major version number.
    :param minor: minor version number.
    :param micro: micro (patch) version number.
    :param pre_type: pre-release label (``"a"``, ``"b"``, or ``"rc"``), or ``None``.
    :param pre_num: pre-release sequence number, or ``None``.
    :param release: the ``(major, minor, micro)`` tuple.
    """

    version_str: str
    major: int
    minor: int
    micro: int
    pre_type: str | None
    pre_num: int | None
    release: tuple[int, int, int]

    @classmethod
    def from_string(cls, version_str: str) -> SimpleVersion:
        """
        Parse a PEP 440 version string (e.g. ``3.12.1``).

        :param version_str: the version string to parse.
        """
        stripped = version_str.strip()
        if not (match := _VERSION_RE.match(stripped)):
            msg = f"Invalid version: {version_str}"
            raise ValueError(msg)
        major = int(match.group(1))
        minor = int(match.group(2)) if match.group(2) else 0
        micro = int(match.group(3)) if match.group(3) else 0
        return cls(
            version_str=stripped,
            major=major,
            minor=minor,
            micro=micro,
            pre_type=match.group(4),
            pre_num=int(match.group(5)) if match.group(5) else None,
            release=(major, minor, micro),
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        return self.release == other.release and self.pre_type == other.pre_type and self.pre_num == other.pre_num

    def __hash__(self) -> int:
        return hash((self.release, self.pre_type, self.pre_num))

    def __lt__(self, other: object) -> bool:  # ruff:ignore[too-many-return-statements]
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        if self.release != other.release:
            return self.release < other.release
        if self.pre_type is None and other.pre_type is None:
            return False
        if self.pre_type is None:
            return False
        if other.pre_type is None:
            return True
        if _PRE_ORDER[self.pre_type] != _PRE_ORDER[other.pre_type]:
            return _PRE_ORDER[self.pre_type] < _PRE_ORDER[other.pre_type]
        return (self.pre_num or 0) < (other.pre_num or 0)

    def __le__(self, other: object) -> bool:
        return self == other or self < other

    def __gt__(self, other: object) -> bool:
        if not isinstance(other, SimpleVersion):
            return NotImplemented
        return not self <= other

    def __ge__(self, other: object) -> bool:
        return not self < other

    def __str__(self) -> str:
        return self.version_str

    def __repr__(self) -> str:
        return f"SimpleVersion('{self.version_str}')"


@dataclass(**_DC_KW)
class SimpleSpecifier:
    """
    Simple PEP 440-like version specifier using only standard library.

    :param spec_str: the original specifier string (e.g. ``>=3.10``).
    :param operator: the comparison operator (``==``, ``>=``, ``<``, etc.).
    :param version_str: the version portion of the specifier, without the operator.
    :param is_wildcard: ``True`` if the specifier uses a wildcard suffix (``.*``).
    :param wildcard_precision: number of version components before the wildcard, or ``None``.
    :param version: the parsed version, or ``None`` if parsing failed.
    """

    spec_str: str
    operator: str
    version_str: str
    is_wildcard: bool
    wildcard_precision: int | None
    version: SimpleVersion | None

    @classmethod
    def from_string(cls, spec_str: str) -> SimpleSpecifier:
        """
        Parse a single PEP 440 specifier (e.g. ``>=3.10``).

        :param spec_str: the specifier string to parse.
        """
        stripped = spec_str.strip()
        if not (match := _SPECIFIER_RE.match(stripped)):
            msg = f"Invalid specifier: {spec_str}"
            raise ValueError(msg)
        op = match.group(1)
        version_str = match.group(2).strip()
        is_wildcard = version_str.endswith(".*")
        wildcard_precision: int | None = None
        if is_wildcard:
            version_str = version_str[:-2]
            wildcard_precision = len(version_str.split("."))
        try:
            version = SimpleVersion.from_string(version_str)
        except ValueError:
            version = None
        return cls(
            spec_str=stripped,
            operator=op,
            version_str=version_str,
            is_wildcard=is_wildcard,
            wildcard_precision=wildcard_precision,
            version=version,
        )

    def contains(self, version_str: str) -> bool:
        """
        Check if a version string satisfies this specifier.

        :param version_str: the version string to test.
        """
        try:
            candidate = SimpleVersion.from_string(version_str) if isinstance(version_str, str) else version_str
        except ValueError:
            return False
        if self.version is None:
            return False
        if self.is_wildcard:
            return self._check_wildcard(candidate)
        return self._check_standard(candidate)

    def _check_wildcard(self, candidate: SimpleVersion) -> bool:
        if self.version is None:  # pragma: no branch
            return False  # pragma: no cover
        if self.operator == "==":
            return candidate.release[: self.wildcard_precision] == self.version.release[: self.wildcard_precision]
        if self.operator == "!=":
            return candidate.release[: self.wildcard_precision] != self.version.release[: self.wildcard_precision]
        return False

    def _check_standard(self, candidate: SimpleVersion) -> bool:
        if self.version is None:  # pragma: no branch
            return False  # pragma: no cover
        if self.operator == "===":
            return str(candidate) == str(self.version)
        if self.operator == "~=":
            return self._check_compatible_release(candidate)
        cmp_ops = {
            "==": operator.eq,
            "!=": operator.ne,
            "<": operator.lt,
            "<=": operator.le,
            ">": operator.gt,
            ">=": operator.ge,
        }
        if self.operator in cmp_ops:
            return cmp_ops[self.operator](candidate, self.version)
        return False

    def _check_compatible_release(self, candidate: SimpleVersion) -> bool:
        if self.version is None:
            return False
        if candidate < self.version:
            return False
        if len(self.version.release) >= 2:  # ruff:ignore[magic-value-comparison]  # pragma: no branch # SimpleVersion always has 3-part release
            upper_parts = list(self.version.release[:-1])
            upper_parts[-1] += 1
            upper = SimpleVersion.from_string(".".join(str(p) for p in upper_parts))
            return candidate < upper
        return True  # pragma: no cover

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SimpleSpecifier):
            return NotImplemented
        return self.spec_str == other.spec_str

    def __hash__(self) -> int:
        return hash(self.spec_str)

    def __str__(self) -> str:
        return self.spec_str

    def __repr__(self) -> str:
        return f"SimpleSpecifier('{self.spec_str}')"


@dataclass(**_DC_KW)
class SimpleSpecifierSet:
    """
    Simple PEP 440-like specifier set using only standard library.

    :param specifiers_str: the original comma-separated specifier string.
    :param specifiers: the parsed individual specifiers.
    """

    specifiers_str: str
    specifiers: tuple[SimpleSpecifier, ...]

    @classmethod
    def from_string(cls, specifiers_str: str = "") -> SimpleSpecifierSet:
        """
        Parse a comma-separated PEP 440 specifier string (e.g. ``>=3.10,<4``).

        :param specifiers_str: the specifier string to parse.
        """
        stripped = specifiers_str.strip()
        specs: list[SimpleSpecifier] = []
        if stripped:
            for spec_item in stripped.split(","):
                item = spec_item.strip()
                if item:
                    with contextlib.suppress(ValueError):
                        specs.append(SimpleSpecifier.from_string(item))
        return cls(specifiers_str=stripped, specifiers=tuple(specs))

    def contains(self, version_str: str) -> bool:
        """
        Check if a version satisfies all specifiers in the set.

        :param version_str: the version string to test.
        """
        if not self.specifiers:
            return True
        return all(spec.contains(version_str) for spec in self.specifiers)

    def __iter__(self) -> Iterator[SimpleSpecifier]:
        return iter(self.specifiers)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SimpleSpecifierSet):
            return NotImplemented
        return self.specifiers_str == other.specifiers_str

    def __hash__(self) -> int:
        return hash(self.specifiers_str)

    def __str__(self) -> str:
        return self.specifiers_str

    def __repr__(self) -> str:
        return f"SimpleSpecifierSet('{self.specifiers_str}')"


__all__ = [
    "SimpleSpecifier",
    "SimpleSpecifierSet",
    "SimpleVersion",
]
