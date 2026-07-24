"""A Python specification is an abstract requirement definition of an interpreter."""

from __future__ import annotations

import contextlib
import pathlib
import re
from typing import Final

from ._py_info import normalize_isa
from ._specifier import SimpleSpecifier, SimpleSpecifierSet, SimpleVersion

PATTERN = re.compile(
    r"""
    ^
    (?P<impl>[a-zA-Z]+)?            # implementation (e.g. cpython, pypy)
    (?P<version>[0-9.]+)?           # version (e.g. 3.12, 3.12.1)
    (?P<threaded>t)?                # free-threaded flag
    (?P<debug>d|(?:-dbg|-debug)(?=-|$))?  # debug build flag (d, -dbg or -debug)
    (?:-(?P<arch>32|64))?           # architecture bitness
    (?:-(?P<machine>[a-zA-Z0-9_.]+))?  # ISA (e.g. arm64, x86_64, i86pc.64bit)
    $
    """,
    re.VERBOSE,
)
SPECIFIER_PATTERN = re.compile(
    r"""
    ^
    (?:(?P<impl>[A-Za-z]+)\s*)?     # optional implementation prefix
    (?P<spec>(?:===|==|~=|!=|<=|>=|<|>).+)  # PEP 440 version specifier
    $
    """,
    re.VERBOSE,
)

_MAX_VERSION_PARTS: Final[int] = 3
_SINGLE_DIGIT_MAX: Final[int] = 9

KNOWN_IMPLEMENTATIONS: Final[tuple[str, ...]] = ("python", "cpython", "pypy", "graalpy")

SpecifierSet = SimpleSpecifierSet
Version = SimpleVersion
InvalidSpecifier = ValueError
InvalidVersion = ValueError


def _int_or_none(val: str | None) -> int | None:
    return None if val is None else int(val)


def _parse_version_parts(version: str) -> tuple[int | None, int | None, int | None]:
    versions = tuple(int(i) for i in version.split(".") if i)
    if len(versions) > _MAX_VERSION_PARTS:
        msg = "too many version parts"
        raise ValueError(msg)
    if len(versions) == _MAX_VERSION_PARTS:
        return versions[0], versions[1], versions[2]
    if len(versions) == 2:  # ruff:ignore[magic-value-comparison]
        return versions[0], versions[1], None
    version_data = versions[0]
    major = int(str(version_data)[0])
    minor = int(str(version_data)[1:]) if version_data > _SINGLE_DIGIT_MAX else None
    return major, minor, None


def _parse_spec_pattern(string_spec: str) -> PythonSpec | None:
    if not (match := re.match(PATTERN, string_spec)):
        return None
    groups = match.groupdict()
    version = groups["version"]
    major, minor, micro, threaded = None, None, None, None
    debug = True if groups["debug"] else None  # unconstrained unless an explicit d/-dbg/-debug marker is present
    if version is not None:
        try:
            major, minor, micro = _parse_version_parts(version)
        except ValueError:
            return None
        threaded = bool(groups["threaded"])
    impl = groups["impl"]
    if impl in {"py", "python"}:
        impl = None
    if impl == "graalvm":
        impl = "graalpy"
    arch = _int_or_none(groups["arch"])
    machine = groups.get("machine")
    if machine is not None:
        machine = normalize_isa(machine)
    return PythonSpec(
        string_spec, impl, major, minor, micro, arch, None, free_threaded=threaded, machine=machine, debug=debug
    )


def _parse_specifier(string_spec: str) -> PythonSpec | None:
    if not (specifier_match := SPECIFIER_PATTERN.match(string_spec.strip())):
        return None
    if SpecifierSet is None:  # pragma: no cover
        return None
    impl = specifier_match.group("impl")
    spec_text = specifier_match.group("spec").strip()
    try:
        version_specifier = SpecifierSet.from_string(spec_text)
    except InvalidSpecifier:  # pragma: no cover
        return None
    if impl in {"py", "python"}:
        impl = None
    if impl == "graalvm":
        impl = "graalpy"
    return PythonSpec(string_spec, impl, None, None, None, None, None, version_specifier=version_specifier)


class PythonSpec:
    """
    Contains specification about a Python Interpreter.

    :param str_spec: the raw specification string as provided by the caller.
    :param implementation: interpreter implementation name (e.g. ``"cpython"``, ``"pypy"``), or ``None`` for any.
    :param major: required major version, or ``None`` for any.
    :param minor: required minor version, or ``None`` for any.
    :param micro: required micro (patch) version, or ``None`` for any.
    :param architecture: required pointer-size bitness (``32`` or ``64``), or ``None`` for any.
    :param path: filesystem path to a specific interpreter, or ``None``.
    :param free_threaded: whether a free-threaded build is required, or ``None`` for any.
    :param machine: required ISA (e.g. ``"arm64"``), or ``None`` for any.
    :param debug: whether a debug (``Py_DEBUG``) build is required, or ``None`` for any.
    :param version_specifier:
        `version specifier <https://packaging.python.org/en/latest/specifications/version-specifiers/>`_
        constraints, or ``None``.
    """

    def __init__(  # ruff:ignore[too-many-arguments, too-many-positional-arguments]
        self,
        str_spec: str,
        implementation: str | None,
        major: int | None,
        minor: int | None,
        micro: int | None,
        architecture: int | None,
        path: str | None,
        *,
        free_threaded: bool | None = None,
        machine: str | None = None,
        debug: bool | None = None,
        version_specifier: SpecifierSet | None = None,
    ) -> None:
        self.str_spec = str_spec
        self.implementation = implementation
        self.major = major
        self.minor = minor
        self.micro = micro
        self.free_threaded = free_threaded
        self.architecture = architecture
        self.machine = machine
        self.debug = debug
        self.path = path
        self.version_specifier = version_specifier

    @classmethod
    def from_string_spec(cls, string_spec: str) -> PythonSpec:
        """
        Parse a string specification into a :class:`PythonSpec`.

        :param string_spec: an interpreter spec — an absolute path, a version string, an implementation prefix,
            or a `version specifier <https://packaging.python.org/en/latest/specifications/version-specifiers/>`_.
        """
        if pathlib.Path(string_spec).is_absolute():
            return cls(string_spec, None, None, None, None, None, string_spec)
        if result := _parse_spec_pattern(string_spec):
            return result
        if result := _parse_specifier(string_spec):
            return result
        return cls(string_spec, None, None, None, None, None, string_spec)

    def generate_re(self, *, windows: bool, all_implementations: bool = False) -> re.Pattern:
        """
        Generate a regular expression for matching interpreter filenames.

        :param windows: if ``True``, require a ``.exe`` suffix.
        :param all_implementations: when ``True`` and the spec does not constrain the implementation, match every
            filename in :data:`KNOWN_IMPLEMENTATIONS` instead of only ``python``. Used by enumeration APIs.
        """
        version = r"{}(\.{}(\.{})?)?".format(
            *(r"\d+" if v is None else v for v in (self.major, self.minor, self.micro)),
        )
        if self.implementation is not None:
            impl = f"python|{re.escape(self.implementation)}"
        elif all_implementations:
            impl = "|".join(re.escape(i) for i in KNOWN_IMPLEMENTATIONS)
        else:
            impl = "python"
        mod = "t?" if self.free_threaded else ""
        dbg = "(?:d|-dbg|-debug)?" if self.debug else ""
        suffix = r"\.exe" if windows else ""
        version_conditional = "?" if windows or self.major is None else ""
        return re.compile(
            rf"(?P<impl>{impl})(?P<v>{version}{mod}){version_conditional}{dbg}{suffix}$",
            flags=re.IGNORECASE,
        )

    @property
    def is_abs(self) -> bool:
        """``True`` if the spec refers to an absolute filesystem path."""
        return self.path is not None and pathlib.Path(self.path).is_absolute()

    def _check_version_specifier(self, spec: PythonSpec) -> bool:
        """Check if version specifier is satisfied."""
        components: list[int] = []
        for part in (self.major, self.minor, self.micro):
            if part is None:
                break
            components.append(part)
        if not components:
            return True

        version_str = ".".join(str(part) for part in components)
        if spec.version_specifier is None:
            return True
        with contextlib.suppress(InvalidVersion):
            Version.from_string(version_str)
            for item in spec.version_specifier:
                required_precision = self._get_required_precision(item)
                if required_precision is None or len(components) < required_precision:
                    continue
                if not item.contains(version_str):
                    return False
        return True

    @staticmethod
    def _get_required_precision(item: SimpleSpecifier) -> int | None:
        """Get the required precision for a specifier item."""
        if item.version is None:
            return None
        with contextlib.suppress(AttributeError, ValueError):
            return len(item.version.release)
        return None

    def satisfies(self, spec: PythonSpec) -> bool:  # ruff:ignore[too-many-return-statements]
        """
        Check if this spec is compatible with the given *spec* (e.g. PEP-514 on Windows).

        :param spec: the requirement to check against.
        """
        if spec.is_abs and self.is_abs and self.path != spec.path:
            return False
        if (
            spec.implementation is not None
            and self.implementation is not None
            and spec.implementation.lower() != self.implementation.lower()
        ):
            return False
        if spec.architecture is not None and spec.architecture != self.architecture:
            return False
        if spec.machine is not None and self.machine is not None and spec.machine != self.machine:
            return False
        if spec.free_threaded is not None and spec.free_threaded != self.free_threaded:
            return False
        if spec.version_specifier is not None and not self._check_version_specifier(spec):
            return False
        return all(
            req is None or our is None or our == req
            for our, req in zip((self.major, self.minor, self.micro), (spec.major, spec.minor, spec.micro))
        )

    def __repr__(self) -> str:
        name = type(self).__name__
        params = (
            "implementation",
            "major",
            "minor",
            "micro",
            "architecture",
            "machine",
            "path",
            "free_threaded",
            "debug",
            "version_specifier",
        )
        return f"{name}({', '.join(f'{k}={getattr(self, k)}' for k in params if getattr(self, k) is not None)})"


__all__ = [
    "KNOWN_IMPLEMENTATIONS",
    "InvalidSpecifier",
    "InvalidVersion",
    "PythonSpec",
    "SpecifierSet",
    "Version",
]
