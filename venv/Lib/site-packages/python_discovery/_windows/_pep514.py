"""Implement https://www.python.org/dev/peps/pep-0514/ to discover interpreters - Windows only."""

from __future__ import annotations

import logging
import os
import re
import sys
import winreg
from logging import basicConfig, getLogger
from typing import TYPE_CHECKING, Any, Final

if TYPE_CHECKING:
    from collections.abc import Generator

    _RegistrySpec = tuple[str, int | None, int | None, int, bool, str, str | None]

_LOGGER: Final[logging.Logger] = getLogger(__name__)
_ARCH_RE: Final[re.Pattern[str]] = re.compile(
    r"""
    ^
    (\d+)   # bitness number
    bit     # literal suffix
    $
    """,
    re.VERBOSE,
)
_VERSION_RE: Final[re.Pattern[str]] = re.compile(
    r"""
    ^
    (\d+)            # major
    (?:\.(\d+))?     # optional minor
    (?:\.(\d+))?     # optional micro
    $
    """,
    re.VERBOSE,
)
_THREADED_TAG_RE: Final[re.Pattern[str]] = re.compile(
    r"""
    ^
    \d+              # major
    (\.\d+){0,2}     # optional minor/micro
    t                # free-threaded flag
    $
    """,
    re.VERBOSE | re.IGNORECASE,
)


def enum_keys(key: Any) -> Generator[str, None, None]:  # ruff:ignore[any-type]
    at = 0
    while True:
        try:
            yield winreg.EnumKey(key, at)  # ty: ignore[unresolved-attribute]
        except OSError:
            break
        at += 1


def get_value(key: Any, value_name: str | None) -> Any:  # ruff:ignore[any-type]
    try:
        return winreg.QueryValueEx(key, value_name)[0]  # ty: ignore[unresolved-attribute]
    except OSError:
        return None


def discover_pythons() -> Generator[_RegistrySpec, None, None]:
    for hive, hive_name, key, flags, default_arch in [
        (winreg.HKEY_CURRENT_USER, "HKEY_CURRENT_USER", r"Software\Python", 0, 64),  # ty: ignore[unresolved-attribute]
        (winreg.HKEY_LOCAL_MACHINE, "HKEY_LOCAL_MACHINE", r"Software\Python", winreg.KEY_WOW64_64KEY, 64),  # ty: ignore[unresolved-attribute]
        (winreg.HKEY_LOCAL_MACHINE, "HKEY_LOCAL_MACHINE", r"Software\Python", winreg.KEY_WOW64_32KEY, 32),  # ty: ignore[unresolved-attribute]
    ]:
        yield from process_set(hive, hive_name, key, flags, default_arch)


def process_set(
    hive: int,
    hive_name: str,
    key: str,
    flags: int,
    default_arch: int,
) -> Generator[_RegistrySpec, None, None]:
    try:
        with winreg.OpenKeyEx(hive, key, 0, winreg.KEY_READ | flags) as root_key:  # ty: ignore[unresolved-attribute]
            for company in enum_keys(root_key):
                if company == "PyLauncher":  # reserved
                    continue
                yield from process_company(hive_name, company, root_key, default_arch)
    except OSError:
        pass


def process_company(
    hive_name: str,
    company: str,
    root_key: Any,  # ruff:ignore[any-type]
    default_arch: int,
) -> Generator[_RegistrySpec, None, None]:
    with winreg.OpenKeyEx(root_key, company) as company_key:  # ty: ignore[unresolved-attribute]
        for tag in enum_keys(company_key):
            spec = process_tag(hive_name, company, company_key, tag, default_arch)
            if spec is not None:
                yield spec


def process_tag(hive_name: str, company: str, company_key: Any, tag: str, default_arch: int) -> _RegistrySpec | None:  # ruff:ignore[any-type]
    with winreg.OpenKeyEx(company_key, tag) as tag_key:  # ty: ignore[unresolved-attribute]
        version = load_version_data(hive_name, company, tag, tag_key)
        if version is not None:  # if failed to get version bail
            major, minor, _ = version
            arch = load_arch_data(hive_name, company, tag, tag_key, default_arch)
            if arch is not None:
                exe_data = load_exe(hive_name, company, company_key, tag)
                if exe_data is not None:
                    exe, args = exe_data
                    threaded = load_threaded(hive_name, company, tag, tag_key)
                    return company, major, minor, arch, threaded, exe, args
                return None
            return None
        return None


def load_exe(hive_name: str, company: str, company_key: Any, tag: str) -> tuple[str, str | None] | None:  # ruff:ignore[any-type]
    key_path = f"{hive_name}/{company}/{tag}"
    try:
        with winreg.OpenKeyEx(company_key, rf"{tag}\InstallPath") as ip_key, ip_key:  # ty: ignore[unresolved-attribute]
            if (exe := _resolve_exe(ip_key, key_path)) is not None and os.path.exists(exe):
                return exe, get_value(ip_key, "ExecutableArguments")
            msg(key_path, f"could not load exe with value {exe}")
    except OSError:
        msg(f"{key_path}/InstallPath", "missing")
    return None


def _resolve_exe(ip_key: Any, key_path: str) -> str | None:  # ruff:ignore[any-type]
    if (exe := get_value(ip_key, "ExecutablePath")) is not None:
        return exe
    if (ip := get_value(ip_key, None)) is None:
        msg(key_path, "no ExecutablePath or default for it")
        return None
    return os.path.join(ip, "python.exe")


def load_arch_data(hive_name: str, company: str, tag: str, tag_key: Any, default_arch: int) -> int | None:  # ruff:ignore[any-type]
    arch_str = get_value(tag_key, "SysArchitecture")
    if arch_str is not None:
        key_path = f"{hive_name}/{company}/{tag}/SysArchitecture"
        try:
            return parse_arch(arch_str)
        except ValueError as sys_arch:
            msg(key_path, sys_arch)
    return default_arch


def parse_arch(arch_str: Any) -> int:  # ruff:ignore[any-type]
    if isinstance(arch_str, str):
        if match := _ARCH_RE.match(arch_str):
            return int(next(iter(match.groups())))
        error = f"invalid format {arch_str}"
    else:
        error = f"arch is not string: {arch_str!r}"
    raise ValueError(error)


def load_version_data(
    hive_name: str,
    company: str,
    tag: str,
    tag_key: Any,  # ruff:ignore[any-type]
) -> tuple[int | None, int | None, int | None] | None:
    for candidate, key_path in [
        (get_value(tag_key, "SysVersion"), f"{hive_name}/{company}/{tag}/SysVersion"),
        (tag, f"{hive_name}/{company}/{tag}"),
    ]:
        if candidate is not None:
            try:
                return parse_version(candidate)
            except ValueError as sys_version:
                msg(key_path, sys_version)
    return None


def parse_version(version_str: Any) -> tuple[int | None, int | None, int | None]:  # ruff:ignore[any-type]
    if isinstance(version_str, str):
        if match := _VERSION_RE.match(version_str):
            g1, g2, g3 = match.groups()
            return (
                int(g1) if g1 is not None else None,
                int(g2) if g2 is not None else None,
                int(g3) if g3 is not None else None,
            )
        error = f"invalid format {version_str}"
    else:
        error = f"version is not string: {version_str!r}"
    raise ValueError(error)


def load_threaded(hive_name: str, company: str, tag: str, tag_key: Any) -> bool:  # ruff:ignore[any-type]
    display_name = get_value(tag_key, "DisplayName")
    if display_name is not None:
        if isinstance(display_name, str):
            if "freethreaded" in display_name.lower():
                return True
        else:
            key_path = f"{hive_name}/{company}/{tag}/DisplayName"
            msg(key_path, f"display name is not string: {display_name!r}")
    return bool(_THREADED_TAG_RE.match(tag))


def msg(path: str, what: object) -> None:
    _LOGGER.warning("PEP-514 violation in Windows Registry at %s error: %s", path, what)


def _run() -> None:
    basicConfig()
    interpreters = [repr(spec) for spec in discover_pythons()]
    sys.stdout.write("\n".join(sorted(interpreters)))
    sys.stdout.write("\n")


if __name__ == "__main__":
    _run()
