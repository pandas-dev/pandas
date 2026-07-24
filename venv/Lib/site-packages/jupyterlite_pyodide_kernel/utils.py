"""Utilities used by multiple addons and tools."""

from __future__ import annotations

import json
import re
from fnmatch import fnmatch
from functools import lru_cache
from urllib.parse import urlparse
from pathlib import Path
from typing import Any
from typing import TYPE_CHECKING

from jupyterlite_core.constants import JSON_FMT, UTF8

from .constants import (
    ALL_WHL,
    RE_WHEEL_DIST_NAME,
    PEP_735_DEP_GROUPS,
    PEP_735_INC_GROUP,
    PYPROJECT_TOML,
)


if TYPE_CHECKING:
    from packaging.utils import NormalizedName
    from pkginfo import Distribution
    from collections.abc import Iterator


@lru_cache(100)
def normalize_names(*names: str) -> list[NormalizedName]:
    """Return a normalized set of Python package names."""
    from packaging.utils import canonicalize_name

    return sorted({*map(canonicalize_name, names)})


@lru_cache(1000)
def get_wheel_metadata(filename: str) -> Distribution | None:
    """Try to get cached metadata for a Python distribution."""
    import pkginfo

    return pkginfo.get_metadata(filename)


def get_wheel_name(wheel: Path) -> NormalizedName | None:
    """Get the normalized package name contained in a wheel"""
    from packaging.utils import canonicalize_name

    info = get_wheel_metadata(f"{wheel}")
    if not (info and info.name):
        return None
    return canonicalize_name(info.name)


def wheel_to_pep508(path_or_url: str) -> str | None:
    """Get a PEP-508 direct URL from a path or URL."""
    from packaging.utils import canonicalize_name

    url = urlparse(path_or_url)

    dist_name_match = re.search(RE_WHEEL_DIST_NAME, path_or_url)

    if not dist_name_match:
        return None

    dist_name = canonicalize_name(dist_name_match.groupdict()["name"])
    final_url = path_or_url if url.scheme else Path(path_or_url).resolve().as_uri()

    return f"{dist_name} @ {final_url}"


def is_pyodide_wheel(filename: str, patterns: list[str] | None = None) -> bool:
    """get whether a wheel is a known-good pyodide wheel."""
    return any(fnmatch(filename, f"**/*{p}") for p in patterns or ALL_WHL)


def list_wheels(
    *wheel_dirs: Path,
    patterns: list[str] | None = None,
    recursive: bool = False,
) -> list[Path]:
    """get all wheels we know how to handle in a directory"""
    wheels = []
    for wheel_dir in wheel_dirs:
        if not wheel_dir.is_dir():
            continue
        wheels += sorted((wheel_dir.rglob if recursive else wheel_dir.glob)("*.whl"))
    return [w for w in wheels if is_pyodide_wheel(w, patterns)]


def patch_dict(old: dict[str, Any], new: dict[str, Any]) -> dict[str, Any]:
    """Recursively update a dict in-place with new values."""
    for key, value in new.items():
        if value is None:
            old.pop(key, None)
        elif isinstance(value, dict):
            old[key] = patch_dict(old.get(key, {}), value)
        else:
            old[key] = value
    return old


def patch_json_path(old_path: Path, patch: dict[str, Any]) -> None:
    """Update an on-disk JSON file with a patch."""
    old = patch_dict(json.loads(old_path.read_text(**UTF8)), patch)
    old_path.write_text(json.dumps(old, **JSON_FMT) + "\n", **UTF8)


def iter_pep508_specs(
    specs: list[str], base_path: Path, seen: set[str] | None = None
) -> Iterator[str]:
    """Parse a set of PEP-508 specs, with relative file references."""
    seen = seen if seen is not None else set()
    for line in specs:
        yield from _iter_one_pep508_spec(line, base_path, seen)


def _iter_one_pep508_spec(line: str, base_path: Path, seen: set[str]) -> Iterator[str]:
    """Parse a single ``requirements.txt`` line."""
    line = line.split("#")[0].strip()
    if not line:
        return
    for pattern, _handler in _PEP_508_FILE_HANDLERS.items():
        match = re.match(pattern, line)
        if match is None:
            continue
        m = match.groupdict()
        m["path"] = (base_path / (m.get("path") or "")).resolve()
        yield from _handler(**m, seen=seen)  # type: ignore[operator]
        return
    yield line


def _iter_one_pep508_reqs_path(path: Path, seen: set[str]) -> Iterator[str]:
    """Parse a ``requirements.txt``-style file as PEP-508 specs."""
    if not path.is_file():  # pragma: no cover
        msg = f"The requirements/constraints file could not be found: {path}"
        raise FileNotFoundError(msg)

    uri = path.as_uri()
    if uri in seen:
        return

    seen |= {uri}

    lines = path.read_text(encoding="utf-8").splitlines()
    yield from iter_pep508_specs(lines, path.parent, seen)


def _iter_one_pep735_group_path(
    path: Path, group: str, seen: set[str], ppt: dict[str, Any] | None = None
) -> Iterator[str]:
    """Parse a named ``dependency-group`` from ``pyproject.toml`` as PEP-508 specs."""
    path = path / PYPROJECT_TOML if path.is_dir() else path
    uri = f"{path.as_uri()}#{group}"
    if uri in seen:
        return

    seen |= {uri}

    if ppt is None:
        try:  # pragma: no cover
            import tomllib
        except ImportError:  # pragma: no cover
            import tomli as tomllib

        ppt = tomllib.loads(path.read_text(**UTF8))

    for spec in ppt[PEP_735_DEP_GROUPS][group]:
        if isinstance(spec, str):
            yield from _iter_one_pep508_spec(spec, path.parent, seen)
        elif isinstance(spec, dict) and PEP_735_INC_GROUP in spec:
            yield from _iter_one_pep735_group_path(
                path, spec[PEP_735_INC_GROUP], seen, ppt
            )
        else:  # pragma: no cover
            msg = f"Can't parse requested spec from {group} in {path}: {spec}"
            raise ValueError(msg)


_PEP_508_FILE_HANDLERS = {
    #: a reference to a ``requirements.txt``-style file
    r"^(-r|--requirements)\s*=?\s*(?P<path>.+)$": _iter_one_pep508_reqs_path,
    #: a reference to a ``dependency-groups`` in a ``pyproject.toml``
    r"^(-g|--group)\s*=?\s*((?P<path>.+):)?(?P<group>[^:/\\]+)$": _iter_one_pep735_group_path,
}
