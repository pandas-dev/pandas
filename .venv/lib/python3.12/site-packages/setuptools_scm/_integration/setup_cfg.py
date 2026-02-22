from __future__ import annotations

import os

from dataclasses import dataclass
from pathlib import Path

import setuptools


@dataclass
class SetuptoolsBasicData:
    path: Path
    name: str | None
    version: str | None


def read_setup_cfg(input: str | os.PathLike[str] = "setup.cfg") -> SetuptoolsBasicData:
    """Parse setup.cfg and return unified data. Does not raise if file is missing."""
    import configparser

    path = Path(input)
    parser = configparser.ConfigParser()
    parser.read([input], encoding="utf-8")

    name = parser.get("metadata", "name", fallback=None)
    version = parser.get("metadata", "version", fallback=None)
    if version is not None and "attr" in version:
        from .deprecation import warn_setup_cfg_dynamic_version

        warn_setup_cfg_dynamic_version(path)
        version = None
    return SetuptoolsBasicData(path=path, name=name, version=version)


def extract_from_legacy(
    dist: setuptools.Distribution,
    *,
    _given_legacy_data: SetuptoolsBasicData | None = None,
) -> SetuptoolsBasicData:
    base = _given_legacy_data if _given_legacy_data is not None else read_setup_cfg()
    if base.name is None:
        base.name = dist.metadata.name
    if base.version is None:
        base.version = dist.metadata.version
    return base
