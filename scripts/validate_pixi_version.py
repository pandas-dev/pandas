#!/usr/bin/env python3
"""
Check the pinned pixi version agrees across the CI config files that pin it.

The version is pinned in three places and nothing links them, so a mismatch
silently produces CI jobs installing under different pixi versions (and can
commit a pixi.lock produced by a different pixi than the one consuming it):

    .github/actions/setup-pixi/action.yml
    .github/workflows/update-pixi-lock.yml

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run validate-pixi-version --all-files
"""

from __future__ import annotations

import pathlib
import re
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

BASE_PATH = pathlib.Path(__file__).parents[1]
ACTION_PATH = BASE_PATH / ".github/actions/setup-pixi/action.yml"
WORKFLOW_PATH = BASE_PATH / ".github/workflows/update-pixi-lock.yml"

PIXI_VERSION_RE = re.compile(r"^\s*pixi-version:\s*(\S+)\s*$")


def check_pixi_versions(paths: Iterable[pathlib.Path]) -> int:
    """Return 0 when every ``pixi-version`` pin agrees, 1 otherwise."""
    locations: dict[str, list[str]] = {}
    for path in paths:
        try:
            f = open(path, encoding="utf-8")
        except OSError as err:
            print(f"Could not read {path}: {err}")
            return 1
        with f:
            for lineno, line in enumerate(f, 1):
                match = PIXI_VERSION_RE.match(line)
                if match:
                    locations.setdefault(match.group(1), []).append(f"{path}:{lineno}")

    if not locations:
        print(
            "No pixi-version pin found. If pixi is no longer used, remove this "
            "check from .pre-commit-config.yaml."
        )
        return 1

    if len(locations) == 1:
        return 0

    print("pixi-version is pinned inconsistently across CI config files:")
    for version, pins in sorted(locations.items()):
        print(f"  {version} found at:")
        for pin in pins:
            print(f"    {pin}")
    print("Please align every pin to the same pixi version.")
    return 1


def main() -> int:
    return check_pixi_versions([ACTION_PATH, WORKFLOW_PATH])


if __name__ == "__main__":
    sys.exit(main())
