#!/usr/bin/env python3
"""
Check pandas required and optional dependencies are synced across:

ci/deps/actions-.*-minimum_versions.yaml
pandas/compat/_optional.py

TODO: doc/source/getting_started/install.rst

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run validate-min-versions-in-sync --all-files
"""
from __future__ import annotations

import pathlib
import sys

DOC_PATH = pathlib.Path("doc/source/getting_started/install.rst").resolve()
CI_PATH = next(
    pathlib.Path("ci/deps").absolute().glob("actions-*-minimum_versions.yaml")
)
CODE_PATH = pathlib.Path("pandas/compat/_optional.py").resolve()
EXCLUDE_DEPS = {"tzdata"}
# pandas package is not available
# in pre-commit environment
sys.path.append("pandas/compat")
sys.path.append("pandas/util")
import _exceptions
import version

sys.modules["pandas.util.version"] = version
sys.modules["pandas.util._exceptions"] = _exceptions
import _optional


def get_versions_from_code() -> dict[str, str]:
    install_map = _optional.INSTALL_MAPPING
    versions = _optional.VERSIONS
    for item in EXCLUDE_DEPS:
        versions.pop(item)
    return {
        install_map.get(k, k).casefold(): v
        for k, v in versions.items()
        if k != "pytest"
    }


def get_versions_from_ci(content: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    # Don't parse with pyyaml because it ignores comments we're looking for
    seen_required = False
    seen_optional = False
    required_deps = {}
    optional_deps = {}
    for line in content:
        if "# required dependencies" in line:
            seen_required = True
        elif "# optional dependencies" in line:
            seen_optional = True
        elif seen_required and line.strip():
            package, version = line.strip().split("=")
            package = package[2:]
            if package in EXCLUDE_DEPS:
                continue
            if not seen_optional:
                required_deps[package] = version
            else:
                optional_deps[package] = version
    return required_deps, optional_deps


def main():
    with open(CI_PATH, encoding="utf-8") as f:
        _, ci_optional = get_versions_from_ci(f.readlines())
    code_optional = get_versions_from_code()
    diff = set(ci_optional.items()).symmetric_difference(code_optional.items())
    if diff:
        sys.stdout.write(
            f"The follow minimum version differences were found between  "
            f"{CI_PATH} and {CODE_PATH}. Please ensure these are aligned: "
            f"{diff}\n"
        )
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
