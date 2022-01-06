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

import ast
import pathlib
import sys

DOC_PATH = pathlib.Path("doc/source/getting_started/install.rst").resolve()
CI_PATH = next(
    pathlib.Path("ci/deps").absolute().glob("actions-*-minimum_versions.yaml")
)
CODE_PATH = pathlib.Path("pandas/compat/_optional.py").resolve()


def get_versions_from_code(content: str) -> dict[str, str]:
    num_dicts = 0
    for node in ast.walk(ast.parse(content)):
        if isinstance(node, ast.Dict):
            if num_dicts == 0:
                version_dict_ast = node
                num_dicts += 1
            elif num_dicts == 1:
                install_map = {k.value: v.value for k, v in zip(node.keys, node.values)}
                return {
                    install_map.get(k.value, k.value).casefold(): v.value
                    for k, v in zip(version_dict_ast.keys, version_dict_ast.values)
                    if k.value != "pytest"
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
            if not seen_optional:
                required_deps[package] = version
            else:
                optional_deps[package] = version
    return required_deps, optional_deps


def main():
    with open(CI_PATH) as f:
        _, ci_optional = get_versions_from_ci(f.readlines())
    with open(CODE_PATH) as f:
        code_optional = get_versions_from_code(f.read())
    diff = set(ci_optional.items()).symmetric_difference(code_optional.items())
    if diff:
        sys.stdout.write(
            f"The follow minimum version differences were found between  "
            f"{CI_PATH} and {CODE_PATH}. Please ensure these are aligned: "
            f"{diff}"
        )
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
