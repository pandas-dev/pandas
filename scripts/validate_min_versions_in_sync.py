#!/usr/bin/env python3
"""
Check pandas required and optional dependencies are synced across:

doc/source/getting_started/install.rst
ci/deps/actions-.*-minimum_versions.yaml
pandas/compat/_optional.py

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
                break
    return {
        install_map.get(k.value, k.value): v.value
        for k, v in zip(version_dict_ast.keys, version_dict_ast.values)
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


def get_versions_from_doc(content: str, ci_verions: dict[str, str]) -> dict[str, str]:
    pass


def main():
    # The CI file is our source of truth since it's what we're testing.
    with open(CI_PATH) as f:
        ci_required, ci_optional = get_versions_from_ci(f.readlines())
    with open(CODE_PATH) as f:
        code_versions = get_versions_from_code(f.read())
    with open(DOC_PATH) as f:
        doc_verions = get_versions_from_doc(f.read(), ci_versions)


if __name__ == "__main__":
    main()
