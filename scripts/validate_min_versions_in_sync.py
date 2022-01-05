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

import yaml

# import re
# import sys

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


def get_versions_from_ci(fle) -> dict[str, str]:
    yml_content = yaml.safe_load(fle)
    yml_version = {}
    for dependency in reversed(yml_content["dependencies"]):
        if "=" not in dependency:
            break
        package, version = dependency.split("=")
        yml_version[package] = version
    return yml_version


if __name__ == "__main__":
    with open(CODE_PATH) as f:
        code_versions = get_versions_from_code(f.read())
    with open(CI_PATH) as f:
        yml_content = get_versions_from_ci(f)
    pass
