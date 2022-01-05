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

# import argparse
import ast
import pathlib

# import re
# import sys

DOC_PATH = pathlib.Path("doc/source/getting_started/install.rst").resolve()
CI_PATH = next(
    pathlib.Path("ci/deps").absolute().glob("actions-*-minimum_versions.yaml")
)
CODE_PATH = pathlib.Path("pandas/compat/_optional.py").resolve()


def get_versions_from_optional(content: str) -> dict[str, str]:
    num_dicts = 0
    for node in ast.walk(ast.parse(content)):
        if isinstance(node, ast.Dict):
            if num_dicts == 0:
                version_dict_ast = node
                num_dicts += 1
            elif num_dicts == 1:
                install_map_ast = node
                break
    install_map = {
        k.value: v.value for k, v in zip(install_map_ast.keys, install_map_ast.values)
    }
    return {
        install_map.get(k.value, k.value): v.value
        for k, v in zip(version_dict_ast.keys, version_dict_ast.values)
    }


if __name__ == "__main__":
    pass
