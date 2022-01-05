#!/usr/bin/env python3
"""
Check pandas required and optional dependencies are synced across:

doc/source/getting_started/install.rst
ci/deps/actions-.*-minimum_versions.yaml
pandas/compat/_optional.py

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run validate-min-versions-in-sync --all-files
"""
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
    for node in ast.walk(ast.parse(content)):
        if isinstance(node, ast.Dict):
            # version_dict = node
            break


if __name__ == "__main__":
    pass
