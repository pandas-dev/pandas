"""
Check that doc/source/reference/testing.rst documents
all exceptions and warnings in pandas/errors/__init__.py.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run pandas-errors-documented --all-files
"""
from __future__ import annotations

import argparse
import ast
import pathlib
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

API_PATH = pathlib.Path("doc/source/reference/testing.rst").resolve()


def get_defined_errors(content: str) -> set[str]:
    errors = set()
    for node in ast.walk(ast.parse(content)):
        if isinstance(node, ast.ClassDef):
            errors.add(node.name)
        elif isinstance(node, ast.ImportFrom) and node.module != "__future__":
            for alias in node.names:
                errors.add(alias.name)
    return errors


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    args = parser.parse_args(argv)
    with open(args.path, encoding="utf-8") as f:
        file_errors = get_defined_errors(f.read())
    with open(API_PATH, encoding="utf-8") as f:
        doc_errors = {
            line.split(".")[1].strip() for line in f.readlines() if "errors" in line
        }
    missing = file_errors.difference(doc_errors)
    if missing:
        sys.stdout.write(
            f"The following exceptions and/or warnings are not documented "
            f"in {API_PATH}: {missing}"
        )
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
