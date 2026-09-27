"""
Check that the test suite reaches the top-level pandas namespace through ``pd``.

Objects in the top-level pandas namespace are used as ``pd.Series``, ``pd.concat``
and so on, rather than imported directly with ``from pandas import Series``.
Anything else is imported from the module which defines it.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run test-imports --all-files
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import (
        Iterator,
        Sequence,
    )


def _public_names() -> frozenset[str]:
    """Return ``pandas.__all__``, read without importing pandas."""
    path = Path(__file__).parent.parent / "pandas" / "__init__.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:
        if (
            isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "__all__"
                for target in node.targets
            )
            and isinstance(node.value, ast.List)
        ):
            return frozenset(
                element.value
                for element in node.value.elts
                if isinstance(element, ast.Constant)
            )
    raise AssertionError(f"could not find __all__ in {path}")


PUBLIC_NAMES = _public_names()


def _check_import(node: ast.Import) -> Iterator[str]:
    for alias in node.names:
        if alias.name == "pandas" and alias.asname != "pd":
            yield "'pandas' should be imported as 'pd'"


def _check_import_from(node: ast.ImportFrom) -> Iterator[str]:
    if node.level != 0 or node.module != "pandas":
        return
    for alias in node.names:
        if alias.name in PUBLIC_NAMES:
            yield (
                f"use 'pd.{alias.name}' rather than importing '{alias.name}' "
                "from the pandas namespace"
            )
        else:
            yield (
                f"import '{alias.name}' from the module which defines it "
                "rather than from the pandas namespace"
            )


def main(content: str, file: str) -> int:
    tree = ast.parse(content)
    ret = 0
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            messages = _check_import(node)
        elif isinstance(node, ast.ImportFrom):
            messages = _check_import_from(node)
        else:
            continue
        for message in messages:
            print(f"{file}:{node.lineno}:{node.col_offset} {message}")
            ret = 1
    return ret


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()

    ret = 0
    for file in args.paths:
        with open(file, encoding="utf-8") as fd:
            content = fd.read()
        ret |= main(content, file)

    sys.exit(ret)
