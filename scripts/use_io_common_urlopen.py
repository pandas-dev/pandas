"""
Check that pandas/core imports pandas.array as pd_array.

This makes it easier to grep for usage of pandas array.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run use-io-common-urlopen --all-files

"""

from __future__ import annotations

import argparse
import ast
import sys
from typing import Sequence

ERROR_MESSAGE = (
    "{path}:{lineno}:{col_offset}: "
    "Don't use urllib.request.urlopen, use pandas.io.common.urlopen instead\n"
)


class Visitor(ast.NodeVisitor):
    def __init__(self, path: str) -> None:
        self.path = path

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        # Check that pandas.io.common.urlopen is used instead of
        # urllib.request.urlopen
        if (
            node.module is not None
            and node.module.startswith("urllib.request")
            and any(i.name == "urlopen" for i in node.names)
        ):
            msg = ERROR_MESSAGE.format(
                path=self.path, lineno=node.lineno, col_offset=node.col_offset
            )
            sys.stdout.write(msg)
            sys.exit(1)
        super().generic_visit(node)


def use_io_common_urlopen(content: str, path: str) -> None:
    tree = ast.parse(content)
    visitor = Visitor(path)
    visitor.visit(tree)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)

    for path in args.paths:
        with open(path, encoding="utf-8") as fd:
            content = fd.read()
        use_io_common_urlopen(content, path)


if __name__ == "__main__":
    main()
