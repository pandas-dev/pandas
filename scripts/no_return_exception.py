"""
Check that functions return isn't use for Exception or its subclasses.

This ensures public methods don't accidentally return exceptions instead of raising.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run no-return-exception --all-files

"""
from __future__ import annotations

import argparse
import ast
import sys
from typing import Sequence


ERROR_MESSAGE = (
    "{path}:{lineno}:{col_offset}: "
    "Don't return an Exception subclass, raise it instead\n"
)


class Visitor(ast.NodeVisitor):
    def __init__(self, path: str) -> None:
        self.path = path

    def visit_FunctionDef(self, node) -> None:
        returns = node.returns
        name = getattr(returns, "id", None)
        exception_ends = ("Exit", "Interrupt", "Exception", "Error", "Iteration")
        if name is not None and any(name.endswith(end) for end in exception_ends):
            msg = ERROR_MESSAGE.format(
                path=self.path, lineno=node.lineno, col_offset=node.col_offset
            )
            sys.stdout.write(msg)
            sys.exit(1)


def no_return_exception(content: str, path: str) -> None:
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
        no_return_exception(content, path)


if __name__ == "__main__":
    main()
