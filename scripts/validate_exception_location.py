"""
Validate that the exceptions and warnings are in appropriate places.

Checks for classes that inherit a python exception and warning and
flags them, unless they are exempted from checking. Exempt meaning
the exception/warning is defined in testing.rst. Testing.rst contains
a list of pandas defined exceptions and warnings. This list is kept
current by other pre-commit hook, pandas_errors_documented.py.
This hook maintains that errors.__init__.py and testing.rst are in-sync.
Therefore, the exception or warning should be defined or imported in
errors.__init__.py. Ideally, the exception or warning is defined unless
there's special reason to import it.

Prints the exception/warning that do not follow this convention.

Usage::

As a pre-commit hook:
    pre-commit run validate-errors-locations --all-files
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
ERROR_MESSAGE = (
    "The following exception(s) and/or warning(s): {errors} exist(s) outside of "
    "pandas/errors/__init__.py. Please either define them in "
    "pandas/errors/__init__.py. Or, if not possible then import them in "
    "pandas/errors/__init__.py.\n"
)


def get_warnings_and_exceptions_from_api_path() -> set[str]:
    with open(API_PATH, encoding="utf-8") as f:
        doc_errors = {
            line.split(".")[1].strip() for line in f.readlines() if "errors" in line
        }
        return doc_errors


class Visitor(ast.NodeVisitor):
    def __init__(self, path: str, exception_set: set[str]) -> None:
        self.path = path
        self.exception_set = exception_set
        self.found_exceptions = set()

    def visit_ClassDef(self, node) -> None:
        def is_an_exception_subclass(base_id: str) -> bool:
            return base_id == "Exception" or base_id.endswith(("Warning", "Error"))

        exception_classes = []

        # Go through the class's bases and check if they are an Exception or Warning.
        for base in node.bases:
            base_id = getattr(base, "id", None)
            if base_id and is_an_exception_subclass(base_id):
                exception_classes.append(base_id)

        # The class subclassed an Exception or Warning so add it to the list.
        if exception_classes:
            self.found_exceptions.add(node.name)


def validate_exception_and_warning_placement(
    file_path: str, file_content: str, errors: set[str]
):
    tree = ast.parse(file_content)
    visitor = Visitor(file_path, errors)
    visitor.visit(tree)

    misplaced_exceptions = visitor.found_exceptions.difference(errors)

    # If misplaced_exceptions isn't an empty list then there exists
    # pandas-defined Exception or Warnings outside of pandas/errors/__init__.py, so
    # we should flag them.
    if misplaced_exceptions:
        msg = ERROR_MESSAGE.format(errors=", ".join(misplaced_exceptions))
        sys.stdout.write(msg)
        sys.exit(1)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)

    error_set = get_warnings_and_exceptions_from_api_path()

    for path in args.paths:
        with open(path, encoding="utf-8") as fd:
            content = fd.read()
        validate_exception_and_warning_placement(path, content, error_set)


if __name__ == "__main__":
    main()
