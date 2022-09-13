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
from typing import Sequence

API_PATH = pathlib.Path("doc/source/reference/testing.rst").resolve()
ERROR_MESSAGE = (
    "The following exception(s) and/or warning(s): {errors} exist(s) outside of "
    "pandas/errors/__init__.py. Please either define them in "
    "pandas/errors/__init__.py. Or, if not possible then import them in "
    "pandas/errors/__init__.py.\n"
)
exception_warning_list = {
    "ArithmeticError",
    "AssertionError",
    "AttributeError",
    "EOFError",
    "Exception",
    "FloatingPointError",
    "GeneratorExit",
    "ImportError",
    "IndentationError",
    "IndexError",
    "KeyboardInterrupt",
    "KeyError",
    "LookupError",
    "MemoryError",
    "NameError",
    "NotImplementedError",
    "OSError",
    "OverflowError",
    "ReferenceError",
    "RuntimeError",
    "StopIteration",
    "SyntaxError",
    "SystemError",
    "SystemExit",
    "TabError",
    "TypeError",
    "UnboundLocalError",
    "UnicodeDecodeError",
    "UnicodeEncodeError",
    "UnicodeError",
    "UnicodeTranslateError",
    "ValueError",
    "ZeroDivisionError",
    "BytesWarning",
    "DeprecationWarning",
    "FutureWarning",
    "ImportWarning",
    "PendingDeprecationWarning",
    "ResourceWarning",
    "RuntimeWarning",
    "SyntaxWarning",
    "UnicodeWarning",
    "UserWarning",
    "Warning",
}


def get_warnings_and_exceptions_from_api_path() -> set[str]:
    with open(API_PATH) as f:
        doc_errors = {
            line.split(".")[1].strip() for line in f.readlines() if "errors" in line
        }
        return doc_errors


class Visitor(ast.NodeVisitor):
    def __init__(self, path: str, exception_set: set[str]) -> None:
        self.path = path
        self.exception_set = exception_set
        self.possible_exceptions = set()

    def visit_ClassDef(self, node):
        classes = {getattr(n, "id", None) for n in node.bases}

        if classes and classes.issubset(exception_warning_list):
            self.possible_exceptions.add(node.name)


def validate_exception_and_warning_placement(
    file_path: str, file_content: str, errors: set[str]
):
    tree = ast.parse(file_content)
    visitor = Visitor(file_path, errors)
    visitor.visit(tree)

    misplaced_exceptions = visitor.possible_exceptions.difference(errors)

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
