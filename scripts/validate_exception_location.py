from __future__ import annotations

import argparse
import ast
import sys
from typing import Sequence

ERROR_MESSAGE = (
    "{path}:{lineno}:{col_offset}: {exception_name}: "
    "Please don't place exceptions or warnings outside of pandas/errors/__init__.py or "
    "pandas/_libs\n"
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

permisable_exception_warning_list = [
    "LossySetitemError",
    "NoBufferPresent",
    "InvalidComparison",
    "NotThisMethod",
    "OptionError",
]


class Visitor(ast.NodeVisitor):
    def __init__(self, path: str) -> None:
        self.path = path

    def visit_ClassDef(self, node):
        classes = {getattr(n, "id", None) for n in node.bases}

        if (
            classes
            and classes.issubset(exception_warning_list)
            and node.name not in permisable_exception_warning_list
        ):
            msg = ERROR_MESSAGE.format(
                path=self.path,
                lineno=node.lineno,
                col_offset=node.col_offset,
                exception_name=node.name,
            )
            sys.stdout.write(msg)
            sys.exit(1)


def validate_exception_and_warning_placement(file_path: str, file_content: str):
    tree = ast.parse(file_content)
    visitor = Visitor(file_path)
    visitor.visit(tree)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)

    for path in args.paths:
        with open(path, encoding="utf-8") as fd:
            content = fd.read()
        validate_exception_and_warning_placement(path, content)


if __name__ == "__main__":
    main()
