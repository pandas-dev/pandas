"""
Ensure that test fixtures have a `name` argument.

This is to ensure that pylint's redefined-outer-name does not
report tests which use fixtures, see PyCQA/pylint#1535.

You can run this in pre-commit::

    pre-commit run rename-test-fixtures --all-files

or by passing filenames manually::

    python -m scripts.rename_test_fixtures <file_1> <file_2> ... <file_n>
"""
from __future__ import annotations

import argparse
import ast
from ast import NodeVisitor
import re
import sys
from typing import NamedTuple


class Offset(NamedTuple):
    lineno: int
    col_offset: int
    name: str
    new_name: str | None


class FixtureVisitor(NodeVisitor):
    def __init__(self, file) -> None:
        self.file = file
        self.fixture_attributes: list[Offset] = []
        self.fixture_calls: list[Offset] = []
        self.func_defs: list[Offset] = []
        self.misnamed_test: list[Offset] = []

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        for decorator in node.decorator_list:
            if (
                isinstance(decorator, ast.Attribute)
                and isinstance(decorator.value, ast.Name)
                and decorator.value.id == "pytest"
                and decorator.attr == "fixture"
            ):
                self.fixture_attributes.append(
                    Offset(decorator.lineno, decorator.col_offset, node.name, None)
                )
                self.func_defs.append(
                    Offset(
                        node.lineno, node.col_offset, node.name, f"fixture_{node.name}"
                    )
                )
            elif (
                isinstance(decorator, ast.Call)
                and isinstance(decorator.func, ast.Attribute)
                and isinstance(decorator.func.value, ast.Name)
                and decorator.func.value.id == "pytest"
                and decorator.func.attr == "fixture"
            ):
                if not any(keyword.arg == "name" for keyword in decorator.keywords):
                    self.fixture_calls.append(
                        Offset(decorator.lineno, decorator.col_offset, node.name, None)
                    )
                    self.func_defs.append(
                        Offset(
                            node.lineno,
                            node.col_offset,
                            node.name,
                            f"fixture_{node.name}",
                        )
                    )
                for keyword in decorator.keywords:
                    if (
                        (keyword.arg == "name")
                        and isinstance(keyword.value, ast.Constant)
                        and (f"fixture_{keyword.value.value}" != node.name)
                    ):
                        self.misnamed_test.append(
                            Offset(
                                node.lineno,
                                node.col_offset,
                                node.name,
                                f"fixture_{keyword.value.value}",
                            )
                        )
        self.generic_visit(node)


def main(content: str, file: str) -> str | None:
    """
    Rename fixture so it can be used with pylint's redefined-outer-name.

    If there is a fixture which doesn't have `name=`, then it is rewritten.
    If no rewrite is necessary, then `None` is returned.
    """
    lines = content.splitlines(keepends=True)
    tree = ast.parse(content)
    visitor = FixtureVisitor(file)
    visitor.visit(tree)

    ret = 0
    for lineno, col_offset, name, new_name in visitor.func_defs:
        line = lines[lineno - 1]
        subbed = re.sub(rf"^def {name}", f"def {new_name}", line[col_offset:])
        lines[lineno - 1] = line[:col_offset] + subbed
        print(f"{file}:{lineno}:{col_offset}: renamed {name} to {new_name}")
        ret |= 1
    for lineno, col_offset, name, _ in visitor.fixture_calls:
        line = lines[lineno - 1]
        # If there are existing arguments, e.g. `pytest.fixture(autouse=True)`
        subbed = re.sub(
            r"^pytest\.fixture\(([^\)])",
            f'pytest.fixture(name="{name}", \\1',
            line[col_offset:],
        )
        # If there are no existing arguments, e.g. `pytest.fixture()`
        subbed = re.sub(
            r"^pytest\.fixture\(\)", f'pytest.fixture(name="{name}")', subbed
        )
        lines[lineno - 1] = line[:col_offset] + subbed
        ret |= 1
    for lineno, col_offset, name, _ in visitor.fixture_attributes:
        line = lines[lineno - 1]
        subbed = re.sub(
            r"^pytest\.fixture", f'pytest.fixture(name="{name}")', line[col_offset:]
        )
        lines[lineno - 1] = line[:col_offset] + subbed
        ret |= 1
    for lineno, col_offset, name, new_name in visitor.misnamed_test:
        line = lines[lineno - 1]
        subbed = re.sub(rf"^def {name}", f"def {new_name}", line[col_offset:])
        lines[lineno - 1] = line[:col_offset] + subbed
        print(f"{file}:{lineno}:{col_offset}: renamed {name} to {new_name}")
        ret |= 1

    if ret:
        new_content = "".join(lines)
        return new_content
    return None


if __name__ == "__main__":  # pragma: no cover
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args()

    ret = 0
    for file in args.paths:
        with open(file, encoding="utf-8") as fd:
            content = fd.read()
        new_content = main(content, file)
        if new_content is not None:
            with open(file, "w", encoding="utf-8") as fd:
                fd.write(new_content)
            ret |= 1

    sys.exit(ret)
