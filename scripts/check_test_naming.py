from __future__ import annotations

import argparse
import ast
import os
from pathlib import Path
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

PRAGMA = "# not a test"


def find_names(node: ast.Module) -> Iterator[str]:
    for _node in ast.walk(node):
        if isinstance(_node, (ast.Name, ast.Attribute)):
            yield _node.id if isinstance(_node, ast.Name) else _node.attr


def is_fixture(node: ast.expr) -> bool:
    if isinstance(node, ast.Call):
        node = node.func
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "fixture"
        and isinstance(node.value, ast.Name)
        and node.value.id == "pytest"
    )


def is_register_dtype(node):
    return isinstance(node, ast.Name) and node.id == "register_extension_dtype"


def is_misnamed_test_func(node: ast.expr | ast.stmt, names: Sequence[str], line: str) -> bool:
    return (
        isinstance(node, ast.FunctionDef)
        and not node.name.startswith("test")
        and names.count(node.name) == 0
        and not any(_is_fixture(decorator) for decorator in node.decorator_list)
        and PRAGMA not in line
        and node.name not in ("teardown_method", "setup_method", "teardown_class", "setup_class")
    )


def is_misnamed_test_class(node: ast.expr | ast.stmt, names: Sequence[str], line: str) -> bool:
    return (
        isinstance(node, ast.ClassDef)
        and not node.name.startswith("Test")
        and names.count(node.name) == 0
        and not any(_is_register_dtype(decorator) for decorator in node.decorator_list)
        and PRAGMA not in line
    )


def main(content: str, file: str) -> int:
    lines = content.splitlines()
    tree = ast.parse(content)
    names = list(find_names(tree))
    ret = 0

    for node in tree.body:
        if is_misnamed_test_func(node, names, lines[node.lineno - 1]):
            print(
                f"{file}:{node.lineno}:{node.col_offset} found test function which does not start with 'test'"
            )
            ret = 1
        elif is_misnamed_test_class(node, names, lines[node.lineno - 1]):
            print(
                f"{file}:{node.lineno}:{node.col_offset} found test class which does not start with 'Test'"
            )
            ret = 1

        if isinstance(node, ast.ClassDef) and names.count(node.name) == 0 and not any(
            _is_register_dtype(decorator) for decorator in node.decorator_list
        ) and PRAGMA not in lines[node.lineno - 1]:
            for _node in node.body:
                if is_misnamed_test_func(_node, names, lines[_node.lineno - 1]):
                    should_continue = False
                    for _file in (Path("pandas") / "tests").rglob("*.py"):
                        with open(os.path.join(_file), encoding="utf-8") as fd:
                            _content = fd.read()
                        if f"self.{_node.name}" in _content:
                            should_continue = True
                            break

                    if should_continue:
                        continue

                    print(
                        f"{file}:{_node.lineno}:{_node.col_offset} found test function which does not start with 'test'"
                    )
                    ret = 1

    return ret


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args()

    ret = 0

    for file in args.paths:
        filename = os.path.basename(file)
        if not (filename.startswith("test") and filename.endswith(".py")):
            continue
        with open(file, encoding="utf-8") as fd:
            content = fd.read()
        ret |= main(content, file)

    sys.exit(ret)
