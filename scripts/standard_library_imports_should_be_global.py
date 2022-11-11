"""
Check that standard library imports appear at the top of modules.

Imports within functions should only be used to prevent circular imports
or for optional dependencies.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run stdlib-imports --all-files

"""
import argparse
import ast
from ast import NodeVisitor
import importlib
import sys


class Visitor(NodeVisitor):
    def __init__(self, file) -> None:
        self.ret = 0
        self.file = file

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        for _node in ast.walk(node):
            if isinstance(_node, ast.ImportFrom):
                try:
                    importlib.import_module(_node.module)
                except Exception as exp:  # noqa: F841
                    pass
                else:
                    print(
                        f"{self.file}:{_node.lineno}:{_node.col_offset} standard "
                        f"library import {_node.module} should be global"
                    )
                    self.ret = 1
            elif isinstance(_node, ast.Import):
                try:
                    _name = _node.names[0].name
                    importlib.import_module(_name)
                except Exception as _:  # noqa: F841
                    pass
                else:
                    print(
                        f"{self.file}:{_node.lineno}:{_node.col_offset} standard "
                        f"library import {_name} should be global"
                    )
                    self.ret = 1
        self.generic_visit(node)

    # def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
    #     if node.col_offset > 0:
    #     print(f"{self.file}:{node.lineno}:{node.col_offset}")
    #     self.generic_visit(node)

    # def visit_Import(self, node: ast.Import) -> None:
    #     breakpoint()
    #     print(f"{self.file}:{node.lineno}:{node.col_offset}")
    #     self.generic_visit(node)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args()
    ret = 0
    for file in args.paths:
        with open(file, encoding="utf-8") as fd:
            content = fd.read()
        tree = ast.parse(content)
        visitor = Visitor(file)
        visitor.visit(tree)
        ret |= visitor.ret
    sys.exit(ret)
