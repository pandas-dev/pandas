"""
Check that standard library imports appear at the top of modules.

Imports within functions should only be used to prevent circular imports
, for optional dependencies, or if an import is slow.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run stdlib-imports --all-files

"""
import argparse
import ast
from ast import NodeVisitor
import importlib
import sys

BLOCKLIST = {
    "bs4",
    "ctypes",
    "gcsfs",
    "html5lib",
    "http",
    "importlib.metadata",
    "ipython",
    "jinja2",
    "hypothesis",
    "lxml",
    "matplotlib",
    "openpyxl",
    "py",
    "pytest",
    "s3fs",
    "scipy",
    "sqlite3",
    "tables",
    "urllib.error",
    "urllib.request",
    "xlrd",
    "xlsxwriter",
    "xml",
}


class Visitor(NodeVisitor):
    def __init__(self, file) -> None:
        self.ret = 0
        self.file = file

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        for _node in ast.walk(node):
            if (
                isinstance(_node, ast.ImportFrom)
                and _node.__module__ != "__main__"
                and _node.module not in BLOCKLIST
                and _node.module.split(".")[0] not in BLOCKLIST
            ):
                try:
                    importlib.import_module(_node.module)
                except Exception as exp:  # noqa: F841
                    pass
                else:
                    print(
                        f"{self.file}:{_node.lineno}:{_node.col_offset} standard "
                        f"library import '{_node.module}' should be at the top of "
                        "the file"
                    )
                    self.ret = 1
            elif isinstance(_node, ast.Import):
                for _name in _node.names:
                    if (
                        _name.name == "__main__"
                        or _name.name in BLOCKLIST
                        or _name.name.split(".")[0] in BLOCKLIST
                    ):
                        continue
                    try:
                        importlib.import_module(_name.name)
                    except Exception as exp:  # noqa: F841
                        pass
                    else:
                        print(
                            f"{self.file}:{_node.lineno}:{_node.col_offset} standard "
                            f"library import '{_name.name}' should be at the top of "
                            "the file"
                        )
                        self.ret = 1
        self.generic_visit(node)


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
