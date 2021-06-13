"""
Check that pandas/core/generic.py doesn't use bool as a type annotation.

There is already the method `bool`, so the alias `bool_t` should be used instead.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run no-bool-in-core-generic --all-files

The function `visit` is adapted from a function by the same name in pyupgrade:
https://github.com/asottile/pyupgrade/blob/5495a248f2165941c5d3b82ac3226ba7ad1fa59d/pyupgrade/_data.py#L70-L113
"""
from __future__ import annotations

import argparse
import ast
import collections
from typing import Sequence


def visit(tree: ast.Module) -> dict[int, list[int]]:
    "Step through tree, recording when nodes are in annotations."
    in_annotation = False
    nodes: list[tuple[bool, ast.AST]] = [(in_annotation, tree)]
    to_replace = collections.defaultdict(list)

    while nodes:
        in_annotation, node = nodes.pop()

        if isinstance(node, ast.Name) and in_annotation and node.id == "bool":
            to_replace[node.lineno].append(node.col_offset)

        for name in reversed(node._fields):
            value = getattr(node, name)
            if name in {"annotation", "returns"}:
                next_in_annotation = True
            else:
                next_in_annotation = in_annotation
            if isinstance(value, ast.AST):
                nodes.append((next_in_annotation, value))
            elif isinstance(value, list):
                for value in reversed(value):
                    if isinstance(value, ast.AST):
                        nodes.append((next_in_annotation, value))

    return to_replace


def replace_bool_with_bool_t(to_replace, content: str) -> str:
    new_lines = []

    for n, line in enumerate(content.splitlines(), start=1):
        if n in to_replace:
            for col_offset in reversed(to_replace[n]):
                line = line[:col_offset] + "bool_t" + line[col_offset + 4 :]
        new_lines.append(line)
    return "\n".join(new_lines)


def check_for_bool_in_generic(content: str) -> tuple[bool, str]:
    tree = ast.parse(content)
    to_replace = visit(tree)

    if not to_replace:
        mutated = False
        return mutated, content

    mutated = True
    return mutated, replace_bool_with_bool_t(to_replace, content)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)

    for path in args.paths:
        with open(path, encoding="utf-8") as fd:
            content = fd.read()
        mutated, new_content = check_for_bool_in_generic(content)
        if mutated:
            with open(path, "w", encoding="utf-8") as fd:
                fd.write(new_content)


if __name__ == "__main__":
    main()
