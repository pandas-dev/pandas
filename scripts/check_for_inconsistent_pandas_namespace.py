"""
Check that test suite file doesn't use the pandas namespace inconsistently.

We check for cases of ``Series`` and ``pd.Series`` appearing in the same file
(likewise for some other common classes).

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run inconsistent-namespace-usage --all-files

To automatically fixup a given file, you can pass `--replace`, e.g.

    python scripts/check_for_inconsistent_pandas_namespace.py test_me.py --replace

Note that you may need to manually fixup some imports.
"""

import argparse
import ast
from typing import MutableMapping, Optional, Sequence, Set

from tokenize_rt import Offset, reversed_enumerate, src_to_tokens, tokens_to_src

ERROR_MESSAGE = "Found both `pd.{name}` and `{name}` in {path}"
EXCLUDE = {
    "np",  # pd.np is deprecated but still tested
    "eval",
}


class Visitor(ast.NodeVisitor):
    def __init__(self) -> None:
        self.pandas_namespace: MutableMapping[Offset, str] = {}
        self.no_namespace: Set[str] = set()

    def visit_Attribute(self, node: ast.Attribute) -> None:
        if (
            isinstance(node.value, ast.Name)
            and node.value.id == "pd"
            and node.attr not in EXCLUDE
        ):
            self.pandas_namespace[Offset(node.lineno, node.col_offset)] = node.attr
        self.generic_visit(node)

    def visit_Name(self, node: ast.Name) -> None:
        if node.id not in EXCLUDE:
            self.no_namespace.add(node.id)
        self.generic_visit(node)


def check_for_inconsistent_pandas_namespace(
    content: str, path: str, *, replace: bool
) -> Optional[str]:
    tree = ast.parse(content)

    visitor = Visitor()
    visitor.visit(tree)

    tokens = src_to_tokens(content)

    if not visitor.no_namespace.intersection(visitor.pandas_namespace.values()):
        # No inconsistent namespace usage, nothing to replace.
        return content

    for n, i in reversed_enumerate(tokens):
        if (
            i.offset in visitor.pandas_namespace
            and visitor.pandas_namespace[i.offset] in visitor.no_namespace
        ):
            if not replace:
                raise RuntimeError(
                    ERROR_MESSAGE.format(
                        name=visitor.pandas_namespace[i.offset], path=path
                    )
                )
            # Replace `pd`
            tokens[n] = i._replace(src="")
            # Replace `.`
            tokens[n + 1] = i._replace(src="")

    new_src: str = tokens_to_src(tokens)
    return new_src


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args(argv)

    for path in args.paths:
        with open(path, encoding="utf-8") as fd:
            content = fd.read()
        new_content = check_for_inconsistent_pandas_namespace(
            content, path, replace=args.replace
        )
        if not args.replace or new_content is None:
            continue
        with open(path, "w", encoding="utf-8") as fd:
            fd.write(new_content)


if __name__ == "__main__":
    main()
