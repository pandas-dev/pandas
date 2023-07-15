"""
Check that test suite file doesn't use the pandas namespace inconsistently.

We check for cases of ``Series`` and ``pd.Series`` appearing in the same file
(likewise for other pandas objects).

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run inconsistent-namespace-usage --all-files

To automatically fixup a given file, you can pass `--replace`, e.g.

    python scripts/check_for_inconsistent_pandas_namespace.py test_me.py --replace

though note that you may need to manually fixup some imports and that you will also
need the additional dependency `tokenize-rt` (which is left out from the pre-commit
hook so that it uses the same virtualenv as the other local ones).

The general structure is similar to that of some plugins from
https://github.com/asottile/pyupgrade .
"""

import argparse
import ast
from collections.abc import (
    MutableMapping,
    Sequence,
)
import sys
from typing import (
    NamedTuple,
    Optional,
)

ERROR_MESSAGE = (
    "{path}:{lineno}:{col_offset}: "
    "Found both '{prefix}.{name}' and '{name}' in {path}"
)


class OffsetWithNamespace(NamedTuple):
    lineno: int
    col_offset: int
    namespace: str


class Visitor(ast.NodeVisitor):
    def __init__(self) -> None:
        self.pandas_namespace: MutableMapping[OffsetWithNamespace, str] = {}
        self.imported_from_pandas: set[str] = set()

    def visit_Attribute(self, node: ast.Attribute) -> None:
        if isinstance(node.value, ast.Name) and node.value.id in {"pandas", "pd"}:
            offset_with_namespace = OffsetWithNamespace(
                node.lineno, node.col_offset, node.value.id
            )
            self.pandas_namespace[offset_with_namespace] = node.attr
        self.generic_visit(node)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        if node.module is not None and "pandas" in node.module:
            self.imported_from_pandas.update(name.name for name in node.names)
        self.generic_visit(node)


def replace_inconsistent_pandas_namespace(visitor: Visitor, content: str) -> str:
    from tokenize_rt import (
        reversed_enumerate,
        src_to_tokens,
        tokens_to_src,
    )

    tokens = src_to_tokens(content)
    for n, i in reversed_enumerate(tokens):
        offset_with_namespace = OffsetWithNamespace(i.offset[0], i.offset[1], i.src)
        if (
            offset_with_namespace in visitor.pandas_namespace
            and visitor.pandas_namespace[offset_with_namespace]
            in visitor.imported_from_pandas
        ):
            # Replace `pd`
            tokens[n] = i._replace(src="")
            # Replace `.`
            tokens[n + 1] = tokens[n + 1]._replace(src="")

    new_src: str = tokens_to_src(tokens)
    return new_src


def check_for_inconsistent_pandas_namespace(
    content: str, path: str, *, replace: bool
) -> Optional[str]:
    tree = ast.parse(content)

    visitor = Visitor()
    visitor.visit(tree)

    inconsistencies = visitor.imported_from_pandas.intersection(
        visitor.pandas_namespace.values()
    )

    if not inconsistencies:
        # No inconsistent namespace usage, nothing to replace.
        return None

    if not replace:
        inconsistency = inconsistencies.pop()
        lineno, col_offset, prefix = next(
            key for key, val in visitor.pandas_namespace.items() if val == inconsistency
        )
        msg = ERROR_MESSAGE.format(
            lineno=lineno,
            col_offset=col_offset,
            prefix=prefix,
            name=inconsistency,
            path=path,
        )
        sys.stdout.write(msg)
        sys.exit(1)

    return replace_inconsistent_pandas_namespace(visitor, content)


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
