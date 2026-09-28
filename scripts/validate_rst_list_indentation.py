"""
Validate that lists in rst files are not indented relative to their context.

An indented list is parsed as a list inside a block quote, which renders with
extra indentation. Top-level lists should not be indented, and nested lists
should align with the text of their parent item.

Usage::

As pre-commit hook (recommended):
    pre-commit run rst-list-indentation --all-files

From the command-line:
    python scripts/validate_rst_list_indentation.py <rst file>
"""

from __future__ import annotations

import argparse
import sys

from docutils import nodes
from docutils.core import publish_doctree

LIST_NODES = (nodes.bullet_list, nodes.enumerated_list)


def find_indented_lists(content: str, source_path: str = "<string>") -> list[int]:
    """
    Find lists that docutils parses as the sole content of a block quote.

    Parameters
    ----------
    content : str
        The rst text to check.
    source_path : str, default "<string>"
        Path used by docutils when resolving relative references.

    Returns
    -------
    list[int]
        Line numbers of the indented lists.
    """
    doctree = publish_doctree(
        content,
        source_path=source_path,
        settings_overrides={
            # Sphinx-only directives and roles are unknown to docutils; ignore
            # the resulting errors, their content is not parsed as rst.
            "report_level": 5,
            "halt_level": 5,
            "file_insertion_enabled": False,
        },
    )
    line_numbers = []
    for block_quote in doctree.findall(nodes.block_quote):
        children = [
            child
            for child in block_quote.children
            if not isinstance(child, nodes.attribution)
        ]
        if children and all(isinstance(child, LIST_NODES) for child in children):
            # list nodes don't carry a line number; use the first descendant's
            line = next(node.line for node in block_quote.findall() if node.line)
            line_numbers.append(line)
    return line_numbers


def main(source_paths: list[str]) -> int:
    """
    Print the location of every indented list.

    Parameters
    ----------
    source_paths : list[str]
        rst files to validate.

    Returns
    -------
    int
        1 if any indented lists were found, 0 otherwise.
    """
    number_of_errors = 0
    for filename in source_paths:
        with open(filename, encoding="utf-8") as file_obj:
            content = file_obj.read()
        for line in find_indented_lists(content, filename):
            print(
                f"{filename}:{line}: List is indented relative to its context "
                "(rendered as a block quote)"
            )
            number_of_errors += 1
    return int(number_of_errors > 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate rst list indentation")
    parser.add_argument("paths", nargs="*", help="rst files to check.")
    args = parser.parse_args()
    sys.exit(main(args.paths))
