"""Sphinx extension that warns about lists indented relative to their context.

An indented list is parsed as a list inside a block quote, which renders with
extra indentation. Top-level lists should not be indented, and nested lists
should align with the text of their parent item.
"""

import re

from docutils import nodes
from sphinx.util import logging

logger = logging.getLogger(__name__)

LIST_NODES = (nodes.bullet_list, nodes.enumerated_list)

# fixing the pre-2.0 release notes would rewrite ~900 rarely-read lines
EXCLUDE_DOCNAMES = re.compile(r"whatsnew/v(0|1|2\.0\.0)")


def find_indented_lists(doctree: nodes.document) -> list[nodes.block_quote]:
    """
    Find lists that are parsed as the sole content of a block quote.

    Parameters
    ----------
    doctree : docutils.nodes.document
        The parsed document to check.

    Returns
    -------
    list[docutils.nodes.block_quote]
        The block quotes wrapping each indented list.
    """
    found = []
    for block_quote in doctree.findall(nodes.block_quote):
        children = [
            child
            for child in block_quote.children
            if not isinstance(child, nodes.attribution)
        ]
        if children and all(isinstance(child, LIST_NODES) for child in children):
            found.append(block_quote)
    return found


def check_list_indentation(app, doctree: nodes.document) -> None:
    if EXCLUDE_DOCNAMES.match(app.env.docname):
        return
    for node in find_indented_lists(doctree):
        logger.warning(
            "List is indented relative to its context (rendered as a block quote)",
            location=node,
            type="pandas",
            subtype="list_indentation",
        )


def setup(app):
    app.connect("doctree-read", check_list_indentation)

    return {"version": "0.1", "parallel_read_safe": True, "parallel_write_safe": True}
