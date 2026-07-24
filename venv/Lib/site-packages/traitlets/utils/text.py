"""
Utilities imported from ipython_genutils
"""
from __future__ import annotations

import re
import textwrap
from textwrap import indent as _indent


def indent(val: str) -> str:
    return _indent(val, "    ")


def _dedent(text: str) -> str:
    """Equivalent of textwrap.dedent that ignores unindented first line."""

    if text.startswith("\n"):
        # text starts with blank line, don't ignore the first line
        return textwrap.dedent(text)

    # split first line
    splits = text.split("\n", 1)
    if len(splits) == 1:
        # only one line
        return textwrap.dedent(text)

    first, rest = splits
    # dedent everything but the first line
    rest = textwrap.dedent(rest)
    return "\n".join([first, rest])


def wrap_paragraphs(text: str, ncols: int = 80) -> list[str]:
    """Wrap multiple paragraphs to fit a specified width.

    This is equivalent to textwrap.wrap, but with support for multiple
    paragraphs, as separated by empty lines.

    Returns
    -------

    list of complete paragraphs, wrapped to fill `ncols` columns.
    """
    paragraph_re = re.compile(r"\n(\s*\n)+", re.MULTILINE)
    text = _dedent(text).strip()
    paragraphs = paragraph_re.split(text)[::2]  # every other entry is space
    out_ps = []
    indent_re = re.compile(r"\n\s+", re.MULTILINE)
    for p in paragraphs:
        # presume indentation that survives dedent is meaningful formatting,
        # so don't fill unless text is flush.
        if indent_re.search(p) is None:
            # wrap paragraph
            p = textwrap.fill(p, ncols)
        out_ps.append(p)
    return out_ps
