"""The basic dict based notebook format.

Authors:

* Brian Granger
"""

# -----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE, distributed as part of this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------
from __future__ import annotations

from nbformat._struct import Struct

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------


class NotebookNode(Struct):
    """A notebook node object."""


def from_dict(d):
    """Create notebook node(s) from an object."""
    if isinstance(d, dict):
        newd = NotebookNode()
        for k, v in d.items():
            newd[k] = from_dict(v)
        return newd
    if isinstance(d, (tuple, list)):
        return [from_dict(i) for i in d]
    return d


def new_code_cell(code=None, prompt_number=None):
    """Create a new code cell with input and output"""
    cell = NotebookNode()
    cell.cell_type = "code"
    if code is not None:
        cell.code = str(code)
    if prompt_number is not None:
        cell.prompt_number = int(prompt_number)
    return cell


def new_text_cell(text=None):
    """Create a new text cell."""
    cell = NotebookNode()
    if text is not None:
        cell.text = str(text)
    cell.cell_type = "text"
    return cell


def new_notebook(cells=None):
    """Create a notebook by name, id and a list of worksheets."""
    nb = NotebookNode()
    if cells is not None:
        nb.cells = cells
    else:
        nb.cells = []
    return nb
