"""Code for converting notebooks to and from the v2 format.

Authors:

* Brian Granger
* Jonathan Frederic
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

from .nbbase import new_code_cell, new_notebook, new_text_cell, new_worksheet

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------


def upgrade(nb, from_version=1):
    """Convert a notebook to the v2 format.

    Parameters
    ----------
    nb : NotebookNode
        The Python representation of the notebook to convert.
    from_version : int
        The version of the notebook to convert from.
    """
    if from_version == 1:
        newnb = new_notebook()
        ws = new_worksheet()
        for cell in nb.cells:
            if cell.cell_type == "code":
                newcell = new_code_cell(
                    input=cell.get("code"), prompt_number=cell.get("prompt_number")
                )
            elif cell.cell_type == "text":
                newcell = new_text_cell("markdown", source=cell.get("text"))
            ws.cells.append(newcell)
        newnb.worksheets.append(ws)
        return newnb

    raise ValueError("Cannot convert a notebook from v%s to v2" % from_version)


def downgrade(nb):
    """Convert a v2 notebook to v1.

    Parameters
    ----------
    nb : NotebookNode
        The Python representation of the notebook to convert.
    """
    msg = "Downgrade from notebook v2 to v1 is not supported."
    raise Exception(msg)
