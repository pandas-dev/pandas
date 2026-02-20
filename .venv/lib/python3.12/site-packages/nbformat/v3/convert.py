"""Code for converting notebooks to and from the v2 format."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from .nbbase import nbformat, nbformat_minor


def _unbytes(obj):
    """There should be no bytes objects in a notebook

    v2 stores png/jpeg as b64 ascii bytes
    """
    if isinstance(obj, dict):
        for k, v in obj.items():
            obj[k] = _unbytes(v)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            obj[i] = _unbytes(v)
    elif isinstance(obj, bytes):
        # only valid bytes are b64-encoded ascii
        obj = obj.decode("ascii")
    return obj


def upgrade(nb, from_version=2, from_minor=0):
    """Convert a notebook to v3.

    Parameters
    ----------
    nb : NotebookNode
        The Python representation of the notebook to convert.
    from_version : int
        The original version of the notebook to convert.
    from_minor : int
        The original minor version of the notebook to convert (only relevant for v >= 3).
    """
    if from_version == 2:
        # Mark the original nbformat so consumers know it has been converted.
        nb.nbformat = nbformat
        nb.nbformat_minor = nbformat_minor

        nb.orig_nbformat = 2
        nb = _unbytes(nb)
        for ws in nb["worksheets"]:
            for cell in ws["cells"]:
                cell.setdefault("metadata", {})
        return nb
    if from_version == 3:
        if from_minor != nbformat_minor:
            nb.orig_nbformat_minor = from_minor
        nb.nbformat_minor = nbformat_minor
        return nb
    msg = (
        "Cannot convert a notebook directly from v%s to v3.  "
        "Try using the nbformat.convert module." % from_version
    )
    raise ValueError(msg)


def heading_to_md(cell):
    """turn heading cell into corresponding markdown"""
    cell.cell_type = "markdown"
    level = cell.pop("level", 1)
    cell.source = "#" * level + " " + cell.source


def raw_to_md(cell):
    """let raw passthrough as markdown"""
    cell.cell_type = "markdown"


def downgrade(nb):
    """Convert a v3 notebook to v2.

    Parameters
    ----------
    nb : NotebookNode
        The Python representation of the notebook to convert.
    """
    if nb.nbformat != 3:
        return nb
    nb.nbformat = 2
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == "heading":
                heading_to_md(cell)
            elif cell.cell_type == "raw":
                raw_to_md(cell)
    return nb
