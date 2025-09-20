"""The main API for the v3 notebook format."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

__all__ = [
    "NotebookNode",
    "new_code_cell",
    "new_text_cell",
    "new_notebook",
    "new_output",
    "new_worksheet",
    "new_metadata",
    "new_author",
    "new_heading_cell",
    "nbformat",
    "nbformat_minor",
    "nbformat_schema",
    "reads_json",
    "writes_json",
    "read_json",
    "write_json",
    "to_notebook_json",
    "reads_py",
    "writes_py",
    "read_py",
    "write_py",
    "to_notebook_py",
    "downgrade",
    "upgrade",
    "parse_filename",
]

import os

from .convert import downgrade, upgrade
from .nbbase import (
    NotebookNode,
    nbformat,
    nbformat_minor,
    nbformat_schema,
    new_author,
    new_code_cell,
    new_heading_cell,
    new_metadata,
    new_notebook,
    new_output,
    new_text_cell,
    new_worksheet,
)
from .nbjson import reads as read_json
from .nbjson import reads as reads_json
from .nbjson import to_notebook as to_notebook_json
from .nbjson import writes as write_json
from .nbjson import writes as writes_json
from .nbpy import reads as read_py
from .nbpy import reads as reads_py
from .nbpy import to_notebook as to_notebook_py
from .nbpy import writes as write_py
from .nbpy import writes as writes_py


def parse_filename(fname):
    """Parse a notebook filename.

    This function takes a notebook filename and returns the notebook
    format (json/py) and the notebook name. This logic can be
    summarized as follows:

    * notebook.ipynb -> (notebook.ipynb, notebook, json)
    * notebook.json  -> (notebook.json, notebook, json)
    * notebook.py    -> (notebook.py, notebook, py)
    * notebook       -> (notebook.ipynb, notebook, json)

    Parameters
    ----------
    fname : unicode
        The notebook filename. The filename can use a specific filename
        extension (.ipynb, .json, .py) or none, in which case .ipynb will
        be assumed.

    Returns
    -------
    (fname, name, format) : (unicode, unicode, unicode)
        The filename, notebook name and format.
    """
    basename, ext = os.path.splitext(fname)  # noqa: PTH122
    if ext in [".ipynb", ".json"]:
        format_ = "json"
    elif ext == ".py":
        format_ = "py"
    else:
        basename = fname
        fname = fname + ".ipynb"
        format_ = "json"
    return fname, basename, format_
