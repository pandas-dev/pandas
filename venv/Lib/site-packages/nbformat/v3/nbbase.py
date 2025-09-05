"""The basic dict based notebook format.

The Python representation of a notebook is a nested structure of
dictionary subclasses that support attribute access.
The functions in this module are merely
helpers to build the structs in the right form.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import warnings

from nbformat._struct import Struct

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

# Change this when incrementing the nbformat version
nbformat = 3
nbformat_minor = 0
nbformat_schema = {(3, 0): "nbformat.v3.schema.json"}


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


def str_passthrough(obj):
    """
    Used to be cast_unicode, add this temporarily to make sure no further breakage.
    """
    if not isinstance(obj, str):
        raise AssertionError
    return obj


def cast_str(obj):
    """Cast an object as a string."""
    if isinstance(obj, bytes):
        # really this should never happened, it should
        # have been base64 encoded before.
        warnings.warn(
            "A notebook got bytes instead of likely base64 encoded values."
            "The content will likely be corrupted.",
            UserWarning,
            stacklevel=3,
        )
        return obj.decode("ascii", "replace")
    if not isinstance(obj, str):
        raise AssertionError
    return obj


def new_output(
    output_type,
    output_text=None,
    output_png=None,
    output_html=None,
    output_svg=None,
    output_latex=None,
    output_json=None,
    output_javascript=None,
    output_jpeg=None,
    prompt_number=None,
    ename=None,
    evalue=None,
    traceback=None,
    stream=None,
    metadata=None,
):
    """Create a new output, to go in the ``cell.outputs`` list of a code cell."""
    output = NotebookNode()
    output.output_type = str(output_type)

    if metadata is None:
        metadata = {}
    if not isinstance(metadata, dict):
        msg = "metadata must be dict"
        raise TypeError(msg)

    if output_type in {"pyout", "display_data"}:
        output.metadata = metadata

    if output_type != "pyerr":
        if output_text is not None:
            output.text = str_passthrough(output_text)
        if output_png is not None:
            output.png = cast_str(output_png)
        if output_jpeg is not None:
            output.jpeg = cast_str(output_jpeg)
        if output_html is not None:
            output.html = str_passthrough(output_html)
        if output_svg is not None:
            output.svg = str_passthrough(output_svg)
        if output_latex is not None:
            output.latex = str_passthrough(output_latex)
        if output_json is not None:
            output.json = str_passthrough(output_json)
        if output_javascript is not None:
            output.javascript = str_passthrough(output_javascript)

    if output_type == "pyout" and prompt_number is not None:
        output.prompt_number = int(prompt_number)

    if output_type == "pyerr":
        if ename is not None:
            output.ename = str_passthrough(ename)
        if evalue is not None:
            output.evalue = str_passthrough(evalue)
        if traceback is not None:
            output.traceback = [str_passthrough(frame) for frame in list(traceback)]

    if output_type == "stream":
        output.stream = "stdout" if stream is None else str_passthrough(stream)

    return output


def new_code_cell(
    input=None,
    prompt_number=None,
    outputs=None,
    language="python",
    collapsed=False,
    metadata=None,
):
    """Create a new code cell with input and output"""
    cell = NotebookNode()
    cell.cell_type = "code"
    if language is not None:
        cell.language = str_passthrough(language)
    if input is not None:
        cell.input = str_passthrough(input)
    if prompt_number is not None:
        cell.prompt_number = int(prompt_number)
    if outputs is None:
        cell.outputs = []
    else:
        cell.outputs = outputs
    if collapsed is not None:
        cell.collapsed = bool(collapsed)
    cell.metadata = NotebookNode(metadata or {})

    return cell


def new_text_cell(cell_type, source=None, rendered=None, metadata=None):
    """Create a new text cell."""
    cell = NotebookNode()
    # VERSIONHACK: plaintext -> raw
    # handle never-released plaintext name for raw cells
    if cell_type == "plaintext":
        cell_type = "raw"
    if source is not None:
        cell.source = str_passthrough(source)
    cell.metadata = NotebookNode(metadata or {})
    cell.cell_type = cell_type
    return cell


def new_heading_cell(source=None, level=1, rendered=None, metadata=None):
    """Create a new section cell with a given integer level."""
    cell = NotebookNode()
    cell.cell_type = "heading"
    if source is not None:
        cell.source = str_passthrough(source)
    cell.level = int(level)
    cell.metadata = NotebookNode(metadata or {})
    return cell


def new_worksheet(name=None, cells=None, metadata=None):
    """Create a worksheet by name with with a list of cells."""
    ws = NotebookNode()
    if cells is None:
        ws.cells = []
    else:
        ws.cells = list(cells)
    ws.metadata = NotebookNode(metadata or {})
    return ws


def new_notebook(name=None, metadata=None, worksheets=None):
    """Create a notebook by name, id and a list of worksheets."""
    nb = NotebookNode()
    nb.nbformat = nbformat
    nb.nbformat_minor = nbformat_minor
    if worksheets is None:
        nb.worksheets = []
    else:
        nb.worksheets = list(worksheets)
    if metadata is None:
        nb.metadata = new_metadata()
    else:
        nb.metadata = NotebookNode(metadata)
    if name is not None:
        nb.metadata.name = str_passthrough(name)
    return nb


def new_metadata(
    name=None,
    authors=None,
    license=None,
    created=None,
    modified=None,
    gistid=None,
):
    """Create a new metadata node."""
    metadata = NotebookNode()
    if name is not None:
        metadata.name = str_passthrough(name)
    if authors is not None:
        metadata.authors = list(authors)
    if created is not None:
        metadata.created = str_passthrough(created)
    if modified is not None:
        metadata.modified = str_passthrough(modified)
    if license is not None:
        metadata.license = str_passthrough(license)
    if gistid is not None:
        metadata.gistid = str_passthrough(gistid)
    return metadata


def new_author(name=None, email=None, affiliation=None, url=None):
    """Create a new author."""
    author = NotebookNode()
    if name is not None:
        author.name = str_passthrough(name)
    if email is not None:
        author.email = str_passthrough(email)
    if affiliation is not None:
        author.affiliation = str_passthrough(affiliation)
    if url is not None:
        author.url = str_passthrough(url)
    return author
