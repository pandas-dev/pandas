"""Deprecated API for working with notebooks

- use nbformat for read/write/validate public API
- use nbformat.vX directly for Python API for composing notebooks
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import re
import warnings

from traitlets.log import get_logger

from nbformat import v3 as _v_latest
from nbformat.v3 import (
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
    parse_filename,
    to_notebook_json,
)

from . import versions
from .converter import convert
from .reader import reads as reader_reads
from .validator import ValidationError, validate

warnings.warn(
    """nbformat.current is deprecated since before nbformat 3.0

- use nbformat for read/write/validate public API
- use nbformat.vX directly to composing notebooks of a particular version
""",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "NotebookNode",
    "new_code_cell",
    "new_text_cell",
    "new_notebook",
    "new_output",
    "new_worksheet",
    "parse_filename",
    "new_metadata",
    "new_author",
    "new_heading_cell",
    "nbformat",
    "nbformat_minor",
    "nbformat_schema",
    "to_notebook_json",
    "convert",
    "validate",
    "NBFormatError",
    "parse_py",
    "reads_json",
    "writes_json",
    "reads_py",
    "writes_py",
    "reads",
    "writes",
    "read",
    "write",
]

current_nbformat = nbformat
current_nbformat_minor = nbformat_minor
current_nbformat_module = _v_latest.__name__


class NBFormatError(ValueError):
    """An error raised for an nbformat error."""


def _warn_format():
    warnings.warn(
        """Non-JSON file support in nbformat is deprecated since nbformat 1.0.
    Use nbconvert to create files of other formats.""",
        stacklevel=2,
    )


def parse_py(s, **kwargs):
    """Parse a string into a (nbformat, string) tuple."""
    nbf = current_nbformat
    nbm = current_nbformat_minor

    pattern = r"# <nbformat>(?P<nbformat>\d+[\.\d+]*)</nbformat>"
    m = re.search(pattern, s)
    if m is not None:
        digits = m.group("nbformat").split(".")
        nbf = int(digits[0])
        if len(digits) > 1:
            nbm = int(digits[1])

    return nbf, nbm, s


def reads_json(nbjson, **kwargs):
    """DEPRECATED, use reads"""
    warnings.warn(
        "reads_json is deprecated since nbformat 3.0, use reads",
        DeprecationWarning,
        stacklevel=2,
    )
    return reads(nbjson)


def writes_json(nb, **kwargs):
    """DEPRECATED, use writes"""
    warnings.warn(
        "writes_json is deprecated since nbformat 3.0, use writes",
        DeprecationWarning,
        stacklevel=2,
    )
    return writes(nb, **kwargs)


def reads_py(s, **kwargs):
    """DEPRECATED: use nbconvert"""
    _warn_format()
    nbf, nbm, s = parse_py(s, **kwargs)
    if nbf in (2, 3):
        nb = versions[nbf].to_notebook_py(s, **kwargs)
    else:
        raise NBFormatError("Unsupported PY nbformat version: %i" % nbf)
    return nb


def writes_py(nb, **kwargs):
    """DEPRECATED: use nbconvert"""
    _warn_format()
    return versions[3].writes_py(nb, **kwargs)


# High level API


def reads(s, format="DEPRECATED", version=current_nbformat, **kwargs):
    """Read a notebook from a string and return the NotebookNode object.

    This function properly handles notebooks of any version. The notebook
    returned will always be in the current version's format.

    Parameters
    ----------
    s : unicode
        The raw unicode string to read the notebook from.

    Returns
    -------
    nb : NotebookNode
        The notebook that was read.
    """
    if format not in {"DEPRECATED", "json"}:
        _warn_format()
    nb = reader_reads(s, **kwargs)
    nb = convert(nb, version)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            validate(nb, repair_duplicate_cell_ids=False)
    except ValidationError as e:
        get_logger().error("Notebook JSON is invalid: %s", e)
    return nb


def writes(nb, format="DEPRECATED", version=current_nbformat, **kwargs):
    """Write a notebook to a string in a given format in the current nbformat version.

    This function always writes the notebook in the current nbformat version.

    Parameters
    ----------
    nb : NotebookNode
        The notebook to write.
    version : int
        The nbformat version to write.
        Used for downgrading notebooks.

    Returns
    -------
    s : unicode
        The notebook string.
    """
    if format not in {"DEPRECATED", "json"}:
        _warn_format()
    nb = convert(nb, version)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            validate(nb, repair_duplicate_cell_ids=False)
    except ValidationError as e:
        get_logger().error("Notebook JSON is invalid: %s", e)
    return versions[version].writes_json(nb, **kwargs)


def read(fp, format="DEPRECATED", **kwargs):
    """Read a notebook from a file and return the NotebookNode object.

    This function properly handles notebooks of any version. The notebook
    returned will always be in the current version's format.

    Parameters
    ----------
    fp : file
        Any file-like object with a read method.

    Returns
    -------
    nb : NotebookNode
        The notebook that was read.
    """
    return reads(fp.read(), **kwargs)


def write(nb, fp, format="DEPRECATED", **kwargs):
    """Write a notebook to a file in a given format in the current nbformat version.

    This function always writes the notebook in the current nbformat version.

    Parameters
    ----------
    nb : NotebookNode
        The notebook to write.
    fp : file
        Any file-like object with a write method.
    """
    s = writes(nb, **kwargs)
    if isinstance(s, bytes):
        s = s.decode("utf8")
    return fp.write(s)
