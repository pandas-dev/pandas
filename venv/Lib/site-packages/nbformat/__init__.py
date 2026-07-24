"""The Jupyter notebook format

Use this module to read or write notebook files as particular nbformat versions.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from pathlib import Path

from traitlets.log import get_logger

from . import v1, v2, v3, v4
from ._version import __version__, version_info
from .sentinel import Sentinel

__all__ = [
    "versions",
    "validate",
    "ValidationError",
    "convert",
    "from_dict",
    "NotebookNode",
    "current_nbformat",
    "current_nbformat_minor",
    "NBFormatError",
    "NO_CONVERT",
    "reads",
    "read",
    "writes",
    "write",
    "version_info",
    "__version__",
    "Sentinel",
]

versions = {
    1: v1,
    2: v2,
    3: v3,
    4: v4,
}

from . import reader  # noqa: E402
from .converter import convert  # noqa: E402
from .notebooknode import NotebookNode, from_dict  # noqa: E402
from .v4 import nbformat as current_nbformat  # noqa: E402
from .v4 import nbformat_minor as current_nbformat_minor  # noqa: E402
from .validator import ValidationError, validate  # noqa: E402


class NBFormatError(ValueError):
    pass


# no-conversion singleton
NO_CONVERT = Sentinel(
    "NO_CONVERT",
    __name__,
    """Value to prevent nbformat to convert notebooks to most recent version.
    """,
)


def reads(s, as_version, capture_validation_error=None, **kwargs):
    """Read a notebook from a string and return the NotebookNode object as the given version.

    The string can contain a notebook of any version.
    The notebook will be returned `as_version`, converting, if necessary.

    Notebook format errors will be logged.

    Parameters
    ----------
    s : unicode
        The raw unicode string to read the notebook from.
    as_version : int
        The version of the notebook format to return.
        The notebook will be converted, if necessary.
        Pass nbformat.NO_CONVERT to prevent conversion.
    capture_validation_error : dict, optional
        If provided, a key of "ValidationError" with a
        value of the ValidationError instance will be added
        to the dictionary.

    Returns
    -------
    nb : NotebookNode
        The notebook that was read.
    """
    nb = reader.reads(s, **kwargs)
    if as_version is not NO_CONVERT:
        nb = convert(nb, as_version)
    try:
        validate(nb)
    except ValidationError as e:
        get_logger().error("Notebook JSON is invalid: %s", e)
        if isinstance(capture_validation_error, dict):
            capture_validation_error["ValidationError"] = e
    return nb


def writes(nb, version=NO_CONVERT, capture_validation_error=None, **kwargs):
    """Write a notebook to a string in a given format in the given nbformat version.

    Any notebook format errors will be logged.

    Parameters
    ----------
    nb : NotebookNode
        The notebook to write.
    version : int, optional
        The nbformat version to write.
        If unspecified, or specified as nbformat.NO_CONVERT,
        the notebook's own version will be used and no conversion performed.
    capture_validation_error : dict, optional
        If provided, a key of "ValidationError" with a
        value of the ValidationError instance will be added
        to the dictionary.

    Returns
    -------
    s : unicode
        The notebook as a JSON string.
    """
    if version is not NO_CONVERT:
        nb = convert(nb, version)
    else:
        version, _ = reader.get_version(nb)
    try:
        validate(nb)
    except ValidationError as e:
        get_logger().error("Notebook JSON is invalid: %s", e)
        if isinstance(capture_validation_error, dict):
            capture_validation_error["ValidationError"] = e
    return versions[version].writes_json(nb, **kwargs)


def read(fp, as_version, capture_validation_error=None, **kwargs):
    """Read a notebook from a file as a NotebookNode of the given version.

    The string can contain a notebook of any version.
    The notebook will be returned `as_version`, converting, if necessary.

    Notebook format errors will be logged.

    Parameters
    ----------
    fp : file or str
        A file-like object with a read method that returns unicode (use
        ``io.open()`` in Python 2), or a path to a file.
    as_version : int
        The version of the notebook format to return.
        The notebook will be converted, if necessary.
        Pass nbformat.NO_CONVERT to prevent conversion.
    capture_validation_error : dict, optional
        If provided, a key of "ValidationError" with a
        value of the ValidationError instance will be added
        to the dictionary.

    Returns
    -------
    nb : NotebookNode
        The notebook that was read.
    """

    try:
        buf = fp.read()
    except AttributeError:
        with open(fp, encoding="utf8") as f:  # noqa: PTH123
            return reads(f.read(), as_version, capture_validation_error, **kwargs)

    return reads(buf, as_version, capture_validation_error, **kwargs)


def write(nb, fp, version=NO_CONVERT, capture_validation_error=None, **kwargs):
    """Write a notebook to a file in a given nbformat version.

    The file-like object must accept unicode input.

    Parameters
    ----------
    nb : NotebookNode
        The notebook to write.
    fp : file or str
        Any file-like object with a write method that accepts unicode, or
        a path to write a file.
    version : int, optional
        The nbformat version to write.
        If nb is not this version, it will be converted.
        If unspecified, or specified as nbformat.NO_CONVERT,
        the notebook's own version will be used and no conversion performed.
    capture_validation_error : dict, optional
        If provided, a key of "ValidationError" with a
        value of the ValidationError instance will be added
        to the dictionary.
    """
    s = writes(nb, version, capture_validation_error, **kwargs)
    if isinstance(s, bytes):
        s = s.decode("utf8")

    try:
        fp.write(s)
        if not s.endswith("\n"):
            fp.write("\n")
    except AttributeError:
        with Path(fp).open("w", encoding="utf8") as f:
            f.write(s)
            if not s.endswith("\n"):
                f.write("\n")
