"""Code for converting notebooks to and from v3."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import json
import re

from traitlets.log import get_logger

from nbformat import v3, validator
from nbformat.corpus.words import generate_corpus_id as random_cell_id
from nbformat.notebooknode import NotebookNode

from .nbbase import nbformat, nbformat_minor


def _warn_if_invalid(nb, version):
    """Log validation errors, if there are any."""
    from nbformat import ValidationError, validate

    try:
        validate(nb, version=version)
    except ValidationError as e:
        get_logger().error("Notebook JSON is not valid v%i: %s", version, e)


def upgrade(nb, from_version=None, from_minor=None):
    """Convert a notebook to latest v4.

    Parameters
    ----------
    nb : NotebookNode
        The Python representation of the notebook to convert.
    from_version : int
        The original version of the notebook to convert.
    from_minor : int
        The original minor version of the notebook to convert (only relevant for v >= 3).
    """
    if not from_version:
        from_version = nb["nbformat"]
    if not from_minor:
        if "nbformat_minor" not in nb:
            if from_version == 4:
                msg = "The v4 notebook does not include the nbformat minor, which is needed."
                raise validator.ValidationError(msg)
            from_minor = 0
        else:
            from_minor = nb["nbformat_minor"]

    if from_version == 3:
        # Validate the notebook before conversion
        _warn_if_invalid(nb, from_version)

        # Mark the original nbformat so consumers know it has been converted
        orig_nbformat = nb.pop("orig_nbformat", None)
        orig_nbformat_minor = nb.pop("orig_nbformat_minor", None)
        nb.metadata.orig_nbformat = orig_nbformat or 3
        nb.metadata.orig_nbformat_minor = orig_nbformat_minor or 0

        # Mark the new format
        nb.nbformat = nbformat
        nb.nbformat_minor = nbformat_minor

        # remove worksheet(s)
        nb["cells"] = cells = []
        # In the unlikely event of multiple worksheets,
        # they will be flattened
        for ws in nb.pop("worksheets", []):
            # upgrade each cell
            for cell in ws["cells"]:
                cells.append(upgrade_cell(cell))
        # upgrade metadata
        nb.metadata.pop("name", "")
        nb.metadata.pop("signature", "")
        # Validate the converted notebook before returning it
        _warn_if_invalid(nb, nbformat)
        return nb
    if from_version == 4:
        if from_minor == nbformat_minor:
            return nb

        # other versions migration code e.g.
        # if from_minor < 3:
        # if from_minor < 4:

        if from_minor < 5:
            for cell in nb.cells:
                cell.id = random_cell_id()

        nb.metadata.orig_nbformat_minor = from_minor
        nb.nbformat_minor = nbformat_minor

        return nb
    raise ValueError(
        "Cannot convert a notebook directly from v%s to v4.  "
        "Try using the nbformat.convert module." % from_version
    )


def upgrade_cell(cell):
    """upgrade a cell from v3 to v4

    heading cell:
        - -> markdown heading
    code cell:
        - remove language metadata
        - cell.input -> cell.source
        - cell.prompt_number -> cell.execution_count
        - update outputs
    """
    cell.setdefault("metadata", NotebookNode())
    cell.id = random_cell_id()
    if cell.cell_type == "code":
        cell.pop("language", "")
        if "collapsed" in cell:
            cell.metadata["collapsed"] = cell.pop("collapsed")
        cell.source = cell.pop("input", "")
        cell.execution_count = cell.pop("prompt_number", None)
        cell.outputs = upgrade_outputs(cell.outputs)
    elif cell.cell_type == "heading":
        cell.cell_type = "markdown"
        level = cell.pop("level", 1)
        cell.source = "{hashes} {single_line}".format(
            hashes="#" * level,
            single_line=" ".join(cell.get("source", "").splitlines()),
        )
    elif cell.cell_type == "html":
        # Technically, this exists. It will never happen in practice.
        cell.cell_type = "markdown"
    return cell


def downgrade_cell(cell):
    """downgrade a cell from v4 to v3

    code cell:
        - set cell.language
        - cell.input <- cell.source
        - cell.prompt_number <- cell.execution_count
        - update outputs
    markdown cell:
        - single-line heading -> heading cell
    """
    if cell.cell_type == "code":
        cell.language = "python"
        cell.input = cell.pop("source", "")
        cell.prompt_number = cell.pop("execution_count", None)
        cell.collapsed = cell.metadata.pop("collapsed", False)
        cell.outputs = downgrade_outputs(cell.outputs)
    elif cell.cell_type == "markdown":
        source = cell.get("source", "")
        if "\n" not in source and source.startswith("#"):
            match = re.match(r"(#+)\s*(.*)", source)
            assert match is not None
            prefix, text = match.groups()
            cell.cell_type = "heading"
            cell.source = text
            cell.level = len(prefix)
    cell.pop("id", None)
    cell.pop("attachments", None)
    return cell


_mime_map = {
    "text": "text/plain",
    "html": "text/html",
    "svg": "image/svg+xml",
    "png": "image/png",
    "jpeg": "image/jpeg",
    "latex": "text/latex",
    "json": "application/json",
    "javascript": "application/javascript",
}


def to_mime_key(d):
    """convert dict with v3 aliases to plain mime-type keys"""
    for alias, mime in _mime_map.items():
        if alias in d:
            d[mime] = d.pop(alias)
    return d


def from_mime_key(d):
    """convert dict with mime-type keys to v3 aliases"""
    d2 = {}
    for alias, mime in _mime_map.items():
        if mime in d:
            d2[alias] = d[mime]
    return d2


def upgrade_output(output):
    """upgrade a single code cell output from v3 to v4

    - pyout -> execute_result
    - pyerr -> error
    - output.type -> output.data.mime/type
    - mime-type keys
    - stream.stream -> stream.name
    """
    if output["output_type"] in {"pyout", "display_data"}:
        output.setdefault("metadata", NotebookNode())
        if output["output_type"] == "pyout":
            output["output_type"] = "execute_result"
            output["execution_count"] = output.pop("prompt_number", None)

        # move output data into data sub-dict
        data = {}
        for key in list(output):
            if key in {"output_type", "execution_count", "metadata"}:
                continue
            data[key] = output.pop(key)
        to_mime_key(data)
        output["data"] = data
        to_mime_key(output.metadata)
        if "application/json" in data:
            data["application/json"] = json.loads(data["application/json"])
        # promote ascii bytes (from v2) to unicode
        for key in ("image/png", "image/jpeg"):
            if key in data and isinstance(data[key], bytes):
                data[key] = data[key].decode("ascii")
    elif output["output_type"] == "pyerr":
        output["output_type"] = "error"
    elif output["output_type"] == "stream":
        output["name"] = output.pop("stream", "stdout")
    return output


def downgrade_output(output):
    """downgrade a single code cell output to v3 from v4

    - pyout <- execute_result
    - pyerr <- error
    - output.data.mime/type -> output.type
    - un-mime-type keys
    - stream.stream <- stream.name
    """
    if output["output_type"] in {"execute_result", "display_data"}:
        if output["output_type"] == "execute_result":
            output["output_type"] = "pyout"
            output["prompt_number"] = output.pop("execution_count", None)

        # promote data dict to top-level output namespace
        data = output.pop("data", {})
        if "application/json" in data:
            data["application/json"] = json.dumps(data["application/json"])
        data = from_mime_key(data)
        output.update(data)
        from_mime_key(output.get("metadata", {}))
    elif output["output_type"] == "error":
        output["output_type"] = "pyerr"
    elif output["output_type"] == "stream":
        output["stream"] = output.pop("name")
    return output


def upgrade_outputs(outputs):
    """upgrade outputs of a code cell from v3 to v4"""
    return [upgrade_output(op) for op in outputs]


def downgrade_outputs(outputs):
    """downgrade outputs of a code cell to v3 from v4"""
    return [downgrade_output(op) for op in outputs]


def downgrade(nb):
    """Convert a v4 notebook to v3.

    Parameters
    ----------
    nb : NotebookNode
        The Python representation of the notebook to convert.
    """
    if nb.nbformat != nbformat:
        return nb

    # Validate the notebook before conversion
    _warn_if_invalid(nb, nbformat)

    nb.nbformat = v3.nbformat
    nb.nbformat_minor = v3.nbformat_minor
    cells = [downgrade_cell(cell) for cell in nb.pop("cells")]
    nb.worksheets = [v3.new_worksheet(cells=cells)]
    nb.metadata.setdefault("name", "")

    # Validate the converted notebook before returning it
    _warn_if_invalid(nb, v3.nbformat)

    nb.orig_nbformat = nb.metadata.pop("orig_nbformat", nbformat)
    nb.orig_nbformat_minor = nb.metadata.pop("orig_nbformat_minor", nbformat_minor)

    return nb
