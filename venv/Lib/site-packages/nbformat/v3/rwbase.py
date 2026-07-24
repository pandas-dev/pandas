"""Base classes and utilities for readers and writers."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from base64 import decodebytes, encodebytes


def restore_bytes(nb):
    """Restore bytes of image data from unicode-only formats.

    Base64 encoding is handled elsewhere.  Bytes objects in the notebook are
    always b64-encoded. We DO NOT encode/decode around file formats.

    Note: this is never used
    """
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == "code":
                for output in cell.outputs:
                    if "png" in output:
                        output.png = output.png.encode("ascii", "replace")
                    if "jpeg" in output:
                        output.jpeg = output.jpeg.encode("ascii", "replace")
    return nb


# output keys that are likely to have multiline values
_multiline_outputs = ["text", "html", "svg", "latex", "javascript", "json"]


# FIXME: workaround for old splitlines()
def _join_lines(lines):
    """join lines that have been written by splitlines()

    Has logic to protect against `splitlines()`, which
    should have been `splitlines(True)`
    """
    if lines and lines[0].endswith(("\n", "\r")):
        # created by splitlines(True)
        return "".join(lines)
    # created by splitlines()
    return "\n".join(lines)


def rejoin_lines(nb):
    """rejoin multiline text into strings

    For reversing effects of ``split_lines(nb)``.

    This only rejoins lines that have been split, so if text objects were not split
    they will pass through unchanged.

    Used when reading JSON files that may have been passed through split_lines.
    """
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == "code":
                if "input" in cell and isinstance(cell.input, list):
                    cell.input = _join_lines(cell.input)
                for output in cell.outputs:
                    for key in _multiline_outputs:
                        item = output.get(key, None)
                        if isinstance(item, list):
                            output[key] = _join_lines(item)
            else:  # text, heading cell
                for key in ["source", "rendered"]:
                    item = cell.get(key, None)
                    if isinstance(item, list):
                        cell[key] = _join_lines(item)
    return nb


def split_lines(nb):
    """split likely multiline text into lists of strings

    For file output more friendly to line-based VCS. ``rejoin_lines(nb)`` will
    reverse the effects of ``split_lines(nb)``.

    Used when writing JSON files.
    """
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == "code":
                if "input" in cell and isinstance(cell.input, str):
                    cell.input = cell.input.splitlines(True)
                for output in cell.outputs:
                    for key in _multiline_outputs:
                        item = output.get(key, None)
                        if isinstance(item, str):
                            output[key] = item.splitlines(True)
            else:  # text, heading cell
                for key in ["source", "rendered"]:
                    item = cell.get(key, None)
                    if isinstance(item, str):
                        cell[key] = item.splitlines(True)
    return nb


# b64 encode/decode are never actually used, because all bytes objects in
# the notebook are already b64-encoded, and we don't need/want to double-encode


def base64_decode(nb):
    """Restore all bytes objects in the notebook from base64-encoded strings.

    Note: This is never used
    """
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == "code":
                for output in cell.outputs:
                    if "png" in output:
                        if isinstance(output.png, str):
                            output.png = output.png.encode("ascii")
                        output.png = decodebytes(output.png)
                    if "jpeg" in output:
                        if isinstance(output.jpeg, str):
                            output.jpeg = output.jpeg.encode("ascii")
                        output.jpeg = decodebytes(output.jpeg)
    return nb


def base64_encode(nb):
    """Base64 encode all bytes objects in the notebook.

    These will be b64-encoded unicode strings

    Note: This is never used
    """
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type == "code":
                for output in cell.outputs:
                    if "png" in output:
                        output.png = encodebytes(output.png).decode("ascii")
                    if "jpeg" in output:
                        output.jpeg = encodebytes(output.jpeg).decode("ascii")
    return nb


def strip_transient(nb):
    """Strip transient values that shouldn't be stored in files.

    This should be called in *both* read and write.
    """
    nb.pop("orig_nbformat", None)
    nb.pop("orig_nbformat_minor", None)
    for ws in nb["worksheets"]:
        for cell in ws["cells"]:
            cell.get("metadata", {}).pop("trusted", None)
            # strip cell.trusted even though it shouldn't be used,
            # since it's where the transient value used to be stored.
            cell.pop("trusted", None)
    return nb


class NotebookReader:
    """A class for reading notebooks."""

    def reads(self, s, **kwargs):
        """Read a notebook from a string."""
        msg = "loads must be implemented in a subclass"
        raise NotImplementedError(msg)

    def read(self, fp, **kwargs):
        """Read a notebook from a file like object"""
        nbs = fp.read()
        return self.reads(nbs, **kwargs)


class NotebookWriter:
    """A class for writing notebooks."""

    def writes(self, nb, **kwargs):
        """Write a notebook to a string."""
        msg = "loads must be implemented in a subclass"
        raise NotImplementedError(msg)

    def write(self, nb, fp, **kwargs):
        """Write a notebook to a file like object"""
        nbs = self.writes(nb, **kwargs)
        return fp.write(nbs)
