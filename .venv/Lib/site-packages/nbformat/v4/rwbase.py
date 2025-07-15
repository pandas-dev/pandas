"""Base classes and utilities for readers and writers."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations


def _is_json_mime(mime):
    """Is a key a JSON mime-type that should be left alone?"""
    return mime == "application/json" or (
        mime.startswith("application/") and mime.endswith("+json")
    )


def _rejoin_mimebundle(data):
    """Rejoin the multi-line string fields in a mimebundle (in-place)"""
    for key, value in list(data.items()):
        if (
            not _is_json_mime(key)
            and isinstance(value, list)
            and all(isinstance(line, str) for line in value)
        ):
            data[key] = "".join(value)
    return data


def rejoin_lines(nb):
    """rejoin multiline text into strings

    For reversing effects of ``split_lines(nb)``.

    This only rejoins lines that have been split, so if text objects were not split
    they will pass through unchanged.

    Used when reading JSON files that may have been passed through split_lines.
    """
    for cell in nb.cells:
        if "source" in cell and isinstance(cell.source, list):
            cell.source = "".join(cell.source)

        attachments = cell.get("attachments", {})
        for _, attachment in attachments.items():
            _rejoin_mimebundle(attachment)

        if cell.get("cell_type", None) == "code":
            for output in cell.get("outputs", []):
                output_type = output.get("output_type", "")
                if output_type in {"execute_result", "display_data"}:
                    _rejoin_mimebundle(output.get("data", {}))
                elif output_type and isinstance(output.get("text", ""), list):
                    output.text = "".join(output.text)
    return nb


_non_text_split_mimes = {
    "application/javascript",
    "image/svg+xml",
}


def _split_mimebundle(data):
    """Split multi-line string fields in a mimebundle (in-place)"""
    for key, value in list(data.items()):
        if isinstance(value, str) and (key.startswith("text/") or key in _non_text_split_mimes):
            data[key] = value.splitlines(True)
    return data


def split_lines(nb):
    """split likely multiline text into lists of strings

    For file output more friendly to line-based VCS. ``rejoin_lines(nb)`` will
    reverse the effects of ``split_lines(nb)``.

    Used when writing JSON files.
    """
    for cell in nb.cells:
        source = cell.get("source", None)
        if isinstance(source, str):
            cell["source"] = source.splitlines(True)

        attachments = cell.get("attachments", {})
        for _, attachment in attachments.items():
            _split_mimebundle(attachment)

        if cell.cell_type == "code":
            for output in cell.outputs:
                if output.output_type in {"execute_result", "display_data"}:
                    _split_mimebundle(output.get("data", {}))
                elif output.output_type == "stream" and isinstance(output.text, str):
                    output.text = output.text.splitlines(True)
    return nb


def strip_transient(nb):
    """Strip transient values that shouldn't be stored in files.

    This should be called in *both* read and write.
    """
    nb.metadata.pop("orig_nbformat", None)
    nb.metadata.pop("orig_nbformat_minor", None)
    nb.metadata.pop("signature", None)
    for cell in nb.cells:
        cell.metadata.pop("trusted", None)
    return nb


class NotebookReader:
    """A class for reading notebooks."""

    def reads(self, s, **kwargs):
        """Read a notebook from a string."""
        msg = "reads must be implemented in a subclass"
        raise NotImplementedError(msg)

    def read(self, fp, **kwargs):
        """Read a notebook from a file like object"""
        nbs = fp.read()
        return self.reads(nbs, **kwargs)


class NotebookWriter:
    """A class for writing notebooks."""

    def writes(self, nb, **kwargs):
        """Write a notebook to a string."""
        msg = "writes must be implemented in a subclass"
        raise NotImplementedError(msg)

    def write(self, nb, fp, **kwargs):
        """Write a notebook to a file like object"""
        nbs = self.writes(nb, **kwargs)
        return fp.write(nbs)
