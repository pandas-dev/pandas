"""The main API for the v4 notebook format."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

__all__ = [
    "nbformat",
    "nbformat_minor",
    "nbformat_schema",
    "new_code_cell",
    "new_markdown_cell",
    "new_raw_cell",
    "new_notebook",
    "new_output",
    "output_from_msg",
    "reads",
    "writes",
    "to_notebook",
    "downgrade",
    "upgrade",
]

from .convert import downgrade, upgrade
from .nbbase import (
    nbformat,
    nbformat_minor,
    nbformat_schema,
    new_code_cell,
    new_markdown_cell,
    new_notebook,
    new_output,
    new_raw_cell,
    output_from_msg,
)
from .nbjson import reads, to_notebook, writes

reads_json = reads
writes_json = writes
to_notebook_json = to_notebook
