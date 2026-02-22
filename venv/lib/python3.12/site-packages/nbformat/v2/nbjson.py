"""Read and write notebooks in JSON format.

Authors:

* Brian Granger
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

import copy
import json

from .nbbase import from_dict
from .rwbase import NotebookReader, NotebookWriter, rejoin_lines, restore_bytes, split_lines

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------


class BytesEncoder(json.JSONEncoder):
    """A JSON encoder that accepts b64 (and other *ascii*) bytestrings."""

    def default(self, obj):
        """The default value of an object."""
        if isinstance(obj, bytes):
            return obj.decode("ascii")
        return json.JSONEncoder.default(self, obj)


class JSONReader(NotebookReader):
    """A JSON notebook reader."""

    def reads(self, s, **kwargs):
        """Convert a string to a notebook."""
        nb = json.loads(s, **kwargs)
        nb = self.to_notebook(nb, **kwargs)
        return nb  # noqa: RET504

    def to_notebook(self, d, **kwargs):
        """Convert a string to a notebook."""
        return restore_bytes(rejoin_lines(from_dict(d)))


class JSONWriter(NotebookWriter):
    """A JSON notebook writer."""

    def writes(self, nb, **kwargs):
        """Convert a notebook object to a string."""
        kwargs["cls"] = BytesEncoder
        kwargs["indent"] = 1
        kwargs["sort_keys"] = True
        if kwargs.pop("split_lines", True):
            nb = split_lines(copy.deepcopy(nb))
        return json.dumps(nb, **kwargs)


_reader = JSONReader()
_writer = JSONWriter()

reads = _reader.reads
read = _reader.read
to_notebook = _reader.to_notebook
write = _writer.write
writes = _writer.writes
