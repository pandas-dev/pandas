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

import json

from .nbbase import from_dict
from .rwbase import NotebookReader, NotebookWriter

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------


class JSONReader(NotebookReader):
    """A JSON notebook reader."""

    def reads(self, s, **kwargs):
        """Convert a string to a notebook object."""
        nb = json.loads(s, **kwargs)
        return self.to_notebook(nb, **kwargs)

    def to_notebook(self, d, **kwargs):
        """Convert from a raw JSON dict to a nested NotebookNode structure."""
        return from_dict(d)


class JSONWriter(NotebookWriter):
    """A JSON notebook writer."""

    def writes(self, nb, **kwargs):
        """Convert a notebook object to a string."""
        kwargs["indent"] = 4
        return json.dumps(nb, **kwargs)


_reader = JSONReader()
_writer = JSONWriter()

reads = _reader.reads
read = _reader.read
to_notebook = _reader.to_notebook
write = _writer.write
writes = _writer.writes
