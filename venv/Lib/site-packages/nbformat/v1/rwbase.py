"""Base classes and function for readers and writers.

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

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------
from __future__ import annotations


class NotebookReader:
    """The base notebook reader."""

    def reads(self, s, **kwargs):
        """Read a notebook from a string."""
        msg = "loads must be implemented in a subclass"
        raise NotImplementedError(msg)

    def read(self, fp, **kwargs):
        """Read a notebook from a file like object"""
        return self.reads(fp.read(), **kwargs)


class NotebookWriter:
    """The base notebook writer."""

    def writes(self, nb, **kwargs):
        """Write a notebook to a string."""
        msg = "loads must be implemented in a subclass"
        raise NotImplementedError(msg)

    def write(self, nb, fp, **kwargs):
        """Write a notebook to a file like object"""
        return fp.write(self.writes(nb, **kwargs))
