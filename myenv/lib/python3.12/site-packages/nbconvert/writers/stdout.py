"""
Contains Stdout writer
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from nbconvert.utils import io

from .base import WriterBase


class StdoutWriter(WriterBase):
    """Consumes output from nbconvert export...() methods and writes to the
    stdout stream."""

    def write(self, output, resources, **kw):
        """
        Consume and write Jinja output.

        See base for more...
        """
        stream = io.unicode_std_stream()
        stream.write(output)
