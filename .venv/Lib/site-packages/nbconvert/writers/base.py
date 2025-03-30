"""
Contains writer base class.
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from traitlets import List, Unicode

from nbconvert.utils.base import NbConvertBase


class WriterBase(NbConvertBase):
    """Consumes output from nbconvert export...() methods and writes to a
    useful location."""

    files = List(
        Unicode(),
        help="""
        List of the files that the notebook references.  Files will be
        included with written output.""",
    ).tag(config=True)

    def __init__(self, config=None, **kw):
        """
        Constructor
        """
        super().__init__(config=config, **kw)

    def write(self, output, resources, **kw):
        """
        Consume and write Jinja output.

        Parameters
        ----------
        output : string
            Conversion results.  This string contains the file contents of the
            converted file.
        resources : dict
            Resources created and filled by the nbconvert conversion process.
            Includes output from preprocessors, such as the extract figure
            preprocessor.
        """

        raise NotImplementedError()
