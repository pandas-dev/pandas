"""
Basic post processor
"""
# -----------------------------------------------------------------------------
# Copyright (c) 2013, the IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------

from nbconvert.utils.base import NbConvertBase


# -----------------------------------------------------------------------------
# Classes
# -----------------------------------------------------------------------------
class PostProcessorBase(NbConvertBase):
    """The base class for post processors."""

    def __call__(self, input_):
        """
        See def postprocess() ...
        """
        self.postprocess(input_)

    def postprocess(self, input_):
        """
        Post-process output from a writer.
        """
        msg = "postprocess"
        raise NotImplementedError(msg)
