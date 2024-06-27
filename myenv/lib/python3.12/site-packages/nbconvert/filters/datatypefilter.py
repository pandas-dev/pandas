"""Filter used to select the first preferred output format available.

The filter contained in the file allows the converter templates to select
the output format that is most valuable to the active export format.  The
value of the different formats is set via
NbConvertBase.display_data_priority
"""
# -----------------------------------------------------------------------------
# Copyright (c) 2013, the IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Classes and functions
# -----------------------------------------------------------------------------

from warnings import warn

from nbconvert.utils.base import NbConvertBase

__all__ = ["DataTypeFilter"]


class DataTypeFilter(NbConvertBase):
    """Returns the preferred display format"""

    def __call__(self, output):
        """Return the first available format in the priority.

        Produces a UserWarning if no compatible mimetype is found.

        `output` is dict with structure {mimetype-of-element: value-of-element}

        """
        for fmt in self.display_data_priority:
            if fmt in output:
                return [fmt]
        warn(
            f"Your element with mimetype(s) {output.keys()} is not able to be represented.",
            stacklevel=2,
        )

        return []
