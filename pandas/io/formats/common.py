# -*- coding: utf-8 -*-
"""
Common helper methods used in different submodules of pandas.io.formats
"""

from pandas.core.dtypes.common import is_integer
from pandas import compat

def get_level_lengths(levels, sentinel=''):
    """For each index in each level the function returns lengths of indexes.

    Parameters
    ----------
    levels : list of lists
        List of values on for level.
    sentinel : string, optional
        Value which states that no new index starts on there.

    Returns
    ----------
    Returns list of maps. For each level returns map of indexes (key is index
    in row and value is length of index).
    """
    if len(levels) == 0:
        return []

    control = [True for x in levels[0]]

    result = []
    for level in levels:
        last_index = 0

        lengths = {}
        for i, key in enumerate(level):
            if control[i] and key == sentinel:
                pass
            else:
                control[i] = False
                lengths[last_index] = i - last_index
                last_index = i

        lengths[last_index] = len(level) - last_index

        result.append(lengths)

    return result

def _put_lines(buf, lines):
    if any(isinstance(x, compat.text_type) for x in lines):
        lines = [compat.text_type(x) for x in lines]
    buf.write('\n'.join(lines))


def _binify(cols, line_width):
    adjoin_width = 1
    bins = []
    curr_width = 0
    i_last_column = len(cols) - 1
    for i, w in enumerate(cols):
        w_adjoined = w + adjoin_width
        curr_width += w_adjoined
        if i_last_column == i:
            wrap = curr_width + 1 > line_width and i > 0
        else:
            wrap = curr_width + 2 > line_width and i > 0
        if wrap:
            bins.append(i)
            curr_width = w_adjoined

    bins.append(len(cols))
    return bins

class TableFormatter(object):
    is_truncated = False
    show_dimensions = None

    @property
    def should_show_dimensions(self):
        return (self.show_dimensions is True or
                (self.show_dimensions == 'truncate' and self.is_truncated))

    def _get_formatter(self, i):
        if isinstance(self.formatters, (list, tuple)):
            if is_integer(i):
                return self.formatters[i]
            else:
                return None
        else:
            if is_integer(i) and i not in self.columns:
                i = self.columns[i]
            return self.formatters.get(i, None)
