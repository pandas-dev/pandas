import collections.abc
import itertools

import numpy as np

import matplotlib.cbook as cbook
import matplotlib.colors as mcolors
import matplotlib.lines as mlines


def check_non_empty(key, value):
    """Raise a TypeError if an empty sequence is passed"""
    if (not cbook.is_scalar_or_string(value) and
            isinstance(value, collections.abc.Sized) and len(value) == 0):
        raise TypeError(f'{key} must not be an empty sequence')


def style_generator(kw):
    """
    Helper for handling style sequences (e.g. facecolor=['r', 'b', 'k']) within plotting
    methods that repeatedly call other plotting methods (e.g. hist, stackplot).  Remove
    style keywords from the given dictionary.  Return the reduced dictionary together
    with a generator which provides a series of dictionaries to be used in each call to
    the wrapped function.
    """
    kw_iterators = {}
    remaining_kw = {}
    for key, value in kw.items():
        if key in ['facecolor', 'edgecolor']:
            if value is None or cbook._str_lower_equal(value, 'none'):
                kw_iterators[key] = itertools.repeat(value)
            else:
                check_non_empty(key, value)
                kw_iterators[key] = itertools.cycle(mcolors.to_rgba_array(value))

        elif key in ['hatch', 'linewidth']:
            check_non_empty(key, value)
            kw_iterators[key] = itertools.cycle(np.atleast_1d(value))

        elif key == 'linestyle':
            check_non_empty(key, value)
            kw_iterators[key] = itertools.cycle(mlines._get_dash_patterns(value))

        else:
            remaining_kw[key] = value

    def style_gen():
        while True:
            yield {key: next(val) for key, val in kw_iterators.items()}

    return remaining_kw, style_gen()
