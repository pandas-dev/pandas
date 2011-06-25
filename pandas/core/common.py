"""
Misc tools for implementing data structures
"""

from cStringIO import StringIO

from numpy.lib.format import read_array, write_array
import numpy as np

import pandas.lib.tseries as tseries

# XXX: HACK for NumPy 1.5.1 to suppress warnings
try:
    np.seterr(all='ignore')
except Exception: # pragma: no cover
    pass

def isnull(input):
    '''
    Replacement for numpy.isnan / -numpy.isfinite which is suitable
    for use on object arrays.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    '''
    if isinstance(input, np.ndarray):
        if input.dtype.kind in ('O', 'S'):
            # Working around NumPy ticket 1542
            result = input.copy().astype(bool)
            result[:] = tseries.isnullobj(input)
        else:
            result = -np.isfinite(input)
    else:
        result = tseries.checknull(input)

    return result

def notnull(input):
    '''
    Replacement for numpy.isfinite / -numpy.isnan which is suitable
    for use on object arrays.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    '''
    if isinstance(input, np.ndarray):
        return -isnull(input)
    else:
        return not tseries.checknull(input)

def _pickle_array(arr):
    arr = arr.view(np.ndarray)

    buf = StringIO()
    write_array(buf, arr)

    return buf.getvalue()

def _unpickle_array(bytes):
    arr = read_array(StringIO(bytes))
    return arr

def null_out_axis(arr, mask, axis):
    if axis == 0:
        arr[mask] = np.NaN
    else:
        indexer = [slice(None)] * arr.ndim
        indexer[axis] = mask

        arr[tuple(indexer)] = np.NaN

#-------------------------------------------------------------------------------
# Lots of little utilities

def ensure_float(arr):
    if issubclass(arr.dtype.type, np.integer):
        arr = arr.astype(float)

    return arr

def _mut_exclusive(arg1, arg2):
    if arg1 is not None and arg2 is not None:
        raise Exception('mutually exclusive arguments')
    elif arg1 is not None:
        return arg1
    else:
        return arg2

def _ensure_index(index_like):
    from pandas.core.index import Index
    if not isinstance(index_like, Index):
        index_like = Index(index_like)

    return index_like

def _any_none(*args):
    for arg in args:
        if arg is None:
            return True
    return False

def _all_not_none(*args):
    for arg in args:
        if arg is None:
            return False
    return True

def _try_sort(iterable):
    listed = list(iterable)
    try:
        return sorted(listed)
    except Exception:
        return listed


def _pfixed(s, space, nanRep=None, float_format=None):
    if isinstance(s, float):
        if nanRep is not None and isnull(s):
            if np.isnan(s):
                s = nanRep
            return (' %s' % s).ljust(space)

        if float_format:
            formatted = float_format(s)
        else:
            is_neg = s < 0
            formatted = '%.4g' % np.abs(s)

            if is_neg:
                formatted = '-' + formatted
            else:
                formatted = ' ' + formatted

        return formatted.ljust(space)
    else:
        return (' %s' % s)[:space].ljust(space)

