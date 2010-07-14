from cStringIO import StringIO

from numpy.lib.format import read_array, write_array
import numpy as np

import pandas.lib.tseries as tseries

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

def _pfixed(s, space, nanRep=None, float_format=None):
    if isinstance(s, float):
        if nanRep is not None and isnull(s):
            if np.isnan(s):
                return nanRep.ljust(space)
            else:
                return str(s).ljust(space)

        if float_format:
            formatted = float_format(s)
        else:
            formatted = '%.4g' % s

        return formatted.ljust(space)
    else:
        return str(s)[:space-4].ljust(space)



def get_indexer(source, target, fill_method):
    if fill_method:
        fill_method = fill_method.upper()

    indexer, mask = tseries.getFillVec(source, target, source.indexMap,
                                       target.indexMap, fill_method)

    return indexer, mask

def null_out_axis(arr, mask, axis):
    if axis == 0:
        arr[mask] = np.NaN
    else:
        indexer = [slice(None)] * arr.ndim
        indexer[axis] = mask

        arr[tuple(indexer)] = np.NaN
