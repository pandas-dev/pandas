from cStringIO import StringIO

from numpy.lib.format import read_array, write_array
import numpy as np

from pandas.lib.tseries import isnullobj, checknull

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
    from pandas.core.series import Series

    if isinstance(input, np.ndarray):
        if input.dtype.kind in ('O', 'S'):
            result = isnullobj(input).astype(bool)
            if isinstance(input, Series):
                result = Series(result, index=input.index)
        else:
            result = -np.isfinite(input)
    else:
        result = checknull(input)

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
        return not checknull(input)

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
            return nanRep.ljust(space)

        if float_format:
            formatted = float_format(s)
        else:
            formatted = '%.4g' % s

        return formatted.ljust(space)
    else:
        return str(s)[:space-4].ljust(space)

#-------------------------------------------------------------------------------
# Functions needed from scipy
