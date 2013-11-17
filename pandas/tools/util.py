from pandas.compat import reduce
from pandas.core.index import Index
import numpy as np
from pandas import algos
import pandas.core.common as com


def match(needles, haystack):
    haystack = Index(haystack)
    needles = Index(needles)
    return haystack.get_indexer(needles)


def cartesian_product(X):
    '''
    Numpy version of itertools.product or pandas.compat.product.
    Sometimes faster (for large inputs)...

    Examples
    --------
    >>> cartesian_product([list('ABC'), [1, 2]])
    [array(['A', 'A', 'B', 'B', 'C', 'C'], dtype='|S1'),
    array([1, 2, 1, 2, 1, 2])]

    '''

    lenX = np.fromiter((len(x) for x in X), dtype=int)
    cumprodX = np.cumproduct(lenX)

    a = np.roll(cumprodX, 1)
    a[0] = 1

    b = cumprodX[-1] / cumprodX

    return [np.tile(np.repeat(np.asarray(x), b[i]),
                    np.product(a[i]))
               for i, x in enumerate(X)]


def _compose2(f, g):
    """Compose 2 callables"""
    return lambda *args, **kwargs: f(g(*args, **kwargs))


def compose(*funcs):
    """Compose 2 or more callables"""
    assert len(funcs) > 1, 'At least 2 callables must be passed to compose'
    return reduce(_compose2, funcs)


def nsmallest(arr, n=5, take_last=False):
    '''
    Find the indices of the n smallest values of a numpy array.

    Note: Fails silently with NaN.

    '''
    if n <= 0:
        return np.array([])  # empty
    elif n >= len(arr):
        n = len(arr)

    if arr.dtype == object:
        try:
            arr = arr.astype(float)
        except:
            raise TypeError("An object array must convert to float.")

    if com.needs_i8_conversion(arr):
        dtype = 'i8'
        kth_s = algos.kth_smallest_int64
    elif arr.dtype in ['int64']:
        dtype = 'int64'
        kth_s = algos.kth_smallest_int64
    elif arr.dtype in ['float64']:
        dtype = 'float64'
        kth_s = algos.kth_smallest_float64
    else:
        raise NotImplementedError("Not implemented for %s dtype, "
                                  "perhaps convert to int64 or float64, "
                                  "or use .order().head(n)") % arr.dtype

    if take_last:
        arr = arr.view(dtype)[::-1]
    else:
        arr = arr.view(dtype)

    kth_val = kth_s(arr.copy(), n - 1)

    ns = np.nonzero(arr <= kth_val)[0]
    inds = ns[arr[ns].argsort(kind='mergesort')][:n]

    if take_last:
        # reverse indices
        return len(arr) - 1 - inds
    else:
        return inds


def nlargest(arr, n=5, take_last=False):
    '''
    Find the indices of the n largest values of a numpy array.

    Note: Fails silently with NaN.

    '''
    if n <= 0:
        return np.array([])  # empty
    elif n >= len(arr):
        n = len(arr)

    if arr.dtype == object:
        try:
            arr = arr.astype(float)
        except:
            raise TypeError("An object array must convert to float.")

    if com.needs_i8_conversion(arr):
        arr = -arr.view('i8')
    else:
        arr = -arr

    return nsmallest(arr, n, take_last=take_last)
