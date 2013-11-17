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


_dtype_map = {'datetime64[ns]': 'int64', 'int64': 'int64',
              'float64': 'float64'}


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
        # just sort and take n
        return arr.argsort(kind='mergesort')[:n]

    try:
        dtype = _dtype_map[str(arr.dtype)]
    except KeyError:
        raise NotImplementedError("Not implemented for %s dtype, "
                                  "perhaps convert to int64 or float64, "
                                  "or use .order().head(n)") % arr.dtype

    arr = arr.view(dtype)

    if take_last:
        arr = arr[::-1]

    kth_val = algos.kth_smallest(arr.copy(), n - 1)

    ns, = np.nonzero(arr <= kth_val)
    inds = ns[arr[ns].argsort(kind='mergesort')][:n]

    if take_last:
        # reverse indices
        return len(arr) - 1 - inds
    return inds


def nlargest(arr, n=5, take_last=False):
    '''
    Find the indices of the n largest values of a numpy array.

    Note: Fails silently with NaN.

    '''
    if n <= 0:
        return np.array([])  # empty

    n = min(n, len(arr))

    if arr.dtype == object:
        try:
            arr = arr.astype(float)
        except:
            raise TypeError("An object array must convert to float.")

    arr = -arr.view(_dtype_map[str(arr.dtype)])
    return nsmallest(arr, n, take_last=take_last)
