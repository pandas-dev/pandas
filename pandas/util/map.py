import numpy as np
from pandas import _tseries as lib
from pandas import notnull, Series
from functools import wraps

class repeat(object):
    def __init__(self, obj):
        self.obj = obj
    
    def __getitem__(self, i):
        return self.obj

class azip(object):
    def __init__(self, *args):
        self.cols = []
        for a in args:
            if np.isscalar(a):
                self.cols.append(repeat(a))
            else:
                self.cols.append(a)
    
    def __getitem__(self, i):
        return [col[i] for col in self.cols]

def map_iter_args(arr, f, otherargs, n_otherargs, required, n_results):
    '''
    Substitute for np.vectorize with pandas-friendly dtype inference

    Parameters
    ----------
    arr : ndarray
    f : function

    Returns
    -------
    mapped : ndarray
    '''
    n = len(arr)
    result = np.empty((n, n_results), dtype=object)
    for i, val in enumerate(arr):
        args = otherargs[i]
        if notnull(val) and all(notnull(args[r]) for r in required):
            result[i] = f(val, *args)
        else:
            result[i] = [np.nan] * n_results

    return [lib.maybe_convert_objects(col, try_float=0) for col in result.T]
    
def auto_map(arr, f, otherargs, n_results=1, required='all'):
    if all(np.isscalar(a) for a in otherargs):
        res = lib.map_infer(arr, lambda v: f(v, *otherargs))
        return Series(res, index=arr.index, copy=False)
    
    n_otherargs = len(otherargs)
    if required == 'all':
        required = list(range(n_otherargs))
    res = map_iter_args(arr, f, azip(*otherargs), n_otherargs, required, n_results)
    res = [Series(col, index=arr.index, copy=False) for col in res]
    if n_results == 1:
        return res[0]
    return res

def mapwrap(f, n_results_default=1, required='all'):
    @wraps(f)
    def wrapped(arr, n_results=None, *otherargs):
        n_results = n_results or n_results_default
        return auto_map(arr, f, otherargs, n_results, required)
    
    return wrapped
