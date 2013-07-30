""" support pre 0.12 series pickle compatibility """

import sys
import pickle
import numpy as np
import pandas
from pandas import compat
from pandas.core.series import Series
from pandas.sparse.series import SparseSeries

def load_reduce(self):
    stack = self.stack
    args = stack.pop()
    func = stack[-1]
    if type(args[0]) is type:
        n = args[0].__name__
        if n == 'DeprecatedSeries':
            stack[-1] = object.__new__(Series)
            return
        elif n == 'DeprecatedSparseSeries':
            stack[-1] = object.__new__(SparseSeries)
            return

    try:
        value = func(*args)
    except:
        print(sys.exc_info())
        print(func, args)
        raise

    stack[-1] = value

if compat.PY3:
    class Unpickler(pickle._Unpickler):
        pass
else:
    class Unpickler(pickle.Unpickler):
        pass

Unpickler.dispatch[pickle.REDUCE[0]] = load_reduce

def load(file):
    # try to load a compatibility pickle
    # fake the old class hierarchy
    # if it works, then return the new type objects

    try:
        pandas.core.series.Series = DeprecatedSeries
        pandas.sparse.series.SparseSeries = DeprecatedSparseSeries
        with open(file,'rb') as fh:
            return Unpickler(fh).load()
    except:
        raise
    finally:
        pandas.core.series.Series = Series
        pandas.sparse.series.SparseSeries = SparseSeries

class DeprecatedSeries(Series, np.ndarray):
    pass

class DeprecatedSparseSeries(DeprecatedSeries):
    pass
