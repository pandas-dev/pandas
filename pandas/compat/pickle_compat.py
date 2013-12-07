""" support pre 0.12 series pickle compatibility """

import sys
import numpy as np
import pandas
import copy
import pickle as pkl
from pandas import compat
from pandas.compat import u, string_types
from pandas.core.series import Series, TimeSeries
from pandas.sparse.series import SparseSeries, SparseTimeSeries


def load_reduce(self):
    stack = self.stack
    args = stack.pop()
    func = stack[-1]
    if type(args[0]) is type:
        n = args[0].__name__
        if n == u('DeprecatedSeries') or n == u('DeprecatedTimeSeries'):
            stack[-1] = object.__new__(Series)
            return
        elif (n == u('DeprecatedSparseSeries') or
              n == u('DeprecatedSparseTimeSeries')):
            stack[-1] = object.__new__(SparseSeries)
            return

    try:
        value = func(*args)
    except:

        # try to reencode the arguments
        if getattr(self,'encoding',None) is not None:
            args = tuple([arg.encode(self.encoding)
                          if isinstance(arg, string_types)
                          else arg for arg in args])
            try:
                stack[-1] = func(*args)
                return
            except:
                pass

        if getattr(self,'is_verbose',None):
            print(sys.exc_info())
            print(func, args)
        raise

    stack[-1] = value

if compat.PY3:
    class Unpickler(pkl._Unpickler):
        pass
else:
    class Unpickler(pkl.Unpickler):
        pass

Unpickler.dispatch = copy.copy(Unpickler.dispatch)
Unpickler.dispatch[pkl.REDUCE[0]] = load_reduce


def load(fh, encoding=None, compat=False, is_verbose=False):
    """load a pickle, with a provided encoding

    if compat is True:
       fake the old class hierarchy
       if it works, then return the new type objects

    Parameters
    ----------
    fh: a filelike object
    encoding: an optional encoding
    compat: provide Series compatibility mode, boolean, default False
    is_verbose: show exception output
    """

    try:
        if compat:
            pandas.core.series.Series = DeprecatedSeries
            pandas.core.series.TimeSeries = DeprecatedTimeSeries
            pandas.sparse.series.SparseSeries = DeprecatedSparseSeries
            pandas.sparse.series.SparseTimeSeries = DeprecatedSparseTimeSeries
        fh.seek(0)
        if encoding is not None:
            up = Unpickler(fh, encoding=encoding)
        else:
            up = Unpickler(fh)
        up.is_verbose = is_verbose

        return up.load()
    except:
        raise
    finally:
        if compat:
            pandas.core.series.Series = Series
            pandas.core.series.Series = TimeSeries
            pandas.sparse.series.SparseSeries = SparseSeries
            pandas.sparse.series.SparseTimeSeries = SparseTimeSeries


class DeprecatedSeries(np.ndarray, Series):
    pass


class DeprecatedTimeSeries(DeprecatedSeries):
    pass


class DeprecatedSparseSeries(DeprecatedSeries):
    pass


class DeprecatedSparseTimeSeries(DeprecatedSparseSeries):
    pass
