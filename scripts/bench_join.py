import numpy as np
import pandas._tseries as lib
from pandas import *
from copy import deepcopy

a = np.arange(100000, dtype=np.int64)
b = np.arange(20000, 120000, dtype=np.int64)

dr1 = DateRange('1/1/2000', periods=100000, offset=datetools.Minute())
dr2 = DateRange(dr1[20000], periods=100000, offset=datetools.Minute(2))

aobj = a.astype(object)
bobj = b.astype(object)

av = np.random.randn(100000)
bv = np.random.randn(100000)

K = 5

avf = np.random.randn(100000, K)
bvf = np.random.randn(100000, K)

a_series = Series(av, index=a)
b_series = Series(bv, index=b)

a_frame = DataFrame(avf, index=dr1, columns=range(K))
b_frame = DataFrame(bvf, index=dr2, columns=range(K, 2 * K))

def do_left_join(a, b, av, bv):
    indexer, mask = lib.ordered_left_join_int64(a, b)
    result = bv.take(indexer)
    np.putmask(result, mask, np.nan)
    return result

from line_profiler import LineProfiler
prof = LineProfiler()

from pandas.util.testing import set_trace

def do_left_join_multi(a, b, av, bv):
    indexer, mask = lib.ordered_left_join_int64(a, b)

    n, ak = av.shape
    _, bk = bv.shape
    result_width = ak + bk

    result = np.empty((result_width, n), dtype=np.float64)
    result[:ak] = av.T

    bchunk = result[ak:]
    _take_multi(bv.T, indexer, bchunk)
    np.putmask(bchunk, np.tile(mask, bk), np.nan)
    return result

def _take_multi(data, indexer, out):
    if not data.flags.c_contiguous:
        data = data.copy()
    for i in xrange(data.shape[0]):
        data[i].take(indexer, out=out[i])

def do_left_join_multi_put(a, b, av, bv):
    n, ak = av.shape
    _, bk = bv.shape
    result = np.empty((n, ak + bk), dtype=np.float64)
    lib.ordered_left_join_put(a, b, av, bv, result)
    return result

def do_left_join_multi_v2(a, b, av, bv):
    indexer, mask = lib.ordered_left_join_int64(a, b)
    bv_taken = bv.take(indexer, axis=0)
    np.putmask(bv_taken, mask.repeat(bv.shape[1]), np.nan)
    return np.concatenate((av, bv_taken), axis=1)

def do_left_join_series(a, b):
    return b.reindex(a.index)

def do_left_join_frame(a, b):
    a.index._indexMap = None
    b.index._indexMap = None
    return a.join(b, how='left')
