from pandas.compat import range, lrange
import numpy as np
import pandas.lib as lib
from pandas import *
from copy import deepcopy
import time

n = 1000000
K = 1
pct_overlap = 0.2

a = np.arange(n, dtype=np.int64)
b = np.arange(n * pct_overlap, n * (1 + pct_overlap), dtype=np.int64)

dr1 = DatetimeIndex('1/1/2000', periods=n, offset=offsets.Minute())
dr2 = DatetimeIndex(
    dr1[int(pct_overlap * n)], periods=n, offset=offsets.Minute(2))

aobj = a.astype(object)
bobj = b.astype(object)

av = np.random.randn(n)
bv = np.random.randn(n)

avf = np.random.randn(n, K)
bvf = np.random.randn(n, K)

a_series = Series(av, index=a)
b_series = Series(bv, index=b)

a_frame = DataFrame(avf, index=a, columns=lrange(K))
b_frame = DataFrame(bvf, index=b, columns=lrange(K, 2 * K))


def do_left_join(a, b, av, bv):
    out = np.empty((len(a), 2))
    lib.left_join_1d(a, b, av, bv, out)
    return out


def do_outer_join(a, b, av, bv):
    result_index, aindexer, bindexer = lib.outer_join_indexer(a, b)
    result = np.empty((2, len(result_index)))
    lib.take_1d(av, aindexer, result[0])
    lib.take_1d(bv, bindexer, result[1])
    return result_index, result


def do_inner_join(a, b, av, bv):
    result_index, aindexer, bindexer = lib.inner_join_indexer(a, b)
    result = np.empty((2, len(result_index)))
    lib.take_1d(av, aindexer, result[0])
    lib.take_1d(bv, bindexer, result[1])
    return result_index, result

from line_profiler import LineProfiler
prof = LineProfiler()

from pandas.util.testing import set_trace


def do_left_join_python(a, b, av, bv):
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
    for i in range(data.shape[0]):
        data[i].take(indexer, out=out[i])


def do_left_join_multi(a, b, av, bv):
    n, ak = av.shape
    _, bk = bv.shape
    result = np.empty((n, ak + bk), dtype=np.float64)
    lib.left_join_2d(a, b, av, bv, result)
    return result


def do_outer_join_multi(a, b, av, bv):
    n, ak = av.shape
    _, bk = bv.shape
    result_index, rindexer, lindexer = lib.outer_join_indexer(a, b)
    result = np.empty((len(result_index), ak + bk), dtype=np.float64)
    lib.take_join_contiguous(av, bv, lindexer, rindexer, result)
    # result = np.empty((ak + bk, len(result_index)), dtype=np.float64)
    # lib.take_axis0(av, rindexer, out=result[:ak].T)
    # lib.take_axis0(bv, lindexer, out=result[ak:].T)
    return result_index, result


def do_inner_join_multi(a, b, av, bv):
    n, ak = av.shape
    _, bk = bv.shape
    result_index, rindexer, lindexer = lib.inner_join_indexer(a, b)
    result = np.empty((len(result_index), ak + bk), dtype=np.float64)
    lib.take_join_contiguous(av, bv, lindexer, rindexer, result)
    # result = np.empty((ak + bk, len(result_index)), dtype=np.float64)
    # lib.take_axis0(av, rindexer, out=result[:ak].T)
    # lib.take_axis0(bv, lindexer, out=result[ak:].T)
    return result_index, result


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


# a = np.array([1, 2, 3, 4, 5], dtype=np.int64)
# b = np.array([0, 3, 5, 7, 9], dtype=np.int64)
# print(lib.inner_join_indexer(a, b))

out = np.empty((10, 120000))


def join(a, b, av, bv, how="left"):
    func_dict = {'left': do_left_join_multi,
                 'outer': do_outer_join_multi,
                 'inner': do_inner_join_multi}

    f = func_dict[how]
    return f(a, b, av, bv)


def bench_python(n=100000, pct_overlap=0.20, K=1):
    import gc
    ns = [2, 3, 4, 5, 6]
    iterations = 200
    pct_overlap = 0.2
    kinds = ['outer', 'left', 'inner']

    all_results = {}
    for logn in ns:
        n = 10 ** logn
        a = np.arange(n, dtype=np.int64)
        b = np.arange(n * pct_overlap, n * pct_overlap + n, dtype=np.int64)

        avf = np.random.randn(n, K)
        bvf = np.random.randn(n, K)

        a_frame = DataFrame(avf, index=a, columns=lrange(K))
        b_frame = DataFrame(bvf, index=b, columns=lrange(K, 2 * K))

        all_results[logn] = result = {}

        for kind in kinds:
            gc.disable()
            elapsed = 0
            _s = time.clock()
            for i in range(iterations):
                if i % 10 == 0:
                    elapsed += time.clock() - _s
                    gc.collect()
                    _s = time.clock()
                a_frame.join(b_frame, how=kind)
                # join(a, b, avf, bvf, how=kind)
            elapsed += time.clock() - _s
            gc.enable()
            result[kind] = (elapsed / iterations) * 1000

    return DataFrame(all_results, index=kinds)


def bench_xts(n=100000, pct_overlap=0.20):
    from pandas.rpy.common import r
    r('a <- 5')

    xrng = '1:%d' % n

    start = n * pct_overlap + 1
    end = n + start - 1
    yrng = '%d:%d' % (start, end)

    r('library(xts)')

    iterations = 500

    kinds = ['left', 'outer', 'inner']
    result = {}
    for kind in kinds:
        r('x <- xts(rnorm(%d), as.POSIXct(Sys.Date()) + %s)' % (n, xrng))
        r('y <- xts(rnorm(%d), as.POSIXct(Sys.Date()) + %s)' % (n, yrng))
        stmt = 'for (i in 1:%d) merge(x, y, join="%s")' % (iterations, kind)
        elapsed = r('as.list(system.time(%s, gcFirst=F))$elapsed' % stmt)[0]
        result[kind] = (elapsed / iterations) * 1000
    return Series(result)
