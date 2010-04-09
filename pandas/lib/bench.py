import time

import numpy as np

from pandas import Series, Index, isnull
import pandas.lib.tseries as tseries
from pandas.util.testing import assert_almost_equal, assert_dict_equal

def _timeit(f, n=10):
    _s = time.clock()
    for i in xrange(n):
        f()

    return (time.clock() - _s) / n

def bench_reindex():
    K = 100000
    index = Index(np.arange(K))
    values = np.arange(float(K))
    obj_vals = values.astype(object)

    new_index = np.arange(K)
    np.random.shuffle(new_index)
    new_index = Index(new_index)

    f = lambda: tseries.reindex(new_index, values, index.indexMap)
    print 'tseries.reindex: %.2f ms per iteration' % (_timeit(f, n=50) * 1000)

    def _test():
        filler, mask = tseries.getMergeVec(new_index, index.indexMap)
        result = values.take(filler)
        np.putmask(result, -mask, np.NaN)

        return result

    timing = _timeit(_test, n=50) * 1000
    print 'getMergeVec method: %.2f ms per iteration' % timing

    f2 = lambda: tseries.reindexObj(new_index, values, index.indexMap)
    print ('tseries.reindexObj with floats: %.2f ms per iteration'
           % (_timeit(f2, n=50) * 1000))

    f3 = lambda: tseries.reindexObj(new_index, obj_vals, index.indexMap)
    print ('tseries.reindexObj with objects: %.2f ms per iteration'
           % (_timeit(f3, n=50) * 1000))

    f4 = lambda: tseries.reindexObject(new_index, obj_vals, index.indexMap)
    print ('tseries.reindexObject buffers: %.2f ms per iteration'
           % (_timeit(f4, n=50) * 1000))

    def _test2():
        filler, mask = tseries.getMergeVec(new_index, index.indexMap)
        result = obj_vals.take(filler)
        np.putmask(result, -mask, np.NaN)

        return result

    timing = _timeit(_test2, n=50) * 1000
    print 'getMergeVec method: %.2f ms per iteration' % timing

    assert_almost_equal(_test(), f())
    assert_almost_equal(f2(), f3())
    assert_almost_equal(f3(), f4())
    assert_almost_equal(f2(), f4())
    assert_almost_equal(f2(), _test2())


def _isnan(obj):
    return obj != obj

def test_groupby():
    mapping = Series({
        1 : 2.,
        2 : 2.,
        3 : np.NaN,
        4 : np.NaN,
        5 : 3.,
        6 : 3.,
        7 : np.NaN
    })

    index = Index([1, 2, 3, 4, 5, 6, 7])

    expected = {
        2 : [1, 2],
        3 : [5, 6],
        np.NaN : [3, 4, 7]
    }

    def compare_with_null(d1, d2):
        d1_nulls = None
        d2_nulls = None
        for k, v in d1.iteritems():
            if _isnan(k):
                d1_nulls = v
            else:
                assert(k in d2)
                assert(np.array_equal(v, d2[k]))

        for k, v in d2.iteritems():
            if _isnan(k):
                d2_nulls = v
            else:
                assert(k in d1)

        if d1_nulls is not None or d2_nulls is not None:
            assert(np.array_equal(d1_nulls, d2_nulls))

    grouped = tseries.groupby(index, mapping.get)
    compare_with_null(grouped, expected)

def groupby_nocython(index, mapper, output=None):
    if output is None:
        result = {}
    else:
        result = output

    index = np.asarray(index)
    mapped_index = np.array([mapper(x) for x in index])

    # A little hack here
    if issubclass(mapped_index.dtype.type, basestring):
        mapped_index = mapped_index.astype(object)

    mask = isnull(mapped_index)
    nullkeys = index[mask]

    if nullkeys is not None and len(nullkeys) > 0:
        result[np.NaN] = nullkeys

    notmask = -mask
    index = index[notmask]
    mapped_index = mapped_index[notmask]

    for idx, key in zip(index, mapped_index):
        result.setdefault(key, []).append(idx)

    return result

def bench_groupby():
    N = 200

    arr = np.arange(10000).astype(object)
    values = np.random.randn(10000)
    keys = arr // 10
    d = dict(zip(arr, keys))

    f = lambda: groupby_nocython(arr, d.get)
    print 'no cython: %.2f ms per iteration' % (_timeit(f, n=N) * 1000)

    f = lambda: tseries.arrmap(arr, d.get)
    timing = _timeit(f, n=N) * 1000
    print 'arrmap: %.2f ms per iteration' % timing

    f = lambda: isnull(tseries.arrmap(arr, d.get))
    print 'isnull: %.2f ms per iteration' % (_timeit(f, n=N) * 1000 - timing)

    f = lambda: tseries.groupby(arr, d.get)
    print 'groupby: %.2f ms per iteration' % (_timeit(f, n=N) * 1000)

    f = lambda: tseries.groupby_indices(arr, d.get)
    print 'groupby_inds: %.2f ms per iteration' % (_timeit(f, n=N) * 1000)

    def _test():
        groups = tseries.groupby_indices(arr, d.get)

        result = {}
        for k, v in groups.iteritems():
            result[k] = np.mean(values.take(v))

        return result

    print 'test: %.2f ms per iteration' % (_timeit(_test, n=N) * 1000)

def bench_map_indices():
    pass
