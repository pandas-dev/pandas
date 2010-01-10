# pylint: disable-msg=W0402

from datetime import datetime
import random
import string

from numpy.random import randn
import numpy as np

from pandas.core.api import (DateRange, Index, Series, DataFrame,
                             DataMatrix, WidePanel)

N = 30
K = 4

def rands(n):
    choices = string.letters + string.digits
    return ''.join([random.choice(choices) for _ in xrange(n)])

def equalContents(arr1, arr2):
    """Checks if the set of unique elements of arr1 and arr2 are equivalent.
    """
    return frozenset(arr1) == frozenset(arr2)

def isiterable(obj):
    return hasattr(obj, '__iter__')

def assert_almost_equal(a, b):
    if isiterable(a):
        np.testing.assert_(isiterable(b))
        np.testing.assert_equal(len(a), len(b))
        for i in xrange(len(a)):
            assert_almost_equal(a[i], b[i])
        return True

    err_msg = lambda a, b: 'expected %.5f but got %.5f' % (a, b)

    if np.isnan(a):
        np.testing.assert_(np.isnan(b))
        return

    # case for zero
    if abs(a) < 1e-5:
        np.testing.assert_almost_equal(
            a, b, decimal=5, err_msg=err_msg(a, b), verbose=False)
    else:
        np.testing.assert_almost_equal(
            1, a/b, decimal=5, err_msg=err_msg(a, b), verbose=False)

def is_sorted(seq):
    return assert_almost_equal(seq, np.sort(np.array(seq)))

def assert_dict_equal(a, b, compare_keys=True):
    a_keys = frozenset(a.keys())
    b_keys = frozenset(b.keys())

    if compare_keys:
        assert(a_keys == b_keys)

    for k in a_keys:
        assert_almost_equal(a[k], b[k])

def assert_series_equal(left, right):
    assert_almost_equal(left, right)
    assert(np.array_equal(left.index, right.index))

def assert_frame_equal(left, right):
    for col, series in left.iteritems():
        assert(col in right)
        assert_series_equal(series, right[col])

    for col in right:
        assert(col in left)

def assert_contains_all(iterable, dic):
    for k in iterable:
        assert(k in dic)

def getCols(k):
    return string.ascii_uppercase[:k]

def makeStringIndex(k):
    return Index([rands(10) for _ in xrange(k)])

def makeIntIndex(k):
    return Index(np.arange(k))

def makeDateIndex(k):
    dates = list(DateRange(datetime(2000, 1, 1), periods=k))
    return Index(dates)

def makeFloatSeries():
    index = makeStringIndex(N)
    return Series(randn(N), index=index)

def makeStringSeries():
    index = makeStringIndex(N)
    return Series(randn(N), index=index)

def makeObjectSeries():
    dateIndex = makeDateIndex(N)
    index = makeStringIndex(N)
    return Series(dateIndex, index=index)

def makeTimeSeries():
    return Series(randn(N), index=makeDateIndex(N))

def getArangeMat():
    return np.arange(N * K).reshape((N, K))

def getSeriesData():
    index = makeStringIndex(N)

    return dict((c, Series(randn(N), index=index)) for c in getCols(K))

def getTimeSeriesData():
    return dict((c, makeTimeSeries()) for c in getCols(K))

def getMixedTypeDict():
    index = Index(['a', 'b', 'c', 'd', 'e'])

    data = {
        'A' : [0., 1., 2., 3., 4.],
        'B' : [0., 1., 0., 1., 0.],
        'C' : ['foo1', 'foo2', 'foo3', 'foo4', 'foo5'],
        'D' : DateRange('1/1/2009', periods=5)
    }

    return index, data

def makeDataFrame():
    data = getSeriesData()
    return DataFrame(data)

def makeTimeDataFrame():
    data = getTimeSeriesData()
    return DataFrame(data)

def makeDataMatrix():
    data = getSeriesData()
    return DataMatrix(data)

def makeTimeDataMatrix():
    data = getTimeSeriesData()
    return DataMatrix(data)

def makeWidePanel():
    cols = ['Item' + c for c in string.ascii_uppercase[:K - 1]]
    data = dict((c, makeTimeDataMatrix()) for c in cols)
    return WidePanel.fromDict(data)

def add_nans(panel):
    I, J, N = panel.dims
    for i, item in enumerate(panel.items):
        dm = panel[item]
        for j, col in enumerate(dm.columns):
            dm[col][:i + j] = np.NaN

def makeLongPanel():
    wp = makeWidePanel()
    add_nans(wp)

    return wp.toLong()
