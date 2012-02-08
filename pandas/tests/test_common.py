from datetime import datetime

import unittest

from pandas import Series, DataFrame
from pandas.core.common import notnull, isnull
import pandas.core.common as com
import pandas.util.testing as tm

import numpy as np

def test_notnull():
    assert notnull(1.)
    assert not notnull(None)
    assert not notnull(np.NaN)
    assert not notnull(np.inf)
    assert not notnull(-np.inf)

    float_series = Series(np.random.randn(5))
    obj_series = Series(np.random.randn(5), dtype=object)
    assert(isinstance(notnull(float_series), Series))
    assert(isinstance(notnull(obj_series), Series))

def test_isnull():
    assert not isnull(1.)
    assert isnull(None)
    assert isnull(np.NaN)
    assert isnull(np.inf)
    assert isnull(-np.inf)

    float_series = Series(np.random.randn(5))
    obj_series = Series(np.random.randn(5), dtype=object)
    assert(isinstance(isnull(float_series), Series))
    assert(isinstance(isnull(obj_series), Series))

    # call on DataFrame
    df = DataFrame(np.random.randn(10, 5))
    df['foo'] = 'bar'
    result = isnull(df)
    expected = result.apply(isnull)
    tm.assert_frame_equal(result, expected)

def test_isnull_datetime():
    assert (not isnull(datetime.now()))
    assert notnull(datetime.now())

def test_any_none():
    assert(com._any_none(1, 2, 3, None))
    assert(not com._any_none(1, 2, 3, 4))

def test_all_not_none():
    assert(com._all_not_none(1, 2, 3, 4))
    assert(not com._all_not_none(1, 2, 3, None))
    assert(not com._all_not_none(None, None, None, None))

def test_rands():
    r = com.rands(10)
    assert(len(r) == 10)

def test_adjoin():
    data = [['a', 'b', 'c'],
            ['dd', 'ee', 'ff'],
            ['ggg', 'hhh', 'iii']]
    expected = 'a  dd  ggg\nb  ee  hhh\nc  ff  iii'

    adjoined = com.adjoin(2, *data)

    assert(adjoined == expected)

def test_iterpairs():
    data = [1, 2, 3, 4]
    expected = [(1, 2),
                (2, 3),
                (3, 4)]

    result = list(com.iterpairs(data))

    assert(result == expected)

def test_indent():
    s = 'a b c\nd e f'
    result = com.indent(s, spaces=6)

    assert(result == '      a b c\n      d e f')

def test_banner():
    ban = com.banner('hi')
    assert(ban == ('%s\nhi\n%s' % ('=' * 80, '=' * 80)))

def test_map_indices_py():
    data = [4, 3, 2, 1]
    expected = {4 : 0, 3 : 1, 2 : 2, 1 : 3}

    result = com.map_indices_py(data)

    assert(result == expected)

def test_union():
    a = [1, 2, 3]
    b = [4, 5, 6]

    union = sorted(com.union(a, b))

    assert((a + b) == union)

def test_difference():
    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5, 6]

    inter = sorted(com.difference(b, a))

    assert([4, 5, 6] == inter)

def test_intersection():
    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5, 6]

    inter = sorted(com.intersection(a, b))

    assert(a == inter)

def test_groupby():
    values = ['foo', 'bar', 'baz', 'baz2', 'qux', 'foo3']
    expected = {'f' : ['foo', 'foo3'],
                'b' : ['bar', 'baz', 'baz2'],
                'q' : ['qux']}

    grouped = com.groupby(values, lambda x: x[0])

    for k, v in grouped:
        assert v == expected[k]

def test_ensure_int32():
    values = np.arange(10, dtype=np.int32)
    result = com._ensure_int32(values)
    assert(result.dtype == np.int32)

    values = np.arange(10, dtype=np.int64)
    result = com._ensure_int32(values)
    assert(result.dtype == np.int32)


class TestTake(unittest.TestCase):

    def test_1d_with_out(self):
        def _test_dtype(dtype):
            out = np.empty(5, dtype=dtype)
            arr = np.random.randn(10).astype(dtype)
            indexer = [0, 2, 4, 7, 1]

            arr.take(indexer, out=out)
            expected = arr.take(indexer)
            tm.assert_almost_equal(out, expected)

        _test_dtype(np.float64)
        _test_dtype(np.float32)
        _test_dtype(np.int32)
        _test_dtype(np.int64)
        _test_dtype(np.object_)
        _test_dtype(np.bool)

    def test_1d_upcast_with_out(self):
        def _test_dtype(dtype):
            out = np.empty(4, dtype=dtype)
            data = np.random.randint(0, 2, 5).astype(dtype)

            indexer = [2, 1, 0, -1]
            self.assertRaises(Exception, com.take_1d, data,
                              indexer, out=out)

        _test_dtype(np.int64)
        _test_dtype(np.int32)
        _test_dtype(np.int16)
        _test_dtype(np.int8)
        _test_dtype(np.bool)

    def test_2d_upcast_with_out(self):
        def _test_dtype(dtype):
            out0 = np.empty((4, 3), dtype=dtype)
            out1 = np.empty((5, 4), dtype=dtype)

            data = np.random.randint(0, 2, (5, 3)).astype(dtype)

            indexer = [2, 1, 0, -1]
            self.assertRaises(Exception, com.take_2d, data,
                              indexer, out=out0, axis=0)
            self.assertRaises(Exception, com.take_2d, data,
                              indexer, out=out1, axis=1)

            # no exception o/w
            data.take(indexer, out=out0, axis=0)
            data.take(indexer, out=out1, axis=1)

        _test_dtype(np.int64)
        _test_dtype(np.int32)
        _test_dtype(np.int16)
        _test_dtype(np.int8)
        _test_dtype(np.bool)

    def test_1d_other_dtypes(self):
        arr = np.random.randn(10).astype(np.float32)

        indexer = [1, 2, 3, -1]
        result = com.take_1d(arr, indexer)
        expected = arr.take(indexer)
        expected[-1] = np.nan
        tm.assert_almost_equal(result, expected)

    def test_2d_other_dtypes(self):
        arr = np.random.randn(10, 5).astype(np.float32)

        indexer = [1, 2, 3, -1]

        # axis=0
        result = com.take_2d(arr, indexer, axis=0)
        expected = arr.take(indexer, axis=0)
        expected[-1] = np.nan
        tm.assert_almost_equal(result, expected)

        # axis=1
        result = com.take_2d(arr, indexer, axis=1)
        expected = arr.take(indexer, axis=1)
        expected[:, -1] = np.nan
        tm.assert_almost_equal(result, expected)

    def test_1d_bool(self):
        arr = np.array([0, 1, 0], dtype=bool)

        result = com.take_1d(arr, [0, 2, 2, 1])
        expected = arr.take([0, 2, 2, 1])
        self.assert_(np.array_equal(result, expected))

        result = com.take_1d(arr, [0, 2, -1])
        self.assert_(result.dtype == np.object_)

    def test_2d_bool(self):
        arr = np.array([[0, 1, 0],
                        [1, 0, 1],
                        [0, 1, 1]], dtype=bool)

        result = com.take_2d(arr, [0, 2, 2, 1])
        expected = arr.take([0, 2, 2, 1], axis=0)
        self.assert_(np.array_equal(result, expected))

        result = com.take_2d(arr, [0, 2, 2, 1], axis=1)
        expected = arr.take([0, 2, 2, 1], axis=1)
        self.assert_(np.array_equal(result, expected))

        result = com.take_2d(arr, [0, 2, -1])
        self.assert_(result.dtype == np.object_)

    def test_2d_float32(self):
        arr = np.random.randn(4, 3).astype(np.float32)
        indexer = [0, 2, -1, 1, -1]

        # axis=0
        result = com.take_2d(arr, indexer)
        result2 = np.empty_like(result)
        com.take_2d(arr, indexer, out=result2)
        tm.assert_almost_equal(result, result)

        expected = arr.take(indexer, axis=0)
        expected[[2, 4]] = np.nan
        tm.assert_almost_equal(result, expected)

        # test with float64 out buffer
        out = np.empty((len(indexer), arr.shape[1]), dtype='f8')
        com.take_2d(arr, indexer, out=out) # it works!

        # axis=1
        result = com.take_2d(arr, indexer, axis=1)
        result2 = np.empty_like(result)
        com.take_2d(arr, indexer, axis=1, out=result2)
        tm.assert_almost_equal(result, result)

        expected = arr.take(indexer, axis=1)
        expected[:, [2, 4]] = np.nan
        tm.assert_almost_equal(result, expected)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

