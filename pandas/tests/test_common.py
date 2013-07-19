from datetime import datetime
import sys
import re

import nose
import unittest

from pandas import Series, DataFrame, date_range, DatetimeIndex, Timestamp
from pandas.core.common import notnull, isnull
import pandas.core.common as com
import pandas.util.testing as tm
import pandas.core.config as cf

import numpy as np

from pandas.tslib import iNaT
from pandas.util import py3compat

_multiprocess_can_split_ = True


def test_is_sequence():
    is_seq = com._is_sequence
    assert(is_seq((1, 2)))
    assert(is_seq([1, 2]))
    assert(not is_seq("abcd"))
    assert(not is_seq(u"abcd"))
    assert(not is_seq(np.int64))

    class A(object):
        def __getitem__(self):
            return 1

    assert(not is_seq(A()))

def test_notnull():
    assert notnull(1.)
    assert not notnull(None)
    assert not notnull(np.NaN)

    with cf.option_context("mode.use_inf_as_null", False):
        assert notnull(np.inf)
        assert notnull(-np.inf)

        arr = np.array([1.5, np.inf, 3.5, -np.inf])
        result = notnull(arr)
        assert result.all()

    with cf.option_context("mode.use_inf_as_null", True):
        assert not notnull(np.inf)
        assert not notnull(-np.inf)

        arr = np.array([1.5, np.inf, 3.5, -np.inf])
        result = notnull(arr)
        assert result.sum() == 2

    with cf.option_context("mode.use_inf_as_null", False):
        float_series = Series(np.random.randn(5))
        obj_series = Series(np.random.randn(5), dtype=object)
        assert(isinstance(notnull(float_series), Series))
        assert(isinstance(notnull(obj_series), Series))


def test_isnull():
    assert not isnull(1.)
    assert isnull(None)
    assert isnull(np.NaN)
    assert not isnull(np.inf)
    assert not isnull(-np.inf)

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


def test_isnull_lists():
    result = isnull([[False]])
    exp = np.array([[False]])
    assert(np.array_equal(result, exp))

    result = isnull([[1], [2]])
    exp = np.array([[False], [False]])
    assert(np.array_equal(result, exp))

    # list of strings / unicode
    result = isnull(['foo', 'bar'])
    assert(not result.any())

    result = isnull([u'foo', u'bar'])
    assert(not result.any())


def test_isnull_datetime():
    assert (not isnull(datetime.now()))
    assert notnull(datetime.now())

    idx = date_range('1/1/1990', periods=20)
    assert(notnull(idx).all())

    idx = np.asarray(idx)
    idx[0] = iNaT
    idx = DatetimeIndex(idx)
    mask = isnull(idx)
    assert(mask[0])
    assert(not mask[1:].any())

def test_datetimeindex_from_empty_datetime64_array():
    for unit in [ 'ms', 'us', 'ns' ]:
        idx = DatetimeIndex(np.array([], dtype='datetime64[%s]' % unit))
        assert(len(idx) == 0)

def test_nan_to_nat_conversions():

    df = DataFrame(dict({
        'A' : np.asarray(range(10),dtype='float64'),
        'B' : Timestamp('20010101') }))
    df.iloc[3:6,:] = np.nan
    result = df.loc[4,'B'].value
    assert(result == iNaT)

    values = df['B'].values
    result, changed = com._maybe_upcast_indexer(values,tuple([slice(8,9)]),np.nan)
    assert(isnull(result[8]))

    # numpy < 1.7.0 is wrong
    from distutils.version import LooseVersion
    if LooseVersion(np.__version__) >= '1.7.0':
        assert(result[8] == np.datetime64('NaT'))

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


def test_split_ranges():
    def _bin(x, width):
        "return int(x) as a base2 string of given width"
        return ''.join(str((x >> i) & 1) for i in xrange(width - 1, -1, -1))

    def test_locs(mask):
        nfalse = sum(np.array(mask) == 0)

        remaining = 0
        for s, e in com.split_ranges(mask):
            remaining += e - s

            assert 0 not in mask[s:e]

        # make sure the total items covered by the ranges are a complete cover
        assert remaining + nfalse == len(mask)

    # exhaustively test all possible mask sequences of length 8
    ncols = 8
    for i in range(2 ** ncols):
        cols = map(int, list(_bin(i, ncols)))  # count up in base2
        mask = [cols[i] == 1 for i in range(len(cols))]
        test_locs(mask)

    # base cases
    test_locs([])
    test_locs([0])
    test_locs([1])


def test_indent():
    s = 'a b c\nd e f'
    result = com.indent(s, spaces=6)

    assert(result == '      a b c\n      d e f')


def test_banner():
    ban = com.banner('hi')
    assert(ban == ('%s\nhi\n%s' % ('=' * 80, '=' * 80)))


def test_map_indices_py():
    data = [4, 3, 2, 1]
    expected = {4: 0, 3: 1, 2: 2, 1: 3}

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
    expected = {'f': ['foo', 'foo3'],
                'b': ['bar', 'baz', 'baz2'],
                'q': ['qux']}

    grouped = com.groupby(values, lambda x: x[0])

    for k, v in grouped:
        assert v == expected[k]


def test_is_list_like():
    passes = ([], [1], (1,), (1, 2), {'a': 1}, set([1, 'a']), Series([1]),
              Series([]), Series(['a']).str)
    fails = (1, '2', object())

    for p in passes:
        assert com.is_list_like(p)

    for f in fails:
        assert not com.is_list_like(f)


def test_ensure_int32():
    values = np.arange(10, dtype=np.int32)
    result = com._ensure_int32(values)
    assert(result.dtype == np.int32)

    values = np.arange(10, dtype=np.int64)
    result = com._ensure_int32(values)
    assert(result.dtype == np.int32)

def test_ensure_platform_int():

    # verify that when we create certain types of indices
    # they remain the correct type under platform conversions
    from pandas.core.index import Int64Index

    # int64
    x = Int64Index([1, 2, 3], dtype='int64')
    assert(x.dtype == np.int64)

    pi = com._ensure_platform_int(x)
    assert(pi.dtype == np.int_)

    # int32
    x = Int64Index([1, 2, 3], dtype='int32')
    assert(x.dtype == np.int32)

    pi = com._ensure_platform_int(x)
    assert(pi.dtype == np.int_)

# TODO: fix this broken test

# def test_console_encode():
#     """
#     On Python 2, if sys.stdin.encoding is None (IPython with zmq frontend)
#     common.console_encode should encode things as utf-8.
#     """
#     if py3compat.PY3:
#         raise nose.SkipTest

#     with tm.stdin_encoding(encoding=None):
#         result = com.console_encode(u"\u05d0")
#         expected = u"\u05d0".encode('utf-8')
#         assert (result == expected)


def test_is_re():
    passes = re.compile('ad'),
    fails = 'x', 2, 3, object()

    for p in passes:
        assert com.is_re(p)

    for f in fails:
        assert not com.is_re(f)


def test_is_recompilable():
    passes = (r'a', u'x', r'asdf', re.compile('adsf'), ur'\u2233\s*',
              re.compile(r''))
    fails = 1, [], object()

    for p in passes:
        assert com.is_re_compilable(p)

    for f in fails:
        assert not com.is_re_compilable(f)


class TestTake(unittest.TestCase):

    _multiprocess_can_split_ = True

    def test_1d_with_out(self):
        def _test_dtype(dtype, can_hold_na):
            data = np.random.randint(0, 2, 4).astype(dtype)

            indexer = [2, 1, 0, 1]
            out = np.empty(4, dtype=dtype)
            com.take_1d(data, indexer, out=out)
            expected = data.take(indexer)
            tm.assert_almost_equal(out, expected)

            indexer = [2, 1, 0, -1]
            out = np.empty(4, dtype=dtype)
            if can_hold_na:
                com.take_1d(data, indexer, out=out)
                expected = data.take(indexer)
                expected[3] = np.nan
                tm.assert_almost_equal(out, expected)
            else:
                self.assertRaises(Exception, com.take_1d, data,
                                  indexer, out=out)
                # no exception o/w
                data.take(indexer, out=out)

        _test_dtype(np.float64, True)
        _test_dtype(np.float32, True)
        _test_dtype(np.uint64, False)
        _test_dtype(np.uint32, False)
        _test_dtype(np.uint16, False)
        _test_dtype(np.uint8, False)
        _test_dtype(np.int64, False)
        _test_dtype(np.int32, False)
        _test_dtype(np.int16, False)
        _test_dtype(np.int8, False)
        _test_dtype(np.object_, True)
        _test_dtype(np.bool, False)

    def test_1d_fill_nonna(self):
        def _test_dtype(dtype, fill_value, out_dtype):
            data = np.random.randint(0, 2, 4).astype(dtype)

            indexer = [2, 1, 0, -1]

            result = com.take_1d(data, indexer, fill_value=fill_value)
            assert((result[[0, 1, 2]] == data[[2, 1, 0]]).all())
            assert(result[3] == fill_value)
            assert(result.dtype == out_dtype)

            indexer = [2, 1, 0, 1]

            result = com.take_1d(data, indexer, fill_value=fill_value)
            assert((result[[0, 1, 2, 3]] == data[indexer]).all())
            assert(result.dtype == dtype)

        _test_dtype(np.int8, np.int16(127), np.int8)
        _test_dtype(np.int8, np.int16(128), np.int16)
        _test_dtype(np.int32, 1, np.int32)
        _test_dtype(np.int32, 2.0, np.float64)
        _test_dtype(np.int32, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.int32, True, np.object_)
        _test_dtype(np.int32, '', np.object_)
        _test_dtype(np.float64, 1, np.float64)
        _test_dtype(np.float64, 2.0, np.float64)
        _test_dtype(np.float64, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.float64, True, np.object_)
        _test_dtype(np.float64, '', np.object_)
        _test_dtype(np.complex128, 1, np.complex128)
        _test_dtype(np.complex128, 2.0, np.complex128)
        _test_dtype(np.complex128, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.complex128, True, np.object_)
        _test_dtype(np.complex128, '', np.object_)
        _test_dtype(np.bool_, 1, np.object_)
        _test_dtype(np.bool_, 2.0, np.object_)
        _test_dtype(np.bool_, 3.0 + 4.0j, np.object_)
        _test_dtype(np.bool_, True, np.bool_)
        _test_dtype(np.bool_, '', np.object_)

    def test_2d_with_out(self):
        def _test_dtype(dtype, can_hold_na):
            data = np.random.randint(0, 2, (5, 3)).astype(dtype)

            indexer = [2, 1, 0, 1]
            out0 = np.empty((4, 3), dtype=dtype)
            out1 = np.empty((5, 4), dtype=dtype)
            com.take_nd(data, indexer, out=out0, axis=0)
            com.take_nd(data, indexer, out=out1, axis=1)
            expected0 = data.take(indexer, axis=0)
            expected1 = data.take(indexer, axis=1)
            tm.assert_almost_equal(out0, expected0)
            tm.assert_almost_equal(out1, expected1)

            indexer = [2, 1, 0, -1]
            out0 = np.empty((4, 3), dtype=dtype)
            out1 = np.empty((5, 4), dtype=dtype)
            if can_hold_na:
                com.take_nd(data, indexer, out=out0, axis=0)
                com.take_nd(data, indexer, out=out1, axis=1)
                expected0 = data.take(indexer, axis=0)
                expected1 = data.take(indexer, axis=1)
                expected0[3, :] = np.nan
                expected1[:, 3] = np.nan
                tm.assert_almost_equal(out0, expected0)
                tm.assert_almost_equal(out1, expected1)
            else:
                self.assertRaises(Exception, com.take_nd, data,
                                  indexer, out=out0, axis=0)
                self.assertRaises(Exception, com.take_nd, data,
                                  indexer, out=out1, axis=1)
                # no exception o/w
                data.take(indexer, out=out0, axis=0)
                data.take(indexer, out=out1, axis=1)

        _test_dtype(np.float64, True)
        _test_dtype(np.float32, True)
        _test_dtype(np.uint64, False)
        _test_dtype(np.uint32, False)
        _test_dtype(np.uint16, False)
        _test_dtype(np.uint8, False)
        _test_dtype(np.int64, False)
        _test_dtype(np.int32, False)
        _test_dtype(np.int16, False)
        _test_dtype(np.int8, False)
        _test_dtype(np.object_, True)
        _test_dtype(np.bool, False)

    def test_2d_fill_nonna(self):
        def _test_dtype(dtype, fill_value, out_dtype):
            data = np.random.randint(0, 2, (5, 3)).astype(dtype)

            indexer = [2, 1, 0, -1]

            result = com.take_nd(data, indexer, axis=0, fill_value=fill_value)
            assert((result[[0, 1, 2], :] == data[[2, 1, 0], :]).all())
            assert((result[3, :] == fill_value).all())
            assert(result.dtype == out_dtype)

            result = com.take_nd(data, indexer, axis=1, fill_value=fill_value)
            assert((result[:, [0, 1, 2]] == data[:, [2, 1, 0]]).all())
            assert((result[:, 3] == fill_value).all())
            assert(result.dtype == out_dtype)

            indexer = [2, 1, 0, 1]

            result = com.take_nd(data, indexer, axis=0, fill_value=fill_value)
            assert((result[[0, 1, 2, 3], :] == data[indexer, :]).all())
            assert(result.dtype == dtype)

            result = com.take_nd(data, indexer, axis=1, fill_value=fill_value)
            assert((result[:, [0, 1, 2, 3]] == data[:, indexer]).all())
            assert(result.dtype == dtype)

        _test_dtype(np.int8, np.int16(127), np.int8)
        _test_dtype(np.int8, np.int16(128), np.int16)
        _test_dtype(np.int32, 1, np.int32)
        _test_dtype(np.int32, 2.0, np.float64)
        _test_dtype(np.int32, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.int32, True, np.object_)
        _test_dtype(np.int32, '', np.object_)
        _test_dtype(np.float64, 1, np.float64)
        _test_dtype(np.float64, 2.0, np.float64)
        _test_dtype(np.float64, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.float64, True, np.object_)
        _test_dtype(np.float64, '', np.object_)
        _test_dtype(np.complex128, 1, np.complex128)
        _test_dtype(np.complex128, 2.0, np.complex128)
        _test_dtype(np.complex128, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.complex128, True, np.object_)
        _test_dtype(np.complex128, '', np.object_)
        _test_dtype(np.bool_, 1, np.object_)
        _test_dtype(np.bool_, 2.0, np.object_)
        _test_dtype(np.bool_, 3.0 + 4.0j, np.object_)
        _test_dtype(np.bool_, True, np.bool_)
        _test_dtype(np.bool_, '', np.object_)

    def test_3d_with_out(self):
        def _test_dtype(dtype, can_hold_na):
            data = np.random.randint(0, 2, (5, 4, 3)).astype(dtype)

            indexer = [2, 1, 0, 1]
            out0 = np.empty((4, 4, 3), dtype=dtype)
            out1 = np.empty((5, 4, 3), dtype=dtype)
            out2 = np.empty((5, 4, 4), dtype=dtype)
            com.take_nd(data, indexer, out=out0, axis=0)
            com.take_nd(data, indexer, out=out1, axis=1)
            com.take_nd(data, indexer, out=out2, axis=2)
            expected0 = data.take(indexer, axis=0)
            expected1 = data.take(indexer, axis=1)
            expected2 = data.take(indexer, axis=2)
            tm.assert_almost_equal(out0, expected0)
            tm.assert_almost_equal(out1, expected1)
            tm.assert_almost_equal(out2, expected2)

            indexer = [2, 1, 0, -1]
            out0 = np.empty((4, 4, 3), dtype=dtype)
            out1 = np.empty((5, 4, 3), dtype=dtype)
            out2 = np.empty((5, 4, 4), dtype=dtype)
            if can_hold_na:
                com.take_nd(data, indexer, out=out0, axis=0)
                com.take_nd(data, indexer, out=out1, axis=1)
                com.take_nd(data, indexer, out=out2, axis=2)
                expected0 = data.take(indexer, axis=0)
                expected1 = data.take(indexer, axis=1)
                expected2 = data.take(indexer, axis=2)
                expected0[3, :, :] = np.nan
                expected1[:, 3, :] = np.nan
                expected2[:, :, 3] = np.nan
                tm.assert_almost_equal(out0, expected0)
                tm.assert_almost_equal(out1, expected1)
                tm.assert_almost_equal(out2, expected2)
            else:
                self.assertRaises(Exception, com.take_nd, data,
                                  indexer, out=out0, axis=0)
                self.assertRaises(Exception, com.take_nd, data,
                                  indexer, out=out1, axis=1)
                self.assertRaises(Exception, com.take_nd, data,
                                  indexer, out=out2, axis=2)
                # no exception o/w
                data.take(indexer, out=out0, axis=0)
                data.take(indexer, out=out1, axis=1)
                data.take(indexer, out=out2, axis=2)

        _test_dtype(np.float64, True)
        _test_dtype(np.float32, True)
        _test_dtype(np.uint64, False)
        _test_dtype(np.uint32, False)
        _test_dtype(np.uint16, False)
        _test_dtype(np.uint8, False)
        _test_dtype(np.int64, False)
        _test_dtype(np.int32, False)
        _test_dtype(np.int16, False)
        _test_dtype(np.int8, False)
        _test_dtype(np.object_, True)
        _test_dtype(np.bool, False)

    def test_3d_fill_nonna(self):
        def _test_dtype(dtype, fill_value, out_dtype):
            data = np.random.randint(0, 2, (5, 4, 3)).astype(dtype)

            indexer = [2, 1, 0, -1]

            result = com.take_nd(data, indexer, axis=0, fill_value=fill_value)
            assert((result[[0, 1, 2], :, :] == data[[2, 1, 0], :, :]).all())
            assert((result[3, :, :] == fill_value).all())
            assert(result.dtype == out_dtype)

            result = com.take_nd(data, indexer, axis=1, fill_value=fill_value)
            assert((result[:, [0, 1, 2], :] == data[:, [2, 1, 0], :]).all())
            assert((result[:, 3, :] == fill_value).all())
            assert(result.dtype == out_dtype)

            result = com.take_nd(data, indexer, axis=2, fill_value=fill_value)
            assert((result[:, :, [0, 1, 2]] == data[:, :, [2, 1, 0]]).all())
            assert((result[:, :, 3] == fill_value).all())
            assert(result.dtype == out_dtype)

            indexer = [2, 1, 0, 1]

            result = com.take_nd(data, indexer, axis=0, fill_value=fill_value)
            assert((result[[0, 1, 2, 3], :, :] == data[indexer, :, :]).all())
            assert(result.dtype == dtype)

            result = com.take_nd(data, indexer, axis=1, fill_value=fill_value)
            assert((result[:, [0, 1, 2, 3], :] == data[:, indexer, :]).all())
            assert(result.dtype == dtype)

            result = com.take_nd(data, indexer, axis=2, fill_value=fill_value)
            assert((result[:, :, [0, 1, 2, 3]] == data[:, :, indexer]).all())
            assert(result.dtype == dtype)

        _test_dtype(np.int8, np.int16(127), np.int8)
        _test_dtype(np.int8, np.int16(128), np.int16)
        _test_dtype(np.int32, 1, np.int32)
        _test_dtype(np.int32, 2.0, np.float64)
        _test_dtype(np.int32, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.int32, True, np.object_)
        _test_dtype(np.int32, '', np.object_)
        _test_dtype(np.float64, 1, np.float64)
        _test_dtype(np.float64, 2.0, np.float64)
        _test_dtype(np.float64, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.float64, True, np.object_)
        _test_dtype(np.float64, '', np.object_)
        _test_dtype(np.complex128, 1, np.complex128)
        _test_dtype(np.complex128, 2.0, np.complex128)
        _test_dtype(np.complex128, 3.0 + 4.0j, np.complex128)
        _test_dtype(np.complex128, True, np.object_)
        _test_dtype(np.complex128, '', np.object_)
        _test_dtype(np.bool_, 1, np.object_)
        _test_dtype(np.bool_, 2.0, np.object_)
        _test_dtype(np.bool_, 3.0 + 4.0j, np.object_)
        _test_dtype(np.bool_, True, np.bool_)
        _test_dtype(np.bool_, '', np.object_)

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
        result = com.take_nd(arr, indexer, axis=0)
        expected = arr.take(indexer, axis=0)
        expected[-1] = np.nan
        tm.assert_almost_equal(result, expected)

        # axis=1
        result = com.take_nd(arr, indexer, axis=1)
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

        result = com.take_nd(arr, [0, 2, 2, 1])
        expected = arr.take([0, 2, 2, 1], axis=0)
        self.assert_(np.array_equal(result, expected))

        result = com.take_nd(arr, [0, 2, 2, 1], axis=1)
        expected = arr.take([0, 2, 2, 1], axis=1)
        self.assert_(np.array_equal(result, expected))

        result = com.take_nd(arr, [0, 2, -1])
        self.assert_(result.dtype == np.object_)

    def test_2d_float32(self):
        arr = np.random.randn(4, 3).astype(np.float32)
        indexer = [0, 2, -1, 1, -1]

        # axis=0
        result = com.take_nd(arr, indexer, axis=0)
        result2 = np.empty_like(result)
        com.take_nd(arr, indexer, axis=0, out=result2)
        tm.assert_almost_equal(result, result2)

        expected = arr.take(indexer, axis=0)
        expected[[2, 4], :] = np.nan
        tm.assert_almost_equal(result, expected)

        #### this now accepts a float32! # test with float64 out buffer
        out = np.empty((len(indexer), arr.shape[1]), dtype='float32')
        com.take_nd(arr, indexer, out=out)  # it works!

        # axis=1
        result = com.take_nd(arr, indexer, axis=1)
        result2 = np.empty_like(result)
        com.take_nd(arr, indexer, axis=1, out=result2)
        tm.assert_almost_equal(result, result2)

        expected = arr.take(indexer, axis=1)
        expected[:, [2, 4]] = np.nan
        tm.assert_almost_equal(result, expected)

    def test_2d_datetime64(self):
        # 2005/01/01 - 2006/01/01
        arr = np.random.randint(11045376L, 11360736L, (5,3))*100000000000
        arr = arr.view(dtype='datetime64[ns]')
        indexer = [0, 2, -1, 1, -1]

        # axis=0
        result = com.take_nd(arr, indexer, axis=0)
        result2 = np.empty_like(result)
        com.take_nd(arr, indexer, axis=0, out=result2)
        tm.assert_almost_equal(result, result2)

        expected = arr.take(indexer, axis=0)
        expected.view(np.int64)[[2, 4], :] = iNaT
        tm.assert_almost_equal(result, expected)

        result = com.take_nd(arr, indexer, axis=0,
                             fill_value=datetime(2007, 1, 1))
        result2 = np.empty_like(result)
        com.take_nd(arr, indexer, out=result2, axis=0,
                    fill_value=datetime(2007, 1, 1))
        tm.assert_almost_equal(result, result2)

        expected = arr.take(indexer, axis=0)
        expected[[2, 4], :] = datetime(2007, 1, 1)
        tm.assert_almost_equal(result, expected)

        # axis=1
        result = com.take_nd(arr, indexer, axis=1)
        result2 = np.empty_like(result)
        com.take_nd(arr, indexer, axis=1, out=result2)
        tm.assert_almost_equal(result, result2)

        expected = arr.take(indexer, axis=1)
        expected.view(np.int64)[:, [2, 4]] = iNaT
        tm.assert_almost_equal(result, expected)

        result = com.take_nd(arr, indexer, axis=1,
                             fill_value=datetime(2007, 1, 1))
        result2 = np.empty_like(result)
        com.take_nd(arr, indexer, out=result2, axis=1,
                    fill_value=datetime(2007, 1, 1))
        tm.assert_almost_equal(result, result2)

        expected = arr.take(indexer, axis=1)
        expected[:, [2, 4]] = datetime(2007, 1, 1)
        tm.assert_almost_equal(result, expected)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
