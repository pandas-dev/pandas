# -*- coding: utf-8 -*-
import collections
from datetime import datetime, timedelta
import re

import nose
import numpy as np
import pandas as pd
from pandas.tslib import iNaT, NaT
from pandas import (Series, DataFrame, date_range, DatetimeIndex,
                    TimedeltaIndex, Timestamp, Float64Index)
from pandas import compat
from pandas.compat import range, lrange, lmap, u
from pandas.core.common import notnull, isnull, array_equivalent
import pandas.core.common as com
import pandas.core.convert as convert
import pandas.util.testing as tm
import pandas.core.config as cf

_multiprocess_can_split_ = True


def test_mut_exclusive():
    msg = "mutually exclusive arguments: '[ab]' and '[ab]'"
    with tm.assertRaisesRegexp(TypeError, msg):
        com._mut_exclusive(a=1, b=2)
    assert com._mut_exclusive(a=1, b=None) == 1
    assert com._mut_exclusive(major=None, major_axis=None) is None


def test_is_sequence():
    is_seq = com.is_sequence
    assert (is_seq((1, 2)))
    assert (is_seq([1, 2]))
    assert (not is_seq("abcd"))
    assert (not is_seq(u("abcd")))
    assert (not is_seq(np.int64))

    class A(object):

        def __getitem__(self):
            return 1

    assert (not is_seq(A()))


def test_get_callable_name():
    from functools import partial
    getname = com._get_callable_name

    def fn(x):
        return x

    lambda_ = lambda x: x
    part1 = partial(fn)
    part2 = partial(part1)

    class somecall(object):

        def __call__(self):
            return x  # noqa

    assert getname(fn) == 'fn'
    assert getname(lambda_)
    assert getname(part1) == 'fn'
    assert getname(part2) == 'fn'
    assert getname(somecall()) == 'somecall'
    assert getname(1) is None


class TestInferDtype(tm.TestCase):

    def test_infer_dtype_from_scalar(self):
        # Test that _infer_dtype_from_scalar is returning correct dtype for int
        # and float.

        for dtypec in [np.uint8, np.int8, np.uint16, np.int16, np.uint32,
                       np.int32, np.uint64, np.int64]:
            data = dtypec(12)
            dtype, val = com._infer_dtype_from_scalar(data)
            self.assertEqual(dtype, type(data))

        data = 12
        dtype, val = com._infer_dtype_from_scalar(data)
        self.assertEqual(dtype, np.int64)

        for dtypec in [np.float16, np.float32, np.float64]:
            data = dtypec(12)
            dtype, val = com._infer_dtype_from_scalar(data)
            self.assertEqual(dtype, dtypec)

        data = np.float(12)
        dtype, val = com._infer_dtype_from_scalar(data)
        self.assertEqual(dtype, np.float64)

        for data in [True, False]:
            dtype, val = com._infer_dtype_from_scalar(data)
            self.assertEqual(dtype, np.bool_)

        for data in [np.complex64(1), np.complex128(1)]:
            dtype, val = com._infer_dtype_from_scalar(data)
            self.assertEqual(dtype, np.complex_)

        import datetime
        for data in [np.datetime64(1, 'ns'), pd.Timestamp(1),
                     datetime.datetime(2000, 1, 1, 0, 0)]:
            dtype, val = com._infer_dtype_from_scalar(data)
            self.assertEqual(dtype, 'M8[ns]')

        for data in [np.timedelta64(1, 'ns'), pd.Timedelta(1),
                     datetime.timedelta(1)]:
            dtype, val = com._infer_dtype_from_scalar(data)
            self.assertEqual(dtype, 'm8[ns]')

        for data in [datetime.date(2000, 1, 1),
                     pd.Timestamp(1, tz='US/Eastern'), 'foo']:
            dtype, val = com._infer_dtype_from_scalar(data)
            self.assertEqual(dtype, np.object_)


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
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(),
                  tm.makeObjectSeries(), tm.makeTimeSeries(),
                  tm.makePeriodSeries()]:
            assert (isinstance(isnull(s), Series))


def test_isnull():
    assert not isnull(1.)
    assert isnull(None)
    assert isnull(np.NaN)
    assert not isnull(np.inf)
    assert not isnull(-np.inf)

    # series
    for s in [tm.makeFloatSeries(), tm.makeStringSeries(),
              tm.makeObjectSeries(), tm.makeTimeSeries(),
              tm.makePeriodSeries()]:
        assert (isinstance(isnull(s), Series))

    # frame
    for df in [tm.makeTimeDataFrame(), tm.makePeriodFrame(),
               tm.makeMixedDataFrame()]:
        result = isnull(df)
        expected = df.apply(isnull)
        tm.assert_frame_equal(result, expected)

    # panel
    for p in [tm.makePanel(), tm.makePeriodPanel(), tm.add_nans(tm.makePanel())
              ]:
        result = isnull(p)
        expected = p.apply(isnull)
        tm.assert_panel_equal(result, expected)

    # panel 4d
    for p in [tm.makePanel4D(), tm.add_nans_panel4d(tm.makePanel4D())]:
        result = isnull(p)
        expected = p.apply(isnull)
        tm.assert_panel4d_equal(result, expected)


def test_isnull_lists():
    result = isnull([[False]])
    exp = np.array([[False]])
    assert (np.array_equal(result, exp))

    result = isnull([[1], [2]])
    exp = np.array([[False], [False]])
    assert (np.array_equal(result, exp))

    # list of strings / unicode
    result = isnull(['foo', 'bar'])
    assert (not result.any())

    result = isnull([u('foo'), u('bar')])
    assert (not result.any())


def test_isnull_nat():
    result = isnull([NaT])
    exp = np.array([True])
    assert (np.array_equal(result, exp))

    result = isnull(np.array([NaT], dtype=object))
    exp = np.array([True])
    assert (np.array_equal(result, exp))


def test_isnull_numpy_nat():
    arr = np.array([NaT, np.datetime64('NaT'), np.timedelta64('NaT'),
                    np.datetime64('NaT', 's')])
    result = isnull(arr)
    expected = np.array([True] * 4)
    tm.assert_numpy_array_equal(result, expected)


def test_isnull_datetime():
    assert (not isnull(datetime.now()))
    assert notnull(datetime.now())

    idx = date_range('1/1/1990', periods=20)
    assert (notnull(idx).all())

    idx = np.asarray(idx)
    idx[0] = iNaT
    idx = DatetimeIndex(idx)
    mask = isnull(idx)
    assert (mask[0])
    assert (not mask[1:].any())

    # GH 9129
    pidx = idx.to_period(freq='M')
    mask = isnull(pidx)
    assert (mask[0])
    assert (not mask[1:].any())

    mask = isnull(pidx[1:])
    assert (not mask.any())


class TestIsNull(tm.TestCase):

    def test_0d_array(self):
        self.assertTrue(isnull(np.array(np.nan)))
        self.assertFalse(isnull(np.array(0.0)))
        self.assertFalse(isnull(np.array(0)))
        # test object dtype
        self.assertTrue(isnull(np.array(np.nan, dtype=object)))
        self.assertFalse(isnull(np.array(0.0, dtype=object)))
        self.assertFalse(isnull(np.array(0, dtype=object)))


class TestNumberScalar(tm.TestCase):

    def test_is_number(self):

        self.assertTrue(com.is_number(True))
        self.assertTrue(com.is_number(1))
        self.assertTrue(com.is_number(1.1))
        self.assertTrue(com.is_number(1 + 3j))
        self.assertTrue(com.is_number(np.bool(False)))
        self.assertTrue(com.is_number(np.int64(1)))
        self.assertTrue(com.is_number(np.float64(1.1)))
        self.assertTrue(com.is_number(np.complex128(1 + 3j)))
        self.assertTrue(com.is_number(np.nan))

        self.assertFalse(com.is_number(None))
        self.assertFalse(com.is_number('x'))
        self.assertFalse(com.is_number(datetime(2011, 1, 1)))
        self.assertFalse(com.is_number(np.datetime64('2011-01-01')))
        self.assertFalse(com.is_number(pd.Timestamp('2011-01-01')))
        self.assertFalse(com.is_number(pd.Timestamp('2011-01-01',
                                                    tz='US/Eastern')))
        self.assertFalse(com.is_number(timedelta(1000)))
        self.assertFalse(com.is_number(pd.Timedelta('1 days')))

        # questionable
        self.assertFalse(com.is_number(np.bool_(False)))
        self.assertTrue(com.is_number(np.timedelta64(1, 'D')))

    def test_is_bool(self):
        self.assertTrue(com.is_bool(True))
        self.assertTrue(com.is_bool(np.bool(False)))
        self.assertTrue(com.is_bool(np.bool_(False)))

        self.assertFalse(com.is_bool(1))
        self.assertFalse(com.is_bool(1.1))
        self.assertFalse(com.is_bool(1 + 3j))
        self.assertFalse(com.is_bool(np.int64(1)))
        self.assertFalse(com.is_bool(np.float64(1.1)))
        self.assertFalse(com.is_bool(np.complex128(1 + 3j)))
        self.assertFalse(com.is_bool(np.nan))
        self.assertFalse(com.is_bool(None))
        self.assertFalse(com.is_bool('x'))
        self.assertFalse(com.is_bool(datetime(2011, 1, 1)))
        self.assertFalse(com.is_bool(np.datetime64('2011-01-01')))
        self.assertFalse(com.is_bool(pd.Timestamp('2011-01-01')))
        self.assertFalse(com.is_bool(pd.Timestamp('2011-01-01',
                                                  tz='US/Eastern')))
        self.assertFalse(com.is_bool(timedelta(1000)))
        self.assertFalse(com.is_bool(np.timedelta64(1, 'D')))
        self.assertFalse(com.is_bool(pd.Timedelta('1 days')))

    def test_is_integer(self):
        self.assertTrue(com.is_integer(1))
        self.assertTrue(com.is_integer(np.int64(1)))

        self.assertFalse(com.is_integer(True))
        self.assertFalse(com.is_integer(1.1))
        self.assertFalse(com.is_integer(1 + 3j))
        self.assertFalse(com.is_integer(np.bool(False)))
        self.assertFalse(com.is_integer(np.bool_(False)))
        self.assertFalse(com.is_integer(np.float64(1.1)))
        self.assertFalse(com.is_integer(np.complex128(1 + 3j)))
        self.assertFalse(com.is_integer(np.nan))
        self.assertFalse(com.is_integer(None))
        self.assertFalse(com.is_integer('x'))
        self.assertFalse(com.is_integer(datetime(2011, 1, 1)))
        self.assertFalse(com.is_integer(np.datetime64('2011-01-01')))
        self.assertFalse(com.is_integer(pd.Timestamp('2011-01-01')))
        self.assertFalse(com.is_integer(pd.Timestamp('2011-01-01',
                                                     tz='US/Eastern')))
        self.assertFalse(com.is_integer(timedelta(1000)))
        self.assertFalse(com.is_integer(pd.Timedelta('1 days')))

        # questionable
        self.assertTrue(com.is_integer(np.timedelta64(1, 'D')))

    def test_is_float(self):
        self.assertTrue(com.is_float(1.1))
        self.assertTrue(com.is_float(np.float64(1.1)))
        self.assertTrue(com.is_float(np.nan))

        self.assertFalse(com.is_float(True))
        self.assertFalse(com.is_float(1))
        self.assertFalse(com.is_float(1 + 3j))
        self.assertFalse(com.is_float(np.bool(False)))
        self.assertFalse(com.is_float(np.bool_(False)))
        self.assertFalse(com.is_float(np.int64(1)))
        self.assertFalse(com.is_float(np.complex128(1 + 3j)))
        self.assertFalse(com.is_float(None))
        self.assertFalse(com.is_float('x'))
        self.assertFalse(com.is_float(datetime(2011, 1, 1)))
        self.assertFalse(com.is_float(np.datetime64('2011-01-01')))
        self.assertFalse(com.is_float(pd.Timestamp('2011-01-01')))
        self.assertFalse(com.is_float(pd.Timestamp('2011-01-01',
                                                   tz='US/Eastern')))
        self.assertFalse(com.is_float(timedelta(1000)))
        self.assertFalse(com.is_float(np.timedelta64(1, 'D')))
        self.assertFalse(com.is_float(pd.Timedelta('1 days')))


def test_downcast_conv():
    # test downcasting

    arr = np.array([8.5, 8.6, 8.7, 8.8, 8.9999999999995])
    result = com._possibly_downcast_to_dtype(arr, 'infer')
    assert (np.array_equal(result, arr))

    arr = np.array([8., 8., 8., 8., 8.9999999999995])
    result = com._possibly_downcast_to_dtype(arr, 'infer')
    expected = np.array([8, 8, 8, 8, 9])
    assert (np.array_equal(result, expected))

    arr = np.array([8., 8., 8., 8., 9.0000000000005])
    result = com._possibly_downcast_to_dtype(arr, 'infer')
    expected = np.array([8, 8, 8, 8, 9])
    assert (np.array_equal(result, expected))

    # conversions

    expected = np.array([1, 2])
    for dtype in [np.float64, object, np.int64]:
        arr = np.array([1.0, 2.0], dtype=dtype)
        result = com._possibly_downcast_to_dtype(arr, 'infer')
        tm.assert_almost_equal(result, expected)

    expected = np.array([1.0, 2.0, np.nan])
    for dtype in [np.float64, object]:
        arr = np.array([1.0, 2.0, np.nan], dtype=dtype)
        result = com._possibly_downcast_to_dtype(arr, 'infer')
        tm.assert_almost_equal(result, expected)

    # empties
    for dtype in [np.int32, np.float64, np.float32, np.bool_, np.int64, object
                  ]:
        arr = np.array([], dtype=dtype)
        result = com._possibly_downcast_to_dtype(arr, 'int64')
        tm.assert_almost_equal(result, np.array([], dtype=np.int64))
        assert result.dtype == np.int64


def test_array_equivalent():
    assert array_equivalent(np.array([np.nan, np.nan]),
                            np.array([np.nan, np.nan]))
    assert array_equivalent(np.array([np.nan, 1, np.nan]),
                            np.array([np.nan, 1, np.nan]))
    assert array_equivalent(np.array([np.nan, None], dtype='object'),
                            np.array([np.nan, None], dtype='object'))
    assert array_equivalent(np.array([np.nan, 1 + 1j], dtype='complex'),
                            np.array([np.nan, 1 + 1j], dtype='complex'))
    assert not array_equivalent(
        np.array([np.nan, 1 + 1j], dtype='complex'), np.array(
            [np.nan, 1 + 2j], dtype='complex'))
    assert not array_equivalent(
        np.array([np.nan, 1, np.nan]), np.array([np.nan, 2, np.nan]))
    assert not array_equivalent(
        np.array(['a', 'b', 'c', 'd']), np.array(['e', 'e']))
    assert array_equivalent(Float64Index([0, np.nan]),
                            Float64Index([0, np.nan]))
    assert not array_equivalent(
        Float64Index([0, np.nan]), Float64Index([1, np.nan]))
    assert array_equivalent(DatetimeIndex([0, np.nan]),
                            DatetimeIndex([0, np.nan]))
    assert not array_equivalent(
        DatetimeIndex([0, np.nan]), DatetimeIndex([1, np.nan]))
    assert array_equivalent(TimedeltaIndex([0, np.nan]),
                            TimedeltaIndex([0, np.nan]))
    assert not array_equivalent(
        TimedeltaIndex([0, np.nan]), TimedeltaIndex([1, np.nan]))
    assert array_equivalent(DatetimeIndex([0, np.nan], tz='US/Eastern'),
                            DatetimeIndex([0, np.nan], tz='US/Eastern'))
    assert not array_equivalent(
        DatetimeIndex([0, np.nan], tz='US/Eastern'), DatetimeIndex(
            [1, np.nan], tz='US/Eastern'))
    assert not array_equivalent(
        DatetimeIndex([0, np.nan]), DatetimeIndex(
            [0, np.nan], tz='US/Eastern'))
    assert not array_equivalent(
        DatetimeIndex([0, np.nan], tz='CET'), DatetimeIndex(
            [0, np.nan], tz='US/Eastern'))
    assert not array_equivalent(
        DatetimeIndex([0, np.nan]), TimedeltaIndex([0, np.nan]))


def test_datetimeindex_from_empty_datetime64_array():
    for unit in ['ms', 'us', 'ns']:
        idx = DatetimeIndex(np.array([], dtype='datetime64[%s]' % unit))
        assert (len(idx) == 0)


def test_nan_to_nat_conversions():

    df = DataFrame(dict({
        'A': np.asarray(
            lrange(10), dtype='float64'),
        'B': Timestamp('20010101')
    }))
    df.iloc[3:6, :] = np.nan
    result = df.loc[4, 'B'].value
    assert (result == iNaT)

    s = df['B'].copy()
    s._data = s._data.setitem(indexer=tuple([slice(8, 9)]), value=np.nan)
    assert (isnull(s[8]))

    # numpy < 1.7.0 is wrong
    from distutils.version import LooseVersion
    if LooseVersion(np.__version__) >= '1.7.0':
        assert (s[8].value == np.datetime64('NaT').astype(np.int64))


def test_any_none():
    assert (com._any_none(1, 2, 3, None))
    assert (not com._any_none(1, 2, 3, 4))


def test_all_not_none():
    assert (com._all_not_none(1, 2, 3, 4))
    assert (not com._all_not_none(1, 2, 3, None))
    assert (not com._all_not_none(None, None, None, None))


def test_iterpairs():
    data = [1, 2, 3, 4]
    expected = [(1, 2), (2, 3), (3, 4)]

    result = list(com.iterpairs(data))

    assert (result == expected)


def test_split_ranges():
    def _bin(x, width):
        "return int(x) as a base2 string of given width"
        return ''.join(str((x >> i) & 1) for i in range(width - 1, -1, -1))

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
        cols = lmap(int, list(_bin(i, ncols)))  # count up in base2
        mask = [cols[i] == 1 for i in range(len(cols))]
        test_locs(mask)

    # base cases
    test_locs([])
    test_locs([0])
    test_locs([1])


def test_map_indices_py():
    data = [4, 3, 2, 1]
    expected = {4: 0, 3: 1, 2: 2, 1: 3}

    result = com.map_indices_py(data)

    assert (result == expected)


def test_union():
    a = [1, 2, 3]
    b = [4, 5, 6]

    union = sorted(com.union(a, b))

    assert ((a + b) == union)


def test_difference():
    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5, 6]

    inter = sorted(com.difference(b, a))

    assert ([4, 5, 6] == inter)


def test_intersection():
    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5, 6]

    inter = sorted(com.intersection(a, b))

    assert (a == inter)


def test_groupby():
    values = ['foo', 'bar', 'baz', 'baz2', 'qux', 'foo3']
    expected = {'f': ['foo', 'foo3'],
                'b': ['bar', 'baz', 'baz2'],
                'q': ['qux']}

    grouped = com.groupby(values, lambda x: x[0])

    for k, v in grouped:
        assert v == expected[k]


def test_is_list_like():
    passes = ([], [1], (1, ), (1, 2), {'a': 1}, set([1, 'a']), Series([1]),
              Series([]), Series(['a']).str)
    fails = (1, '2', object())

    for p in passes:
        assert com.is_list_like(p)

    for f in fails:
        assert not com.is_list_like(f)


def test_is_dict_like():
    passes = [{}, {'A': 1}, pd.Series([1])]
    fails = ['1', 1, [1, 2], (1, 2), range(2), pd.Index([1])]

    for p in passes:
        assert com.is_dict_like(p)

    for f in fails:
        assert not com.is_dict_like(f)


def test_is_named_tuple():
    passes = (collections.namedtuple('Test', list('abc'))(1, 2, 3), )
    fails = ((1, 2, 3), 'a', Series({'pi': 3.14}))

    for p in passes:
        assert com.is_named_tuple(p)

    for f in fails:
        assert not com.is_named_tuple(f)


def test_is_hashable():

    # all new-style classes are hashable by default
    class HashableClass(object):
        pass

    class UnhashableClass1(object):
        __hash__ = None

    class UnhashableClass2(object):

        def __hash__(self):
            raise TypeError("Not hashable")

    hashable = (1,
                3.14,
                np.float64(3.14),
                'a',
                tuple(),
                (1, ),
                HashableClass(), )
    not_hashable = ([], UnhashableClass1(), )
    abc_hashable_not_really_hashable = (([], ), UnhashableClass2(), )

    for i in hashable:
        assert com.is_hashable(i)
    for i in not_hashable:
        assert not com.is_hashable(i)
    for i in abc_hashable_not_really_hashable:
        assert not com.is_hashable(i)

    # numpy.array is no longer collections.Hashable as of
    # https://github.com/numpy/numpy/pull/5326, just test
    # pandas.common.is_hashable()
    assert not com.is_hashable(np.array([]))

    # old-style classes in Python 2 don't appear hashable to
    # collections.Hashable but also seem to support hash() by default
    if compat.PY2:

        class OldStyleClass():
            pass

        c = OldStyleClass()
        assert not isinstance(c, collections.Hashable)
        assert com.is_hashable(c)
        hash(c)  # this will not raise


def test_ensure_int32():
    values = np.arange(10, dtype=np.int32)
    result = com._ensure_int32(values)
    assert (result.dtype == np.int32)

    values = np.arange(10, dtype=np.int64)
    result = com._ensure_int32(values)
    assert (result.dtype == np.int32)


def test_is_re():
    passes = re.compile('ad'),
    fails = 'x', 2, 3, object()

    for p in passes:
        assert com.is_re(p)

    for f in fails:
        assert not com.is_re(f)


def test_is_recompilable():
    passes = (r'a', u('x'), r'asdf', re.compile('adsf'), u(r'\u2233\s*'),
              re.compile(r''))
    fails = 1, [], object()

    for p in passes:
        assert com.is_re_compilable(p)

    for f in fails:
        assert not com.is_re_compilable(f)


def test_random_state():
    import numpy.random as npr
    # Check with seed
    state = com._random_state(5)
    tm.assert_equal(state.uniform(), npr.RandomState(5).uniform())

    # Check with random state object
    state2 = npr.RandomState(10)
    tm.assert_equal(
        com._random_state(state2).uniform(), npr.RandomState(10).uniform())

    # check with no arg random state
    assert isinstance(com._random_state(), npr.RandomState)

    # Error for floats or strings
    with tm.assertRaises(ValueError):
        com._random_state('test')

    with tm.assertRaises(ValueError):
        com._random_state(5.5)


def test_maybe_match_name():

    matched = com._maybe_match_name(
        Series([1], name='x'), Series(
            [2], name='x'))
    assert (matched == 'x')

    matched = com._maybe_match_name(
        Series([1], name='x'), Series(
            [2], name='y'))
    assert (matched is None)

    matched = com._maybe_match_name(Series([1]), Series([2], name='x'))
    assert (matched is None)

    matched = com._maybe_match_name(Series([1], name='x'), Series([2]))
    assert (matched is None)

    matched = com._maybe_match_name(Series([1], name='x'), [2])
    assert (matched == 'x')

    matched = com._maybe_match_name([1], Series([2], name='y'))
    assert (matched == 'y')


class TestMaybe(tm.TestCase):

    def test_maybe_convert_string_to_array(self):
        result = com._maybe_convert_string_to_object('x')
        tm.assert_numpy_array_equal(result, np.array(['x'], dtype=object))
        self.assertTrue(result.dtype == object)

        result = com._maybe_convert_string_to_object(1)
        self.assertEqual(result, 1)

        arr = np.array(['x', 'y'], dtype=str)
        result = com._maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 'y'], dtype=object))
        self.assertTrue(result.dtype == object)

        # unicode
        arr = np.array(['x', 'y']).astype('U')
        result = com._maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 'y'], dtype=object))
        self.assertTrue(result.dtype == object)

        # object
        arr = np.array(['x', 2], dtype=object)
        result = com._maybe_convert_string_to_object(arr)
        tm.assert_numpy_array_equal(result, np.array(['x', 2], dtype=object))
        self.assertTrue(result.dtype == object)

    def test_maybe_convert_scalar(self):

        # pass thru
        result = com._maybe_convert_scalar('x')
        self.assertEqual(result, 'x')
        result = com._maybe_convert_scalar(np.array([1]))
        self.assertEqual(result, np.array([1]))

        # leave scalar dtype
        result = com._maybe_convert_scalar(np.int64(1))
        self.assertEqual(result, np.int64(1))
        result = com._maybe_convert_scalar(np.int32(1))
        self.assertEqual(result, np.int32(1))
        result = com._maybe_convert_scalar(np.float32(1))
        self.assertEqual(result, np.float32(1))
        result = com._maybe_convert_scalar(np.int64(1))
        self.assertEqual(result, np.float64(1))

        # coerce
        result = com._maybe_convert_scalar(1)
        self.assertEqual(result, np.int64(1))
        result = com._maybe_convert_scalar(1.0)
        self.assertEqual(result, np.float64(1))
        result = com._maybe_convert_scalar(pd.Timestamp('20130101'))
        self.assertEqual(result, pd.Timestamp('20130101').value)
        result = com._maybe_convert_scalar(datetime(2013, 1, 1))
        self.assertEqual(result, pd.Timestamp('20130101').value)
        result = com._maybe_convert_scalar(pd.Timedelta('1 day 1 min'))
        self.assertEqual(result, pd.Timedelta('1 day 1 min').value)


class TestConvert(tm.TestCase):

    def test_possibly_convert_objects_copy(self):
        values = np.array([1, 2])

        out = convert._possibly_convert_objects(values, copy=False)
        self.assertTrue(values is out)

        out = convert._possibly_convert_objects(values, copy=True)
        self.assertTrue(values is not out)

        values = np.array(['apply', 'banana'])
        out = convert._possibly_convert_objects(values, copy=False)
        self.assertTrue(values is out)

        out = convert._possibly_convert_objects(values, copy=True)
        self.assertTrue(values is not out)


def test_dict_compat():
    data_datetime64 = {np.datetime64('1990-03-15'): 1,
                       np.datetime64('2015-03-15'): 2}
    data_unchanged = {1: 2, 3: 4, 5: 6}
    expected = {Timestamp('1990-3-15'): 1, Timestamp('2015-03-15'): 2}
    assert (com._dict_compat(data_datetime64) == expected)
    assert (com._dict_compat(expected) == expected)
    assert (com._dict_compat(data_unchanged) == data_unchanged)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
