# -*- coding: utf-8 -*-
import collections
from datetime import datetime
import re

import nose
from nose.tools import assert_equal, assert_true
import numpy as np
import pandas as pd
from pandas.tslib import iNaT, NaT
from pandas import Series, DataFrame, date_range, DatetimeIndex, Timestamp, Float64Index
from pandas import compat
from pandas.compat import range, long, lrange, lmap, u
from pandas.core.common import notnull, isnull, array_equivalent
import pandas.core.common as com
import pandas.core.convert as convert
import pandas.core.format as fmt
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
    assert(is_seq((1, 2)))
    assert(is_seq([1, 2]))
    assert(not is_seq("abcd"))
    assert(not is_seq(u("abcd")))
    assert(not is_seq(np.int64))

    class A(object):
        def __getitem__(self):
            return 1

    assert(not is_seq(A()))


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
            return x

    assert getname(fn) == 'fn'
    assert getname(lambda_)
    assert getname(part1) == 'fn'
    assert getname(part2) == 'fn'
    assert getname(somecall()) == 'somecall'
    assert getname(1) is None

#Issue 10859
class TestABCClasses(tm.TestCase):
    tuples = [[1, 2, 2], ['red', 'blue', 'red']]
    multi_index = pd.MultiIndex.from_arrays(tuples, names=('number', 'color'))
    datetime_index = pd.to_datetime(['2000/1/1', '2010/1/1'])
    timedelta_index = pd.to_timedelta(np.arange(5), unit='s')
    period_index = pd.period_range('2000/1/1', '2010/1/1/', freq='M')
    categorical = pd.Categorical([1, 2, 3], categories=[2, 3, 1])
    categorical_df = pd.DataFrame({"values": [1, 2, 3]}, index=categorical)
    df = pd.DataFrame({'names': ['a', 'b', 'c']}, index=multi_index)
    sparse_series = pd.Series([1, 2, 3]).to_sparse()
    sparse_array = pd.SparseArray(np.random.randn(10))

    def test_abc_types(self):
        self.assertIsInstance(pd.Index(['a', 'b', 'c']), com.ABCIndex)
        self.assertIsInstance(pd.Int64Index([1, 2, 3]), com.ABCInt64Index)
        self.assertIsInstance(pd.Float64Index([1, 2, 3]), com.ABCFloat64Index)
        self.assertIsInstance(self.multi_index, com.ABCMultiIndex)
        self.assertIsInstance(self.datetime_index, com.ABCDatetimeIndex)
        self.assertIsInstance(self.timedelta_index, com.ABCTimedeltaIndex)
        self.assertIsInstance(self.period_index, com.ABCPeriodIndex)
        self.assertIsInstance(self.categorical_df.index, com.ABCCategoricalIndex)
        self.assertIsInstance(pd.Index(['a', 'b', 'c']), com.ABCIndexClass)
        self.assertIsInstance(pd.Int64Index([1, 2, 3]), com.ABCIndexClass)
        self.assertIsInstance(pd.Series([1, 2, 3]), com.ABCSeries)
        self.assertIsInstance(self.df, com.ABCDataFrame)
        self.assertIsInstance(self.df.to_panel(), com.ABCPanel)
        self.assertIsInstance(self.sparse_series, com.ABCSparseSeries)
        self.assertIsInstance(self.sparse_array, com.ABCSparseArray)
        self.assertIsInstance(self.categorical, com.ABCCategorical)
        self.assertIsInstance(pd.Period('2012', freq='A-DEC'), com.ABCPeriod)


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
        for s in [tm.makeFloatSeries(),tm.makeStringSeries(),
                  tm.makeObjectSeries(),tm.makeTimeSeries(),tm.makePeriodSeries()]:
            assert(isinstance(isnull(s), Series))

def test_isnull():
    assert not isnull(1.)
    assert isnull(None)
    assert isnull(np.NaN)
    assert not isnull(np.inf)
    assert not isnull(-np.inf)

    # series
    for s in [tm.makeFloatSeries(),tm.makeStringSeries(),
              tm.makeObjectSeries(),tm.makeTimeSeries(),tm.makePeriodSeries()]:
        assert(isinstance(isnull(s), Series))

    # frame
    for df in [tm.makeTimeDataFrame(),tm.makePeriodFrame(),tm.makeMixedDataFrame()]:
        result = isnull(df)
        expected = df.apply(isnull)
        tm.assert_frame_equal(result, expected)

    # panel
    for p in [ tm.makePanel(), tm.makePeriodPanel(), tm.add_nans(tm.makePanel()) ]:
        result = isnull(p)
        expected = p.apply(isnull)
        tm.assert_panel_equal(result, expected)

    # panel 4d
    for p in [ tm.makePanel4D(), tm.add_nans_panel4d(tm.makePanel4D()) ]:
        result = isnull(p)
        expected = p.apply(isnull)
        tm.assert_panel4d_equal(result, expected)

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

    result = isnull([u('foo'), u('bar')])
    assert(not result.any())

def test_isnull_nat():
    result = isnull([NaT])
    exp = np.array([True])
    assert(np.array_equal(result, exp))

    result = isnull(np.array([NaT], dtype=object))
    exp = np.array([True])
    assert(np.array_equal(result, exp))

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

    # GH 9129
    pidx = idx.to_period(freq='M')
    mask = isnull(pidx)
    assert(mask[0])
    assert(not mask[1:].any())

    mask = isnull(pidx[1:])
    assert(not mask.any())


class TestIsNull(tm.TestCase):
    def test_0d_array(self):
        self.assertTrue(isnull(np.array(np.nan)))
        self.assertFalse(isnull(np.array(0.0)))
        self.assertFalse(isnull(np.array(0)))
        # test object dtype
        self.assertTrue(isnull(np.array(np.nan, dtype=object)))
        self.assertFalse(isnull(np.array(0.0, dtype=object)))
        self.assertFalse(isnull(np.array(0, dtype=object)))


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

    expected = np.array([1,2])
    for dtype in [np.float64,object,np.int64]:
        arr = np.array([1.0,2.0],dtype=dtype)
        result = com._possibly_downcast_to_dtype(arr,'infer')
        tm.assert_almost_equal(result, expected)

    expected = np.array([1.0,2.0,np.nan])
    for dtype in [np.float64,object]:
        arr = np.array([1.0,2.0,np.nan],dtype=dtype)
        result = com._possibly_downcast_to_dtype(arr,'infer')
        tm.assert_almost_equal(result, expected)

    # empties
    for dtype in [np.int32,np.float64,np.float32,np.bool_,np.int64,object]:
        arr = np.array([],dtype=dtype)
        result = com._possibly_downcast_to_dtype(arr,'int64')
        tm.assert_almost_equal(result, np.array([],dtype=np.int64))
        assert result.dtype == np.int64

def test_array_equivalent():
    assert array_equivalent(np.array([np.nan, np.nan]),
                            np.array([np.nan, np.nan]))
    assert array_equivalent(np.array([np.nan, 1, np.nan]),
                            np.array([np.nan, 1, np.nan]))
    assert array_equivalent(np.array([np.nan, None], dtype='object'),
                            np.array([np.nan, None], dtype='object'))
    assert array_equivalent(np.array([np.nan, 1+1j], dtype='complex'),
                            np.array([np.nan, 1+1j], dtype='complex'))
    assert not array_equivalent(np.array([np.nan, 1+1j], dtype='complex'),
                                np.array([np.nan, 1+2j], dtype='complex'))
    assert not array_equivalent(np.array([np.nan, 1, np.nan]),
                                np.array([np.nan, 2, np.nan]))
    assert not array_equivalent(np.array(['a', 'b', 'c', 'd']), np.array(['e', 'e']))
    assert array_equivalent(Float64Index([0, np.nan]), Float64Index([0, np.nan]))
    assert not array_equivalent(Float64Index([0, np.nan]), Float64Index([1, np.nan]))
    assert array_equivalent(DatetimeIndex([0, np.nan]), DatetimeIndex([0, np.nan]))
    assert not array_equivalent(DatetimeIndex([0, np.nan]), DatetimeIndex([1, np.nan]))

def test_datetimeindex_from_empty_datetime64_array():
    for unit in [ 'ms', 'us', 'ns' ]:
        idx = DatetimeIndex(np.array([], dtype='datetime64[%s]' % unit))
        assert(len(idx) == 0)


def test_nan_to_nat_conversions():

    df = DataFrame(dict({
        'A' : np.asarray(lrange(10),dtype='float64'),
        'B' : Timestamp('20010101') }))
    df.iloc[3:6,:] = np.nan
    result = df.loc[4,'B'].value
    assert(result == iNaT)

    s = df['B'].copy()
    s._data = s._data.setitem(indexer=tuple([slice(8,9)]),value=np.nan)
    assert(isnull(s[8]))

    # numpy < 1.7.0 is wrong
    from distutils.version import LooseVersion
    if LooseVersion(np.__version__) >= '1.7.0':
        assert(s[8].value == np.datetime64('NaT').astype(np.int64))


def test_any_none():
    assert(com._any_none(1, 2, 3, None))
    assert(not com._any_none(1, 2, 3, 4))


def test_all_not_none():
    assert(com._all_not_none(1, 2, 3, 4))
    assert(not com._all_not_none(1, 2, 3, None))
    assert(not com._all_not_none(None, None, None, None))


def test_repr_binary_type():
    import string
    letters = string.ascii_letters
    btype = compat.binary_type
    try:
        raw = btype(letters, encoding=cf.get_option('display.encoding'))
    except TypeError:
        raw = btype(letters)
    b = compat.text_type(compat.bytes_to_str(raw))
    res = com.pprint_thing(b, quote_strings=True)
    assert_equal(res, repr(b))
    res = com.pprint_thing(b, quote_strings=False)
    assert_equal(res, b)


def test_adjoin():
    data = [['a', 'b', 'c'],
            ['dd', 'ee', 'ff'],
            ['ggg', 'hhh', 'iii']]
    expected = 'a  dd  ggg\nb  ee  hhh\nc  ff  iii'

    adjoined = com.adjoin(2, *data)

    assert(adjoined == expected)



class TestFormattBase(tm.TestCase):

    def test_adjoin(self):
        data = [['a', 'b', 'c'],
                ['dd', 'ee', 'ff'],
                ['ggg', 'hhh', 'iii']]
        expected = 'a  dd  ggg\nb  ee  hhh\nc  ff  iii'

        adjoined = com.adjoin(2, *data)

        self.assertEqual(adjoined, expected)

    def test_adjoin_unicode(self):
        data = [[u'あ', 'b', 'c'],
                ['dd', u'ええ', 'ff'],
                ['ggg', 'hhh', u'いいい']]
        expected = u'あ  dd  ggg\nb  ええ  hhh\nc  ff  いいい'
        adjoined = com.adjoin(2, *data)
        self.assertEqual(adjoined, expected)

        adj = fmt.EastAsianTextAdjustment()

        expected = u"""あ  dd    ggg
b   ええ  hhh
c   ff    いいい"""
        adjoined = adj.adjoin(2, *data)
        self.assertEqual(adjoined, expected)
        cols = adjoined.split('\n')
        self.assertEqual(adj.len(cols[0]), 13)
        self.assertEqual(adj.len(cols[1]), 13)
        self.assertEqual(adj.len(cols[2]), 16)

        expected = u"""あ       dd         ggg
b        ええ       hhh
c        ff         いいい"""
        adjoined = adj.adjoin(7, *data)
        self.assertEqual(adjoined, expected)
        cols = adjoined.split('\n')
        self.assertEqual(adj.len(cols[0]), 23)
        self.assertEqual(adj.len(cols[1]), 23)
        self.assertEqual(adj.len(cols[2]), 26)

    def test_justify(self):
        adj = fmt.EastAsianTextAdjustment()

        def just(x, *args, **kwargs):
            # wrapper to test single str
            return adj.justify([x], *args, **kwargs)[0]

        self.assertEqual(just('abc', 5, mode='left'), 'abc  ')
        self.assertEqual(just('abc', 5, mode='center'), ' abc ')
        self.assertEqual(just('abc', 5, mode='right'), '  abc')
        self.assertEqual(just(u'abc', 5, mode='left'), 'abc  ')
        self.assertEqual(just(u'abc', 5, mode='center'), ' abc ')
        self.assertEqual(just(u'abc', 5, mode='right'), '  abc')

        self.assertEqual(just(u'パンダ', 5, mode='left'), u'パンダ')
        self.assertEqual(just(u'パンダ', 5, mode='center'), u'パンダ')
        self.assertEqual(just(u'パンダ', 5, mode='right'), u'パンダ')

        self.assertEqual(just(u'パンダ', 10, mode='left'), u'パンダ    ')
        self.assertEqual(just(u'パンダ', 10, mode='center'), u'  パンダ  ')
        self.assertEqual(just(u'パンダ', 10, mode='right'), u'    パンダ')

    def test_east_asian_len(self):
        adj = fmt.EastAsianTextAdjustment()

        self.assertEqual(adj.len('abc'), 3)
        self.assertEqual(adj.len(u'abc'), 3)

        self.assertEqual(adj.len(u'パンダ'), 6)
        self.assertEqual(adj.len(u'ﾊﾟﾝﾀﾞ'), 5)
        self.assertEqual(adj.len(u'パンダpanda'), 11)
        self.assertEqual(adj.len(u'ﾊﾟﾝﾀﾞpanda'), 10)


    def test_ambiguous_width(self):
        adj = fmt.EastAsianTextAdjustment()
        self.assertEqual(adj.len(u'¡¡ab'), 4)

        with cf.option_context('display.unicode.ambiguous_as_wide', True):
            adj = fmt.EastAsianTextAdjustment()
            self.assertEqual(adj.len(u'¡¡ab'), 6)

        data = [[u'あ', 'b', 'c'],
                ['dd', u'ええ', 'ff'],
                ['ggg', u'¡¡ab', u'いいい']]
        expected = u'あ  dd    ggg \nb   ええ  ¡¡ab\nc   ff    いいい'
        adjoined = adj.adjoin(2, *data)
        self.assertEqual(adjoined, expected)


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


def test_is_hashable():

    # all new-style classes are hashable by default
    class HashableClass(object):
        pass

    class UnhashableClass1(object):
        __hash__ = None

    class UnhashableClass2(object):
        def __hash__(self):
            raise TypeError("Not hashable")

    hashable = (
        1, 3.14, np.float64(3.14), 'a', tuple(), (1,), HashableClass(),
    )
    not_hashable = (
        [], UnhashableClass1(),
    )
    abc_hashable_not_really_hashable = (
        ([],), UnhashableClass2(),
    )

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
#     if compat.PY3:
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
    passes = (r'a', u('x'), r'asdf', re.compile('adsf'),
              u(r'\u2233\s*'), re.compile(r''))
    fails = 1, [], object()

    for p in passes:
        assert com.is_re_compilable(p)

    for f in fails:
        assert not com.is_re_compilable(f)

def test_random_state():
    import numpy.random as npr
    # Check with seed
    state = com._random_state(5)
    assert_equal(state.uniform(), npr.RandomState(5).uniform())

    # Check with random state object
    state2 = npr.RandomState(10)
    assert_equal(com._random_state(state2).uniform(), npr.RandomState(10).uniform())

    # check with no arg random state
    assert isinstance(com._random_state(), npr.RandomState)

    # Error for floats or strings
    with tm.assertRaises(ValueError):
        com._random_state('test')

    with tm.assertRaises(ValueError):
        com._random_state(5.5)


def test_maybe_match_name():

    matched = com._maybe_match_name(Series([1], name='x'), Series([2], name='x'))
    assert(matched == 'x')

    matched = com._maybe_match_name(Series([1], name='x'), Series([2], name='y'))
    assert(matched is None)

    matched = com._maybe_match_name(Series([1]), Series([2], name='x'))
    assert(matched is None)

    matched = com._maybe_match_name(Series([1], name='x'), Series([2]))
    assert(matched is None)

    matched = com._maybe_match_name(Series([1], name='x'), [2])
    assert(matched == 'x')

    matched = com._maybe_match_name([1], Series([2], name='y'))
    assert(matched == 'y')


class TestTake(tm.TestCase):
    # standard incompatible fill error
    fill_error = re.compile("Incompatible type for fill_value")

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
                with tm.assertRaisesRegexp(TypeError, self.fill_error):
                    com.take_1d(data, indexer, out=out)
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
        def _test_dtype(dtype, can_hold_na, writeable=True):
            data = np.random.randint(0, 2, (5, 3)).astype(dtype)
            data.flags.writeable = writeable

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
                for i, out in enumerate([out0, out1]):
                    with tm.assertRaisesRegexp(TypeError, self.fill_error):
                        com.take_nd(data, indexer, out=out, axis=i)
                    # no exception o/w
                    data.take(indexer, out=out, axis=i)

        for writeable in [True, False]:
            # Check that take_nd works both with writeable arrays (in which
            # case fast typed memoryviews implementation) and read-only
            # arrays alike.
            _test_dtype(np.float64, True, writeable=writeable)
            _test_dtype(np.float32, True, writeable=writeable)
            _test_dtype(np.uint64, False, writeable=writeable)
            _test_dtype(np.uint32, False, writeable=writeable)
            _test_dtype(np.uint16, False, writeable=writeable)
            _test_dtype(np.uint8, False, writeable=writeable)
            _test_dtype(np.int64, False, writeable=writeable)
            _test_dtype(np.int32, False, writeable=writeable)
            _test_dtype(np.int16, False, writeable=writeable)
            _test_dtype(np.int8, False, writeable=writeable)
            _test_dtype(np.object_, True, writeable=writeable)
            _test_dtype(np.bool, False, writeable=writeable)

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
                for i, out in enumerate([out0, out1, out2]):
                    with tm.assertRaisesRegexp(TypeError, self.fill_error):
                        com.take_nd(data, indexer, out=out, axis=i)
                    # no exception o/w
                    data.take(indexer, out=out, axis=i)

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
        self.assert_numpy_array_equal(result, expected)

        result = com.take_1d(arr, [0, 2, -1])
        self.assertEqual(result.dtype, np.object_)

    def test_2d_bool(self):
        arr = np.array([[0, 1, 0],
                        [1, 0, 1],
                        [0, 1, 1]], dtype=bool)

        result = com.take_nd(arr, [0, 2, 2, 1])
        expected = arr.take([0, 2, 2, 1], axis=0)
        self.assert_numpy_array_equal(result, expected)

        result = com.take_nd(arr, [0, 2, 2, 1], axis=1)
        expected = arr.take([0, 2, 2, 1], axis=1)
        self.assert_numpy_array_equal(result, expected)

        result = com.take_nd(arr, [0, 2, -1])
        self.assertEqual(result.dtype, np.object_)

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
        arr = np.random.randint(long(11045376), long(11360736), (5, 3))*100000000000
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

def test_possibly_convert_objects_copy():
    values = np.array([1, 2])

    out = convert._possibly_convert_objects(values, copy=False)
    assert_true(values is out)

    out = convert._possibly_convert_objects(values, copy=True)
    assert_true(values is not out)

    values = np.array(['apply','banana'])
    out = convert._possibly_convert_objects(values, copy=False)
    assert_true(values is out)

    out = convert._possibly_convert_objects(values, copy=True)
    assert_true(values is not out)


def test_dict_compat():
    data_datetime64 = {np.datetime64('1990-03-15'): 1,
                       np.datetime64('2015-03-15'): 2}
    data_unchanged = {1: 2, 3: 4, 5: 6}
    expected = {Timestamp('1990-3-15'): 1, Timestamp('2015-03-15'): 2}
    assert(com._dict_compat(data_datetime64) == expected)
    assert(com._dict_compat(expected) == expected)
    assert(com._dict_compat(data_unchanged) == data_unchanged)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
