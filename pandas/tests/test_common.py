from datetime import datetime
import re

import nose
from nose.tools import assert_equal
import numpy as np
from pandas.tslib import iNaT, NaT
from pandas import Series, DataFrame, date_range, DatetimeIndex, Timestamp
from pandas import compat
from pandas.compat import range, long, lrange, lmap, u
from pandas.core.common import notnull, isnull, array_equivalent
import pandas.core.common as com
import pandas.util.testing as tm
import pandas.core.config as cf
from pandas.core import nanops

_multiprocess_can_split_ = True


def test_mut_exclusive():
    msg = "mutually exclusive arguments: '[ab]' and '[ab]'"
    with tm.assertRaisesRegexp(TypeError, msg):
        com._mut_exclusive(a=1, b=2)
    assert com._mut_exclusive(a=1, b=None) == 1
    assert com._mut_exclusive(major=None, major_axis=None) is None


def test_is_sequence():
    is_seq = com._is_sequence
    assert(is_seq((1, 2)))
    assert(is_seq([1, 2]))
    assert(not is_seq("abcd"))
    assert(not is_seq(u("abcd")))
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
    s._data = s._data.setitem(tuple([slice(8,9)]),np.nan)
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


def test_ensure_int32():
    values = np.arange(10, dtype=np.int32)
    result = com._ensure_int32(values)
    assert(result.dtype == np.int32)

    values = np.arange(10, dtype=np.int64)
    result = com._ensure_int32(values)
    assert(result.dtype == np.int32)


class TestEnsureNumeric(tm.TestCase):
    def test_numeric_values(self):
        # Test integer
        self.assertEqual(nanops._ensure_numeric(1), 1, 'Failed for int')
        # Test float
        self.assertEqual(nanops._ensure_numeric(1.1), 1.1, 'Failed for float')
        # Test complex
        self.assertEqual(nanops._ensure_numeric(1 + 2j), 1 + 2j,
                         'Failed for complex')

    def test_ndarray(self):
        # Test numeric ndarray
        values = np.array([1, 2, 3])
        self.assertTrue(np.allclose(nanops._ensure_numeric(values), values),
                        'Failed for numeric ndarray')

        # Test object ndarray
        o_values = values.astype(object)
        self.assertTrue(np.allclose(nanops._ensure_numeric(o_values), values),
                        'Failed for object ndarray')

        # Test convertible string ndarray
        s_values = np.array(['1', '2', '3'], dtype=object)
        self.assertTrue(np.allclose(nanops._ensure_numeric(s_values), values),
                        'Failed for convertible string ndarray')

        # Test non-convertible string ndarray
        s_values = np.array(['foo', 'bar', 'baz'], dtype=object)
        self.assertRaises(ValueError,
                          lambda: nanops._ensure_numeric(s_values))

    def test_convertable_values(self):
        self.assertTrue(np.allclose(nanops._ensure_numeric('1'), 1.0),
                        'Failed for convertible integer string')
        self.assertTrue(np.allclose(nanops._ensure_numeric('1.1'), 1.1),
                        'Failed for convertible float string')
        self.assertTrue(np.allclose(nanops._ensure_numeric('1+1j'), 1 + 1j),
                        'Failed for convertible complex string')

    def test_non_convertable_values(self):
        self.assertRaises(TypeError,
                          lambda: nanops._ensure_numeric('foo'))
        self.assertRaises(TypeError,
                          lambda: nanops._ensure_numeric({}))
        self.assertRaises(TypeError,
                          lambda: nanops._ensure_numeric([]))


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
                for i, out in enumerate([out0, out1]):
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
        arr = np.random.randint(long(11045376), long(11360736), (5,3))*100000000000
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
