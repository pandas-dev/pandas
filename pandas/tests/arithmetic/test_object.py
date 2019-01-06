# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
# Specifically for object dtype
import operator

import numpy as np
import pytest

import pandas as pd
from pandas import Series, Timestamp
from pandas.core import ops
import pandas.util.testing as tm

# ------------------------------------------------------------------
# Comparisons


class TestObjectComparisons(object):

    def test_comparison_object_numeric_nas(self):
        ser = Series(np.random.randn(10), dtype=object)
        shifted = ser.shift(2)

        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            func = getattr(operator, op)

            result = func(ser, shifted)
            expected = func(ser.astype(float), shifted.astype(float))
            tm.assert_series_equal(result, expected)

    def test_object_comparisons(self):
        ser = Series(['a', 'b', np.nan, 'c', 'a'])

        result = ser == 'a'
        expected = Series([True, False, False, False, True])
        tm.assert_series_equal(result, expected)

        result = ser < 'a'
        expected = Series([False, False, False, False, False])
        tm.assert_series_equal(result, expected)

        result = ser != 'a'
        expected = -(ser == 'a')
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_more_na_comparisons(self, dtype):
        left = Series(['a', np.nan, 'c'], dtype=dtype)
        right = Series(['a', np.nan, 'd'], dtype=dtype)

        result = left == right
        expected = Series([True, False, False])
        tm.assert_series_equal(result, expected)

        result = left != right
        expected = Series([False, True, True])
        tm.assert_series_equal(result, expected)

        result = left == np.nan
        expected = Series([False, False, False])
        tm.assert_series_equal(result, expected)

        result = left != np.nan
        expected = Series([True, True, True])
        tm.assert_series_equal(result, expected)


# ------------------------------------------------------------------
# Arithmetic

class TestTimedeltaNaTArithmetic(object):
    # Tests for arithmetic with np.timedelta64('NaT') which has some tough
    #  corner cases

    def test_tdarr_rfloordiv_nat(self):
        # TODO: test belongs elsewhere, mostly just putting this here because
        #  it is the only TDA method patched for the proof of concept
        td = np.timedelta64('NaT')

        arr = np.arange(3) * 10**9
        tda = pd.TimedeltaIndex(arr)._data

        result = td // tda

        expected = np.array([np.nan, np.nan, np.nan])
        tm.assert_numpy_array_equal(result, expected)

    def test_numeric_with_timedelta_nat(self, box):
        arr = np.array([1, 2, 3, 4], dtype=np.int64)
        obj = tm.box_expected(arr, box)

        td = np.timedelta64('NaT')

        expected = np.array([td] * 4)
        expected = tm.box_expected(expected, box)

        # TODO: RangeIndex
        for dtype in [np.int64, np.float64, np.uint64, object]:
            dobj = obj.astype(dtype)
            for op in [operator.add, operator.sub, ops.radd, ops.rsub]:
                with pytest.raises(TypeError):
                    op(td, dobj)

            if type(dobj) is pd.Index:
                # FIXME: implement these on pd.Index
                continue

            result = dobj * td
            tm.assert_equal(result, expected)
            result = td * dobj
            tm.assert_equal(result, expected)
            result = td / dobj
            tm.assert_equal(result, expected)
            result = td // dobj
            tm.assert_equal(result, expected)
            result = td % dobj
            tm.assert_equal(result, expected)

            # ops that are invalid with tdnat on the right
            with pytest.raises(TypeError):
                dobj / td
            with pytest.raises(TypeError):
                dobj // td
            with pytest.raises(TypeError):
                dobj % td
            with pytest.raises(TypeError):
                divmod(dobj, td)

    @pytest.mark.xfail(reason="I haven't fixed it yet...")
    def test_object_with_timedelta_nat(self, box):
        td = np.timedelta64('NaT')

        arr = np.array([
            pd.offsets.Minute(2),
            pd.Timedelta(hours=2),
            pd.Timedelta(seconds=1).to_pytimedelta(),
            pd.Timedelta(days=3).to_timedelta64()])
        obj = tm.box_expected(arr, box)

        # FIXME: obj + td raises incorrectly
        result = obj + td

        expected = tm.box_expected(np.array([td, td, td, td]))
        tm.assert_equal(result, expected)


class TestArithmetic(object):

    # TODO: parametrize
    def test_pow_ops_object(self):
        # GH#22922
        # pow is weird with masking & 1, so testing here
        a = Series([1, np.nan, 1, np.nan], dtype=object)
        b = Series([1, np.nan, np.nan, 1], dtype=object)
        result = a ** b
        expected = Series(a.values ** b.values, dtype=object)
        tm.assert_series_equal(result, expected)

        result = b ** a
        expected = Series(b.values ** a.values, dtype=object)

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    @pytest.mark.parametrize("other", ["category", "Int64"])
    def test_add_extension_scalar(self, other, box, op):
        # GH#22378
        # Check that scalars satisfying is_extension_array_dtype(obj)
        # do not incorrectly try to dispatch to an ExtensionArray operation

        arr = pd.Series(['a', 'b', 'c'])
        expected = pd.Series([op(x, other) for x in arr])

        arr = tm.box_expected(arr, box)
        expected = tm.box_expected(expected, box)

        result = op(arr, other)
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pytest.param(pd.Index,
                     marks=pytest.mark.xfail(reason="Does not mask nulls",
                                             raises=TypeError)),
        pd.Series,
        pd.DataFrame
    ], ids=lambda x: x.__name__)
    def test_objarr_add_str(self, box):
        ser = pd.Series(['x', np.nan, 'x'])
        expected = pd.Series(['xa', np.nan, 'xa'])

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

        result = ser + 'a'
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pytest.param(pd.Index,
                     marks=pytest.mark.xfail(reason="Does not mask nulls",
                                             raises=TypeError)),
        pd.Series,
        pd.DataFrame
    ], ids=lambda x: x.__name__)
    def test_objarr_radd_str(self, box):
        ser = pd.Series(['x', np.nan, 'x'])
        expected = pd.Series(['ax', np.nan, 'ax'])

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

        result = 'a' + ser
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('data', [
        [1, 2, 3],
        [1.1, 2.2, 3.3],
        [Timestamp('2011-01-01'), Timestamp('2011-01-02'), pd.NaT],
        ['x', 'y', 1]])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_objarr_radd_str_invalid(self, dtype, data, box):
        ser = Series(data, dtype=dtype)

        ser = tm.box_expected(ser, box)
        with pytest.raises(TypeError):
            'foo_' + ser

    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub])
    def test_objarr_add_invalid(self, op, box):
        # invalid ops

        obj_ser = tm.makeObjectSeries()
        obj_ser.name = 'objects'

        obj_ser = tm.box_expected(obj_ser, box)
        with pytest.raises(Exception):
            op(obj_ser, 1)
        with pytest.raises(Exception):
            op(obj_ser, np.array(1, dtype=np.int64))

    # TODO: Moved from tests.series.test_operators; needs cleanup
    def test_operators_na_handling(self):
        ser = Series(['foo', 'bar', 'baz', np.nan])
        result = 'prefix_' + ser
        expected = pd.Series(['prefix_foo', 'prefix_bar',
                              'prefix_baz', np.nan])
        tm.assert_series_equal(result, expected)

        result = ser + '_suffix'
        expected = pd.Series(['foo_suffix', 'bar_suffix',
                              'baz_suffix', np.nan])
        tm.assert_series_equal(result, expected)

    # TODO: parametrize over box
    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_timedelta(self, dtype):
        # note this test is _not_ aimed at timedelta64-dtyped Series
        ser = pd.Series([pd.Timedelta('1 days'), pd.Timedelta('2 days'),
                         pd.Timedelta('3 days')], dtype=dtype)
        expected = pd.Series([pd.Timedelta('4 days'), pd.Timedelta('5 days'),
                              pd.Timedelta('6 days')])

        result = pd.Timedelta('3 days') + ser
        tm.assert_series_equal(result, expected)

        result = ser + pd.Timedelta('3 days')
        tm.assert_series_equal(result, expected)

    # TODO: cleanup & parametrize over box
    def test_mixed_timezone_series_ops_object(self):
        # GH#13043
        ser = pd.Series([pd.Timestamp('2015-01-01', tz='US/Eastern'),
                         pd.Timestamp('2015-01-01', tz='Asia/Tokyo')],
                        name='xxx')
        assert ser.dtype == object

        exp = pd.Series([pd.Timestamp('2015-01-02', tz='US/Eastern'),
                         pd.Timestamp('2015-01-02', tz='Asia/Tokyo')],
                        name='xxx')
        tm.assert_series_equal(ser + pd.Timedelta('1 days'), exp)
        tm.assert_series_equal(pd.Timedelta('1 days') + ser, exp)

        # object series & object series
        ser2 = pd.Series([pd.Timestamp('2015-01-03', tz='US/Eastern'),
                          pd.Timestamp('2015-01-05', tz='Asia/Tokyo')],
                         name='xxx')
        assert ser2.dtype == object
        exp = pd.Series([pd.Timedelta('2 days'), pd.Timedelta('4 days')],
                        name='xxx')
        tm.assert_series_equal(ser2 - ser, exp)
        tm.assert_series_equal(ser - ser2, -exp)

        ser = pd.Series([pd.Timedelta('01:00:00'), pd.Timedelta('02:00:00')],
                        name='xxx', dtype=object)
        assert ser.dtype == object

        exp = pd.Series([pd.Timedelta('01:30:00'), pd.Timedelta('02:30:00')],
                        name='xxx')
        tm.assert_series_equal(ser + pd.Timedelta('00:30:00'), exp)
        tm.assert_series_equal(pd.Timedelta('00:30:00') + ser, exp)
