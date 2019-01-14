# -*- coding: utf-8 -*-
from datetime import datetime, timedelta
import warnings

import numpy as np
import pytest

import pandas.util._test_decorators as td
from pandas.compat import lrange

import pandas as pd
from pandas import (
    Categorical, DataFrame, Index, MultiIndex,
    PeriodIndex, Series, compat, date_range, isna)
from pandas.core import nanops
import pandas.util.testing as tm


# TODO: copied from tests.frame.test_analytics; belongs in pd.util.testing?
def assert_stat_op_calc(opname, alternative, frame, has_skipna=True,
                        check_dtype=True, check_dates=False,
                        check_less_precise=False, skipna_alternative=None):
    """
    Check that operator opname works as advertised on frame

    Parameters
    ----------
    opname : string
        Name of the operator to test on frame
    alternative : function
        Function that opname is tested against; i.e. "frame.opname()" should
        equal "alternative(frame)".
    frame : DataFrame
        The object that the tests are executed on
    has_skipna : bool, default True
        Whether the method "opname" has the kwarg "skip_na"
    check_dtype : bool, default True
        Whether the dtypes of the result of "frame.opname()" and
        "alternative(frame)" should be checked.
    check_dates : bool, default false
        Whether opname should be tested on a Datetime Series
    check_less_precise : bool, default False
        Whether results should only be compared approximately;
        passed on to tm.assert_series_equal
    skipna_alternative : function, default None
        NaN-safe version of alternative
    """

    f = getattr(frame, opname)

    if check_dates:
        df = DataFrame({'b': date_range('1/1/2001', periods=2)})
        result = getattr(df, opname)()
        assert isinstance(result, Series)

        df['a'] = lrange(len(df))
        result = getattr(df, opname)()
        assert isinstance(result, Series)
        assert len(result)

    if has_skipna:
        def wrapper(x):
            return alternative(x.values)

        skipna_wrapper = tm._make_skipna_wrapper(alternative,
                                                 skipna_alternative)
        result0 = f(axis=0, skipna=False)
        result1 = f(axis=1, skipna=False)
        tm.assert_series_equal(result0, frame.apply(wrapper),
                               check_dtype=check_dtype,
                               check_less_precise=check_less_precise)
        # HACK: win32
        tm.assert_series_equal(result1, frame.apply(wrapper, axis=1),
                               check_dtype=False,
                               check_less_precise=check_less_precise)
    else:
        skipna_wrapper = alternative

    result0 = f(axis=0)
    result1 = f(axis=1)
    tm.assert_series_equal(result0, frame.apply(skipna_wrapper),
                           check_dtype=check_dtype,
                           check_less_precise=check_less_precise)

    if opname in ['sum', 'prod']:
        expected = frame.apply(skipna_wrapper, axis=1)
        tm.assert_series_equal(result1, expected, check_dtype=False,
                               check_less_precise=check_less_precise)

    # check dtypes
    if check_dtype:
        lcd_dtype = frame.values.dtype
        assert lcd_dtype == result0.dtype
        assert lcd_dtype == result1.dtype

    # bad axis
    with pytest.raises(ValueError, match='No axis named 2'):
        f(axis=2)

    # all NA case
    if has_skipna:
        all_na = frame * np.NaN
        r0 = getattr(all_na, opname)(axis=0)
        r1 = getattr(all_na, opname)(axis=1)
        if opname in ['sum', 'prod']:
            unit = 1 if opname == 'prod' else 0  # result for empty sum/prod
            expected = pd.Series(unit, index=r0.index, dtype=r0.dtype)
            tm.assert_series_equal(r0, expected)
            expected = pd.Series(unit, index=r1.index, dtype=r1.dtype)
            tm.assert_series_equal(r1, expected)


# TODO: copied from tests.frame.test_analytics; belongs in pd.util.testing?
def assert_stat_op_api(opname, float_frame, float_string_frame,
                       has_numeric_only=False):
    """
    Check that API for operator opname works as advertised on frame

    Parameters
    ----------
    opname : string
        Name of the operator to test on frame
    float_frame : DataFrame
        DataFrame with columns of type float
    float_string_frame : DataFrame
        DataFrame with both float and string columns
    has_numeric_only : bool, default False
        Whether the method "opname" has the kwarg "numeric_only"
    """

    # make sure works on mixed-type frame
    getattr(float_string_frame, opname)(axis=0)
    getattr(float_string_frame, opname)(axis=1)

    if has_numeric_only:
        getattr(float_string_frame, opname)(axis=0, numeric_only=True)
        getattr(float_string_frame, opname)(axis=1, numeric_only=True)
        getattr(float_frame, opname)(axis=0, numeric_only=False)
        getattr(float_frame, opname)(axis=1, numeric_only=False)


def assert_bool_op_calc(opname, alternative, frame, has_skipna=True):
    """
    Check that bool operator opname works as advertised on frame

    Parameters
    ----------
    opname : string
        Name of the operator to test on frame
    alternative : function
        Function that opname is tested against; i.e. "frame.opname()" should
        equal "alternative(frame)".
    frame : DataFrame
        The object that the tests are executed on
    has_skipna : bool, default True
        Whether the method "opname" has the kwarg "skip_na"
    """

    f = getattr(frame, opname)

    if has_skipna:
        def skipna_wrapper(x):
            nona = x.dropna().values
            return alternative(nona)

        def wrapper(x):
            return alternative(x.values)

        result0 = f(axis=0, skipna=False)
        result1 = f(axis=1, skipna=False)

        tm.assert_series_equal(result0, frame.apply(wrapper))
        tm.assert_series_equal(result1, frame.apply(wrapper, axis=1),
                               check_dtype=False)  # HACK: win32
    else:
        skipna_wrapper = alternative
        wrapper = alternative

    result0 = f(axis=0)
    result1 = f(axis=1)

    tm.assert_series_equal(result0, frame.apply(skipna_wrapper))
    tm.assert_series_equal(result1, frame.apply(skipna_wrapper, axis=1),
                           check_dtype=False)

    # bad axis
    with pytest.raises(ValueError, match='No axis named 2'):
        f(axis=2)

    # all NA case
    if has_skipna:
        all_na = frame * np.NaN
        r0 = getattr(all_na, opname)(axis=0)
        r1 = getattr(all_na, opname)(axis=1)
        if opname == 'any':
            assert not r0.any()
            assert not r1.any()
        else:
            assert r0.all()
            assert r1.all()


def assert_bool_op_api(opname, bool_frame_with_na, float_string_frame,
                       has_bool_only=False):
    """
    Check that API for boolean operator opname works as advertised on frame

    Parameters
    ----------
    opname : string
        Name of the operator to test on frame
    float_frame : DataFrame
        DataFrame with columns of type float
    float_string_frame : DataFrame
        DataFrame with both float and string columns
    has_bool_only : bool, default False
        Whether the method "opname" has the kwarg "bool_only"
    """
    # make sure op works on mixed-type frame
    mixed = float_string_frame
    mixed['_bool_'] = np.random.randn(len(mixed)) > 0.5
    getattr(mixed, opname)(axis=0)
    getattr(mixed, opname)(axis=1)

    if has_bool_only:
        getattr(mixed, opname)(axis=0, bool_only=True)
        getattr(mixed, opname)(axis=1, bool_only=True)
        getattr(bool_frame_with_na, opname)(axis=0, bool_only=False)
        getattr(bool_frame_with_na, opname)(axis=1, bool_only=False)


def get_float_frame_with_na():
    """
    Fixture for DataFrame of floats with index of unique strings

    Columns are ['A', 'B', 'C', 'D']; some entries are missing
    """
    df = DataFrame(tm.getSeriesData())
    # set some NAs
    df.loc[5:10] = np.nan
    df.loc[15:20, -2:] = np.nan
    return df


def get_float_string_frame():
    """
    Fixture for DataFrame of floats and strings with index of unique strings

    Columns are ['A', 'B', 'C', 'D', 'foo'].
    """
    df = DataFrame(tm.getSeriesData())
    df['foo'] = 'bar'
    return df


def get_int_frame():
    """
    Fixture for DataFrame of ints with index of unique strings

    Columns are ['A', 'B', 'C', 'D']
    """
    df = DataFrame({k: v.astype(int)
                   for k, v in compat.iteritems(tm.getSeriesData())})
    # force these all to int64 to avoid platform testing issues
    return DataFrame({c: s for c, s in compat.iteritems(df)}, dtype=np.int64)


def get_bool_frame_with_na():
    """
    Fixture for DataFrame of booleans with index of unique strings

    Columns are ['A', 'B', 'C', 'D']; some entries are missing
    """
    df = DataFrame(tm.getSeriesData()) > 0
    df = df.astype(object)
    # set some NAs
    df.loc[5:10] = np.nan
    df.loc[15:20, -2:] = np.nan
    return df


def get_objs():
    indexes = [
        tm.makeBoolIndex(10, name='a'),
        tm.makeIntIndex(10, name='a'),
        tm.makeFloatIndex(10, name='a'),
        tm.makeDateIndex(10, name='a'),
        tm.makeDateIndex(10, name='a').tz_localize(tz='US/Eastern'),
        tm.makePeriodIndex(10, name='a'),
        tm.makeStringIndex(10, name='a'),
        tm.makeUnicodeIndex(10, name='a')
    ]

    arr = np.random.randn(10)
    series = [Series(arr, index=idx, name='a') for idx in indexes]

    objs = indexes + series
    return objs


objs = get_objs()


class TestReductions(object):

    @pytest.mark.parametrize('opname', ['max', 'min'])
    @pytest.mark.parametrize('obj', objs)
    def test_ops(self, opname, obj):
        result = getattr(obj, opname)()
        if not isinstance(obj, PeriodIndex):
            expected = getattr(obj.values, opname)()
        else:
            expected = pd.Period(
                ordinal=getattr(obj._ndarray_values, opname)(),
                freq=obj.freq)
        try:
            assert result == expected
        except TypeError:
            # comparing tz-aware series with np.array results in
            # TypeError
            expected = expected.astype('M8[ns]').astype('int64')
            assert result.value == expected

    def test_nanops(self):
        # GH#7261
        for opname in ['max', 'min']:
            for klass in [Index, Series]:
                arg_op = 'arg' + opname if klass is Index else 'idx' + opname

                obj = klass([np.nan, 2.0])
                assert getattr(obj, opname)() == 2.0

                obj = klass([np.nan])
                assert pd.isna(getattr(obj, opname)())
                assert pd.isna(getattr(obj, opname)(skipna=False))

                obj = klass([])
                assert pd.isna(getattr(obj, opname)())
                assert pd.isna(getattr(obj, opname)(skipna=False))

                obj = klass([pd.NaT, datetime(2011, 11, 1)])
                # check DatetimeIndex monotonic path
                assert getattr(obj, opname)() == datetime(2011, 11, 1)
                assert getattr(obj, opname)(skipna=False) is pd.NaT

                assert getattr(obj, arg_op)() == 1
                result = getattr(obj, arg_op)(skipna=False)
                if klass is Series:
                    assert np.isnan(result)
                else:
                    assert result == -1

                obj = klass([pd.NaT, datetime(2011, 11, 1), pd.NaT])
                # check DatetimeIndex non-monotonic path
                assert getattr(obj, opname)(), datetime(2011, 11, 1)
                assert getattr(obj, opname)(skipna=False) is pd.NaT

                assert getattr(obj, arg_op)() == 1
                result = getattr(obj, arg_op)(skipna=False)
                if klass is Series:
                    assert np.isnan(result)
                else:
                    assert result == -1

                for dtype in ["M8[ns]", "datetime64[ns, UTC]"]:
                    # cases with empty Series/DatetimeIndex
                    obj = klass([], dtype=dtype)

                    assert getattr(obj, opname)() is pd.NaT
                    assert getattr(obj, opname)(skipna=False) is pd.NaT

                    with pytest.raises(ValueError, match="empty sequence"):
                        getattr(obj, arg_op)()
                    with pytest.raises(ValueError, match="empty sequence"):
                        getattr(obj, arg_op)(skipna=False)

        # argmin/max
        obj = Index(np.arange(5, dtype='int64'))
        assert obj.argmin() == 0
        assert obj.argmax() == 4

        obj = Index([np.nan, 1, np.nan, 2])
        assert obj.argmin() == 1
        assert obj.argmax() == 3
        assert obj.argmin(skipna=False) == -1
        assert obj.argmax(skipna=False) == -1

        obj = Index([np.nan])
        assert obj.argmin() == -1
        assert obj.argmax() == -1
        assert obj.argmin(skipna=False) == -1
        assert obj.argmax(skipna=False) == -1

        obj = Index([pd.NaT, datetime(2011, 11, 1), datetime(2011, 11, 2),
                     pd.NaT])
        assert obj.argmin() == 1
        assert obj.argmax() == 2
        assert obj.argmin(skipna=False) == -1
        assert obj.argmax(skipna=False) == -1

        obj = Index([pd.NaT])
        assert obj.argmin() == -1
        assert obj.argmax() == -1
        assert obj.argmin(skipna=False) == -1
        assert obj.argmax(skipna=False) == -1


class TestSeriesReductions(object):
    # Note: the name TestSeriesReductions indicates these tests
    #  were moved from a series-specific test file, _not_ that these tests are
    #  intended long-term to be series-specific

    def test_sum_inf(self):
        s = Series(np.random.randn(10))
        s2 = s.copy()

        s[5:8] = np.inf
        s2[5:8] = np.nan

        assert np.isinf(s.sum())

        arr = np.random.randn(100, 100).astype('f4')
        arr[:, 2] = np.inf

        with pd.option_context("mode.use_inf_as_na", True):
            tm.assert_almost_equal(s.sum(), s2.sum())

        res = nanops.nansum(arr, axis=1)
        assert np.isinf(res).all()

    @pytest.mark.parametrize("use_bottleneck", [True, False])
    @pytest.mark.parametrize("method, unit", [
        ("sum", 0.0),
        ("prod", 1.0)
    ])
    def test_empty(self, method, unit, use_bottleneck):
        with pd.option_context("use_bottleneck", use_bottleneck):
            # GH#9422 / GH#18921
            # Entirely empty
            s = Series([])
            # NA by default
            result = getattr(s, method)()
            assert result == unit

            # Explicit
            result = getattr(s, method)(min_count=0)
            assert result == unit

            result = getattr(s, method)(min_count=1)
            assert pd.isna(result)

            # Skipna, default
            result = getattr(s, method)(skipna=True)
            result == unit

            # Skipna, explicit
            result = getattr(s, method)(skipna=True, min_count=0)
            assert result == unit

            result = getattr(s, method)(skipna=True, min_count=1)
            assert pd.isna(result)

            # All-NA
            s = Series([np.nan])
            # NA by default
            result = getattr(s, method)()
            assert result == unit

            # Explicit
            result = getattr(s, method)(min_count=0)
            assert result == unit

            result = getattr(s, method)(min_count=1)
            assert pd.isna(result)

            # Skipna, default
            result = getattr(s, method)(skipna=True)
            result == unit

            # skipna, explicit
            result = getattr(s, method)(skipna=True, min_count=0)
            assert result == unit

            result = getattr(s, method)(skipna=True, min_count=1)
            assert pd.isna(result)

            # Mix of valid, empty
            s = Series([np.nan, 1])
            # Default
            result = getattr(s, method)()
            assert result == 1.0

            # Explicit
            result = getattr(s, method)(min_count=0)
            assert result == 1.0

            result = getattr(s, method)(min_count=1)
            assert result == 1.0

            # Skipna
            result = getattr(s, method)(skipna=True)
            assert result == 1.0

            result = getattr(s, method)(skipna=True, min_count=0)
            assert result == 1.0

            result = getattr(s, method)(skipna=True, min_count=1)
            assert result == 1.0

            # GH#844 (changed in GH#9422)
            df = DataFrame(np.empty((10, 0)))
            assert (getattr(df, method)(1) == unit).all()

            s = pd.Series([1])
            result = getattr(s, method)(min_count=2)
            assert pd.isna(result)

            s = pd.Series([np.nan])
            result = getattr(s, method)(min_count=2)
            assert pd.isna(result)

            s = pd.Series([np.nan, 1])
            result = getattr(s, method)(min_count=2)
            assert pd.isna(result)

    @pytest.mark.parametrize('method, unit', [
        ('sum', 0.0),
        ('prod', 1.0),
    ])
    def test_empty_multi(self, method, unit):
        s = pd.Series([1, np.nan, np.nan, np.nan],
                      index=MultiIndex.from_product([('a', 'b'), (0, 1)]))
        # 1 / 0 by default
        result = getattr(s, method)(level=0)
        expected = pd.Series([1, unit], index=['a', 'b'])
        tm.assert_series_equal(result, expected)

        # min_count=0
        result = getattr(s, method)(level=0, min_count=0)
        expected = pd.Series([1, unit], index=['a', 'b'])
        tm.assert_series_equal(result, expected)

        # min_count=1
        result = getattr(s, method)(level=0, min_count=1)
        expected = pd.Series([1, np.nan], index=['a', 'b'])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "method", ['mean', 'median', 'std', 'var'])
    def test_ops_consistency_on_empty(self, method):

        # GH#7869
        # consistency on empty

        # float
        result = getattr(Series(dtype=float), method)()
        assert pd.isna(result)

        # timedelta64[ns]
        result = getattr(Series(dtype='m8[ns]'), method)()
        assert result is pd.NaT

    def test_nansum_buglet(self):
        ser = Series([1.0, np.nan], index=[0, 1])
        result = np.nansum(ser)
        tm.assert_almost_equal(result, 1)

    @pytest.mark.parametrize("use_bottleneck", [True, False])
    def test_sum_overflow(self, use_bottleneck):

        with pd.option_context('use_bottleneck', use_bottleneck):
            # GH#6915
            # overflowing on the smaller int dtypes
            for dtype in ['int32', 'int64']:
                v = np.arange(5000000, dtype=dtype)
                s = Series(v)

                result = s.sum(skipna=False)
                assert int(result) == v.sum(dtype='int64')
                result = s.min(skipna=False)
                assert int(result) == 0
                result = s.max(skipna=False)
                assert int(result) == v[-1]

            for dtype in ['float32', 'float64']:
                v = np.arange(5000000, dtype=dtype)
                s = Series(v)

                result = s.sum(skipna=False)
                assert result == v.sum(dtype=dtype)
                result = s.min(skipna=False)
                assert np.allclose(float(result), 0.0)
                result = s.max(skipna=False)
                assert np.allclose(float(result), v[-1])

    def test_empty_timeseries_reductions_return_nat(self):
        # covers GH#11245
        for dtype in ('m8[ns]', 'm8[ns]', 'M8[ns]', 'M8[ns, UTC]'):
            assert Series([], dtype=dtype).min() is pd.NaT
            assert Series([], dtype=dtype).max() is pd.NaT
            assert Series([], dtype=dtype).min(skipna=False) is pd.NaT
            assert Series([], dtype=dtype).max(skipna=False) is pd.NaT

    def test_numpy_argmin_deprecated(self):
        # See GH#16830
        data = np.arange(1, 11)

        s = Series(data, index=data)
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            # The deprecation of Series.argmin also causes a deprecation
            # warning when calling np.argmin. This behavior is temporary
            # until the implementation of Series.argmin is corrected.
            result = np.argmin(s)

        assert result == 1

        with tm.assert_produces_warning(FutureWarning):
            # argmin is aliased to idxmin
            result = s.argmin()

        assert result == 1

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            msg = "the 'out' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                np.argmin(s, out=data)

    def test_numpy_argmax_deprecated(self):
        # See GH#16830
        data = np.arange(1, 11)

        s = Series(data, index=data)
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            # The deprecation of Series.argmax also causes a deprecation
            # warning when calling np.argmax. This behavior is temporary
            # until the implementation of Series.argmax is corrected.
            result = np.argmax(s)
        assert result == 10

        with tm.assert_produces_warning(FutureWarning):
            # argmax is aliased to idxmax
            result = s.argmax()

        assert result == 10

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            msg = "the 'out' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                np.argmax(s, out=data)

    def test_idxmin(self):
        # test idxmin
        # _check_stat_op approach can not be used here because of isna check.
        string_series = tm.makeStringSeries().rename('series')

        # add some NaNs
        string_series[5:15] = np.NaN

        # skipna or no
        assert string_series[string_series.idxmin()] == string_series.min()
        assert pd.isna(string_series.idxmin(skipna=False))

        # no NaNs
        nona = string_series.dropna()
        assert nona[nona.idxmin()] == nona.min()
        assert (nona.index.values.tolist().index(nona.idxmin()) ==
                nona.values.argmin())

        # all NaNs
        allna = string_series * np.nan
        assert pd.isna(allna.idxmin())

        # datetime64[ns]
        s = Series(date_range('20130102', periods=6))
        result = s.idxmin()
        assert result == 0

        s[0] = np.nan
        result = s.idxmin()
        assert result == 1

    def test_idxmax(self):
        # test idxmax
        # _check_stat_op approach can not be used here because of isna check.
        string_series = tm.makeStringSeries().rename('series')

        # add some NaNs
        string_series[5:15] = np.NaN

        # skipna or no
        assert string_series[string_series.idxmax()] == string_series.max()
        assert pd.isna(string_series.idxmax(skipna=False))

        # no NaNs
        nona = string_series.dropna()
        assert nona[nona.idxmax()] == nona.max()
        assert (nona.index.values.tolist().index(nona.idxmax()) ==
                nona.values.argmax())

        # all NaNs
        allna = string_series * np.nan
        assert pd.isna(allna.idxmax())

        s = Series(date_range('20130102', periods=6))
        result = s.idxmax()
        assert result == 5

        s[5] = np.nan
        result = s.idxmax()
        assert result == 4

        # Float64Index
        # GH#5914
        s = pd.Series([1, 2, 3], [1.1, 2.1, 3.1])
        result = s.idxmax()
        assert result == 3.1
        result = s.idxmin()
        assert result == 1.1

        s = pd.Series(s.index, s.index)
        result = s.idxmax()
        assert result == 3.1
        result = s.idxmin()
        assert result == 1.1

    def test_all_any(self):
        ts = tm.makeTimeSeries()
        bool_series = ts > 0
        assert not bool_series.all()
        assert bool_series.any()

        # Alternative types, with implicit 'object' dtype.
        s = Series(['abc', True])
        assert 'abc' == s.any()  # 'abc' || True => 'abc'

    def test_all_any_params(self):
        # Check skipna, with implicit 'object' dtype.
        s1 = Series([np.nan, True])
        s2 = Series([np.nan, False])
        assert s1.all(skipna=False)  # nan && True => True
        assert s1.all(skipna=True)
        assert np.isnan(s2.any(skipna=False))  # nan || False => nan
        assert not s2.any(skipna=True)

        # Check level.
        s = pd.Series([False, False, True, True, False, True],
                      index=[0, 0, 1, 1, 2, 2])
        tm.assert_series_equal(s.all(level=0), Series([False, True, False]))
        tm.assert_series_equal(s.any(level=0), Series([False, True, True]))

        # bool_only is not implemented with level option.
        with pytest.raises(NotImplementedError):
            s.any(bool_only=True, level=0)
        with pytest.raises(NotImplementedError):
            s.all(bool_only=True, level=0)

        # bool_only is not implemented alone.
        with pytest.raises(NotImplementedError):
            s.any(bool_only=True,)
        with pytest.raises(NotImplementedError):
            s.all(bool_only=True)

    def test_timedelta64_analytics(self):

        # index min/max
        dti = date_range('2012-1-1', periods=3, freq='D')
        td = Series(dti) - pd.Timestamp('20120101')

        result = td.idxmin()
        assert result == 0

        result = td.idxmax()
        assert result == 2

        # GH#2982
        # with NaT
        td[0] = np.nan

        result = td.idxmin()
        assert result == 1

        result = td.idxmax()
        assert result == 2

        # abs
        s1 = Series(date_range('20120101', periods=3))
        s2 = Series(date_range('20120102', periods=3))
        expected = Series(s2 - s1)

        # FIXME: don't leave commented-out code
        # this fails as numpy returns timedelta64[us]
        # result = np.abs(s1-s2)
        # assert_frame_equal(result,expected)

        result = (s1 - s2).abs()
        tm.assert_series_equal(result, expected)

        # max/min
        result = td.max()
        expected = pd.Timedelta('2 days')
        assert result == expected

        result = td.min()
        expected = pd.Timedelta('1 days')
        assert result == expected

    @pytest.mark.parametrize(
        "test_input,error_type",
        [
            (pd.Series([]), ValueError),

            # For strings, or any Series with dtype 'O'
            (pd.Series(['foo', 'bar', 'baz']), TypeError),
            (pd.Series([(1,), (2,)]), TypeError),

            # For mixed data types
            (
                pd.Series(['foo', 'foo', 'bar', 'bar', None, np.nan, 'baz']),
                TypeError
            ),
        ]
    )
    def test_assert_idxminmax_raises(self, test_input, error_type):
        """
        Cases where ``Series.argmax`` and related should raise an exception
        """
        with pytest.raises(error_type):
            test_input.idxmin()
        with pytest.raises(error_type):
            test_input.idxmin(skipna=False)
        with pytest.raises(error_type):
            test_input.idxmax()
        with pytest.raises(error_type):
            test_input.idxmax(skipna=False)

    def test_idxminmax_with_inf(self):
        # For numeric data with NA and Inf (GH #13595)
        s = pd.Series([0, -np.inf, np.inf, np.nan])

        assert s.idxmin() == 1
        assert np.isnan(s.idxmin(skipna=False))

        assert s.idxmax() == 2
        assert np.isnan(s.idxmax(skipna=False))

        # Using old-style behavior that treats floating point nan, -inf, and
        # +inf as missing
        with pd.option_context('mode.use_inf_as_na', True):
            assert s.idxmin() == 0
            assert np.isnan(s.idxmin(skipna=False))
            assert s.idxmax() == 0
            np.isnan(s.idxmax(skipna=False))


class TestDatetime64SeriesReductions(object):
    # Note: the name TestDatetime64SeriesReductions indicates these tests
    #  were moved from a series-specific test file, _not_ that these tests are
    #  intended long-term to be series-specific

    @pytest.mark.parametrize('nat_ser', [
        Series([pd.NaT, pd.NaT]),
        Series([pd.NaT, pd.Timedelta('nat')]),
        Series([pd.Timedelta('nat'), pd.Timedelta('nat')])])
    def test_minmax_nat_series(self, nat_ser):
        # GH#23282
        assert nat_ser.min() is pd.NaT
        assert nat_ser.max() is pd.NaT
        assert nat_ser.min(skipna=False) is pd.NaT
        assert nat_ser.max(skipna=False) is pd.NaT

    @pytest.mark.parametrize('nat_df', [
        pd.DataFrame([pd.NaT, pd.NaT]),
        pd.DataFrame([pd.NaT, pd.Timedelta('nat')]),
        pd.DataFrame([pd.Timedelta('nat'), pd.Timedelta('nat')])])
    def test_minmax_nat_dataframe(self, nat_df):
        # GH#23282
        assert nat_df.min()[0] is pd.NaT
        assert nat_df.max()[0] is pd.NaT
        assert nat_df.min(skipna=False)[0] is pd.NaT
        assert nat_df.max(skipna=False)[0] is pd.NaT

    def test_min_max(self):
        rng = date_range('1/1/2000', '12/31/2000')
        rng2 = rng.take(np.random.permutation(len(rng)))

        the_min = rng2.min()
        the_max = rng2.max()
        assert isinstance(the_min, pd.Timestamp)
        assert isinstance(the_max, pd.Timestamp)
        assert the_min == rng[0]
        assert the_max == rng[-1]

        assert rng.min() == rng[0]
        assert rng.max() == rng[-1]

    def test_min_max_series(self):
        rng = date_range('1/1/2000', periods=10, freq='4h')
        lvls = ['A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C', 'C']
        df = DataFrame({'TS': rng, 'V': np.random.randn(len(rng)), 'L': lvls})

        result = df.TS.max()
        exp = pd.Timestamp(df.TS.iat[-1])
        assert isinstance(result, pd.Timestamp)
        assert result == exp

        result = df.TS.min()
        exp = pd.Timestamp(df.TS.iat[0])
        assert isinstance(result, pd.Timestamp)
        assert result == exp


class TestCategoricalSeriesReductions(object):
    # Note: the name TestCategoricalSeriesReductions indicates these tests
    #  were moved from a series-specific test file, _not_ that these tests are
    #  intended long-term to be series-specific

    def test_min_max(self):
        # unordered cats have no min/max
        cat = Series(Categorical(["a", "b", "c", "d"], ordered=False))
        with pytest.raises(TypeError):
            cat.min()
        with pytest.raises(TypeError):
            cat.max()

        cat = Series(Categorical(["a", "b", "c", "d"], ordered=True))
        _min = cat.min()
        _max = cat.max()
        assert _min == "a"
        assert _max == "d"

        cat = Series(Categorical(["a", "b", "c", "d"], categories=[
                     'd', 'c', 'b', 'a'], ordered=True))
        _min = cat.min()
        _max = cat.max()
        assert _min == "d"
        assert _max == "a"

        cat = Series(Categorical(
            [np.nan, "b", "c", np.nan], categories=['d', 'c', 'b', 'a'
                                                    ], ordered=True))
        _min = cat.min()
        _max = cat.max()
        assert np.isnan(_min)
        assert _max == "b"

        cat = Series(Categorical(
            [np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1], ordered=True))
        _min = cat.min()
        _max = cat.max()
        assert np.isnan(_min)
        assert _max == 1


class TestSeriesMode(object):
    # Note: the name TestSeriesMode indicates these tests
    #  were moved from a series-specific test file, _not_ that these tests are
    #  intended long-term to be series-specific

    @pytest.mark.parametrize('dropna, expected', [
        (True, Series([], dtype=np.float64)),
        (False, Series([], dtype=np.float64))
    ])
    def test_mode_empty(self, dropna, expected):
        s = Series([], dtype=np.float64)
        result = s.mode(dropna)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dropna, data, expected', [
        (True, [1, 1, 1, 2], [1]),
        (True, [1, 1, 1, 2, 3, 3, 3], [1, 3]),
        (False, [1, 1, 1, 2], [1]),
        (False, [1, 1, 1, 2, 3, 3, 3], [1, 3]),
    ])
    @pytest.mark.parametrize(
        'dt',
        list(np.typecodes['AllInteger'] + np.typecodes['Float'])
    )
    def test_mode_numerical(self, dropna, data, expected, dt):
        s = Series(data, dtype=dt)
        result = s.mode(dropna)
        expected = Series(expected, dtype=dt)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dropna, expected', [
        (True, [1.0]),
        (False, [1, np.nan]),
    ])
    def test_mode_numerical_nan(self, dropna, expected):
        s = Series([1, 1, 2, np.nan, np.nan])
        result = s.mode(dropna)
        expected = Series(expected)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dropna, expected1, expected2, expected3', [
        (True, ['b'], ['bar'], ['nan']),
        (False, ['b'], [np.nan], ['nan'])
    ])
    def test_mode_str_obj(self, dropna, expected1, expected2, expected3):
        # Test string and object types.
        data = ['a'] * 2 + ['b'] * 3

        s = Series(data, dtype='c')
        result = s.mode(dropna)
        expected1 = Series(expected1, dtype='c')
        tm.assert_series_equal(result, expected1)

        data = ['foo', 'bar', 'bar', np.nan, np.nan, np.nan]

        s = Series(data, dtype=object)
        result = s.mode(dropna)
        expected2 = Series(expected2, dtype=object)
        tm.assert_series_equal(result, expected2)

        data = ['foo', 'bar', 'bar', np.nan, np.nan, np.nan]

        s = Series(data, dtype=object).astype(str)
        result = s.mode(dropna)
        expected3 = Series(expected3, dtype=str)
        tm.assert_series_equal(result, expected3)

    @pytest.mark.parametrize('dropna, expected1, expected2', [
        (True, ['foo'], ['foo']),
        (False, ['foo'], [np.nan])
    ])
    def test_mode_mixeddtype(self, dropna, expected1, expected2):
        s = Series([1, 'foo', 'foo'])
        result = s.mode(dropna)
        expected = Series(expected1)
        tm.assert_series_equal(result, expected)

        s = Series([1, 'foo', 'foo', np.nan, np.nan, np.nan])
        result = s.mode(dropna)
        expected = Series(expected2, dtype=object)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dropna, expected1, expected2', [
        (True, ['1900-05-03', '2011-01-03', '2013-01-02'],
               ['2011-01-03', '2013-01-02']),
        (False, [np.nan], [np.nan, '2011-01-03', '2013-01-02']),
    ])
    def test_mode_datetime(self, dropna, expected1, expected2):
        s = Series(['2011-01-03', '2013-01-02',
                    '1900-05-03', 'nan', 'nan'], dtype='M8[ns]')
        result = s.mode(dropna)
        expected1 = Series(expected1, dtype='M8[ns]')
        tm.assert_series_equal(result, expected1)

        s = Series(['2011-01-03', '2013-01-02', '1900-05-03',
                    '2011-01-03', '2013-01-02', 'nan', 'nan'],
                   dtype='M8[ns]')
        result = s.mode(dropna)
        expected2 = Series(expected2, dtype='M8[ns]')
        tm.assert_series_equal(result, expected2)

    @pytest.mark.parametrize('dropna, expected1, expected2', [
        (True, ['-1 days', '0 days', '1 days'], ['2 min', '1 day']),
        (False, [np.nan], [np.nan, '2 min', '1 day']),
    ])
    def test_mode_timedelta(self, dropna, expected1, expected2):
        # gh-5986: Test timedelta types.

        s = Series(['1 days', '-1 days', '0 days', 'nan', 'nan'],
                   dtype='timedelta64[ns]')
        result = s.mode(dropna)
        expected1 = Series(expected1, dtype='timedelta64[ns]')
        tm.assert_series_equal(result, expected1)

        s = Series(['1 day', '1 day', '-1 day', '-1 day 2 min',
                    '2 min', '2 min', 'nan', 'nan'],
                   dtype='timedelta64[ns]')
        result = s.mode(dropna)
        expected2 = Series(expected2, dtype='timedelta64[ns]')
        tm.assert_series_equal(result, expected2)

    @pytest.mark.parametrize('dropna, expected1, expected2, expected3', [
        (True, Categorical([1, 2], categories=[1, 2]),
         Categorical(['a'], categories=[1, 'a']),
         Categorical([3, 1], categories=[3, 2, 1], ordered=True)),
        (False, Categorical([np.nan], categories=[1, 2]),
         Categorical([np.nan, 'a'], categories=[1, 'a']),
         Categorical([np.nan, 3, 1], categories=[3, 2, 1], ordered=True)),
    ])
    def test_mode_category(self, dropna, expected1, expected2, expected3):
        s = Series(Categorical([1, 2, np.nan, np.nan]))
        result = s.mode(dropna)
        expected1 = Series(expected1, dtype='category')
        tm.assert_series_equal(result, expected1)

        s = Series(Categorical([1, 'a', 'a', np.nan, np.nan]))
        result = s.mode(dropna)
        expected2 = Series(expected2, dtype='category')
        tm.assert_series_equal(result, expected2)

        s = Series(Categorical([1, 1, 2, 3, 3, np.nan, np.nan],
                               categories=[3, 2, 1], ordered=True))
        result = s.mode(dropna)
        expected3 = Series(expected3, dtype='category')
        tm.assert_series_equal(result, expected3)

    @pytest.mark.parametrize('dropna, expected1, expected2', [
        (True, [2**63], [1, 2**63]),
        (False, [2**63], [1, 2**63])
    ])
    def test_mode_intoverflow(self, dropna, expected1, expected2):
        # Test for uint64 overflow.
        s = Series([1, 2**63, 2**63], dtype=np.uint64)
        result = s.mode(dropna)
        expected1 = Series(expected1, dtype=np.uint64)
        tm.assert_series_equal(result, expected1)

        s = Series([1, 2**63], dtype=np.uint64)
        result = s.mode(dropna)
        expected2 = Series(expected2, dtype=np.uint64)
        tm.assert_series_equal(result, expected2)

    @pytest.mark.skipif(not compat.PY3, reason="only PY3")
    def test_mode_sortwarning(self):
        # Check for the warning that is raised when the mode
        # results cannot be sorted

        expected = Series(['foo', np.nan])
        s = Series([1, 'foo', 'foo', np.nan, np.nan])

        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            result = s.mode(dropna=False)
            result = result.sort_values().reset_index(drop=True)

        tm.assert_series_equal(result, expected)


class TestDataFrameReductions(object):
    # Note: the name TestDataFrameReductions indicates these tests
    #  were moved from a DataFrame-specific test file, _not_ that these
    #  tests are intended long-term to be DataFrame-specific

    def test_reduce_mixed_frame(self):
        # GH#6806
        df = DataFrame({
            'bool_data': [True, True, False, False, False],
            'int_data': [10, 20, 30, 40, 50],
            'string_data': ['a', 'b', 'c', 'd', 'e'],
        })
        df.reindex(columns=['bool_data', 'int_data', 'string_data'])
        test = df.sum(axis=0)
        tm.assert_numpy_array_equal(test.values,
                                    np.array([2, 150, 'abcde'], dtype=object))
        tm.assert_series_equal(test, df.T.sum(axis=1))

    # ----------------------------------------------------------------
    # Min/Max

    def test_max(self):

        int_frame = get_int_frame()
        float_frame = DataFrame(tm.getSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore", RuntimeWarning)
            assert_stat_op_calc('max', np.max, float_frame_with_na,
                                check_dates=True)
        assert_stat_op_calc('max', np.max, int_frame)
        assert_stat_op_api('max', float_frame, float_string_frame)

    def test_min(self):

        int_frame = get_int_frame()
        float_frame = DataFrame(tm.getSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore", RuntimeWarning)
            assert_stat_op_calc('min', np.min, float_frame_with_na,
                                check_dates=True)
        assert_stat_op_calc('min', np.min, int_frame)
        assert_stat_op_api('min', float_frame, float_string_frame)

    # ----------------------------------------------------------------
    # Any/All

    @pytest.mark.parametrize('opname', ['any', 'all'])
    def test_any_all(self, opname):
        float_string_frame = get_float_string_frame()
        bool_frame_with_na = get_bool_frame_with_na()

        assert_bool_op_calc(opname, getattr(np, opname), bool_frame_with_na,
                            has_skipna=True)
        assert_bool_op_api(opname, bool_frame_with_na, float_string_frame,
                           has_bool_only=True)

    def test_any_all_extra(self):
        df = DataFrame({
            'A': [True, False, False],
            'B': [True, True, False],
            'C': [True, True, True],
        }, index=['a', 'b', 'c'])
        result = df[['A', 'B']].any(1)
        expected = Series([True, True, False], index=['a', 'b', 'c'])
        tm.assert_series_equal(result, expected)

        result = df[['A', 'B']].any(1, bool_only=True)
        tm.assert_series_equal(result, expected)

        result = df.all(1)
        expected = Series([True, False, False], index=['a', 'b', 'c'])
        tm.assert_series_equal(result, expected)

        result = df.all(1, bool_only=True)
        tm.assert_series_equal(result, expected)

        # Axis is None
        result = df.all(axis=None).item()
        assert result is False

        result = df.any(axis=None).item()
        assert result is True

        result = df[['C']].all(axis=None).item()
        assert result is True

    def test_any_datetime(self):
        # GH#23070
        float_data = [1, np.nan, 3, np.nan]
        datetime_data = [pd.Timestamp('1960-02-15'),
                         pd.Timestamp('1960-02-16'),
                         pd.NaT,
                         pd.NaT]
        df = DataFrame({
            "A": float_data,
            "B": datetime_data
        })

        result = df.any(1)
        expected = Series([True, True, True, False])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('func, data, expected', [
        (np.any, {}, False),
        (np.all, {}, True),
        (np.any, {'A': []}, False),
        (np.all, {'A': []}, True),
        (np.any, {'A': [False, False]}, False),
        (np.all, {'A': [False, False]}, False),
        (np.any, {'A': [True, False]}, True),
        (np.all, {'A': [True, False]}, False),
        (np.any, {'A': [True, True]}, True),
        (np.all, {'A': [True, True]}, True),

        (np.any, {'A': [False], 'B': [False]}, False),
        (np.all, {'A': [False], 'B': [False]}, False),

        (np.any, {'A': [False, False], 'B': [False, True]}, True),
        (np.all, {'A': [False, False], 'B': [False, True]}, False),

        # other types
        (np.all, {'A': pd.Series([0.0, 1.0], dtype='float')}, False),
        (np.any, {'A': pd.Series([0.0, 1.0], dtype='float')}, True),
        (np.all, {'A': pd.Series([0, 1], dtype=int)}, False),
        (np.any, {'A': pd.Series([0, 1], dtype=int)}, True),
        pytest.param(np.all, {'A': pd.Series([0, 1], dtype='M8[ns]')}, False,
                     marks=[td.skip_if_np_lt_115]),
        pytest.param(np.any, {'A': pd.Series([0, 1], dtype='M8[ns]')}, True,
                     marks=[td.skip_if_np_lt_115]),
        pytest.param(np.all, {'A': pd.Series([1, 2], dtype='M8[ns]')}, True,
                     marks=[td.skip_if_np_lt_115]),
        pytest.param(np.any, {'A': pd.Series([1, 2], dtype='M8[ns]')}, True,
                     marks=[td.skip_if_np_lt_115]),
        pytest.param(np.all, {'A': pd.Series([0, 1], dtype='m8[ns]')}, False,
                     marks=[td.skip_if_np_lt_115]),
        pytest.param(np.any, {'A': pd.Series([0, 1], dtype='m8[ns]')}, True,
                     marks=[td.skip_if_np_lt_115]),
        pytest.param(np.all, {'A': pd.Series([1, 2], dtype='m8[ns]')}, True,
                     marks=[td.skip_if_np_lt_115]),
        pytest.param(np.any, {'A': pd.Series([1, 2], dtype='m8[ns]')}, True,
                     marks=[td.skip_if_np_lt_115]),
        (np.all, {'A': pd.Series([0, 1], dtype='category')}, False),
        (np.any, {'A': pd.Series([0, 1], dtype='category')}, True),
        (np.all, {'A': pd.Series([1, 2], dtype='category')}, True),
        (np.any, {'A': pd.Series([1, 2], dtype='category')}, True),

        # # Mix
        # GH#21484
        # (np.all, {'A': pd.Series([10, 20], dtype='M8[ns]'),
        #           'B': pd.Series([10, 20], dtype='m8[ns]')}, True),
    ])
    def test_any_all_np_func(self, func, data, expected):
        # GH#19976
        data = DataFrame(data)
        result = func(data)
        assert isinstance(result, np.bool_)
        assert result.item() is expected

        # method version
        result = getattr(DataFrame(data), func.__name__)(axis=None)
        assert isinstance(result, np.bool_)
        assert result.item() is expected

    def test_any_all_object(self):
        # GH#19976
        result = np.all(DataFrame(columns=['a', 'b'])).item()
        assert result is True

        result = np.any(DataFrame(columns=['a', 'b'])).item()
        assert result is False

    @pytest.mark.parametrize('method', ['any', 'all'])
    def test_any_all_level_axis_none_raises(self, method):
        df = DataFrame(
            {"A": 1},
            index=MultiIndex.from_product([['A', 'B'], ['a', 'b']],
                                          names=['out', 'in'])
        )
        xpr = "Must specify 'axis' when aggregating by level."
        with pytest.raises(ValueError, match=xpr):
            getattr(df, method)(axis=None, level='out')

    # ----------------------------------------------------------------
    # Statistical Reductions
    # TODO: belongs in test_stat_reductions

    # TODO: Ensure warning isn't emitted in the first place
    @pytest.mark.filterwarnings("ignore:All-NaN:RuntimeWarning")
    def test_median(self):

        float_frame = DataFrame(tm.getSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        def wrapper(x):
            if isna(x).any():
                return np.nan
            return np.median(x)

        assert_stat_op_calc('median', wrapper, float_frame_with_na,
                            check_dates=True)
        assert_stat_op_api('median', float_frame, float_string_frame)

    def test_mean(self):
        float_frame = DataFrame(tm.getSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        assert_stat_op_calc('mean', np.mean, float_frame_with_na,
                            check_dates=True)
        assert_stat_op_api('mean', float_frame, float_string_frame)

    def test_mean_corner(self):
        # unit test when have object data
        float_string_frame = pd.DataFrame(tm.getSeriesData())
        float_string_frame['foo'] = 'bar'

        the_mean = float_string_frame.mean(axis=0)
        the_sum = float_string_frame.sum(axis=0, numeric_only=True)
        tm.assert_index_equal(the_sum.index, the_mean.index)
        assert len(the_mean.index) < len(float_string_frame.columns)

        # xs sum mixed type, just want to know it works...
        the_mean = float_string_frame.mean(axis=1)
        the_sum = float_string_frame.sum(axis=1, numeric_only=True)
        tm.assert_index_equal(the_sum.index, the_mean.index)

        float_frame = pd.DataFrame(tm.getSeriesData())

        # take mean of boolean column
        float_frame['bool'] = float_frame['A'] > 0
        means = float_frame.mean(0)
        assert means['bool'] == float_frame['bool'].values.mean()

    # TODO: Ensure warning isn't emitted in the first place
    @pytest.mark.filterwarnings("ignore:All-NaN:RuntimeWarning")
    def test_median_corner(self, int_frame, float_frame, float_string_frame):
        def wrapper(x):
            if isna(x).any():
                return np.nan
            return np.median(x)

        assert_stat_op_calc('median', wrapper, int_frame, check_dtype=False,
                            check_dates=True)
        assert_stat_op_api('median', float_frame, float_string_frame)

    def test_stats_mixed_type(self):
        # don't blow up
        float_string_frame = pd.DataFrame(tm.getSeriesData())
        float_string_frame['foo'] = 'bar'

        float_string_frame.std(1)
        float_string_frame.var(1)
        float_string_frame.mean(1)
        float_string_frame.skew(1)

    def test_mad(self):
        float_frame = DataFrame(tm.getSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        f = lambda x: np.abs(x - x.mean()).mean()
        assert_stat_op_calc('mad', f, float_frame_with_na)
        assert_stat_op_api('mad', float_frame, float_string_frame)

    def test_var_std(self):

        float_frame = DataFrame(tm.getSeriesData())
        datetime_frame = DataFrame(tm.getTimeSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        alt = lambda x: np.var(x, ddof=1)
        assert_stat_op_calc('var', alt, float_frame_with_na)
        assert_stat_op_api('var', float_frame, float_string_frame)

        alt = lambda x: np.std(x, ddof=1)
        assert_stat_op_calc('std', alt, float_frame_with_na)
        assert_stat_op_api('std', float_frame, float_string_frame)

        result = datetime_frame.std(ddof=4)
        expected = datetime_frame.apply(lambda x: x.std(ddof=4))
        tm.assert_almost_equal(result, expected)

        result = datetime_frame.var(ddof=4)
        expected = datetime_frame.apply(lambda x: x.var(ddof=4))
        tm.assert_almost_equal(result, expected)

        arr = np.repeat(np.random.random((1, 1000)), 1000, 0)
        result = nanops.nanvar(arr, axis=0)
        assert not (result < 0).any()

        with pd.option_context('use_bottleneck', False):
            result = nanops.nanvar(arr, axis=0)
            assert not (result < 0).any()

    @pytest.mark.parametrize('op', ['mean', 'std', 'var',
                                    'skew', 'kurt', 'sem'])
    def test_mixed_ops(self, op):
        # GH#16116
        df = DataFrame({'int': [1, 2, 3, 4],
                        'float': [1., 2., 3., 4.],
                        'str': ['a', 'b', 'c', 'd']})

        result = getattr(df, op)()
        assert len(result) == 2

        with pd.option_context('use_bottleneck', False):
            result = getattr(df, op)()
            assert len(result) == 2

    @pytest.mark.parametrize("meth", ['sem', 'var', 'std'])
    def test_numeric_only_flag(self, meth):
        # GH#9201
        df1 = DataFrame(np.random.randn(5, 3), columns=['foo', 'bar', 'baz'])
        # set one entry to a number in str format
        df1.loc[0, 'foo'] = '100'

        df2 = DataFrame(np.random.randn(5, 3), columns=['foo', 'bar', 'baz'])
        # set one entry to a non-number str
        df2.loc[0, 'foo'] = 'a'

        result = getattr(df1, meth)(axis=1, numeric_only=True)
        expected = getattr(df1[['bar', 'baz']], meth)(axis=1)
        tm.assert_series_equal(expected, result)

        result = getattr(df2, meth)(axis=1, numeric_only=True)
        expected = getattr(df2[['bar', 'baz']], meth)(axis=1)
        tm.assert_series_equal(expected, result)

        # df1 has all numbers, df2 has a letter inside
        with pytest.raises(TypeError):
            getattr(df1, meth)(axis=1, numeric_only=False)
        with pytest.raises(TypeError):
            getattr(df2, meth)(axis=1, numeric_only=False)

    def test_sem(self):

        float_frame = DataFrame(tm.getSeriesData())
        datetime_frame = DataFrame(tm.getTimeSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        alt = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
        assert_stat_op_calc('sem', alt, float_frame_with_na)
        assert_stat_op_api('sem', float_frame, float_string_frame)

        result = datetime_frame.sem(ddof=4)
        expected = datetime_frame.apply(
            lambda x: x.std(ddof=4) / np.sqrt(len(x)))
        tm.assert_almost_equal(result, expected)

        arr = np.repeat(np.random.random((1, 1000)), 1000, 0)
        result = nanops.nansem(arr, axis=0)
        assert not (result < 0).any()

        with pd.option_context('use_bottleneck', False):
            result = nanops.nansem(arr, axis=0)
            assert not (result < 0).any()

    @td.skip_if_no_scipy
    def test_skew(self):
        from scipy.stats import skew

        float_frame = DataFrame(tm.getSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        def alt(x):
            if len(x) < 3:
                return np.nan
            return skew(x, bias=False)

        assert_stat_op_calc('skew', alt, float_frame_with_na)
        assert_stat_op_api('skew', float_frame, float_string_frame)

    @td.skip_if_no_scipy
    def test_kurt(self):
        from scipy.stats import kurtosis

        float_frame = DataFrame(tm.getSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        def alt(x):
            if len(x) < 4:
                return np.nan
            return kurtosis(x, bias=False)

        assert_stat_op_calc('kurt', alt, float_frame_with_na)
        assert_stat_op_api('kurt', float_frame, float_string_frame)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           codes=[[0, 0, 0, 0, 0, 0],
                                  [0, 1, 2, 0, 1, 2],
                                  [0, 1, 0, 1, 0, 1]])
        df = DataFrame(np.random.randn(6, 3), index=index)

        kurt = df.kurt()
        kurt2 = df.kurt(level=0).xs('bar')
        tm.assert_series_equal(kurt, kurt2, check_names=False)
        assert kurt.name is None
        assert kurt2.name == 'bar'

    @pytest.mark.skipif(not compat.PY3, reason="only PY3")
    def test_mode_sortwarning(self):
        # Check for the warning that is raised when the mode
        # results cannot be sorted

        df = DataFrame({"A": [np.nan, np.nan, 'a', 'a']})
        expected = DataFrame({'A': ['a', np.nan]})

        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            result = df.mode(dropna=False)
            result = result.sort_values(by='A').reset_index(drop=True)

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("dropna, expected", [
        (True, {'A': [12],
                'B': [10.0],
                'C': [1.0],
                'D': ['a'],
                'E': Categorical(['a'], categories=['a']),
                'F': pd.to_datetime(['2000-1-2']),
                'G': pd.to_timedelta(['1 days'])}),
        (False, {'A': [12],
                 'B': [10.0],
                 'C': [np.nan],
                 'D': np.array([np.nan], dtype=object),
                 'E': Categorical([np.nan], categories=['a']),
                 'F': [pd.NaT],
                 'G': pd.to_timedelta([pd.NaT])}),
        (True, {'H': [8, 9, np.nan, np.nan],
                'I': [8, 9, np.nan, np.nan],
                'J': [1, np.nan, np.nan, np.nan],
                'K': Categorical(['a', np.nan, np.nan, np.nan],
                                 categories=['a']),
                'L': pd.to_datetime(['2000-1-2', 'NaT', 'NaT', 'NaT']),
                'M': pd.to_timedelta(['1 days', 'nan', 'nan', 'nan']),
                'N': [0, 1, 2, 3]}),
        (False, {'H': [8, 9, np.nan, np.nan],
                 'I': [8, 9, np.nan, np.nan],
                 'J': [1, np.nan, np.nan, np.nan],
                 'K': Categorical([np.nan, 'a', np.nan, np.nan],
                                  categories=['a']),
                 'L': pd.to_datetime(['NaT', '2000-1-2', 'NaT', 'NaT']),
                 'M': pd.to_timedelta(['nan', '1 days', 'nan', 'nan']),
                 'N': [0, 1, 2, 3]})
    ])
    def test_mode_dropna(self, dropna, expected):

        df = DataFrame({"A": [12, 12, 19, 11],
                        "B": [10, 10, np.nan, 3],
                        "C": [1, np.nan, np.nan, np.nan],
                        "D": [np.nan, np.nan, 'a', np.nan],
                        "E": Categorical([np.nan, np.nan, 'a', np.nan]),
                        "F": pd.to_datetime(['NaT', '2000-1-2', 'NaT', 'NaT']),
                        "G": pd.to_timedelta(['1 days', 'nan', 'nan', 'nan']),
                        "H": [8, 8, 9, 9],
                        "I": [9, 9, 8, 8],
                        "J": [1, 1, np.nan, np.nan],
                        "K": Categorical(['a', np.nan, 'a', np.nan]),
                        "L": pd.to_datetime(['2000-1-2', '2000-1-2',
                                             'NaT', 'NaT']),
                        "M": pd.to_timedelta(['1 days', 'nan',
                                              '1 days', 'nan']),
                        "N": np.arange(4, dtype='int64')})

        result = df[sorted(list(expected.keys()))].mode(dropna=dropna)
        expected = DataFrame(expected)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('method', ['sum', 'mean', 'prod', 'var',
                                        'std', 'skew', 'min', 'max'])
    def test_stat_operators_attempt_obj_array(self, method):
        # GH#676
        data = {
            'a': [-0.00049987540199591344, -0.0016467257772919831,
                  0.00067695870775883013],
            'b': [-0, -0, 0.0],
            'c': [0.00031111847529610595, 0.0014902627951905339,
                  -0.00094099200035979691]
        }
        df1 = DataFrame(data, index=['foo', 'bar', 'baz'], dtype='O')

        df2 = DataFrame({0: [np.nan, 2], 1: [np.nan, 3],
                         2: [np.nan, 4]}, dtype=object)

        for df in [df1, df2]:
            assert df.values.dtype == np.object_
            result = getattr(df, method)(1)
            expected = getattr(df.astype('f8'), method)(1)

            if method in ['sum', 'prod']:
                tm.assert_series_equal(result, expected)

    # ----------------------------------------------------------------
    # Sums

    def test_sum(self):

        float_frame = DataFrame(tm.getSeriesData())
        float_frame_with_na = get_float_frame_with_na()
        float_string_frame = get_float_string_frame()

        mixed_float_frame = DataFrame(tm.getSeriesData())
        mixed_float_frame.A = mixed_float_frame.A.astype('float32')
        mixed_float_frame.B = mixed_float_frame.B.astype('float32')
        mixed_float_frame.C = mixed_float_frame.C.astype('float16')
        mixed_float_frame.D = mixed_float_frame.D.astype('float64')
        # TODO: the only place we use this casts to float32; is this useful?

        assert_stat_op_api('sum', float_frame, float_string_frame,
                           has_numeric_only=True)
        assert_stat_op_calc('sum', np.sum, float_frame_with_na,
                            skipna_alternative=np.nansum)
        # mixed types (with upcasting happening)
        assert_stat_op_calc('sum', np.sum, mixed_float_frame.astype('float32'),
                            check_dtype=False, check_less_precise=True)

    def test_sum_object(self):
        float_frame = DataFrame(tm.getSeriesData())
        values = float_frame.values.astype(int)
        frame = DataFrame(values, index=float_frame.index,
                          columns=float_frame.columns)
        deltas = frame * timedelta(1)
        deltas.sum()

    def test_sum_bool(self):
        float_frame = DataFrame(tm.getSeriesData())

        # ensure this works, bug report
        bools = np.isnan(float_frame)
        bools.sum(1)
        bools.sum(0)

    @pytest.mark.parametrize('method, unit', [
        ('sum', 0),
        ('prod', 1),
    ])
    def test_sum_prod_nanops(self, method, unit):
        idx = ['a', 'b', 'c']
        df = pd.DataFrame({"a": [unit, unit],
                           "b": [unit, np.nan],
                           "c": [np.nan, np.nan]})
        # The default
        result = getattr(df, method)
        expected = pd.Series([unit, unit, unit], index=idx, dtype='float64')

        # min_count=1
        result = getattr(df, method)(min_count=1)
        expected = pd.Series([unit, unit, np.nan], index=idx)
        tm.assert_series_equal(result, expected)

        # min_count=0
        result = getattr(df, method)(min_count=0)
        expected = pd.Series([unit, unit, unit], index=idx, dtype='float64')
        tm.assert_series_equal(result, expected)

        result = getattr(df.iloc[1:], method)(min_count=1)
        expected = pd.Series([unit, np.nan, np.nan], index=idx)
        tm.assert_series_equal(result, expected)

        # min_count > 1
        df = pd.DataFrame({"A": [unit] * 10, "B": [unit] * 5 + [np.nan] * 5})
        result = getattr(df, method)(min_count=5)
        expected = pd.Series(result, index=['A', 'B'])
        tm.assert_series_equal(result, expected)

        result = getattr(df, method)(min_count=6)
        expected = pd.Series(result, index=['A', 'B'])
        tm.assert_series_equal(result, expected)

    def test_sum_nanops_timedelta(self):
        # prod isn't defined on timedeltas
        idx = ['a', 'b', 'c']
        df = pd.DataFrame({"a": [0, 0],
                           "b": [0, np.nan],
                           "c": [np.nan, np.nan]})

        df2 = df.apply(pd.to_timedelta)

        # 0 by default
        result = df2.sum()
        expected = pd.Series([0, 0, 0], dtype='m8[ns]', index=idx)
        tm.assert_series_equal(result, expected)

        # min_count=0
        result = df2.sum(min_count=0)
        tm.assert_series_equal(result, expected)

        # min_count=1
        result = df2.sum(min_count=1)
        expected = pd.Series([0, 0, np.nan], dtype='m8[ns]', index=idx)
        tm.assert_series_equal(result, expected)

    def test_sum_corner(self):
        empty_frame = DataFrame({})

        axis0 = empty_frame.sum(0)
        axis1 = empty_frame.sum(1)
        assert isinstance(axis0, Series)
        assert isinstance(axis1, Series)
        assert len(axis0) == 0
        assert len(axis1) == 0

    def test_sum_bools(self):
        df = pd.DataFrame(index=lrange(1), columns=lrange(10))
        bools = isna(df)
        assert bools.sum(axis=1)[0] == 10

    # ----------------------------------------------------------------

    def test_product(self, float_frame_with_na, float_frame,
                     float_string_frame):
        assert_stat_op_calc('product', np.prod, float_frame_with_na)
        assert_stat_op_api('product', float_frame, float_string_frame)
