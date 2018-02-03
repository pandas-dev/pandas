# -*- coding: utf-8 -*-
import operator

import pytest
import numpy as np

from pandas.compat import PY3

import pandas as pd
import pandas.util.testing as tm

from pandas.tests.frame.common import _check_mixed_float, _check_mixed_int


# -------------------------------------------------------------------
# Comparisons

class TestFrameComparisons(object):
    def test_df_boolean_comparison_error(self):
        # GH#4576
        # boolean comparisons with a tuple/list give unexpected results
        df = pd.DataFrame(np.arange(6).reshape((3, 2)))

        # not shape compatible
        with pytest.raises(ValueError):
            df == (2, 2)
        with pytest.raises(ValueError):
            df == [2, 2]

    def test_df_float_none_comparison(self):
        df = pd.DataFrame(np.random.randn(8, 3), index=range(8),
                          columns=['A', 'B', 'C'])

        with pytest.raises(TypeError):
            df.__eq__(None)

    def test_df_string_comparison(self):
        df = pd.DataFrame([{"a": 1, "b": "foo"}, {"a": 2, "b": "bar"}])
        mask_a = df.a > 1
        tm.assert_frame_equal(df[mask_a], df.loc[1:1, :])
        tm.assert_frame_equal(df[-mask_a], df.loc[0:0, :])

        mask_b = df.b == "foo"
        tm.assert_frame_equal(df[mask_b], df.loc[0:0, :])
        tm.assert_frame_equal(df[-mask_b], df.loc[1:1, :])

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_df_flex_cmp_constant_return_types(self, opname):
        # GH#15077, non-empty DataFrame
        df = pd.DataFrame({'x': [1, 2, 3], 'y': [1., 2., 3.]})
        const = 2

        result = getattr(df, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, pd.Series([2], ['bool']))

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_df_flex_cmp_constant_return_types_empty(self, opname):
        # GH#15077 empty DataFrame
        df = pd.DataFrame({'x': [1, 2, 3], 'y': [1., 2., 3.]})
        const = 2

        empty = df.iloc[:0]
        result = getattr(empty, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, pd.Series([2], ['bool']))


# -------------------------------------------------------------------
# Arithmetic

class TestFrameArithmetic(object):

    @pytest.mark.xfail(reason='GH#7996 datetime64 units not converted to nano')
    def test_df_sub_datetime64_not_ns(self):
        df = pd.DataFrame(pd.date_range('20130101', periods=3))
        dt64 = np.datetime64('2013-01-01')
        assert dt64.dtype == 'datetime64[D]'
        res = df - dt64
        expected = pd.DataFrame([pd.Timedelta(days=0), pd.Timedelta(days=1),
                                 pd.Timedelta(days=2)])
        tm.assert_frame_equal(res, expected)

    @pytest.mark.parametrize('data', [
        [1, 2, 3],
        [1.1, 2.2, 3.3],
        [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02'), pd.NaT],
        ['x', 'y', 1]])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_df_radd_str_invalid(self, dtype, data):
        df = pd.DataFrame(data, dtype=dtype)
        with pytest.raises(TypeError):
            'foo_' + df

    @pytest.mark.parametrize('dtype', [None, object])
    def test_df_with_dtype_radd_int(self, dtype):
        df = pd.DataFrame([1, 2, 3], dtype=dtype)
        expected = pd.DataFrame([2, 3, 4], dtype=dtype)
        result = 1 + df
        tm.assert_frame_equal(result, expected)
        result = df + 1
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_df_with_dtype_radd_nan(self, dtype):
        df = pd.DataFrame([1, 2, 3], dtype=dtype)
        expected = pd.DataFrame([np.nan, np.nan, np.nan], dtype=dtype)
        result = np.nan + df
        tm.assert_frame_equal(result, expected)
        result = df + np.nan
        tm.assert_frame_equal(result, expected)

    def test_df_radd_str(self):
        df = pd.DataFrame(['x', np.nan, 'x'])
        tm.assert_frame_equal('a' + df, pd.DataFrame(['ax', np.nan, 'ax']))
        tm.assert_frame_equal(df + 'a', pd.DataFrame(['xa', np.nan, 'xa']))


class TestPeriodFrameArithmetic(object):

    def test_ops_frame_period(self):
        # GH 13043
        df = pd.DataFrame({'A': [pd.Period('2015-01', freq='M'),
                                 pd.Period('2015-02', freq='M')],
                           'B': [pd.Period('2014-01', freq='M'),
                                 pd.Period('2014-02', freq='M')]})
        assert df['A'].dtype == object
        assert df['B'].dtype == object

        p = pd.Period('2015-03', freq='M')
        # dtype will be object because of original dtype
        exp = pd.DataFrame({'A': np.array([2, 1], dtype=object),
                            'B': np.array([14, 13], dtype=object)})
        tm.assert_frame_equal(p - df, exp)
        tm.assert_frame_equal(df - p, -exp)

        df2 = pd.DataFrame({'A': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')],
                            'B': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')]})
        assert df2['A'].dtype == object
        assert df2['B'].dtype == object

        exp = pd.DataFrame({'A': np.array([4, 4], dtype=object),
                            'B': np.array([16, 16], dtype=object)})
        tm.assert_frame_equal(df2 - df, exp)
        tm.assert_frame_equal(df - df2, -exp)


class TestFrameFlexArithmetic(object):
    """Tests for DataFrame flex operations, mostly focused on alignment
    and broadcasting.
    """
    def test_df_flex_add_int(self):
        frame = pd.DataFrame(tm.getSeriesData())

        result = frame.add(1)
        expected = frame + 1
        tm.assert_frame_equal(result, expected)

    def test_df_flex_add_empty_frame(self):
        # corner cases by convention returns frame * np.nan
        frame = pd.DataFrame(tm.getSeriesData())

        result = frame.add(frame[:0])
        expected = frame * np.nan
        tm.assert_frame_equal(result, expected)

        result = frame[:0].add(frame)
        tm.assert_frame_equal(result, frame * np.nan)

    def test_df_flex_add_fill_value_raises(self):
        # corner case NotImplementeError
        frame = pd.DataFrame(tm.getSeriesData())
        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            frame.add(frame.iloc[0], fill_value=3)

        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            frame.add(frame.iloc[0], axis='index', fill_value=3)

    @pytest.mark.parametrize('opname', ['add', 'sub', 'mul',
                                        'div', 'truediv',
                                        'pow', 'floordiv', 'mod'])
    def test_df_flex_add_higher_dim_array_raises(self, opname):
        frame = pd.DataFrame(tm.getSeriesData())

        alias = 'truediv' if PY3 and opname == 'div' else opname
        dunder_func = getattr(operator, alias)
        flex_method = getattr(frame, opname)

        # ndim >= 3
        ndim_5 = np.ones(frame.shape + (3, 4, 5))
        msg = "Unable to coerce to Series/DataFrame"

        with tm.assert_raises_regex(ValueError, msg):
            dunder_func(frame, ndim_5)

        with tm.assert_raises_regex(ValueError, msg):
            flex_method(ndim_5)

    @pytest.mark.parametrize('opname', ['add', 'sub', 'mul',
                                        'div', 'truediv',
                                        'pow', 'floordiv', 'mod'])
    def test_df_flex_with_df_matches_dunder(self, opname):
        # df.opname(other) is equivalent to df.__opname__(other) if other
        # is also a frame indexed like df.
        frame = pd.DataFrame(tm.getSeriesData())

        alias = 'truediv' if PY3 and opname == 'div' else opname
        dunder_func = getattr(operator, alias)
        result = getattr(frame, opname)(2 * frame)
        expected = dunder_func(frame, 2 * frame)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('opname', ['add', 'sub', 'mul'])
    def test_df_flex_reversed_with_df_matches_dunder(self, opname):
        frame = pd.DataFrame(tm.getSeriesData())

        dunder_func = getattr(operator, opname)
        rfunc = lambda x, y: dunder_func(y, x)

        result = getattr(frame, 'r' + opname)(2 * frame)
        expected = rfunc(frame, 2 * frame)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('opname', ['add', 'sub', 'mul'])
    def test_df_flex_needs_better_name(self, opname):
        # TODO: what is being tested here?  needs better name
        frame = pd.DataFrame(tm.getSeriesData())
        mixed_float = pd.DataFrame({'A': frame['A'].copy().astype('float32'),
                                    'B': frame['B'].copy().astype('float32'),
                                    'C': frame['C'].copy().astype('float16'),
                                    'D': frame['D'].copy().astype('float64')})
        # force these all to int64 to avoid platform testing issues
        intframe = pd.DataFrame({k: v.astype(int)
                                 for k, v in tm.getSeriesData().items()},
                                dtype=np.int64)
        mixed_int = pd.DataFrame({'A': intframe['A'].copy().astype('int32'),
                                  'B': np.ones(len(intframe['B']),
                                               dtype='uint64'),
                                  'C': intframe['C'].copy().astype('uint8'),
                                  'D': intframe['D'].copy().astype('int64')})

        dunder_func = getattr(operator, opname)

        # vs mix int
        result = getattr(mixed_int, opname)(2 + mixed_int)
        expected = dunder_func(mixed_int, 2 + mixed_int)
        tm.assert_frame_equal(result, expected)

        dtype = {'C': None}
        if opname == 'sub':
            dtype['B'] = 'uint64'

        # no overflow in the uint
        _check_mixed_int(result, dtype=dtype)

        # vs mix float
        result = getattr(mixed_float, opname)(2 * mixed_float)
        expected = dunder_func(mixed_float, 2 * mixed_float)
        tm.assert_frame_equal(result, expected)
        _check_mixed_float(result, dtype={'C': None})

        result = getattr(intframe, opname)(2 * intframe)
        expected = dunder_func(intframe, 2 * intframe)
        tm.assert_frame_equal(result, expected)

        # vs mix int
        result = getattr(mixed_int, opname)(2 + mixed_int)
        expected = dunder_func(mixed_int, 2 + mixed_int)
        tm.assert_frame_equal(result, expected)

        # no overflow in the uint
        _check_mixed_int(result, dtype=dtype)

    @pytest.mark.parametrize('opname', ['add', 'sub', 'mul',
                                        'div', 'truediv',
                                        'pow', 'floordiv', 'mod'])
    def test_df_arith_flex_float_dtypes_match_dunder(self, opname):
        # TODO: Is this supposed to be about the dtypes or about the
        # indexed_same?
        frame = pd.DataFrame(tm.getSeriesData())

        mixed_float = pd.DataFrame({'A': frame['A'].copy().astype('float32'),
                                    'B': frame['B'].copy().astype('float32'),
                                    'C': frame['C'].copy().astype('float16'),
                                    'D': frame['D'].copy().astype('float64')})

        alias = 'truediv' if PY3 and opname == 'div' else opname
        try:
            dunder_func = getattr(operator, alias)

            # vs mix float
            result = getattr(mixed_float, opname)(2 * mixed_float)
            expected = dunder_func(mixed_float, 2 * mixed_float)
            tm.assert_frame_equal(result, expected)
            _check_mixed_float(result, dtype={'C': None})

        except:
            import pandas.io.formats.printing as printing
            printing.pprint_thing("Failing operation %r" % opname)
            raise

    def test_df_arith_flex(self):
        pass
        # res_add = frame.add(frame)
        # res_sub = frame.sub(frame)
        # res_mul = frame.mul(frame)
        # res_div = frame.div(2 * frame)

        # tm.assert_frame_equal(res_add, frame + frame)
        # tm.assert_frame_equal(res_sub, frame - frame)
        # tm.assert_frame_equal(res_mul, frame * frame)
        # tm.assert_frame_equal(res_div, frame / (2 * frame))

    def test_df_flex_noaxis_match_dunder(self):
        # Check that df.name(other) without passing an `axis` kwarg matches
        # df.__name__(other)
        df = pd.DataFrame(np.arange(1, 10).reshape(3, 3),
                          columns=['one', 'two', 'three'],
                          index=['a', 'b', 'c'])
        row = df.xs('a')
        # after arithmetic refactor, add truediv here

        result = df.add(row)
        expected = df + row
        tm.assert_frame_equal(result, expected)

        result = df.sub(row)
        expected = df - row
        tm.assert_frame_equal(result, expected)

        result = df.mul(row)
        expected = df * row
        tm.assert_frame_equal(result, expected)

        result = df.mod(row)
        expected = df % row
        tm.assert_frame_equal(result, expected)

    def test_df_flex_arith_methods_axis_match_tranposed_operator_ops(self):
        df = pd.DataFrame(np.arange(1, 10).reshape(3, 3),
                          columns=['one', 'two', 'three'],
                          index=['a', 'b', 'c'])
        row = df.xs('a')
        col = df['two']
        # after arithmetic refactor, add truediv here

        result = df.add(col, axis=0)
        expected = df.T + row
        tm.assert_frame_equal(result, expected.T)

        result = df.sub(col, axis=0)
        expected = df.T - row
        tm.assert_frame_equal(result, expected.T)

        result = df.mul(col, axis=0)
        expected = df.T * row
        tm.assert_frame_equal(result, expected.T)

        result = df.mod(col, axis=0)
        expected = df.T % row
        tm.assert_frame_equal(result, expected.T)

    def test_df_flex_arith_corner_cases(self):
        df = pd.DataFrame(np.arange(1, 10).reshape(3, 3),
                          columns=['one', 'two', 'three'],
                          index=['a', 'b', 'c'])
        row = df.xs('a')
        col = df['two']

        # special case for some reason
        result = df.add(row, axis=None)
        expected = df + row
        tm.assert_frame_equal(result, expected)

        # cases which will be refactored after big arithmetic refactor
        result = df.div(row)
        expected = df / row
        tm.assert_frame_equal(result, expected)

        result = df.div(col, axis=0)
        expected = df.T / col
        tm.assert_frame_equal(result, expected.T)

    def test_df_flex_div_series(self):
        # broadcasting issue in GH#7325
        df = pd.DataFrame(np.arange(3 * 2).reshape((3, 2)), dtype='int64')
        expected = pd.DataFrame([[np.nan, np.inf], [1.0, 1.5], [1.0, 1.25]])

        result = df.div(df[0], axis='index')
        tm.assert_frame_equal(result, expected)

        df = df.astype('float64')
        result = df.div(df[0], axis='index')
        tm.assert_frame_equal(result, expected)
