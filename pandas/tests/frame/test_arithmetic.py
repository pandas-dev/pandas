# -*- coding: utf-8 -*-
import operator

import pytest
import numpy as np

from pandas.compat import range, PY3
import pandas.io.formats.printing as printing

import pandas as pd
import pandas.util.testing as tm

from pandas.tests.frame.common import _check_mixed_float, _check_mixed_int


# -------------------------------------------------------------------
# Comparisons

class TestFrameComparisons(object):
    def test_flex_comparison_nat(self):
        # GH#15697, GH#22163 df.eq(pd.NaT) should behave like df == pd.NaT,
        # and _definitely_ not be NaN
        df = pd.DataFrame([pd.NaT])

        result = df == pd.NaT
        # result.iloc[0, 0] is a np.bool_ object
        assert result.iloc[0, 0].item() is False

        result = df.eq(pd.NaT)
        assert result.iloc[0, 0].item() is False

        result = df != pd.NaT
        assert result.iloc[0, 0].item() is True

        result = df.ne(pd.NaT)
        assert result.iloc[0, 0].item() is True

    def test_mixed_comparison(self):
        # GH#13128, GH#22163 != datetime64 vs non-dt64 should be False,
        # not raise TypeError
        # (this appears to be fixed before #22163, not sure when)
        df = pd.DataFrame([['1989-08-01', 1], ['1989-08-01', 2]])
        other = pd.DataFrame([['a', 'b'], ['c', 'd']])

        result = df == other
        assert not result.any().any()

        result = df != other
        assert result.all().all()

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

        result = df.__eq__(None)
        assert not result.any().any()

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

class TestFrameFlexArithmetic(object):
    def test_df_add_flex_filled_mixed_dtypes(self):
        # GH#19611
        dti = pd.date_range('2016-01-01', periods=3)
        ser = pd.Series(['1 Day', 'NaT', '2 Days'], dtype='timedelta64[ns]')
        df = pd.DataFrame({'A': dti, 'B': ser})
        other = pd.DataFrame({'A': ser, 'B': ser})
        fill = pd.Timedelta(days=1).to_timedelta64()
        result = df.add(other, fill_value=fill)

        expected = pd.DataFrame(
            {'A': pd.Series(['2016-01-02', '2016-01-03', '2016-01-05'],
                            dtype='datetime64[ns]'),
             'B': ser * 2})
        tm.assert_frame_equal(result, expected)

    def test_arith_flex_frame(self):
        seriesd = tm.getSeriesData()
        frame = pd.DataFrame(seriesd).copy()

        mixed_float = pd.DataFrame({'A': frame['A'].copy().astype('float32'),
                                    'B': frame['B'].copy().astype('float32'),
                                    'C': frame['C'].copy().astype('float16'),
                                    'D': frame['D'].copy().astype('float64')})

        intframe = pd.DataFrame({k: v.astype(int)
                                 for k, v in seriesd.items()})
        mixed_int = pd.DataFrame({'A': intframe['A'].copy().astype('int32'),
                                  'B': np.ones(len(intframe), dtype='uint64'),
                                  'C': intframe['C'].copy().astype('uint8'),
                                  'D': intframe['D'].copy().astype('int64')})

        # force these all to int64 to avoid platform testing issues
        intframe = pd.DataFrame({c: s for c, s in intframe.items()},
                                dtype=np.int64)

        ops = ['add', 'sub', 'mul', 'div', 'truediv', 'pow', 'floordiv', 'mod']
        if not PY3:
            aliases = {}
        else:
            aliases = {'div': 'truediv'}

        for op in ops:
            try:
                alias = aliases.get(op, op)
                f = getattr(operator, alias)
                result = getattr(frame, op)(2 * frame)
                exp = f(frame, 2 * frame)
                tm.assert_frame_equal(result, exp)

                # vs mix float
                result = getattr(mixed_float, op)(2 * mixed_float)
                exp = f(mixed_float, 2 * mixed_float)
                tm.assert_frame_equal(result, exp)
                _check_mixed_float(result, dtype=dict(C=None))

                # vs mix int
                if op in ['add', 'sub', 'mul']:
                    result = getattr(mixed_int, op)(2 + mixed_int)
                    exp = f(mixed_int, 2 + mixed_int)

                    # no overflow in the uint
                    dtype = None
                    if op in ['sub']:
                        dtype = dict(B='uint64', C=None)
                    elif op in ['add', 'mul']:
                        dtype = dict(C=None)
                    tm.assert_frame_equal(result, exp)
                    _check_mixed_int(result, dtype=dtype)

                    # rops
                    r_f = lambda x, y: f(y, x)
                    result = getattr(frame, 'r' + op)(2 * frame)
                    exp = r_f(frame, 2 * frame)
                    tm.assert_frame_equal(result, exp)

                    # vs mix float
                    result = getattr(mixed_float, op)(2 * mixed_float)
                    exp = f(mixed_float, 2 * mixed_float)
                    tm.assert_frame_equal(result, exp)
                    _check_mixed_float(result, dtype=dict(C=None))

                    result = getattr(intframe, op)(2 * intframe)
                    exp = f(intframe, 2 * intframe)
                    tm.assert_frame_equal(result, exp)

                    # vs mix int
                    if op in ['add', 'sub', 'mul']:
                        result = getattr(mixed_int, op)(2 + mixed_int)
                        exp = f(mixed_int, 2 + mixed_int)

                        # no overflow in the uint
                        dtype = None
                        if op in ['sub']:
                            dtype = dict(B='uint64', C=None)
                        elif op in ['add', 'mul']:
                            dtype = dict(C=None)
                        tm.assert_frame_equal(result, exp)
                        _check_mixed_int(result, dtype=dtype)
            except:
                printing.pprint_thing("Failing operation %r" % op)
                raise

            # ndim >= 3
            ndim_5 = np.ones(frame.shape + (3, 4, 5))
            msg = "Unable to coerce to Series/DataFrame"
            with tm.assert_raises_regex(ValueError, msg):
                f(frame, ndim_5)

            with tm.assert_raises_regex(ValueError, msg):
                getattr(frame, op)(ndim_5)

        # res_add = frame.add(frame)
        # res_sub = frame.sub(frame)
        # res_mul = frame.mul(frame)
        # res_div = frame.div(2 * frame)

        # tm.assert_frame_equal(res_add, frame + frame)
        # tm.assert_frame_equal(res_sub, frame - frame)
        # tm.assert_frame_equal(res_mul, frame * frame)
        # tm.assert_frame_equal(res_div, frame / (2 * frame))

        const_add = frame.add(1)
        tm.assert_frame_equal(const_add, frame + 1)

        # corner cases
        result = frame.add(frame[:0])
        tm.assert_frame_equal(result, frame * np.nan)

        result = frame[:0].add(frame)
        tm.assert_frame_equal(result, frame * np.nan)
        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            frame.add(frame.iloc[0], fill_value=3)
        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            frame.add(frame.iloc[0], axis='index', fill_value=3)

    def test_arith_flex_series(self):
        arr = np.array([[1., 2., 3.],
                        [4., 5., 6.],
                        [7., 8., 9.]])
        df = pd.DataFrame(arr, columns=['one', 'two', 'three'],
                          index=['a', 'b', 'c'])

        row = df.xs('a')
        col = df['two']
        # after arithmetic refactor, add truediv here
        ops = ['add', 'sub', 'mul', 'mod']
        for op in ops:
            f = getattr(df, op)
            op = getattr(operator, op)
            tm.assert_frame_equal(f(row), op(df, row))
            tm.assert_frame_equal(f(col, axis=0), op(df.T, col).T)

        # special case for some reason
        tm.assert_frame_equal(df.add(row, axis=None), df + row)

        # cases which will be refactored after big arithmetic refactor
        tm.assert_frame_equal(df.div(row), df / row)
        tm.assert_frame_equal(df.div(col, axis=0), (df.T / col).T)

        # broadcasting issue in GH#7325
        df = pd.DataFrame(np.arange(3 * 2).reshape((3, 2)), dtype='int64')
        expected = pd.DataFrame([[np.nan, np.inf], [1.0, 1.5], [1.0, 1.25]])
        result = df.div(df[0], axis='index')
        tm.assert_frame_equal(result, expected)

        df = pd.DataFrame(np.arange(3 * 2).reshape((3, 2)), dtype='float64')
        expected = pd.DataFrame([[np.nan, np.inf], [1.0, 1.5], [1.0, 1.25]])
        result = df.div(df[0], axis='index')
        tm.assert_frame_equal(result, expected)

    def test_arith_flex_zero_len_raises(self):
        # GH#19522 passing fill_value to frame flex arith methods should
        # raise even in the zero-length special cases
        ser_len0 = pd.Series([])
        df_len0 = pd.DataFrame([], columns=['A', 'B'])
        df = pd.DataFrame([[1, 2], [3, 4]], columns=['A', 'B'])

        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            df.add(ser_len0, fill_value='E')

        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            df_len0.sub(df['A'], axis=None, fill_value=3)


class TestFrameArithmetic(object):
    def test_df_bool_mul_int(self):
        # GH#22047, GH#22163 multiplication by 1 should result in int dtype,
        # not object dtype
        df = pd.DataFrame([[False, True], [False, False]])
        result = df * 1

        # On appveyor this comes back as np.int32 instead of np.int64,
        # so we check dtype.kind instead of just dtype
        kinds = result.dtypes.apply(lambda x: x.kind)
        assert (kinds == 'i').all()

        result = 1 * df
        kinds = result.dtypes.apply(lambda x: x.kind)
        assert (kinds == 'i').all()

    def test_td64_df_add_int_frame(self):
        # Check that we don't dispatch to numpy implementation, which treats
        # int64 as m8[ns]
        tdi = pd.timedelta_range('1', periods=3)
        df = tdi.to_frame()
        other = pd.DataFrame([1, 2, 3], index=tdi)  # indexed like `df`
        with pytest.raises(TypeError):
            df + other
        with pytest.raises(TypeError):
            other + df
        with pytest.raises(TypeError):
            df - other
        with pytest.raises(TypeError):
            other - df
