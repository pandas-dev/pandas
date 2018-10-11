# -*- coding: utf-8 -*-
import operator

import pytest
import numpy as np

from pandas.compat import range

import pandas as pd
import pandas.util.testing as tm

from pandas.tests.frame.common import _check_mixed_float, _check_mixed_int


# -------------------------------------------------------------------
# Comparisons

class TestFrameComparisons(object):
    def test_flex_comparison_nat(self):
        # GH 15697, GH 22163 df.eq(pd.NaT) should behave like df == pd.NaT,
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
        # GH 13128, GH 22163 != datetime64 vs non-dt64 should be False,
        # not raise TypeError
        # (this appears to be fixed before #22163, not sure when)
        df = pd.DataFrame([['1989-08-01', 1], ['1989-08-01', 2]])
        other = pd.DataFrame([['a', 'b'], ['c', 'd']])

        result = df == other
        assert not result.any().any()

        result = df != other
        assert result.all().all()

    def test_df_boolean_comparison_error(self):
        # GH 4576
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
        # GH 15077, non-empty DataFrame
        df = pd.DataFrame({'x': [1, 2, 3], 'y': [1., 2., 3.]})
        const = 2

        result = getattr(df, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, pd.Series([2], ['bool']))

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_df_flex_cmp_constant_return_types_empty(self, opname):
        # GH 15077 empty DataFrame
        df = pd.DataFrame({'x': [1, 2, 3], 'y': [1., 2., 3.]})
        const = 2

        empty = df.iloc[:0]
        result = getattr(empty, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, pd.Series([2], ['bool']))


# -------------------------------------------------------------------
# Arithmetic

class TestFrameFlexArithmetic(object):

    def test_df_add_td64_columnwise(self):
        # GH 22534 Check that column-wise addition broadcasts correctly
        dti = pd.date_range('2016-01-01', periods=10)
        tdi = pd.timedelta_range('1', periods=10)
        tser = pd.Series(tdi)
        df = pd.DataFrame({0: dti, 1: tdi})

        result = df.add(tser, axis=0)
        expected = pd.DataFrame({0: dti + tdi,
                                 1: tdi + tdi})
        tm.assert_frame_equal(result, expected)

    def test_df_add_flex_filled_mixed_dtypes(self):
        # GH 19611
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

    def test_arith_flex_frame(self, all_arithmetic_operators, float_frame,
                              mixed_float_frame):
        # one instance of parametrized fixture
        op = all_arithmetic_operators

        def f(x, y):
            # r-versions not in operator-stdlib; get op without "r" and invert
            if op.startswith('__r'):
                return getattr(operator, op.replace('__r', '__'))(y, x)
            return getattr(operator, op)(x, y)

        result = getattr(float_frame, op)(2 * float_frame)
        expected = f(float_frame, 2 * float_frame)
        tm.assert_frame_equal(result, expected)

        # vs mix float
        result = getattr(mixed_float_frame, op)(2 * mixed_float_frame)
        expected = f(mixed_float_frame, 2 * mixed_float_frame)
        tm.assert_frame_equal(result, expected)
        _check_mixed_float(result, dtype=dict(C=None))

    @pytest.mark.parametrize('op', ['__add__', '__sub__', '__mul__'])
    def test_arith_flex_frame_mixed(self, op, int_frame, mixed_int_frame,
                                    mixed_float_frame):
        f = getattr(operator, op)

        # vs mix int
        result = getattr(mixed_int_frame, op)(2 + mixed_int_frame)
        expected = f(mixed_int_frame, 2 + mixed_int_frame)

        # no overflow in the uint
        dtype = None
        if op in ['__sub__']:
            dtype = dict(B='uint64', C=None)
        elif op in ['__add__', '__mul__']:
            dtype = dict(C=None)
        tm.assert_frame_equal(result, expected)
        _check_mixed_int(result, dtype=dtype)

        # vs mix float
        result = getattr(mixed_float_frame, op)(2 * mixed_float_frame)
        expected = f(mixed_float_frame, 2 * mixed_float_frame)
        tm.assert_frame_equal(result, expected)
        _check_mixed_float(result, dtype=dict(C=None))

        # vs plain int
        result = getattr(int_frame, op)(2 * int_frame)
        expected = f(int_frame, 2 * int_frame)
        tm.assert_frame_equal(result, expected)

    def test_arith_flex_frame_raise(self, all_arithmetic_operators,
                                    float_frame):
        # one instance of parametrized fixture
        op = all_arithmetic_operators

        # Check that arrays with dim >= 3 raise
        for dim in range(3, 6):
            arr = np.ones((1,) * dim)
            msg = "Unable to coerce to Series/DataFrame"
            with tm.assert_raises_regex(ValueError, msg):
                getattr(float_frame, op)(arr)

    def test_arith_flex_frame_corner(self, float_frame):

        const_add = float_frame.add(1)
        tm.assert_frame_equal(const_add, float_frame + 1)

        # corner cases
        result = float_frame.add(float_frame[:0])
        tm.assert_frame_equal(result, float_frame * np.nan)

        result = float_frame[:0].add(float_frame)
        tm.assert_frame_equal(result, float_frame * np.nan)

        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            float_frame.add(float_frame.iloc[0], fill_value=3)

        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            float_frame.add(float_frame.iloc[0], axis='index', fill_value=3)

    def test_arith_flex_series(self, simple_frame):
        df = simple_frame

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

        # broadcasting issue in GH 7325
        df = pd.DataFrame(np.arange(3 * 2).reshape((3, 2)), dtype='int64')
        expected = pd.DataFrame([[np.nan, np.inf], [1.0, 1.5], [1.0, 1.25]])
        result = df.div(df[0], axis='index')
        tm.assert_frame_equal(result, expected)

        df = pd.DataFrame(np.arange(3 * 2).reshape((3, 2)), dtype='float64')
        expected = pd.DataFrame([[np.nan, np.inf], [1.0, 1.5], [1.0, 1.25]])
        result = df.div(df[0], axis='index')
        tm.assert_frame_equal(result, expected)

    def test_arith_flex_zero_len_raises(self):
        # GH 19522 passing fill_value to frame flex arith methods should
        # raise even in the zero-length special cases
        ser_len0 = pd.Series([])
        df_len0 = pd.DataFrame([], columns=['A', 'B'])
        df = pd.DataFrame([[1, 2], [3, 4]], columns=['A', 'B'])

        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            df.add(ser_len0, fill_value='E')

        with tm.assert_raises_regex(NotImplementedError, 'fill_value'):
            df_len0.sub(df['A'], axis=None, fill_value=3)


class TestFrameArithmetic(object):
    def test_df_add_2d_array_rowlike_broadcasts(self):
        # GH#23000
        arr = np.arange(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=[True, False], index=['A', 'B', 'C'])

        rowlike = arr[[1], :]  # shape --> (1, ncols)
        assert rowlike.shape == (1, df.shape[1])

        expected = pd.DataFrame([[2, 4],
                                 [4, 6],
                                 [6, 8]],
                                columns=df.columns, index=df.index,
                                # specify dtype explicitly to avoid failing
                                # on 32bit builds
                                dtype=arr.dtype)
        result = df + rowlike
        tm.assert_frame_equal(result, expected)
        result = rowlike + df
        tm.assert_frame_equal(result, expected)

    def test_df_add_2d_array_collike_broadcasts(self):
        # GH#23000
        arr = np.arange(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=[True, False], index=['A', 'B', 'C'])

        collike = arr[:, [1]]  # shape --> (nrows, 1)
        assert collike.shape == (df.shape[0], 1)

        expected = pd.DataFrame([[1, 2],
                                 [5, 6],
                                 [9, 10]],
                                columns=df.columns, index=df.index,
                                # specify dtype explicitly to avoid failing
                                # on 32bit builds
                                dtype=arr.dtype)
        result = df + collike
        tm.assert_frame_equal(result, expected)
        result = collike + df
        tm.assert_frame_equal(result, expected)

    def test_df_arith_2d_array_rowlike_broadcasts(self,
                                                  all_arithmetic_operators):
        # GH#23000
        opname = all_arithmetic_operators

        arr = np.arange(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=[True, False], index=['A', 'B', 'C'])

        rowlike = arr[[1], :]  # shape --> (1, ncols)
        assert rowlike.shape == (1, df.shape[1])

        exvals = [getattr(df.loc['A'], opname)(rowlike.squeeze()),
                  getattr(df.loc['B'], opname)(rowlike.squeeze()),
                  getattr(df.loc['C'], opname)(rowlike.squeeze())]

        expected = pd.DataFrame(exvals, columns=df.columns, index=df.index)

        if opname in ['__rmod__', '__rfloordiv__']:
            # exvals will have dtypes [f8, i8, i8] so expected will be
            #   all-f8, but the DataFrame operation will return mixed dtypes
            # use exvals[-1].dtype instead of "i8" for compat with 32-bit
            # systems/pythons
            expected[False] = expected[False].astype(exvals[-1].dtype)

        result = getattr(df, opname)(rowlike)
        tm.assert_frame_equal(result, expected)

    def test_df_arith_2d_array_collike_broadcasts(self,
                                                  all_arithmetic_operators):
        # GH#23000
        opname = all_arithmetic_operators

        arr = np.arange(6).reshape(3, 2)
        df = pd.DataFrame(arr, columns=[True, False], index=['A', 'B', 'C'])

        collike = arr[:, [1]]  # shape --> (nrows, 1)
        assert collike.shape == (df.shape[0], 1)

        exvals = {True: getattr(df[True], opname)(collike.squeeze()),
                  False: getattr(df[False], opname)(collike.squeeze())}

        dtype = None
        if opname in ['__rmod__', '__rfloordiv__']:
            # Series ops may return mixed int/float dtypes in cases where
            #   DataFrame op will return all-float.  So we upcast `expected`
            dtype = np.common_type(*[x.values for x in exvals.values()])

        expected = pd.DataFrame(exvals, columns=df.columns, index=df.index,
                                dtype=dtype)

        result = getattr(df, opname)(collike)
        tm.assert_frame_equal(result, expected)

    def test_df_bool_mul_int(self):
        # GH 22047, GH 22163 multiplication by 1 should result in int dtype,
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
        # GH#22696 Check that we don't dispatch to numpy implementation,
        # which treats int64 as m8[ns]
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
