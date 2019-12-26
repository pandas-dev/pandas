import operator

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import Categorical, DataFrame, MultiIndex, Series, Timestamp, isna
import pandas.util.testing as tm


class TestSeriesAnalytics:
    def test_argsort(self, datetime_series):
        self._check_accum_op("argsort", datetime_series, check_dtype=False)
        argsorted = datetime_series.argsort()
        assert issubclass(argsorted.dtype.type, np.integer)

        # GH 2967 (introduced bug in 0.11-dev I think)
        s = Series([Timestamp("201301{i:02d}".format(i=i)) for i in range(1, 6)])
        assert s.dtype == "datetime64[ns]"
        shifted = s.shift(-1)
        assert shifted.dtype == "datetime64[ns]"
        assert isna(shifted[4])

        result = s.argsort()
        expected = Series(range(5), dtype="int64")
        tm.assert_series_equal(result, expected)

        result = shifted.argsort()
        expected = Series(list(range(4)) + [-1], dtype="int64")
        tm.assert_series_equal(result, expected)

    def test_argsort_stable(self):
        s = Series(np.random.randint(0, 100, size=10000))
        mindexer = s.argsort(kind="mergesort")
        qindexer = s.argsort()

        mexpected = np.argsort(s.values, kind="mergesort")
        qexpected = np.argsort(s.values, kind="quicksort")

        tm.assert_series_equal(mindexer, Series(mexpected), check_dtype=False)
        tm.assert_series_equal(qindexer, Series(qexpected), check_dtype=False)
        msg = (
            r"ndarray Expected type <class 'numpy\.ndarray'>,"
            r" found <class 'pandas\.core\.series\.Series'> instead"
        )
        with pytest.raises(AssertionError, match=msg):
            tm.assert_numpy_array_equal(qindexer, mindexer)

    def _check_accum_op(self, name, datetime_series_, check_dtype=True):
        func = getattr(np, name)
        tm.assert_numpy_array_equal(
            func(datetime_series_).values,
            func(np.array(datetime_series_)),
            check_dtype=check_dtype,
        )

        # with missing values
        ts = datetime_series_.copy()
        ts[::2] = np.NaN

        result = func(ts)[1::2]
        expected = func(np.array(ts.dropna()))

        tm.assert_numpy_array_equal(result.values, expected, check_dtype=False)

    def test_compress(self):
        cond = [True, False, True, False, False]
        s = Series([1, -1, 5, 8, 7], index=list("abcde"), name="foo")
        expected = Series(s.values.compress(cond), index=list("ac"), name="foo")
        with tm.assert_produces_warning(FutureWarning):
            result = s.compress(cond)
        tm.assert_series_equal(result, expected)

    def test_numpy_compress(self):
        cond = [True, False, True, False, False]
        s = Series([1, -1, 5, 8, 7], index=list("abcde"), name="foo")
        expected = Series(s.values.compress(cond), index=list("ac"), name="foo")
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            tm.assert_series_equal(np.compress(cond, s), expected)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            msg = "the 'axis' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                np.compress(cond, s, axis=1)

            msg = "the 'out' parameter is not supported"
            with pytest.raises(ValueError, match=msg):
                np.compress(cond, s, out=s)

    def test_prod_numpy16_bug(self):
        s = Series([1.0, 1.0, 1.0], index=range(3))
        result = s.prod()

        assert not isinstance(result, Series)

    def test_dot(self):
        a = Series(np.random.randn(4), index=["p", "q", "r", "s"])
        b = DataFrame(
            np.random.randn(3, 4), index=["1", "2", "3"], columns=["p", "q", "r", "s"]
        ).T

        result = a.dot(b)
        expected = Series(np.dot(a.values, b.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        # Check index alignment
        b2 = b.reindex(index=reversed(b.index))
        result = a.dot(b)
        tm.assert_series_equal(result, expected)

        # Check ndarray argument
        result = a.dot(b.values)
        assert np.all(result == expected.values)
        tm.assert_almost_equal(a.dot(b["2"].values), expected["2"])

        # Check series argument
        tm.assert_almost_equal(a.dot(b["1"]), expected["1"])
        tm.assert_almost_equal(a.dot(b2["1"]), expected["1"])

        msg = r"Dot product shape mismatch, \(4,\) vs \(3,\)"
        # exception raised is of type Exception
        with pytest.raises(Exception, match=msg):
            a.dot(a.values[:3])
        msg = "matrices are not aligned"
        with pytest.raises(ValueError, match=msg):
            a.dot(b.T)

    def test_matmul(self):
        # matmul test is for GH #10259
        a = Series(np.random.randn(4), index=["p", "q", "r", "s"])
        b = DataFrame(
            np.random.randn(3, 4), index=["1", "2", "3"], columns=["p", "q", "r", "s"]
        ).T

        # Series @ DataFrame -> Series
        result = operator.matmul(a, b)
        expected = Series(np.dot(a.values, b.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        # DataFrame @ Series -> Series
        result = operator.matmul(b.T, a)
        expected = Series(np.dot(b.T.values, a.T.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        # Series @ Series -> scalar
        result = operator.matmul(a, a)
        expected = np.dot(a.values, a.values)
        tm.assert_almost_equal(result, expected)

        # GH 21530
        # vector (1D np.array) @ Series (__rmatmul__)
        result = operator.matmul(a.values, a)
        expected = np.dot(a.values, a.values)
        tm.assert_almost_equal(result, expected)

        # GH 21530
        # vector (1D list) @ Series (__rmatmul__)
        result = operator.matmul(a.values.tolist(), a)
        expected = np.dot(a.values, a.values)
        tm.assert_almost_equal(result, expected)

        # GH 21530
        # matrix (2D np.array) @ Series (__rmatmul__)
        result = operator.matmul(b.T.values, a)
        expected = np.dot(b.T.values, a.values)
        tm.assert_almost_equal(result, expected)

        # GH 21530
        # matrix (2D nested lists) @ Series (__rmatmul__)
        result = operator.matmul(b.T.values.tolist(), a)
        expected = np.dot(b.T.values, a.values)
        tm.assert_almost_equal(result, expected)

        # mixed dtype DataFrame @ Series
        a["p"] = int(a.p)
        result = operator.matmul(b.T, a)
        expected = Series(np.dot(b.T.values, a.T.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        # different dtypes DataFrame @ Series
        a = a.astype(int)
        result = operator.matmul(b.T, a)
        expected = Series(np.dot(b.T.values, a.T.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        msg = r"Dot product shape mismatch, \(4,\) vs \(3,\)"
        # exception raised is of type Exception
        with pytest.raises(Exception, match=msg):
            a.dot(a.values[:3])
        msg = "matrices are not aligned"
        with pytest.raises(ValueError, match=msg):
            a.dot(b.T)

    def test_ptp(self):
        # GH21614
        N = 1000
        arr = np.random.randn(N)
        ser = Series(arr)
        assert np.ptp(ser) == np.ptp(arr)

    def test_repeat(self):
        s = Series(np.random.randn(3), index=["a", "b", "c"])

        reps = s.repeat(5)
        exp = Series(s.values.repeat(5), index=s.index.values.repeat(5))
        tm.assert_series_equal(reps, exp)

        to_rep = [2, 3, 4]
        reps = s.repeat(to_rep)
        exp = Series(s.values.repeat(to_rep), index=s.index.values.repeat(to_rep))
        tm.assert_series_equal(reps, exp)

    def test_numpy_repeat(self):
        s = Series(np.arange(3), name="x")
        expected = Series(s.values.repeat(2), name="x", index=s.index.values.repeat(2))
        tm.assert_series_equal(np.repeat(s, 2), expected)

        msg = "the 'axis' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.repeat(s, 2, axis=0)

    def test_is_monotonic(self):

        s = Series(np.random.randint(0, 10, size=1000))
        assert not s.is_monotonic
        s = Series(np.arange(1000))
        assert s.is_monotonic is True
        assert s.is_monotonic_increasing is True
        s = Series(np.arange(1000, 0, -1))
        assert s.is_monotonic_decreasing is True

        s = Series(pd.date_range("20130101", periods=10))
        assert s.is_monotonic is True
        assert s.is_monotonic_increasing is True
        s = Series(list(reversed(s.tolist())))
        assert s.is_monotonic is False
        assert s.is_monotonic_decreasing is True

    def test_apply_categorical(self):
        values = pd.Categorical(list("ABBABCD"), categories=list("DCBA"), ordered=True)
        s = pd.Series(values, name="XX", index=list("abcdefg"))
        result = s.apply(lambda x: x.lower())

        # should be categorical dtype when the number of categories are
        # the same
        values = pd.Categorical(list("abbabcd"), categories=list("dcba"), ordered=True)
        exp = pd.Series(values, name="XX", index=list("abcdefg"))
        tm.assert_series_equal(result, exp)
        tm.assert_categorical_equal(result.values, exp.values)

        result = s.apply(lambda x: "A")
        exp = pd.Series(["A"] * 7, name="XX", index=list("abcdefg"))
        tm.assert_series_equal(result, exp)
        assert result.dtype == np.object

    def test_unstack(self):

        index = MultiIndex(
            levels=[["bar", "foo"], ["one", "three", "two"]],
            codes=[[1, 1, 0, 0], [0, 1, 0, 2]],
        )

        s = Series(np.arange(4.0), index=index)
        unstacked = s.unstack()

        expected = DataFrame(
            [[2.0, np.nan, 3.0], [0.0, 1.0, np.nan]],
            index=["bar", "foo"],
            columns=["one", "three", "two"],
        )

        tm.assert_frame_equal(unstacked, expected)

        unstacked = s.unstack(level=0)
        tm.assert_frame_equal(unstacked, expected.T)

        index = MultiIndex(
            levels=[["bar"], ["one", "two", "three"], [0, 1]],
            codes=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1]],
        )
        s = Series(np.random.randn(6), index=index)
        exp_index = MultiIndex(
            levels=[["one", "two", "three"], [0, 1]],
            codes=[[0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1]],
        )
        expected = DataFrame({"bar": s.values}, index=exp_index).sort_index(level=0)
        unstacked = s.unstack(0).sort_index()
        tm.assert_frame_equal(unstacked, expected)

        # GH5873
        idx = pd.MultiIndex.from_arrays([[101, 102], [3.5, np.nan]])
        ts = pd.Series([1, 2], index=idx)
        left = ts.unstack()
        right = DataFrame(
            [[np.nan, 1], [2, np.nan]], index=[101, 102], columns=[np.nan, 3.5]
        )
        tm.assert_frame_equal(left, right)

        idx = pd.MultiIndex.from_arrays(
            [
                ["cat", "cat", "cat", "dog", "dog"],
                ["a", "a", "b", "a", "b"],
                [1, 2, 1, 1, np.nan],
            ]
        )
        ts = pd.Series([1.0, 1.1, 1.2, 1.3, 1.4], index=idx)
        right = DataFrame(
            [[1.0, 1.3], [1.1, np.nan], [np.nan, 1.4], [1.2, np.nan]],
            columns=["cat", "dog"],
        )
        tpls = [("a", 1), ("a", 2), ("b", np.nan), ("b", 1)]
        right.index = pd.MultiIndex.from_tuples(tpls)
        tm.assert_frame_equal(ts.unstack(level=0), right)

    @pytest.mark.parametrize("func", [np.any, np.all])
    @pytest.mark.parametrize("kwargs", [dict(keepdims=True), dict(out=object())])
    @td.skip_if_np_lt("1.15")
    def test_validate_any_all_out_keepdims_raises(self, kwargs, func):
        s = pd.Series([1, 2])
        param = list(kwargs)[0]
        name = func.__name__

        msg = (
            r"the '{arg}' parameter is not "
            r"supported in the pandas "
            r"implementation of {fname}\(\)"
        ).format(arg=param, fname=name)
        with pytest.raises(ValueError, match=msg):
            func(s, **kwargs)

    @td.skip_if_np_lt("1.15")
    def test_validate_sum_initial(self):
        s = pd.Series([1, 2])
        msg = (
            r"the 'initial' parameter is not "
            r"supported in the pandas "
            r"implementation of sum\(\)"
        )
        with pytest.raises(ValueError, match=msg):
            np.sum(s, initial=10)

    def test_validate_median_initial(self):
        s = pd.Series([1, 2])
        msg = (
            r"the 'overwrite_input' parameter is not "
            r"supported in the pandas "
            r"implementation of median\(\)"
        )
        with pytest.raises(ValueError, match=msg):
            # It seems like np.median doesn't dispatch, so we use the
            # method instead of the ufunc.
            s.median(overwrite_input=True)

    @td.skip_if_np_lt("1.15")
    def test_validate_stat_keepdims(self):
        s = pd.Series([1, 2])
        msg = (
            r"the 'keepdims' parameter is not "
            r"supported in the pandas "
            r"implementation of sum\(\)"
        )
        with pytest.raises(ValueError, match=msg):
            np.sum(s, keepdims=True)


class TestCategoricalSeriesAnalytics:
    @pytest.mark.parametrize(
        "dtype",
        ["int_", "uint", "float_", "unicode_", "timedelta64[h]", "datetime64[D]"],
    )
    def test_drop_duplicates_categorical_non_bool(self, dtype, ordered_fixture):
        cat_array = np.array([1, 2, 3, 4, 5], dtype=np.dtype(dtype))

        # Test case 1
        input1 = np.array([1, 2, 3, 3], dtype=np.dtype(dtype))
        tc1 = Series(Categorical(input1, categories=cat_array, ordered=ordered_fixture))
        if dtype == "datetime64[D]":
            # pre-empty flaky xfail, tc1 values are seemingly-random
            if not (np.array(tc1) == input1).all():
                pytest.xfail(reason="GH#7996")

        expected = Series([False, False, False, True])
        tm.assert_series_equal(tc1.duplicated(), expected)
        tm.assert_series_equal(tc1.drop_duplicates(), tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        expected = Series([False, False, True, False])
        tm.assert_series_equal(tc1.duplicated(keep="last"), expected)
        tm.assert_series_equal(tc1.drop_duplicates(keep="last"), tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(keep="last", inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        expected = Series([False, False, True, True])
        tm.assert_series_equal(tc1.duplicated(keep=False), expected)
        tm.assert_series_equal(tc1.drop_duplicates(keep=False), tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        # Test case 2
        input2 = np.array([1, 2, 3, 5, 3, 2, 4], dtype=np.dtype(dtype))
        tc2 = Series(Categorical(input2, categories=cat_array, ordered=ordered_fixture))
        if dtype == "datetime64[D]":
            # pre-empty flaky xfail, tc2 values are seemingly-random
            if not (np.array(tc2) == input2).all():
                pytest.xfail(reason="GH#7996")

        expected = Series([False, False, False, False, True, True, False])
        tm.assert_series_equal(tc2.duplicated(), expected)
        tm.assert_series_equal(tc2.drop_duplicates(), tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

        expected = Series([False, True, True, False, False, False, False])
        tm.assert_series_equal(tc2.duplicated(keep="last"), expected)
        tm.assert_series_equal(tc2.drop_duplicates(keep="last"), tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(keep="last", inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

        expected = Series([False, True, True, False, True, True, False])
        tm.assert_series_equal(tc2.duplicated(keep=False), expected)
        tm.assert_series_equal(tc2.drop_duplicates(keep=False), tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

    def test_drop_duplicates_categorical_bool(self, ordered_fixture):
        tc = Series(
            Categorical(
                [True, False, True, False],
                categories=[True, False],
                ordered=ordered_fixture,
            )
        )

        expected = Series([False, False, True, True])
        tm.assert_series_equal(tc.duplicated(), expected)
        tm.assert_series_equal(tc.drop_duplicates(), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc[~expected])

        expected = Series([True, True, False, False])
        tm.assert_series_equal(tc.duplicated(keep="last"), expected)
        tm.assert_series_equal(tc.drop_duplicates(keep="last"), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(keep="last", inplace=True)
        tm.assert_series_equal(sc, tc[~expected])

        expected = Series([True, True, True, True])
        tm.assert_series_equal(tc.duplicated(keep=False), expected)
        tm.assert_series_equal(tc.drop_duplicates(keep=False), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc[~expected])
