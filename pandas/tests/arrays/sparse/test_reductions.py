import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.sparse import SparseArray


class TestReductions:
    @pytest.mark.parametrize(
        "data,pos,neg",
        [
            ([True, True, True], True, False),
            ([1, 2, 1], 1, 0),
            ([1.0, 2.0, 1.0], 1.0, 0.0),
        ],
    )
    def test_all(self, data, pos, neg):
        # GH#17570
        out = SparseArray(data).all()
        assert out

        out = SparseArray(data, fill_value=pos).all()
        assert out

        data[1] = neg
        out = SparseArray(data).all()
        assert not out

        out = SparseArray(data, fill_value=pos).all()
        assert not out

    @pytest.mark.parametrize(
        "data,pos,neg",
        [
            ([True, True, True], True, False),
            ([1, 2, 1], 1, 0),
            ([1.0, 2.0, 1.0], 1.0, 0.0),
        ],
    )
    def test_numpy_all(self, data, pos, neg):
        # GH#17570
        out = np.all(SparseArray(data))
        assert out

        out = np.all(SparseArray(data, fill_value=pos))
        assert out

        data[1] = neg
        out = np.all(SparseArray(data))
        assert not out

        out = np.all(SparseArray(data, fill_value=pos))
        assert not out

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.all(SparseArray(data), out=np.array([]))

    @pytest.mark.parametrize(
        "data,pos,neg",
        [
            ([False, True, False], True, False),
            ([0, 2, 0], 2, 0),
            ([0.0, 2.0, 0.0], 2.0, 0.0),
        ],
    )
    def test_any(self, data, pos, neg):
        # GH#17570
        out = SparseArray(data).any()
        assert out

        out = SparseArray(data, fill_value=pos).any()
        assert out

        data[1] = neg
        out = SparseArray(data).any()
        assert not out

        out = SparseArray(data, fill_value=pos).any()
        assert not out

    @pytest.mark.parametrize(
        "data,pos,neg",
        [
            ([False, True, False], True, False),
            ([0, 2, 0], 2, 0),
            ([0.0, 2.0, 0.0], 2.0, 0.0),
        ],
    )
    def test_numpy_any(self, data, pos, neg):
        # GH#17570
        out = np.any(SparseArray(data))
        assert out

        out = np.any(SparseArray(data, fill_value=pos))
        assert out

        data[1] = neg
        out = np.any(SparseArray(data))
        assert not out

        out = np.any(SparseArray(data, fill_value=pos))
        assert not out

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.any(SparseArray(data), out=out)

    def test_sum(self):
        data = np.arange(10).astype(float)
        out = SparseArray(data).sum()
        assert out == 45.0

        data[5] = np.nan
        out = SparseArray(data, fill_value=2).sum()
        assert out == 40.0

        out = SparseArray(data, fill_value=np.nan).sum()
        assert out == 40.0

    def test_sum_skipna(self):
        # GH#65478 sum honors skipna for SparseArray-backed reductions
        arr = SparseArray([1.0, np.nan, 3.0], fill_value=np.nan)
        assert arr.sum(skipna=True) == 4.0
        assert pd.isna(arr.sum(skipna=False))

        # a non-null fill value is not a missing value, so skipna=False
        # still returns the full sum
        arr = SparseArray([1, 2, 0, 3], fill_value=0)
        assert arr.sum(skipna=False) == 6

    @pytest.mark.parametrize(
        "arr",
        [[0, 1, np.nan, 1], [0, 1, 1]],
    )
    @pytest.mark.parametrize("fill_value", [0, 1, np.nan])
    @pytest.mark.parametrize("min_count, expected", [(3, 2), (4, np.nan)])
    def test_sum_min_count(self, arr, fill_value, min_count, expected):
        # GH#25777
        sparray = SparseArray(np.array(arr), fill_value=fill_value)
        result = sparray.sum(min_count=min_count)
        if np.isnan(expected):
            assert np.isnan(result)
        else:
            assert result == expected

    def test_bool_sum_min_count(self):
        spar_bool = SparseArray([False, True] * 5, dtype=np.bool_, fill_value=True)
        res = spar_bool.sum(min_count=1)
        assert res == 5
        res = spar_bool.sum(min_count=11)
        assert pd.isna(res)

    def test_numpy_sum(self):
        data = np.arange(10).astype(float)
        out = np.sum(SparseArray(data))
        assert out == 45.0

        data[5] = np.nan
        out = np.sum(SparseArray(data, fill_value=2))
        assert out == 40.0

        out = np.sum(SparseArray(data, fill_value=np.nan))
        assert out == 40.0

        msg = "the 'dtype' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.sum(SparseArray(data), dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.sum(SparseArray(data), out=out)

    def test_mean(self):
        data = np.arange(10).astype(float)
        out = SparseArray(data).mean()
        assert out == 4.5

        data[5] = np.nan
        out = SparseArray(data).mean()
        assert out == 40.0 / 9

        out = SparseArray(data).mean(skipna=True)
        assert out == 40.0 / 9

        out = SparseArray(data).mean(skipna=False)
        assert pd.isna(out)

        arr = SparseArray([1.0, np.nan, 3.0], fill_value=np.nan)
        out = arr.mean(skipna=True)
        assert out == 2.0

        out = arr.mean(skipna=False)
        assert pd.isna(out)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_mean_raises_for_unsupported_object_dtype_with_na(self, skipna):
        arr = SparseArray(["a", np.nan], dtype=pd.SparseDtype(object))

        msg = "unsupported operand type"
        with pytest.raises(TypeError, match=msg):
            arr.mean(skipna=skipna)

    def test_numpy_mean(self):
        data = np.arange(10).astype(float)
        out = np.mean(SparseArray(data))
        assert out == 4.5

        data[5] = np.nan
        out = np.mean(SparseArray(data))
        assert out == 40.0 / 9

        msg = "the 'dtype' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.mean(SparseArray(data), dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.mean(SparseArray(data), out=out)


class TestMinMax:
    @pytest.mark.parametrize(
        "raw_data,max_expected,min_expected",
        [
            (np.arange(5.0), [4], [0]),
            (-np.arange(5.0), [0], [-4]),
            (np.array([0, 1, 2, np.nan, 4]), [4], [0]),
            (np.array([np.nan] * 5), [np.nan], [np.nan]),
            (np.array([]), [np.nan], [np.nan]),
        ],
    )
    def test_nan_fill_value(self, raw_data, max_expected, min_expected):
        arr = SparseArray(raw_data)
        max_result = arr.max()
        min_result = arr.min()
        assert max_result in max_expected
        assert min_result in min_expected

        max_result = arr.max(skipna=False)
        min_result = arr.min(skipna=False)
        if np.isnan(raw_data).any():
            assert np.isnan(max_result)
            assert np.isnan(min_result)
        else:
            assert max_result in max_expected
            assert min_result in min_expected

    @pytest.mark.parametrize(
        "fill_value,max_expected,min_expected",
        [
            (100, 100, 0),
            (-100, 1, -100),
        ],
    )
    def test_fill_value(self, fill_value, max_expected, min_expected):
        arr = SparseArray(
            np.array([fill_value, 0, 1]), dtype=pd.SparseDtype("int", fill_value)
        )
        max_result = arr.max()
        assert max_result == max_expected

        min_result = arr.min()
        assert min_result == min_expected

    @pytest.mark.parametrize("method", ["min", "max"])
    @pytest.mark.parametrize("values", [[np.nan, 0, 1, 0], [np.nan, 0, 0]])
    def test_min_max_skipna_false_with_nonnull_fill_value(self, method, values):
        # GH#65478 skipna=False should see explicit NA values, even when
        # non-null fill-value gaps are also present.
        arr = SparseArray(values, fill_value=0)

        result = getattr(arr, method)(skipna=False)

        assert pd.isna(result)

    def test_only_fill_value(self):
        fv = 100
        arr = SparseArray(np.array([fv, fv, fv]), dtype=pd.SparseDtype("int", fv))
        assert len(arr._valid_sp_values) == 0

        assert arr.max() == fv
        assert arr.min() == fv
        assert arr.max(skipna=False) == fv
        assert arr.min(skipna=False) == fv

    @pytest.mark.parametrize("func", ["min", "max"])
    @pytest.mark.parametrize("data", [np.array([]), np.array([np.nan, np.nan])])
    @pytest.mark.parametrize(
        "dtype,expected",
        [
            (pd.SparseDtype(np.float64, np.nan), np.nan),
            (pd.SparseDtype(np.float64, 5.0), np.nan),
            (pd.SparseDtype("datetime64[ns]", pd.NaT), pd.NaT),
            (pd.SparseDtype("datetime64[ns]", pd.Timestamp("2018-05-05")), pd.NaT),
        ],
    )
    def test_na_value_if_no_valid_values(self, func, data, dtype, expected):
        arr = SparseArray(data, dtype=dtype)
        result = getattr(arr, func)()
        if expected is pd.NaT:
            # TODO: pin down whether we wrap datetime64("NaT")
            assert result is pd.NaT or np.isnat(result)
        else:
            assert np.isnan(result)


class TestArgmaxArgmin:
    @pytest.mark.parametrize(
        "arr,argmax_expected,argmin_expected",
        [
            (SparseArray([1, 2, 0, 1, 2]), 1, 2),
            (SparseArray([-1, -2, 0, -1, -2]), 2, 1),
            (SparseArray([np.nan, 1, 0, 0, np.nan, -1]), 1, 5),
            (SparseArray([np.nan, 1, 0, 0, np.nan, 2]), 5, 2),
            (SparseArray([np.nan, 1, 0, 0, np.nan, 2], fill_value=-1), 5, 2),
            (SparseArray([np.nan, 1, 0, 0, np.nan, 2], fill_value=0), 5, 2),
            (SparseArray([np.nan, 1, 0, 0, np.nan, 2], fill_value=1), 5, 2),
            (SparseArray([np.nan, 1, 0, 0, np.nan, 2], fill_value=2), 5, 2),
            (SparseArray([np.nan, 1, 0, 0, np.nan, 2], fill_value=3), 5, 2),
            (SparseArray([0] * 10 + [-1], fill_value=0), 0, 10),
            (SparseArray([0] * 10 + [-1], fill_value=-1), 0, 10),
            (SparseArray([0] * 10 + [-1], fill_value=1), 0, 10),
            (SparseArray([-1] + [0] * 10, fill_value=0), 1, 0),
            (SparseArray([1] + [0] * 10, fill_value=0), 0, 1),
            (SparseArray([-1] + [0] * 10, fill_value=-1), 1, 0),
            (SparseArray([1] + [0] * 10, fill_value=1), 0, 1),
        ],
    )
    def test_argmax_argmin(self, arr, argmax_expected, argmin_expected):
        argmax_result = arr.argmax()
        argmin_result = arr.argmin()
        assert argmax_result == argmax_expected
        assert argmin_result == argmin_expected

    @pytest.mark.parametrize("method", ["argmax", "argmin"])
    def test_empty_array(self, method):
        msg = f"attempt to get {method} of an empty sequence"
        arr = SparseArray([])
        with pytest.raises(ValueError, match=msg):
            getattr(arr, method)()


@pytest.mark.parametrize("values, expected", [([1, 2, 4], 7 / 3), ([1, 2, 3], 2.0)])
def test_frame_mean_int_column_not_truncated(values, expected):
    # GH#55123 the DataFrame path cast the reduction back to the column's
    #  dtype, truncating the mean of an integer column
    df = pd.DataFrame({"a": SparseArray(np.array(values), fill_value=0)})
    result = df.mean()
    expected = pd.Series([expected], index=["a"], dtype=pd.SparseDtype("float64", 0.0))
    tm.assert_series_equal(result, expected)


def test_frame_sum_bool_column_not_cast_back():
    # GH#55123 summing a bool column gave True instead of the count
    df = pd.DataFrame({"a": SparseArray([True, False, False], fill_value=False)})
    result = df.sum()
    expected = pd.Series([1], index=["a"], dtype=pd.SparseDtype(np.int_, 0))
    tm.assert_series_equal(result, expected)


def test_frame_sum_narrow_int_column_does_not_overflow():
    # GH#55123 casting the sum back to a narrow subtype wrapped around
    df = pd.DataFrame({"a": SparseArray(np.array([100, 100, 100], dtype="int8"))})
    result = df.sum()
    expected = pd.Series([300], index=["a"], dtype=pd.SparseDtype(np.int_, 0))
    tm.assert_series_equal(result, expected)


def test_frame_count_sparse_columns():
    # GH#55123 count went through the same cast and collapsed every sparse
    #  column to 1
    df = pd.DataFrame(
        {"a": SparseArray([1.0, 2.0, 3.0]), "b": SparseArray([1.0, np.nan, 3.0])}
    )
    result = df.count()
    expected = pd.Series([3, 2], index=["a", "b"], dtype="int64")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["min", "max"])
@pytest.mark.parametrize("subtype", ["int64", "bool"])
def test_frame_min_max_empty_column(method, subtype):
    # GH#55123 reducing an empty column gives NaN, which these subtypes cannot
    #  hold; casting back gave True for bool and raised for int64
    df = pd.DataFrame({"a": SparseArray(np.array([], dtype=subtype))})
    result = getattr(df, method)()
    expected = pd.Series([np.nan], index=["a"], dtype=pd.SparseDtype("float64", 0.0))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "arr, subtype",
    [
        (SparseArray([True, False, False], fill_value=False), "float64"),
        (SparseArray(np.array([1, 2, 4]), fill_value=0), "float64"),
        # a float32 subtype already holds the missing value, so it is retained
        (SparseArray(np.array([1, 2, 4], dtype="float32"), fill_value=0), "float32"),
    ],
)
def test_frame_sum_min_count_not_cast_to_column_dtype(arr, subtype):
    # GH#55123 the missing value produced by min_count was swallowed by the
    #  cast back to the column's dtype
    df = pd.DataFrame({"a": arr})
    result = df.sum(min_count=4)
    expected = pd.Series([np.nan], index=["a"], dtype=pd.SparseDtype(subtype, 0.0))
    tm.assert_series_equal(result, expected)


def test_frame_reduction_keeps_datetime64_dtype():
    # the widening must not kick in for a datetime64 subtype
    arr = SparseArray(np.array(["2020-01-01", "2020-01-03"], dtype="M8[ns]"))
    result = pd.DataFrame({"a": arr}).max()
    expected = pd.Series([pd.Timestamp("2020-01-03")], index=["a"], dtype=arr.dtype)
    tm.assert_series_equal(result, expected)
