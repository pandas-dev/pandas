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

    @pytest.mark.parametrize("fill_value", [np.nan, 0.0])
    @pytest.mark.parametrize("method", ["argmax", "argmin"])
    def test_skipna_false_with_na_raises(self, method, fill_value):
        # GH#68387 _reduce reduced self.dropna(), renumbering the elements.  Both
        #  fill values hit it: the NA is a gap under one, a stored value under the
        #  other.
        arr = SparseArray([np.nan, 3.0, 1.0], fill_value=fill_value)
        msg = "Encountered an NA value with skipna=False"
        with pytest.raises(ValueError, match=msg):
            getattr(arr, method)(skipna=False)
        with pytest.raises(ValueError, match=msg):
            arr._reduce(method, skipna=False)

    @pytest.mark.parametrize(
        "data,fill_value,argmax_expected,argmin_expected",
        [
            ([np.nan, 0.0, 0.0], 0.0, 1, 1),
            ([0.0, np.nan, 0.0], 0.0, 0, 0),
            ([0.0, 0.0, np.nan], 0.0, 0, 0),
            ([np.nan, 2.0, np.nan], 2.0, 1, 1),
            ([np.nan, 3.0, np.nan, 3.0], 3.0, 1, 1),
        ],
    )
    def test_argmax_argmin_all_stored_values_na(
        self, data, fill_value, argmax_expected, argmin_expected
    ):
        # GH#68462 the fill value is the extremum when every stored value is NA
        arr = SparseArray(data, fill_value=fill_value)
        assert arr.argmax() == argmax_expected
        assert arr.argmin() == argmin_expected

        ser = pd.Series(data)
        assert arr.argmax() == ser.argmax()
        assert arr.argmin() == ser.argmin()

    def test_argmax_argmin_only_fill_value(self):
        # GH#68462 the argmin/argmax analogue of test_only_fill_value: nothing is
        #  stored at all, so the fill value is the answer at position 0
        fv = 100
        arr = SparseArray(np.array([fv, fv, fv]), dtype=pd.SparseDtype("int", fv))
        assert arr.sp_index.npoints == 0

        assert arr.argmax() == 0
        assert arr.argmin() == 0
        assert arr.argmax(skipna=False) == 0
        assert arr.argmin(skipna=False) == 0

    @pytest.mark.parametrize("method", ["argmax", "argmin"])
    @pytest.mark.parametrize("fill_value", [np.nan, 0.0])
    def test_all_na_still_raises(self, method, fill_value):
        # GH#68462 the two ways the fill value fails to be a candidate: it is NA,
        #  or every position is stored so it holds none
        msg = "Encountered all NA values"
        arr = SparseArray([np.nan, np.nan], fill_value=fill_value)
        with pytest.raises(ValueError, match=msg):
            getattr(arr, method)()

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


@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize(
    "name, kwargs",
    [
        ("prod", {}),
        # 4 is satisfied by the 5 valid values, 6 is not
        ("prod", {"min_count": 4}),
        ("prod", {"min_count": 6}),
        ("median", {}),
        ("var", {}),
        ("var", {"ddof": 0}),
        ("std", {}),
        ("std", {"ddof": 0}),
        ("sem", {}),
        ("sem", {"ddof": 0}),
        ("skew", {}),
        ("kurt", {}),
    ],
)
@pytest.mark.parametrize("fill_value", [0, np.nan])
def test_reductions_without_sparse_kernel_match_dense(name, kwargs, skipna, fill_value):
    # GH#68194 these raised "cannot perform <name> with type Sparse[...]"
    values = [0.0, 1.0, np.nan, -2.0, 4.0, 0.0]
    arr = SparseArray(np.array(values), fill_value=fill_value)
    expected = getattr(pd.Series(values), name)(skipna=skipna, **kwargs)

    assert getattr(arr, name)(skipna=skipna, **kwargs) == pytest.approx(
        expected, nan_ok=True
    )
    result = getattr(pd.Series(arr), name)(skipna=skipna, **kwargs)
    assert result == pytest.approx(expected, nan_ok=True)


def test_frame_prod():
    # GH#68194 DataFrame.prod raised for every SparseDtype column
    df = pd.DataFrame({"a": SparseArray(np.array([1, 2, 3]), fill_value=0)})
    result = df.prod()
    expected = pd.Series([6], index=["a"], dtype=pd.SparseDtype(np.int_, 0))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize("name", ["median", "std"])
@pytest.mark.parametrize("unit", ["M8[us]", "m8[us]"])
def test_datetimelike_reductions_match_dense(name, unit, skipna):
    # GH#68194 median and std are the only two of the newly-routed reductions
    #  nanops accepts for a datetimelike subtype
    values = np.array([1, 3, "NaT", 5], dtype=unit)
    arr = SparseArray(values)

    result = getattr(pd.Series(arr), name)(skipna=skipna)
    expected = getattr(pd.Series(values), name)(skipna=skipna)
    if pd.isna(expected):
        assert pd.isna(result)
    else:
        assert result == expected

    frame_result = getattr(pd.DataFrame({"a": arr}), name)(skipna=skipna)
    frame_expected = getattr(pd.DataFrame({"a": values}), name)(skipna=skipna)
    tm.assert_series_equal(
        frame_result, frame_expected.astype(pd.SparseDtype(frame_expected.dtype))
    )


def test_frame_std_datetime64_widens_to_timedelta64():
    # GH#68194 the keepdims widening skips datetimelike subtypes because
    #  min/max/median stay closed over them; std does not
    arr = SparseArray(
        np.array(["2020-01-01", "2020-01-03", "2020-01-05"], dtype="M8[s]")
    )
    result = pd.DataFrame({"a": arr}).std()
    expected = pd.Series(
        [pd.Timedelta("2 days")], index=["a"], dtype=pd.SparseDtype("m8[s]")
    )
    tm.assert_series_equal(result, expected)


# skew/kurt cast to float64 inside nanops, on the dense path too
@pytest.mark.filterwarnings(
    "ignore:Casting complex values:numpy.exceptions.ComplexWarning"
)
@pytest.mark.parametrize("fill_value", [np.nan, 1 + 1j])
@pytest.mark.parametrize("name", ["var", "std", "sem", "skew", "kurt"])
def test_frame_complex_reduction_narrows_to_real(name, fill_value):
    # GH#68194 these map complex to real, which np.result_type widens straight back
    values = np.array([1 + 2j, 3 - 1j, 5 + 0j])
    arr = SparseArray(values, dtype=pd.SparseDtype("complex128", fill_value))

    result = getattr(pd.DataFrame({"a": arr}), name)()
    expected = getattr(pd.DataFrame({"a": values}), name)()
    tm.assert_series_equal(result, expected.astype(pd.SparseDtype(expected.dtype)))


@pytest.mark.parametrize("kwargs", [{"skipna": False}, {"min_count": 10}])
@pytest.mark.parametrize("name", ["sum", "prod", "mean", "min", "max"])
def test_frame_complex_reduction_na_keeps_complex(name, kwargs):
    # GH#68194 these are closed over complex, so the real NaN that min_count and
    #  skipna=False produce must not narrow the column to float64
    if name in ("mean", "min", "max") and "min_count" in kwargs:
        pytest.skip(f"{name} takes no min_count")

    values = np.array([1 + 2j, 3 - 1j, np.nan])
    arr = SparseArray(values, dtype=pd.SparseDtype("complex128", np.nan))

    result = getattr(pd.DataFrame({"a": arr}), name)(**kwargs)
    expected = getattr(pd.DataFrame({"a": values}), name)(**kwargs)
    tm.assert_series_equal(result, expected.astype(pd.SparseDtype(expected.dtype)))


@pytest.mark.parametrize(
    "name, kwargs",
    [("median", {"skipna": False}), ("sem", {"ddof": 5}), ("std", {"ddof": 5})],
)
@pytest.mark.parametrize("subtype", ["float32", "complex128"])
def test_frame_nan_result_dtype_matches_dense(name, kwargs, subtype):
    # GH#68487 the keepdims widening follows whatever scalar nanops hands back,
    #  so a NaN that did not carry the column's own dtype made this disagree
    #  with dense, in either direction
    values = np.array([1, np.nan, 3], dtype=subtype)
    arr = SparseArray(values, dtype=pd.SparseDtype(subtype, np.nan))

    result = getattr(pd.DataFrame({"a": arr}), name)(**kwargs)
    expected = getattr(pd.DataFrame({"a": values}), name)(**kwargs)
    tm.assert_series_equal(result, expected.astype(pd.SparseDtype(expected.dtype)))


def test_describe():
    # GH#68194 describe computes std, so it raised for every SparseDtype Series.
    #  check_dtype=False: describe gives every EA-backed Series a non-numpy dtype
    values = np.array([1.0, 2.0, 3.0, 4.0])
    result = pd.Series(SparseArray(values)).describe()
    tm.assert_series_equal(result, pd.Series(values).describe(), check_dtype=False)

    frame_result = pd.DataFrame({"a": SparseArray(values)}).describe()
    tm.assert_frame_equal(
        frame_result, pd.DataFrame({"a": values}).describe(), check_dtype=False
    )


@pytest.mark.filterwarnings(
    "ignore:Casting complex values:numpy.exceptions.ComplexWarning"
)
@pytest.mark.parametrize(
    "subtype,fill_value,dense_dtype",
    [
        ("int64", np.nan, "float64"),
        ("uint8", np.nan, "float64"),
        ("bool", np.nan, "float64"),
        # pd.NA fits no numpy subtype at all, not even a float one
        ("float64", pd.NA, "float64"),
        ("float32", pd.NA, "float32"),
        ("complex128", pd.NA, "complex128"),
    ],
)
@pytest.mark.parametrize(
    "name", ["prod", "median", "var", "std", "sem", "skew", "kurt"]
)
@pytest.mark.parametrize("skipna", [True, False])
def test_reductions_with_na_fill_value_match_dense(
    subtype, fill_value, dense_dtype, name, skipna
):
    # GH#68194 densifying to the subtype read the NA gaps back as 0/False, or
    #  raised outright for a fill value the subtype cannot hold at all
    arr = SparseArray([1.0, np.nan, 5.0, np.nan]).astype(
        pd.SparseDtype(subtype, fill_value)
    )
    dense = pd.Series(
        [np.nan if pd.isna(value) else value for value in arr], dtype=dense_dtype
    )

    expected = getattr(dense, name)(skipna=skipna)
    assert getattr(arr, name)(skipna=skipna) == pytest.approx(expected, nan_ok=True)
    assert getattr(pd.Series(arr), name)(skipna=skipna) == pytest.approx(
        expected, nan_ok=True
    )


@pytest.mark.parametrize("name", ["median", "std"])
def test_reduction_keeps_sub_microsecond_fill_value(name):
    # GH#68194 np.full routes a Timedelta fill value through the stdlib datetime
    #  protocol, which floors it to microseconds
    values = np.array([1000, 2500, 2500, 4000], dtype="m8[ns]")
    arr = SparseArray(values, fill_value=pd.Timedelta("2500ns"))
    assert arr.sp_index.ngaps == 2
    assert getattr(pd.Series(arr), name)() == getattr(pd.Series(values), name)()


def test_multiply_reduce_includes_fill_value():
    # GH#68194 with no public prod, np.multiply.reduce fell through to
    #  __array_ufunc__, which reduced the stored values and dropped the gaps
    arr = SparseArray([1, 0, 2, 0, 3], fill_value=0)
    assert np.multiply.reduce(arr) == np.multiply.reduce(arr.to_dense())
    assert np.prod(arr) == np.prod(arr.to_dense())


@pytest.mark.parametrize(
    "ufunc, data, fill_value",
    [
        (np.logaddexp, [1.0, 0.0, 2.0, 0.0, 3.0], 0.0),
        (np.bitwise_or, [1, 4, 0, 0], 4),
        (np.gcd, [12, 15, 8], 15),
        # subtypes np.asarray would widen or objectify; see SparseArray._densify
        (np.bitwise_or, np.array([1, 4, 0, 0], dtype="uint64"), 4),
        (np.left_shift, np.array([1, 3, 4, 4], dtype="int8"), 4),
        (np.fmax, np.array([1000, 2500, 2500, 4000], "m8[ns]"), pd.Timedelta("2500ns")),
    ],
)
def test_unaliased_ufunc_reduce_includes_fill_value(ufunc, data, fill_value):
    # GH#68453 a ufunc outside arraylike.REDUCTION_ALIASES has no named method
    #  to dispatch to, and __array_ufunc__ reduced the stored values alone
    arr = SparseArray(data, fill_value=fill_value)
    assert arr.sp_index.ngaps
    result = ufunc.reduce(arr)
    # reduce the original input, not to_dense(), which floors a sub-microsecond
    #  Timedelta fill value; see test_reduction_keeps_sub_microsecond_fill_value
    # dtype= keeps the comparison off the platform default int width
    expected = ufunc.reduce(np.asarray(data, dtype=arr.dtype.subtype))
    assert result == expected
    assert result.dtype == expected.dtype


@pytest.mark.parametrize("name", ["any", "all"])
@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize(
    "data,fill_value",
    [
        ([0.0, np.nan], np.nan),  # NA is the fill value
        ([0.0, np.nan], 0.0),  # NA is a stored value
        ([1.0, np.nan], np.nan),
        ([np.nan, np.nan], np.nan),  # all-NA
        ([0.0, 0.0], 0.0),  # no NA
        ([0, 1], 0),  # int64: a subtype that cannot hold NA, so masking is skipped
        # object keeps a falsy None; the only case where `all` changes
        (np.array([1, None], dtype=object), 0),
    ],
)
def test_any_all_skipna(name, skipna, data, fill_value):
    # GH#68390 the methods took no skipna at all, and counted NA as a truthy
    #  value even under the default skipna=True
    arr = SparseArray(data, fill_value=fill_value)
    expected = getattr(pd.Series(arr.to_dense()), name)(skipna=skipna)

    assert getattr(arr, name)(skipna=skipna) == expected
    assert arr._reduce(name, skipna=skipna) == expected
    assert getattr(pd.Series(arr), name)(skipna=skipna) == expected
    assert getattr(pd.DataFrame({"A": arr}), name)(skipna=skipna)["A"] == expected


def test_any_all_na_fill_value():
    # GH#68390 a nullable string Series sparsifies to an object subtype with a
    #  pd.NA fill, which is what reaches the fill-value gate; dense raises here too
    arr = SparseArray(pd.array(["", None], dtype="string"))
    assert arr.fill_value is pd.NA

    assert not arr.any()
    assert not arr.all()

    msg = "boolean value of NA is ambiguous"
    with pytest.raises(TypeError, match=msg):
        arr.any(skipna=False)
    with pytest.raises(TypeError, match=msg):
        arr.all(skipna=False)


@pytest.mark.parametrize("name", ["any", "all"])
@pytest.mark.parametrize("unit", ["M8[s]", "M8[ns]"])
def test_any_all_datetime64_raises(name, unit):
    # GH#68438 sparse answered for a datetime64 subtype where dense raises
    values = np.array(["2020-01-01", "NaT"], dtype=unit)
    arr = SparseArray(values)
    msg = f"'{name}' with datetime64 dtypes is not supported"

    with pytest.raises(TypeError, match=msg):
        getattr(arr, name)()
    with pytest.raises(TypeError, match=msg):
        arr._reduce(name, keepdims=True)
    with pytest.raises(TypeError, match=msg):
        getattr(pd.Series(arr), name)()
    with pytest.raises(TypeError, match=msg):
        getattr(pd.DataFrame({"a": arr}), name)()


@pytest.mark.parametrize("name,value", [("any", True), ("all", False)])
@pytest.mark.parametrize("subtype", ["m8[s]", "int64", "float64", "bool"])
def test_any_all_keepdims_is_boolean(name, value, subtype):
    # GH#68438 the keepdims wrapper boxed the bool in self.dtype, so a timedelta64
    #  column came back as a 1ns Timedelta rather than True
    arr = SparseArray(np.array([1, 0], dtype=subtype))

    result = arr._reduce(name, keepdims=True)
    expected = SparseArray([value], dtype=pd.SparseDtype(bool))
    tm.assert_sp_array_equal(result, expected)


def test_numpy_any_all_skip_na():
    # GH#68390 these skip NA like np.sum and np.mean. NaN is truthy and None
    #  falsy, so each entry point needs its own array to be decisive
    assert not np.any(SparseArray([0.0, np.nan]))
    assert np.all(SparseArray(np.array([1, None], dtype=object), fill_value=0))


def test_numpy_std_var_skip_na():
    # GH#68194 np.sum and np.mean already skipped NA here; with no public std/var
    #  numpy densified instead and propagated it
    arr = SparseArray(np.array([1.0, 2.0, np.nan, 4.0]))
    dense = np.asarray(arr)
    assert np.std(arr) == pytest.approx(np.nanstd(dense))
    assert np.var(arr) == pytest.approx(np.nanvar(dense))
    assert np.std(arr, ddof=1) == pytest.approx(np.nanstd(dense, ddof=1))


@pytest.mark.parametrize("sparsify_min", [False, True])
@pytest.mark.parametrize("order", [[5, 3, 9], [3, 5, 9]])
@pytest.mark.parametrize("subtype", ["float64", "complex128", "m8[ns]", "M8[ns]"])
def test_frame_idxmin_idxmax(subtype, order, sparsify_min):
    # GH#68387 the keepdims wrapper cast the position back to the column's dtype,
    #  which for a complex column also warned about discarding the imaginary part.
    #  sparsify_min routes _argmin_argmax through its _first_fill_value_loc branch.
    values = np.array(order).astype(subtype)
    arr = SparseArray(values, fill_value=values.min() if sparsify_min else None)
    assert arr._reduce("argmin", keepdims=True).dtype.subtype == np.intp

    index = pd.Index(list("xyz"))
    df = pd.DataFrame({"a": arr}, index=index)
    with tm.assert_produces_warning(None):
        assert df.idxmin()["a"] == index[values.argmin()]
        assert df.idxmax()["a"] == index[values.argmax()]


@pytest.mark.parametrize("fill_value", [np.nan, 0.0])
@pytest.mark.parametrize("method", ["idxmin", "idxmax"])
def test_frame_idxmin_idxmax_skipna_false(method, fill_value):
    # GH#68387 the position returned indexed the column with its NAs dropped
    df = pd.DataFrame({"a": SparseArray([np.nan, 3.0, 1.0], fill_value=fill_value)})
    msg = "Encountered an NA value with skipna=False"
    with pytest.raises(ValueError, match=msg):
        getattr(df, method)(skipna=False)


@pytest.mark.parametrize("method", ["idxmax", "idxmin"])
def test_frame_idxmax_idxmin_all_stored_values_na(method):
    # GH#68462 the fill value is the extremum when no stored value is non-NA
    arr = SparseArray([np.nan, 0.0, 0.0], fill_value=0.0)
    result = getattr(pd.DataFrame({"a": arr}), method)()
    expected = pd.Series([1], index=["a"])
    tm.assert_series_equal(result, expected)
