import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "op, expected",
    [
        ["sum", np.int64(3)],
        ["prod", np.int64(2)],
        ["min", np.int64(1)],
        ["max", np.int64(2)],
        ["mean", np.float64(1.5)],
        ["median", np.float64(1.5)],
        ["var", np.float64(0.5)],
        ["std", np.float64(0.5**0.5)],
        ["skew", pd.NA],
        ["kurt", pd.NA],
        ["any", True],
        ["all", True],
    ],
)
def test_series_reductions(op, expected):
    ser = pd.Series([1, 2], dtype="Int64")
    result = getattr(ser, op)()
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "op, expected",
    [
        ["sum", pd.Series([3], index=["a"], dtype="Int64")],
        ["prod", pd.Series([2], index=["a"], dtype="Int64")],
        ["min", pd.Series([1], index=["a"], dtype="Int64")],
        ["max", pd.Series([2], index=["a"], dtype="Int64")],
        ["mean", pd.Series([1.5], index=["a"], dtype="Float64")],
        ["median", pd.Series([1.5], index=["a"], dtype="Float64")],
        ["var", pd.Series([0.5], index=["a"], dtype="Float64")],
        ["std", pd.Series([0.5**0.5], index=["a"], dtype="Float64")],
        ["skew", pd.Series([pd.NA], index=["a"], dtype="Float64")],
        ["kurt", pd.Series([pd.NA], index=["a"], dtype="Float64")],
        ["any", pd.Series([True], index=["a"], dtype="boolean")],
        ["all", pd.Series([True], index=["a"], dtype="boolean")],
    ],
)
def test_dataframe_reductions(op, expected):
    df = pd.DataFrame({"a": pd.array([1, 2], dtype="Int64")})
    result = getattr(df, op)()
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "op, expected",
    [
        ["sum", pd.array([1, 3], dtype="Int64")],
        ["prod", pd.array([1, 3], dtype="Int64")],
        ["min", pd.array([1, 3], dtype="Int64")],
        ["max", pd.array([1, 3], dtype="Int64")],
        ["mean", pd.array([1, 3], dtype="Float64")],
        ["median", pd.array([1, 3], dtype="Float64")],
        ["var", pd.array([pd.NA], dtype="Float64")],
        ["std", pd.array([pd.NA], dtype="Float64")],
        ["skew", pd.array([pd.NA], dtype="Float64")],
        ["any", pd.array([True, True], dtype="boolean")],
        ["all", pd.array([True, True], dtype="boolean")],
    ],
)
def test_groupby_reductions(op, expected):
    df = pd.DataFrame(
        {
            "A": ["a", "b", "b"],
            "B": pd.array([1, None, 3], dtype="Int64"),
        }
    )
    result = getattr(df.groupby("A"), op)()
    expected = pd.DataFrame(
        expected, index=pd.Index(["a", "b"], name="A"), columns=["B"]
    )

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "op, expected",
    [
        ["sum", pd.Series([4, 4], index=["B", "C"], dtype="Float64")],
        ["prod", pd.Series([3, 3], index=["B", "C"], dtype="Float64")],
        ["min", pd.Series([1, 1], index=["B", "C"], dtype="Float64")],
        ["max", pd.Series([3, 3], index=["B", "C"], dtype="Float64")],
        ["mean", pd.Series([2, 2], index=["B", "C"], dtype="Float64")],
        ["median", pd.Series([2, 2], index=["B", "C"], dtype="Float64")],
        ["var", pd.Series([2, 2], index=["B", "C"], dtype="Float64")],
        ["std", pd.Series([2**0.5, 2**0.5], index=["B", "C"], dtype="Float64")],
        ["skew", pd.Series([np.nan, pd.NA], index=["B", "C"], dtype="Float64")],
        ["kurt", pd.Series([np.nan, pd.NA], index=["B", "C"], dtype="Float64")],
        ["any", pd.Series([True, True, True], index=["A", "B", "C"], dtype="boolean")],
        ["all", pd.Series([True, True, True], index=["A", "B", "C"], dtype="boolean")],
    ],
)
def test_mixed_reductions(op, expected):
    df = pd.DataFrame(
        {
            "A": ["a", "b", "b"],
            "B": [1, None, 3],
            "C": pd.array([1, None, 3], dtype="Int64"),
        }
    )

    # series
    result = getattr(df.C, op)()
    tm.assert_equal(result, expected["C"])

    # frame
    if op in ["any", "all"]:
        result = getattr(df, op)()
    else:
        result = getattr(df, op)(numeric_only=True)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("op", ["skew", "kurt"])
def test_mixed_reduction_nan_not_influenced_by_nullable_column(op, using_nan_is_na):
    # GH#62024 - presence of nullable column C should not change
    # the result for non-nullable column B
    df = pd.DataFrame(
        {
            "B": [1, None, 3],
            "C": pd.array([1, None, 3], dtype="Int64"),
        }
    )
    result_mixed = getattr(df, op)(numeric_only=True)
    result_alone = getattr(df[["B"]], op)()

    if using_nan_is_na:
        # In default mode both are NA (NaN folded in Float64 result)
        assert result_mixed["B"] is pd.NA
    else:
        # In distinguish mode, B gives NaN in both contexts
        assert np.isnan(result_mixed["B"])
        assert np.isnan(result_alone["B"])
