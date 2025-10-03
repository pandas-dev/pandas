from enum import (
    Enum,
    auto,
)

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_compare_axis(axis):
    # GH#30429
    s1 = pd.Series(["a", "b", "c"])
    s2 = pd.Series(["x", "b", "z"])

    result = s1.compare(s2, align_axis=axis)

    if axis in (1, "columns"):
        indices = pd.Index([0, 2])
        columns = pd.Index(["self", "other"])
        expected = pd.DataFrame(
            [["a", "x"], ["c", "z"]], index=indices, columns=columns
        )
        tm.assert_frame_equal(result, expected)
    else:
        indices = pd.MultiIndex.from_product([[0, 2], ["self", "other"]])
        expected = pd.Series(["a", "x", "c", "z"], index=indices)
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "keep_shape, keep_equal",
    [
        (True, False),
        (False, True),
        (True, True),
        # False, False case is already covered in test_compare_axis
    ],
)
def test_compare_various_formats(keep_shape, keep_equal):
    s1 = pd.Series(["a", "b", "c"])
    s2 = pd.Series(["x", "b", "z"])

    result = s1.compare(s2, keep_shape=keep_shape, keep_equal=keep_equal)

    if keep_shape:
        indices = pd.Index([0, 1, 2])
        columns = pd.Index(["self", "other"])
        if keep_equal:
            expected = pd.DataFrame(
                [["a", "x"], ["b", "b"], ["c", "z"]], index=indices, columns=columns
            )
        else:
            expected = pd.DataFrame(
                [["a", "x"], [np.nan, np.nan], ["c", "z"]],
                index=indices,
                columns=columns,
            )
    else:
        indices = pd.Index([0, 2])
        columns = pd.Index(["self", "other"])
        expected = pd.DataFrame(
            [["a", "x"], ["c", "z"]], index=indices, columns=columns
        )
    tm.assert_frame_equal(result, expected)


def test_compare_with_equal_nulls():
    # We want to make sure two NaNs are considered the same
    # and dropped where applicable
    s1 = pd.Series(["a", "b", np.nan])
    s2 = pd.Series(["x", "b", np.nan])

    result = s1.compare(s2)
    expected = pd.DataFrame([["a", "x"]], columns=["self", "other"])
    tm.assert_frame_equal(result, expected)


def test_compare_with_non_equal_nulls():
    # We want to make sure the relevant NaNs do not get dropped
    s1 = pd.Series(["a", "b", "c"])
    s2 = pd.Series(["x", "b", np.nan])

    result = s1.compare(s2, align_axis=0)

    indices = pd.MultiIndex.from_product([[0, 2], ["self", "other"]])
    expected = pd.Series(["a", "x", "c", np.nan], index=indices)
    tm.assert_series_equal(result, expected)


def test_compare_multi_index():
    index = pd.MultiIndex.from_arrays([[0, 0, 1], [0, 1, 2]])
    s1 = pd.Series(["a", "b", "c"], index=index)
    s2 = pd.Series(["x", "b", "z"], index=index)

    result = s1.compare(s2, align_axis=0)

    indices = pd.MultiIndex.from_arrays(
        [[0, 0, 1, 1], [0, 0, 2, 2], ["self", "other", "self", "other"]]
    )
    expected = pd.Series(["a", "x", "c", "z"], index=indices)
    tm.assert_series_equal(result, expected)


def test_compare_different_indices():
    msg = "Can only compare identically-labeled Series objects"
    ser1 = pd.Series([1, 2, 3], index=["a", "b", "c"])
    ser2 = pd.Series([1, 2, 3], index=["a", "b", "d"])
    with pytest.raises(ValueError, match=msg):
        ser1.compare(ser2)


def test_compare_different_lengths():
    msg = "Can only compare identically-labeled Series objects"
    ser1 = pd.Series([1, 2, 3])
    ser2 = pd.Series([1, 2, 3, 4])
    with pytest.raises(ValueError, match=msg):
        ser1.compare(ser2)


def test_compare_datetime64_and_string():
    # Issue https://github.com/pandas-dev/pandas/issues/45506
    # Catch OverflowError when comparing datetime64 and string
    data = [
        {"a": "2015-07-01", "b": "08335394550"},
        {"a": "2015-07-02", "b": "+49 (0) 0345 300033"},
        {"a": "2015-07-03", "b": "+49(0)2598 04457"},
        {"a": "2015-07-04", "b": "0741470003"},
        {"a": "2015-07-05", "b": "04181 83668"},
    ]
    dtypes = {"a": "datetime64[ns]", "b": "string"}
    df = pd.DataFrame(data=data).astype(dtypes)

    result_eq1 = df["a"].eq(df["b"])
    result_eq2 = df["a"] == df["b"]
    result_neq = df["a"] != df["b"]

    expected_eq = pd.Series([False] * 5)  # For .eq and ==
    expected_neq = pd.Series([True] * 5)  # For !=

    tm.assert_series_equal(result_eq1, expected_eq)
    tm.assert_series_equal(result_eq2, expected_eq)
    tm.assert_series_equal(result_neq, expected_neq)


def test_eq_objects():
    # GH#62191 Test eq with Enum and List elements

    class Thing(Enum):
        FIRST = auto()
        SECOND = auto()

    left = pd.Series([Thing.FIRST, Thing.SECOND])
    py_l = [Thing.FIRST, Thing.SECOND]

    result = left.eq(Thing.FIRST)
    expected = pd.Series([True, False])
    tm.assert_series_equal(result, expected)

    result = left.eq(py_l)
    expected = pd.Series([True, True])
    tm.assert_series_equal(result, expected)

    result = left.eq(np.asarray(py_l))
    expected = pd.Series([True, True])
    tm.assert_series_equal(result, expected)

    result = left.eq(pd.Series(py_l))
    expected = pd.Series([True, True])
    tm.assert_series_equal(result, expected)

    result = pd.Series([[1, 2], [3, 4]]).eq([1, 2])
    expected = pd.Series([True, False])
    with pytest.raises(AssertionError):
        tm.assert_series_equal(result, expected)
    expected = pd.Series([False, False])
    tm.assert_series_equal(result, expected)


def test_eq_with_index():
    # GH#62191 Test eq with non-trivial indices
    left = pd.Series([1, 2], index=[1, 0])
    py_l = [1, 2]

    # assuming Python list has the same index as the Series
    result = left.eq(py_l)
    expected = pd.Series([True, True], index=[1, 0])
    tm.assert_series_equal(result, expected)

    # assuming np.ndarray has the same index as the Series
    result = left.eq(np.asarray(py_l))
    expected = pd.Series([True, True], index=[1, 0])
    tm.assert_series_equal(result, expected)

    result = left.eq(pd.Series(py_l))
    expected = pd.Series([False, False])
    tm.assert_series_equal(result, expected)

    result = left.eq(pd.Series([2, 1]))
    expected = pd.Series([True, True])
    tm.assert_series_equal(result, expected)
