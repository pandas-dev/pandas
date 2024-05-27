import numpy as np
import pytest

from pandas._libs import lib

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


@pytest.mark.parametrize(
    "atol, rtol, check_exact, expected_self, expected_other",
    [
        (lib.no_default, lib.no_default, True, [1.0, 2.0, 4], [0.4, 1.6, 3.5]),
        (0, 0, False, [1.0, 2.0, 4], [0.4, 1.6, 3.5]),
        (0.5, 0, False, [1.0], [0.4]),
        (0, 0.5, False, [1.0], [0.4]),
        (0.5, 0.00000001, False, [1.0], [0.4]),
        (0.00000001, 0.5, False, [1.0], [0.4]),
        (lib.no_default, lib.no_default, False, [1.0, 2.0, 4], [0.4, 1.6, 3.5]),
        (0.5, lib.no_default, False, [1.0], [0.4]),
        (lib.no_default, 0.5, False, [1.0], [0.4]),
        ("a", lib.no_default, False, None, None),
    ],
)
def test_compare_tolerance_float(
    atol, rtol, check_exact, expected_self, expected_other
):
    df1 = pd.Series([1.0, 2.0, 4])

    df2 = pd.Series([0.4, 1.6, 3.5])

    if expected_self is None:
        with pytest.raises(TypeError):
            df1.compare(df2, atol=atol, rtol=rtol, check_exact=check_exact)
        return

    result = df1.compare(df2, atol=atol, rtol=rtol, check_exact=check_exact)

    expected_data = {"self": expected_self, "other": expected_other}
    expected = pd.DataFrame(expected_data)

    tm.assert_frame_equal(result, expected)


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
