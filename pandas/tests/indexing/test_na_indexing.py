import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "values, dtype",
    [
        ([1, 2, 3], "int64"),
        ([1.0, 2.0, 3.0], "float64"),
        (["a", "b", "c"], "object"),
        (["a", "b", "c"], "string"),
        ([1, 2, 3], "datetime64[ns]"),
        ([1, 2, 3], "datetime64[ns, CET]"),
        ([1, 2, 3], "timedelta64[ns]"),
        (["2000", "2001", "2002"], "Period[D]"),
        ([1, 0, 3], "Sparse"),
        ([pd.Interval(0, 1), pd.Interval(1, 2), pd.Interval(3, 4)], "interval"),
    ],
)
@pytest.mark.parametrize(
    "mask", [[True, False, False], [True, True, True], [False, False, False]]
)
@pytest.mark.parametrize("indexer_class", [list, pd.array, pd.Index, pd.Series])
@pytest.mark.parametrize("frame", [True, False])
def test_series_mask_boolean(values, dtype, mask, indexer_class, frame):
    ser = pd.Series(values, dtype=dtype, index=["a", "b", "c"])
    if frame:
        ser = ser.to_frame()

    if indexer_class is pd.array:
        mask = pd.array(mask, dtype="boolean")
    elif indexer_class is pd.Series:
        mask = pd.Series(mask, index=ser.index)
    else:
        mask = indexer_class(mask)

    expected = ser[mask]

    result = ser[mask]
    tm.assert_equal(result, expected)

    msg = "iLocation based boolean indexing cannot use an indexable as a mask"
    if indexer_class is pd.Series:
        with pytest.raises(ValueError, match=msg):
            result = ser.iloc[mask]
            tm.assert_equal(result, expected)
    else:
        result = ser.iloc[mask]
        tm.assert_equal(result, expected)

    result = ser.loc[mask]
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "dtype",
    [
        "int64",
        "float64",
        "object",
        "string",
        "datetime64[ns]",
        "datetime64[ns, CET]",
        "timedelta64[ns]",
        "Period[D]",
        "Sparse",
        "interval",
    ],
)
@pytest.mark.parametrize("indexer_class", [list, pd.array, pd.Index, pd.Series])
@pytest.mark.parametrize("frame", [True, False])
def test_series_mask_boolean_empty(dtype, indexer_class, frame):
    if frame:
        ser = pd.DataFrame(dtype=dtype)
    else:
        ser = pd.Series(dtype=dtype)

    mask = []
    if indexer_class is pd.array:
        mask = pd.array(mask, dtype="boolean")
    elif indexer_class is pd.Series:
        mask = pd.Series(mask, index=ser.index, dtype="boolean")
    else:
        mask = indexer_class(mask)

    expected = ser[mask]
    result = ser[mask]
    tm.assert_equal(result, expected)

    msg = "iLocation based boolean indexing cannot use an indexable as a mask"
    if indexer_class is pd.Series:
        with pytest.raises(ValueError, match=msg):
            result = ser.iloc[mask]
            tm.assert_equal(result, expected)
    else:
        result = ser.iloc[mask]
        tm.assert_equal(result, expected)

    result = ser.loc[mask]
    tm.assert_equal(result, expected)


@pytest.mark.parametrize("frame", [True, False])
def test_na_treated_as_false(frame):
    # https://github.com/pandas-dev/pandas/issues/31503
    s = pd.Series([1, 2, 3], name="name")

    if frame:
        s = s.to_frame()

    mask = pd.array([True, False, None], dtype="boolean")

    result = s[mask]
    expected = s[mask.fillna(False)]

    result_loc = s.loc[mask]
    expected_loc = s.loc[mask.fillna(False)]

    result_iloc = s.iloc[mask]
    expected_iloc = s.iloc[mask.fillna(False)]

    if frame:
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result_loc, expected_loc)
        tm.assert_frame_equal(result_iloc, expected_iloc)
    else:
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result_loc, expected_loc)
        tm.assert_series_equal(result_iloc, expected_iloc)
