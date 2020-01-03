import pytest

import pandas as pd
import pandas.util.testing as tm


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
@pytest.mark.parametrize("box_mask", [True, False])
@pytest.mark.parametrize("frame", [True, False])
def test_series_mask_boolean(values, dtype, mask, box_mask, frame):
    ser = pd.Series(values, dtype=dtype, index=["a", "b", "c"])
    if frame:
        ser = ser.to_frame()
    mask = pd.array(mask, dtype="boolean")
    if box_mask:
        mask = pd.Series(mask, index=ser.index)

    expected = ser[mask.astype("bool")]

    result = ser[mask]
    tm.assert_equal(result, expected)

    if not box_mask:
        # Series.iloc[Series[bool]] isn't allowed
        result = ser.iloc[mask]
        tm.assert_equal(result, expected)

    result = ser.loc[mask]
    tm.assert_equal(result, expected)

    # empty
    mask = mask[:0]
    ser = ser.iloc[:0]
    expected = ser[mask.astype("bool")]
    result = ser[mask]
    tm.assert_equal(result, expected)

    if not box_mask:
        # Series.iloc[Series[bool]] isn't allowed
        result = ser.iloc[mask]
        tm.assert_equal(result, expected)

    result = ser.loc[mask]
    tm.assert_equal(result, expected)


@pytest.mark.parametrize("frame", [True, False])
def test_indexing_with_na_raises(frame):
    s = pd.Series([1, 2, 3], name="name")

    if frame:
        s = s.to_frame()
    mask = pd.array([True, False, None], dtype="boolean")
    match = "cannot mask with array containing NA / NaN values"
    with pytest.raises(ValueError, match=match):
        s[mask]

    with pytest.raises(ValueError, match=match):
        s.loc[mask]

    with pytest.raises(ValueError, match=match):
        s.iloc[mask]
