import pytest

from pandas import DatetimeIndex, IndexSlice, Series, Timestamp
from pandas.util.testing import assert_series_equal


def test_indexslice_bad_kwarg_raises():
    with pytest.raises(ValueError, match="invalid option for 'closed': foo"):
        IndexSlice(closed="foo")


@pytest.mark.parametrize(
    "closed, expected_slice",
    [
        (None, slice(0, 2)),  # default
        ("left", slice(0, 1)),
        ("right", slice(1, 2)),
        ("both", slice(0, 2)),
        ("neither", slice(0, 0)),
    ],
)
@pytest.mark.parametrize("left", [Timestamp("2001-01-01 23:50"), "2001-01-01"])
@pytest.mark.parametrize("right", [Timestamp("2001-01-02 00:00"), "2001-01-02"])
@pytest.mark.parametrize(
    "indexer", [(lambda x: x), (lambda x: x.loc)], ids=["getitem", "loc"]
)
def test_series_getitem_closed_kwarg_dates(
    indexer, closed, left, right, expected_slice
):
    # gh-27209
    dates = ["2001-01-01 23:50", "2001-01-02 00:00", "2001-01-03 00:08"]
    ser = Series(range(3), DatetimeIndex(dates))
    expected = ser.iloc[expected_slice]
    idx = IndexSlice(closed=closed)
    result = indexer(ser)[idx[left:right]]
    assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "closed, expected_slice",
    [
        ("left", slice(0, 1)),
        ("right", slice(1, 2)),
        ("both", slice(0, 2)),
        ("neither", slice(0, 0)),
    ],
)
def test_series_getitem_closed_kwarg_int_labels(closed, expected_slice):
    # gh-27209
    int_labels = [50, 70, 80]
    ser = Series(range(3), index=int_labels)
    expected = ser.iloc[expected_slice]
    idx = IndexSlice(closed=closed)
    result = ser.loc[idx[50:70]]
    assert_series_equal(result, expected)
