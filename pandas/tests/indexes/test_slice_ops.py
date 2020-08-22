from datetime import date, datetime

import pytest

from pandas import DatetimeIndex, Index, Timestamp, bdate_range


def test_get_slice_bounds_invalid_kind():
    with pytest.raises(ValueError, match="Invalid value for kind kwarg"):
        DatetimeIndex([]).get_slice_bound(datetime(2020, 1, 1), kind="foo", side="left")


def test_get_slice_bounds_invalid_side():
    with pytest.raises(ValueError, match="Invalid value for side kwarg"):
        DatetimeIndex([]).get_slice_bound(
            datetime(2020, 1, 1), kind=None, side="middle"
        )


@pytest.mark.parametrize("kind", ["getitem", "loc", None])
@pytest.mark.parametrize("side, expected", [("left", 4), ("right", 5)])
@pytest.mark.parametrize("data, bound", [(range(6), 4), (list("abcdef"), "e")])
def test_get_slice_bounds_within(kind, side, expected, data, bound):
    index = Index(data)
    result = index.get_slice_bound(bound, kind=kind, side=side)
    assert result == expected


@pytest.mark.parametrize("kind", ["getitem", "loc", None])
@pytest.mark.parametrize("side", ["left", "right"])
@pytest.mark.parametrize(
    "data, bound, expected",
    [
        (range(6), -1, 0),
        (range(6), 10, 6),
        (list("abcdef"), "x", 6),
        (list("bcdefg"), "a", 0),
    ],
)
def test_get_slice_bounds_outside(kind, side, expected, data, bound):
    index = Index(data)
    result = index.get_slice_bound(bound, kind=kind, side=side)
    assert result == expected


@pytest.mark.parametrize("box", [date, datetime, Timestamp])
@pytest.mark.parametrize("kind", ["getitem", "loc", None])
@pytest.mark.parametrize("side, expected", [("left", 4), ("right", 5)])
def test_get_slice_bounds_datetime_within(box, kind, side, expected, tz_aware_fixture):
    # GH 35690
    index = bdate_range("2000-01-03", "2000-02-11").tz_localize(tz_aware_fixture)
    result = index.get_slice_bound(box(year=2000, month=1, day=7), kind=kind, side=side)
    assert result == expected


@pytest.mark.parametrize("box", [date, datetime, Timestamp])
@pytest.mark.parametrize("kind", ["getitem", "loc", None])
@pytest.mark.parametrize("side", ["left", "right"])
@pytest.mark.parametrize("year, expected", [(1999, 0), (2020, 30)])
def test_get_slice_bounds_datetime_outside(
    box, kind, side, year, expected, tz_aware_fixture
):
    # GH 35690
    index = bdate_range("2000-01-03", "2000-02-11").tz_localize(tz_aware_fixture)
    result = index.get_slice_bound(box(year=year, month=1, day=7), kind=kind, side=side)
    assert result == expected


@pytest.mark.parametrize("box", [date, datetime, Timestamp])
@pytest.mark.parametrize("kind", ["getitem", "loc", None])
def test_slice_datetime_locs(box, kind, tz_aware_fixture):
    # GH 34077
    index = DatetimeIndex(["2010-01-01", "2010-01-03"]).tz_localize(tz_aware_fixture)
    result = index.slice_locs(box(2010, 1, 1), box(2010, 1, 2))
    expected = (0, 1)
    assert result == expected
