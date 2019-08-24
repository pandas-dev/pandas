import numpy as np
import pytest

from pandas.core.dtypes.cast import maybe_upcast_putmask

from pandas import Series
from pandas.util import testing as tm


@pytest.mark.parametrize("result", [Series([10, 11, 12]), [10, 11, 12], (10, 11, 12)])
def test_upcast_error(result):
    # GH23823
    mask = np.array([False, True, False])
    other = np.array([61, 62, 63])
    with pytest.raises(ValueError):
        result, _ = maybe_upcast_putmask(result, mask, other)


@pytest.mark.parametrize(
    "arr, other, exp_changed, expected",
    [
        (np.arange(1, 6), np.array([61, 62, 63]), False, np.array([1, 61, 3, 62, 63])),
        (
            np.arange(1, 6),
            np.array([61.1, 62.2, 63.3]),
            True,
            np.array([1, 61.1, 3, 62.2, 63.3]),
        ),
        (np.arange(1, 6), np.nan, True, np.array([1, np.nan, 3, np.nan, np.nan])),
        (np.arange(10, 15), np.array([61, 62]), False, np.array([10, 61, 12, 62, 61])),
        (
            np.arange(10, 15),
            np.array([61, np.nan]),
            True,
            np.array([10, 61, 12, np.nan, 61]),
        ),
    ],
)
def test_upcast(arr, other, exp_changed, expected):
    # GH23823
    mask = np.array([False, True, False, True, True])
    result, changed = maybe_upcast_putmask(arr, mask, other)

    assert changed == exp_changed
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize(
    "arr, other, exp_changed, expected",
    [
        (
            np.arange("2019-01-01", "2019-01-06", dtype="datetime64[D]"),
            np.arange("2018-01-01", "2018-01-04", dtype="datetime64[D]"),
            False,
            np.array(
                ["2019-01-01", "2018-01-01", "2019-01-03", "2018-01-02", "2018-01-03"],
                dtype="datetime64[D]",
            ),
        ),
        (
            np.arange("2019-01-01", "2019-01-06", dtype="datetime64[D]"),
            np.nan,
            False,
            np.array(
                [
                    "2019-01-01",
                    np.datetime64("NaT"),
                    "2019-01-03",
                    np.datetime64("NaT"),
                    np.datetime64("NaT"),
                ],
                dtype="datetime64[D]",
            ),
        ),
        (
            np.arange("2019-01-01", "2019-01-06", dtype="datetime64[D]"),
            np.arange("2018-01-01", "2018-01-03", dtype="datetime64[D]"),
            False,
            np.array(
                ["2019-01-01", "2018-01-01", "2019-01-03", "2018-01-02", "2018-01-01"],
                dtype="datetime64[D]",
            ),
        ),
    ],
)
def test_upcast_datetime(arr, other, exp_changed, expected):
    # GH23823
    mask = np.array([False, True, False, True, True])
    result, changed = maybe_upcast_putmask(arr, mask, other)

    assert changed == exp_changed
    tm.assert_numpy_array_equal(result, expected)
