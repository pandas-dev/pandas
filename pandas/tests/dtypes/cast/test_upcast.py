import numpy as np
import pytest

from pandas.core.dtypes.cast import maybe_upcast_putmask

from pandas import Series
import pandas._testing as tm


@pytest.mark.parametrize("result", [Series([10, 11, 12]), [10, 11, 12], (10, 11, 12)])
def test_upcast_error(result):
    # GH23823 require result arg to be ndarray
    mask = np.array([False, True, False])
    other = np.array([61, 62, 63])
    with pytest.raises(ValueError, match="The result input must be a ndarray"):
        result, _ = maybe_upcast_putmask(result, mask, other)


@pytest.mark.parametrize(
    "arr, other",
    [
        (np.arange(1, 6), np.array([61, 62, 63])),
        (np.arange(1, 6), np.array([61.1, 62.2, 63.3])),
        (np.arange(10, 15), np.array([61, 62])),
        (np.arange(10, 15), np.array([61, np.nan])),
        (
            np.arange("2019-01-01", "2019-01-06", dtype="datetime64[D]"),
            np.arange("2018-01-01", "2018-01-04", dtype="datetime64[D]"),
        ),
        (
            np.arange("2019-01-01", "2019-01-06", dtype="datetime64[D]"),
            np.arange("2018-01-01", "2018-01-03", dtype="datetime64[D]"),
        ),
    ],
)
def test_upcast_scalar_other(arr, other):
    # for now we do not support non-scalar `other`
    mask = np.array([False, True, False, True, True])
    with pytest.raises(ValueError, match="other must be a scalar"):
        maybe_upcast_putmask(arr, mask, other)


def test_upcast():
    # GH23823
    arr = np.arange(1, 6)
    mask = np.array([False, True, False, True, True])
    result, changed = maybe_upcast_putmask(arr, mask, other=np.nan)

    expected = np.array([1, np.nan, 3, np.nan, np.nan])
    assert changed
    tm.assert_numpy_array_equal(result, expected)


def test_upcast_datetime():
    # GH23823
    arr = np.arange("2019-01-01", "2019-01-06", dtype="datetime64[D]")
    mask = np.array([False, True, False, True, True])
    result, changed = maybe_upcast_putmask(arr, mask, other=np.nan)

    expected = np.array(
        [
            "2019-01-01",
            np.datetime64("NaT"),
            "2019-01-03",
            np.datetime64("NaT"),
            np.datetime64("NaT"),
        ],
        dtype="datetime64[D]",
    )
    assert not changed
    tm.assert_numpy_array_equal(result, expected)
