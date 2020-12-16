import numpy as np
import pytest

from pandas.core.dtypes.cast import maybe_upcast_putmask

from pandas import Series
import pandas._testing as tm


@pytest.mark.parametrize("result", [Series([10, 11, 12]), [10, 11, 12], (10, 11, 12)])
def test_upcast_error(result):
    # GH23823 require result arg to be ndarray
    mask = np.array([False, True, False])
    with pytest.raises(ValueError, match="The result input must be a ndarray"):
        result = maybe_upcast_putmask(result, mask)


def test_upcast():
    # GH23823
    arr = np.arange(1, 6)
    mask = np.array([False, True, False, True, True])
    result = maybe_upcast_putmask(arr, mask)

    expected = np.array([1, np.nan, 3, np.nan, np.nan])
    tm.assert_numpy_array_equal(result, expected)


def test_maybe_upcast_putmask_bool():
    # a case where maybe_upcast_putmask is *not* equivalent to
    #  try: np.putmask(result, mask, np.nan)
    #  except (ValueError, TypeError): result = np.where(mask, result, np.nan)
    arr = np.array([True, False, True, False, True], dtype=bool)
    mask = np.array([False, True, False, True, True])
    result = maybe_upcast_putmask(arr, mask)

    expected = np.array([True, np.nan, True, np.nan, np.nan], dtype=object)
    tm.assert_numpy_array_equal(result, expected)
