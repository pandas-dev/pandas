import os

import pytest

from pandas import compat
import pandas._testing as tm


def test_rands():
    r = tm.rands(10)
    assert len(r) == 10


def test_rands_array_1d():
    arr = tm.rands_array(5, size=10)
    assert arr.shape == (10,)
    assert len(arr[0]) == 5


def test_rands_array_2d():
    arr = tm.rands_array(7, size=(10, 10))
    assert arr.shape == (10, 10)
    assert len(arr[1, 1]) == 7


def test_numpy_err_state_is_default():
    expected = {"over": "warn", "divide": "warn", "invalid": "warn", "under": "ignore"}
    import numpy as np

    # The error state should be unchanged after that import.
    assert np.geterr() == expected


def test_convert_rows_list_to_csv_str():
    rows_list = ["aaa", "bbb", "ccc"]
    ret = tm.convert_rows_list_to_csv_str(rows_list)

    if compat.is_platform_windows():
        expected = "aaa\r\nbbb\r\nccc\r\n"
    else:
        expected = "aaa\nbbb\nccc\n"

    assert ret == expected


@pytest.mark.parametrize("strict_data_files", [True, False])
def test_datapath_missing(datapath):
    with pytest.raises(ValueError, match="Could not find file"):
        datapath("not_a_file")


def test_datapath(datapath):
    args = ("io", "data", "csv", "iris.csv")

    result = datapath(*args)
    expected = os.path.join(os.path.dirname(os.path.dirname(__file__)), *args)

    assert result == expected


def test_external_error_raised():
    with tm.external_error_raised(TypeError):
        raise TypeError("Should not check this error message, so it will pass")
