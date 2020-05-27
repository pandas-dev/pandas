from pandas import Series
import pandas as pd
import numpy as np
from pandas.core.arrays.dates import DateArray
from pandas.core.arrays.datetimes import DatetimeArray
from pandas.core.dtypes.dtypes import DateDtype
from pandas.core.dtypes.common import is_date_dtype
import pytest


def test_init_date_array_from_ints():
    arr = DateArray(np.arange(5))
    np.array_equal(pd.date_range("1970-01-01", periods=5, freq="H").date, arr.date)


def test_init_date_array_from_strings():
    dt_range = pd.date_range("1970-01-01", periods=5, freq="D")
    date_range_as_string = dt_range.astype("string").to_numpy()
    arr = DateArray(date_range_as_string)
    assert np.array_equal(dt_range.date, arr.date)

def test_init_date_array_from_datetimes():
    dt_range = pd.date_range("1970-01-01", periods=5, freq="D")
    date_range_as_string = dt_range.to_numpy()
    arr = DateArray(date_range_as_string)
    assert np.array_equal(dt_range.date, arr.date)


def test_from_int_array():
    arr = pd.array(np.arange(5, dtype=np.int64)) * 3600 * 10 ** 9
    result = DateArray._from_sequence(arr)
    assert np.array_equal(
        pd.date_range("1970-01-01", periods=5, freq="H").date, result.date
    )


def test_from_object_int_array():
    arr = pd.array(np.arange(5, dtype=np.object)) * 3600 * 10 ** 9
    result = DateArray._from_sequence(arr)
    assert np.array_equal(
        pd.date_range("1970-01-01", periods=5, freq="H").date, result.date
    )


def test_from_object_array():
    dt_range = pd.date_range("1970-01-01", periods=5, freq="D")
    obj_array = dt_range.astype("object")
    result = DateArray._from_sequence(obj_array)
    assert np.array_equal(
        dt_range.date, result.date
    )


def test_from_datetime_array():
    np_array = np.array(
        [
            "2001-01-01T12:00",
            "2002-02-03T13:56:03.172",
            "2007-07-13",
            "2006-01-13",
            "2010-08-13",
        ],
        dtype="datetime64",
    )
    pd_array = pd.array(np_array)
    result = DateArray._from_sequence(pd_array)
    assert np.array_equal(pd_array.date, result.date)


def test_to_int_array():
    pass


def test_to_object_array():
    pass


def test_to_string_array():
    pass


def test_to_datetime_array():
    pass


@pytest.fixture
def series():
    series = Series(["2019-01-01", "2020-12-11", "2020-10-11 12:11:12"])
    series.name = "strings"
    return series


@pytest.fixture
def df(series: Series):
    df = series.to_frame()
    df["strings"] = df["strings"].astype("string")
    return df


def test_dtype_name_display(df: pd.DataFrame, series: Series):
    df["dates"] = series.astype("datetime64").astype("date")
    assert df.dtypes[1] == "date"


def test_date_display_format(df: pd.DataFrame, series: Series):
    print(series.astype("datetime64").astype("date"))
    df["dates"] = series.astype("datetime64").astype("date")
    display = str(df["dates"])
    print(display)
    expected = (
        "0   2019-01-01\n"
        "1   2020-12-11\n"
        "2   2020-10-11\n"
        "Name: dates, dtype: date"
    )
    assert display == expected


# def test_read_data_date():
#     df = read_csv()
#     df["dates"] = series.astype("datetime64").astype("date")
#     df["datetimes"] = series.astype("datetime64")
#     print(df["dates"])
#     print(df)
#     print(df.dtypes)
#     for data in df["dates"]:
#         print(data)


def test_is_date_dtype_for_dtype():
    date_dtype = DateDtype()
    assert is_date_dtype(date_dtype)

    assert not is_date_dtype("yes64")
    assert not is_date_dtype("date27")

    assert is_date_dtype("date")
    assert is_date_dtype("date64")


if __name__ == "__main__":
    # series = Series(["2019-01-01", "2020-12-11", "2020-10-11 12:11:12"])
    # series.name = "strings"
    # df = series.to_frame()
    # df["strings"] = df["strings"].astype("string")
    test_from_datetime_array()
    # test_datetime64_to_date()
    # test_read_data_date()
    # Series([0, 1, 0, 1]).astype("Int64").apply(str).astype("object").astype("string")
