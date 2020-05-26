from pandas import Series, read_csv
import pandas as pd
import numpy as np
from pandas.core.arrays.dates import DateArray
from pandas.core.dtypes.dtypes import DateDtype
from pandas.core.dtypes.common import is_date_dtype
import pytest
from pandas.core.arrays.datetimes import DatetimeTZDtype
import pandas.util.testing as tm


def test_from_int_array():
    arr = pd.array(np.arange(5, dtype=np.int64)) * 3600 * 10 ** 9
    result = DateArray._from_sequence(arr)
    np.array_equal(pd.date_range("1970-01-01", periods=5, freq="H").date, result.date)


def test_from_object_array():
    arr = pd.array(np.arange(5, dtype=np.object)) * 3600 * 10 ** 9
    result = DateArray._from_sequence(arr)
    np.array_equal(pd.date_range("1970-01-01", periods=5, freq="H").date, result.date)


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
    arr = pd.array(np_array)
    result = DateArray._from_sequence(arr)
    np.array_equal(pd.array(np_array).date, result.date)


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

def test_df_datetime64_to_date(df: pd.DataFrame, series: Series):
    df["dates"] = series.astype("datetime64").astype("date")
    print(df["dates"])
    print(df)

def test_df_string_to_date():
    df["dates"] = df["strings"].astype("date")

def test_df_object_to_date():
    pass


def test_df_integer_to_date():
    pass


def test_read_data_date():
    df = read_csv()
    df["dates"] = series.astype("datetime64").astype("date")
    df["datetimes"] = series.astype("datetime64")
    print(df["dates"])
    print(df)
    print(df.dtypes)
    for data in df["dates"]:
        print(data)


def test_is_date_dtype_for_dtype():
    date_dtype = DateDtype()
    assert is_date_dtype(date_dtype)

    assert not is_date_dtype("yes64")
    assert not is_date_dtype("date27")

    assert is_date_dtype("date")
    assert is_date_dtype("date64")


def test_iterate_through_date_data():
    pass


if __name__ == "__main__":
    series = Series(["2019-01-01", "2020-12-11", "2020-10-11 12:11:12"])
    series.name = "strings"
    df = series.to_frame()
    df["strings"] = df["strings"].astype("string")
    test_df_datetime64_to_date(df, series)
    # test_datetime64_to_date()
    # test_read_data_date()
    # Series([0, 1, 0, 1]).astype("Int64").apply(str).astype("object").astype("string")
