from pandas import Series, read_csv
import pandas as pd
import numpy as np
from pandas.core.arrays.dates import Date64Array, Date64Dtype
from pandas.core.arrays.datetimes import DatetimeTZDtype
import pandas.util.testing as tm

series = Series(["2019-01-01", "2020-12-11", "2020-10-11 12:11:12"])
series.name = "strings"
df = series.to_frame()
df["strings"] = df["strings"].astype("string")


def test_from_int_array():
    arr = pd.array(np.arange(5, dtype=np.int64)) * 3600 * 10 ** 9
    result = Date64Array._from_sequence(arr)
    np.array_equal(pd.date_range("1970-01-01", periods=5, freq="H").date, result.date)


def test_from_object_array():
    arr = pd.array(np.arange(5, dtype=np.object)) * 3600 * 10 ** 9
    result = Date64Array._from_sequence(arr)
    np.array_equal(pd.date_range("1970-01-01", periods=5, freq="H").date, result.date)


def test_from_datetime_array():
    arr = pd.array(
        np.array(
            [
                "2001-01-01T12:00",
                "2002-02-03T13:56:03.172",
                "2007-07-13",
                "2006-01-13",
                "2010-08-13",
            ],
            dtype="datetime64",
        )
    )
    result = Date64Array._from_sequence(arr)
    print(result)


def test_datetime64_to_date():
    df["dates"] = series.astype("datetime64").astype("date")
    df["datetimes"] = series.astype("datetime64")


def test_read_data_date():
    df = read_csv()
    df["dates"] = series.astype("datetime64").astype("date")
    df["datetimes"] = series.astype("datetime64")
    print(df["dates"])
    print(df)
    print(df.dtypes)
    for data in df["dates"]:
        print(data)


def test_iterate_through_date_data():
    pass


def test_strings_objects():
    series.astype("string").astype("date")


if __name__ == "__main__":
    test_from_int_array()
    test_from_object_array()
    test_from_datetime_array()
    # test_datetime64_to_date()
    # test_read_data_date()
    # Series([0, 1, 0, 1]).astype("Int64").apply(str).astype("object").astype("string")
