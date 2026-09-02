"""
self-contained to write legacy storage pickle files

To use this script. Create an environment where you want
generate pickles, say its for 0.20.3, with your pandas clone
in ~/pandas

. activate pandas_0.20.3
cd ~/pandas/pandas

$ python -m tests.io.generate_legacy_storage_files \
    tests/io/data/legacy_pickle/0.20.3/ pickle

This script generates a storage file for the current arch, system,
and python version
  pandas version: 0.20.3
  output dir    : pandas/pandas/tests/io/data/legacy_pickle/0.20.3/
  storage format: pickle
created pickle file: 0.20.3_x86_64_darwin_3.5.2.pickle

The idea here is you are using the *current* version of the
generate_legacy_storage_files with an *older* version of pandas to
generate a pickle file. We will then check this file into a current
branch, and test using test_pickle.py. This will load the *older*
pickles and test versus the current data that is generated
(with main). These are then compared.

If we have cases where we changed the signature (e.g. we renamed
offset -> freq in Timestamp). Then we have to conditionally execute
in the generate_legacy_storage_files.py to make it
run under the older AND the newer version.

"""

from datetime import timedelta
import os
import pickle
import platform as pl
import sys

# Remove script directory from path, otherwise Python will try to
# import the JSON test directory as the json module
sys.path.pop(0)

import numpy as np

import pandas as pd
from pandas.arrays import SparseArray

from pandas.tseries.offsets import (
    FY5253,
    BusinessDay,
    BusinessHour,
    CustomBusinessDay,
    DateOffset,
    Day,
    Easter,
    Hour,
    LastWeekOfMonth,
    Minute,
    MonthBegin,
    MonthEnd,
    QuarterBegin,
    QuarterEnd,
    SemiMonthBegin,
    SemiMonthEnd,
    Week,
    WeekOfMonth,
    YearBegin,
    YearEnd,
)


def _create_sp_series():
    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=np.float64)
    arr[7:12] = nan
    arr[-1:] = nan

    bseries = pd.Series(SparseArray(arr, kind="block"))
    bseries.name = "bseries"
    return bseries


def _create_sp_tsseries():
    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=np.float64)
    arr[7:12] = nan
    arr[-1:] = nan

    date_index = pd.bdate_range("1/1/2011", periods=len(arr))
    bseries = pd.Series(SparseArray(arr, kind="block"), index=date_index)
    bseries.name = "btsseries"
    return bseries


def _create_sp_frame():
    nan = np.nan

    data = {
        "A": [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
        "B": [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
        "C": np.arange(10).astype(np.int64),
        "D": [0, 1, 2, 3, 4, 5, nan, nan, nan, nan],
    }

    dates = pd.bdate_range("1/1/2011", periods=10)
    return pd.DataFrame(data, index=dates).apply(SparseArray)


def create_pickle_data(test: bool = True):
    """create the pickle data"""
    data = {
        "A": [0.0, 1.0, 2.0, 3.0, np.nan],
        "B": [0, 1, 0, 1, 0],
        "C": ["foo1", "foo2", "foo3", "foo4", "foo5"],
        "D": pd.date_range("1/1/2009", periods=5),
        "E": [0.0, 1, pd.Timestamp("20100101"), "foo", 2.0],
    }

    scalars = {"timestamp": pd.Timestamp("20130101"), "period": pd.Period("2012", "M")}

    index = {
        "int": pd.Index(np.arange(10)),
        "date": pd.date_range("20130101", periods=10),
        "period": pd.period_range("2013-01-01", freq="M", periods=10),
        "float": pd.Index(np.arange(10, dtype=np.float64)),
        "uint": pd.Index(np.arange(10, dtype=np.uint64)),
        "timedelta": pd.timedelta_range("00:00:00", freq="30min", periods=10),
        "string": pd.Index(["foo", "bar", "baz", "qux", "quux"], dtype="string"),
    }

    index["range"] = pd.RangeIndex(10)

    index["interval"] = pd.interval_range(0, periods=10)

    mi = {
        "reg2": pd.MultiIndex.from_tuples(
            tuple(
                zip(
                    *[
                        ["bar", "bar", "baz", "baz", "foo", "foo", "qux", "qux"],
                        ["one", "two", "one", "two", "one", "two", "one", "two"],
                    ],
                    strict=True,
                )
            ),
            names=["first", "second"],
        )
    }

    series = {
        "float": pd.Series(data["A"]),
        "int": pd.Series(data["B"]),
        "mixed": pd.Series(data["E"]),
        "ts": pd.Series(
            np.arange(10).astype(np.int64), index=pd.date_range("20130101", periods=10)
        ),
        "mi": pd.Series(
            np.arange(5).astype(np.float64),
            index=pd.MultiIndex.from_tuples(
                tuple(zip(*[[1, 1, 2, 2, 2], [3, 4, 3, 4, 5]], strict=True)),
                names=["one", "two"],
            ),
        ),
        "dup": pd.Series(
            np.arange(5).astype(np.float64), index=["A", "B", "C", "D", "A"]
        ),
        "cat": pd.Series(pd.Categorical(["foo", "bar", "baz"])),
        "dt": pd.Series(pd.date_range("20130101", periods=5)),
        "dt_tz": pd.Series(pd.date_range("20130101", periods=5, tz="US/Eastern")),
        "period": pd.Series([pd.Period("2000Q1")] * 5),
        "string": pd.Series(["foo", "bar", "baz", "qux", "quux"], dtype="string"),
    }

    mixed_dup_df = pd.DataFrame(data)
    mixed_dup_df.columns = list("ABCDA")
    frame = {
        "float": pd.DataFrame({"A": series["float"], "B": series["float"] + 1}),
        "int": pd.DataFrame({"A": series["int"], "B": series["int"] + 1}),
        "mixed": pd.DataFrame({k: data[k] for k in ["A", "B", "C", "D"]}),
        "mi": pd.DataFrame(
            {"A": np.arange(5).astype(np.float64), "B": np.arange(5).astype(np.int64)},
            index=pd.MultiIndex.from_tuples(
                tuple(
                    zip(
                        *[
                            ["bar", "bar", "baz", "baz", "baz"],
                            ["one", "two", "one", "two", "three"],
                        ],
                        strict=True,
                    )
                ),
                names=["first", "second"],
            ),
        ),
        "dup": pd.DataFrame(
            np.arange(15).reshape(5, 3).astype(np.float64), columns=["A", "B", "A"]
        ),
        "cat_onecol": pd.DataFrame({"A": pd.Categorical(["foo", "bar"])}),
        "cat_and_float": pd.DataFrame(
            {
                "A": pd.Categorical(["foo", "bar", "baz"]),
                "B": np.arange(3).astype(np.int64),
            }
        ),
        "mixed_dup": mixed_dup_df,
        "dt_mixed_tzs": pd.DataFrame(
            {
                "A": pd.Timestamp("20130102", tz="US/Eastern"),
                "B": pd.Timestamp("20130603", tz="CET"),
            },
            index=range(5),
        ),
        "dt_mixed2_tzs": pd.DataFrame(
            {
                "A": pd.Timestamp("20130102", tz="US/Eastern"),
                "B": pd.Timestamp("20130603", tz="CET"),
                "C": pd.Timestamp("20130603", tz="UTC"),
            },
            index=range(5),
        ),
        "string": pd.DataFrame(
            {
                "A": pd.Series(["foo", "bar", "baz", "qux", "quux"], dtype="string"),
                "B": pd.Series(["one", "two", "one", "two", "three"], dtype="string"),
            }
        ),
    }

    cat = {
        "int8": pd.Categorical(list("abcdefg")),
        "int16": pd.Categorical(np.arange(1000)),
        "int32": pd.Categorical(np.arange(10000)),
    }

    timestamp = {
        "normal": pd.Timestamp("2011-01-01"),
        "nat": pd.NaT,
        "tz": pd.Timestamp("2011-01-01", tz="US/Eastern"),
    }
    if test:
        # kept because those are present in the legacy pickles (<= 1.4)
        timestamp["freq"] = pd.Timestamp("2011-01-01")
        timestamp["both"] = pd.Timestamp("2011-01-01", tz="Asia/Tokyo")

    off = {
        "DateOffset": DateOffset(years=1),
        "DateOffset_h_ns": DateOffset(hour=6, nanoseconds=5824),
        "BusinessDay": BusinessDay(offset=timedelta(seconds=9)),
        "BusinessHour": BusinessHour(normalize=True, n=6, end="15:14"),
        "CustomBusinessDay": CustomBusinessDay(weekmask="Mon Fri"),
        "SemiMonthBegin": SemiMonthBegin(day_of_month=9),
        "SemiMonthEnd": SemiMonthEnd(day_of_month=24),
        "MonthBegin": MonthBegin(1),
        "MonthEnd": MonthEnd(1),
        "QuarterBegin": QuarterBegin(1),
        "QuarterEnd": QuarterEnd(1),
        "Day": Day(1),
        "YearBegin": YearBegin(1),
        "YearEnd": YearEnd(1),
        "Week": Week(1),
        "Week_Tues": Week(2, normalize=False, weekday=1),
        "WeekOfMonth": WeekOfMonth(week=3, weekday=4),
        "LastWeekOfMonth": LastWeekOfMonth(n=1, weekday=3),
        "FY5253": FY5253(n=2, weekday=6, startingMonth=7, variation="last"),
        "Easter": Easter(),
        "Hour": Hour(1),
        "Minute": Minute(1),
    }

    return {
        "series": series,
        "frame": frame,
        "index": index,
        "scalars": scalars,
        "mi": mi,
        "sp_series": {"float": _create_sp_series(), "ts": _create_sp_tsseries()},
        "sp_frame": {"float": _create_sp_frame()},
        "cat": cat,
        "timestamp": timestamp,
        "offsets": off,
    }


def create_dataframe_all_types():
    timestamps = pd.Series(
        [
            pd.Timestamp("2013-01-01"),
            pd.NaT,
            pd.Timestamp("2013-01-03"),
            pd.Timestamp("2013-01-04"),
            pd.Timestamp("2013-01-05"),
        ]
    )
    timedeltas = timestamps - timestamps[0]

    data = {
        # "string": Series(
        #     ["a", "b", "c", None, "e"], dtype=StringDtype(na_value=np.nan)
        # ),
        # "object": Series(["a", "b", "c", None, "e"], dtype=object),
        # "object_nan": Series(["a", "b", "c", np.nan, "e"], dtype=object),
        "int": list(range(1, 6)),
        "uint64": np.arange(3, 8).astype("uint64"),
        "float": [0.1, 0.2, 0.3, 0.4, np.nan],
        "float32": pd.Series([0.1, 0.2, 0.3, 0.4, np.nan], dtype="float32"),
        "bool": [True, False, True, False, True],
        "datetime_ns": timestamps.dt.as_unit("ns"),
        "datetime_us": timestamps.dt.as_unit("us"),
        "datetime_ms": timestamps.dt.as_unit("ms"),
        "datetime_s": timestamps.dt.as_unit("s"),
        "datetimetz_ns": timestamps.dt.tz_localize("US/Eastern").dt.as_unit("ns"),
        "datetimetz_us": timestamps.dt.tz_localize("US/Eastern").dt.as_unit("us"),
        "timedelta_ns": timedeltas.dt.as_unit("ns"),
        "timedelta_us": timedeltas.dt.as_unit("us"),
        "timedelta_ms": timedeltas.dt.as_unit("ms"),
        "timedelta_s": timedeltas.dt.as_unit("s"),
        # "categorical": Categorical(
        #     Series(
        #         ["foo", "bar", "baz",np.nan,"foo"],dtype=StringDtype(na_value=np.nan)
        #     )
        # ),
        # "categorical_object": Categorical(
        #     Series(["foo", "bar", "baz", np.nan, "foo"], dtype=object)
        # ),
        "categorical_int": pd.Categorical([1, 2, 3, np.nan, 1]),
    }
    return pd.DataFrame(data)


def platform_name():
    return "_".join(
        [
            str(pd.__version__),
            str(pl.machine()),
            str(pl.system().lower()),
            str(pl.python_version()),
        ]
    )


def write_legacy_pickles(output_dir):
    pth = f"{platform_name()}.pickle"

    with open(os.path.join(output_dir, pth), "wb") as fh:
        pickle.dump(create_pickle_data(test=False), fh, pickle.DEFAULT_PROTOCOL)

    print(f"created pickle file: {pth}")


def write_legacy_hdf(output_dir, format):
    import tables

    pth = f"{platform_name()}_pytables-{tables.__version__}_{format}.h5"

    df = create_dataframe_all_types()
    if format == "fixed":
        # df = df.drop(columns=["categorical", "categorical_object", "categorical_int"])
        df = df.drop(columns=["categorical_int"])
    complevel = 9 if format == "table" else None
    df.to_hdf(
        os.path.join(output_dir, pth),
        key="df_alltypes",
        format=format,
        complevel=complevel,
    )

    print(f"created hdf file: {pth}")


def write_legacy_file():
    # force our cwd to be the first searched
    sys.path.insert(0, "")

    if not 3 <= len(sys.argv) <= 4:
        sys.exit(
            "Specify output directory and storage type: generate_legacy_"
            "storage_files.py <output_dir> <storage_type> "
        )

    output_dir = str(sys.argv[1])
    storage_type = str(sys.argv[2])

    print(
        "This script generates a storage file for the current arch, system, "
        "and python version"
    )
    print(f"  pandas version: {pd.__version__}")
    print(f"  output dir    : {output_dir}")
    print(f"  storage format: {storage_type}")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if storage_type == "pickle":
        write_legacy_pickles(output_dir=output_dir)
    elif storage_type == "hdf":
        write_legacy_hdf(output_dir=output_dir, format="fixed")
        write_legacy_hdf(output_dir=output_dir, format="table")
    else:
        sys.exit("storage_type must be one of {'pickle', 'hdf'}")


if __name__ == "__main__":
    write_legacy_file()
