import builtins
import datetime as dt
from string import ascii_lowercase

import numpy as np
import pytest

from pandas._libs.tslibs import iNaT

from pandas.core.dtypes.common import pandas_dtype
from pandas.core.dtypes.missing import na_value_for_dtype

import pandas as pd
from pandas import (
    DataFrame,
    MultiIndex,
    Series,
    Timestamp,
    date_range,
    isna,
)
import pandas._testing as tm
from pandas.tests.groupby import get_groupby_method_args
from pandas.util import _test_decorators as td


@pytest.mark.parametrize("dtype", ["int64", "int32", "float64", "float32"])
def test_basic_aggregations(dtype):
    data = Series(np.arange(9) // 3, index=np.arange(9), dtype=dtype)

    index = np.arange(9)
    np.random.default_rng(2).shuffle(index)
    data = data.reindex(index)

    grouped = data.groupby(lambda x: x // 3, group_keys=False)

    for k, v in grouped:
        assert len(v) == 3

    agged = grouped.aggregate(np.mean)
    assert agged[1] == 1

    expected = grouped.agg(np.mean)
    tm.assert_series_equal(agged, expected)  # shorthand
    tm.assert_series_equal(agged, grouped.mean())
    result = grouped.sum()
    expected = grouped.agg(np.sum)
    if dtype == "int32":
        # NumPy's sum returns int64
        expected = expected.astype("int32")
    tm.assert_series_equal(result, expected)

    expected = grouped.apply(lambda x: x * x.sum())
    transformed = grouped.transform(lambda x: x * x.sum())
    assert transformed[7] == 12
    tm.assert_series_equal(transformed, expected)

    value_grouped = data.groupby(data)
    result = value_grouped.aggregate(np.mean)
    tm.assert_series_equal(result, agged, check_index_type=False)

    # complex agg
    agged = grouped.aggregate([np.mean, np.std])

    msg = r"nested renamer is not supported"
    with pytest.raises(pd.errors.SpecificationError, match=msg):
        grouped.aggregate({"one": np.mean, "two": np.std})

    # corner cases
    msg = "Must produce aggregated value"
    # exception raised is type Exception
    with pytest.raises(Exception, match=msg):
        grouped.aggregate(lambda x: x * 2)


@pytest.mark.parametrize(
    "vals",
    [
        ["foo", "bar", "baz"],
        ["foo", "", ""],
        ["", "", ""],
        [1, 2, 3],
        [1, 0, 0],
        [0, 0, 0],
        [1.0, 2.0, 3.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [True, True, True],
        [True, False, False],
        [False, False, False],
        [np.nan, np.nan, np.nan],
    ],
)
def test_groupby_bool_aggs(skipna, all_boolean_reductions, vals):
    df = DataFrame({"key": ["a"] * 3 + ["b"] * 3, "val": vals * 2})

    # Figure out expectation using Python builtin
    exp = getattr(builtins, all_boolean_reductions)(vals)

    # edge case for missing data with skipna and 'any'
    if skipna and all(isna(vals)) and all_boolean_reductions == "any":
        exp = False

    expected = DataFrame(
        [exp] * 2, columns=["val"], index=pd.Index(["a", "b"], name="key")
    )
    result = getattr(df.groupby("key"), all_boolean_reductions)(skipna=skipna)
    tm.assert_frame_equal(result, expected)


def test_any():
    df = DataFrame(
        [[1, 2, "foo"], [1, np.nan, "bar"], [3, np.nan, "baz"]],
        columns=["A", "B", "C"],
    )
    expected = DataFrame(
        [[True, True], [False, True]], columns=["B", "C"], index=[1, 3]
    )
    expected.index.name = "A"
    result = df.groupby("A").any()
    tm.assert_frame_equal(result, expected)


def test_bool_aggs_dup_column_labels(all_boolean_reductions):
    # GH#21668
    df = DataFrame([[True, True]], columns=["a", "a"])
    grp_by = df.groupby([0])
    result = getattr(grp_by, all_boolean_reductions)()

    expected = df.set_axis(np.array([0]))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data",
    [
        [False, False, False],
        [True, True, True],
        [pd.NA, pd.NA, pd.NA],
        [False, pd.NA, False],
        [True, pd.NA, True],
        [True, pd.NA, False],
    ],
)
def test_masked_kleene_logic(all_boolean_reductions, skipna, data):
    # GH#37506
    ser = Series(data, dtype="boolean")

    # The result should match aggregating on the whole series. Correctness
    # there is verified in test_reductions.py::test_any_all_boolean_kleene_logic
    expected_data = getattr(ser, all_boolean_reductions)(skipna=skipna)
    expected = Series(expected_data, index=np.array([0]), dtype="boolean")

    result = ser.groupby([0, 0, 0]).agg(all_boolean_reductions, skipna=skipna)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "dtype1,dtype2,exp_col1,exp_col2",
    [
        (
            "float",
            "Float64",
            np.array([True], dtype=bool),
            pd.array([pd.NA], dtype="boolean"),
        ),
        (
            "Int64",
            "float",
            pd.array([pd.NA], dtype="boolean"),
            np.array([True], dtype=bool),
        ),
        (
            "Int64",
            "Int64",
            pd.array([pd.NA], dtype="boolean"),
            pd.array([pd.NA], dtype="boolean"),
        ),
        (
            "Float64",
            "boolean",
            pd.array([pd.NA], dtype="boolean"),
            pd.array([pd.NA], dtype="boolean"),
        ),
    ],
)
def test_masked_mixed_types(dtype1, dtype2, exp_col1, exp_col2):
    # GH#37506
    data = [1.0, np.nan]
    df = DataFrame(
        {"col1": pd.array(data, dtype=dtype1), "col2": pd.array(data, dtype=dtype2)}
    )
    result = df.groupby([1, 1]).agg("all", skipna=False)

    expected = DataFrame({"col1": exp_col1, "col2": exp_col2}, index=np.array([1]))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("dtype", ["Int64", "Float64", "boolean"])
def test_masked_bool_aggs_skipna(
    all_boolean_reductions, dtype, skipna, frame_or_series
):
    # GH#40585
    obj = frame_or_series([pd.NA, 1], dtype=dtype)
    expected_res = True
    if not skipna and all_boolean_reductions == "all":
        expected_res = pd.NA
    expected = frame_or_series([expected_res], index=np.array([1]), dtype="boolean")

    result = obj.groupby([1, 1]).agg(all_boolean_reductions, skipna=skipna)
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "bool_agg_func,data,expected_res",
    [
        ("any", [pd.NA, np.nan], False),
        ("any", [pd.NA, 1, np.nan], True),
        ("all", [pd.NA, pd.NaT], True),
        ("all", [pd.NA, False, pd.NaT], False),
    ],
)
def test_object_type_missing_vals(bool_agg_func, data, expected_res, frame_or_series):
    # GH#37501
    obj = frame_or_series(data, dtype=object)
    result = obj.groupby([1] * len(data)).agg(bool_agg_func)
    expected = frame_or_series([expected_res], index=np.array([1]), dtype="bool")
    tm.assert_equal(result, expected)


def test_object_NA_raises_with_skipna_false(all_boolean_reductions):
    # GH#37501
    ser = Series([pd.NA], dtype=object)
    with pytest.raises(TypeError, match="boolean value of NA is ambiguous"):
        ser.groupby([1]).agg(all_boolean_reductions, skipna=False)


def test_empty(frame_or_series, all_boolean_reductions):
    # GH 45231
    kwargs = {"columns": ["a"]} if frame_or_series is DataFrame else {"name": "a"}
    obj = frame_or_series(**kwargs, dtype=object)
    result = getattr(obj.groupby(obj.index), all_boolean_reductions)()
    expected = frame_or_series(**kwargs, dtype=bool)
    tm.assert_equal(result, expected)


@pytest.mark.parametrize("how", ["idxmin", "idxmax"])
def test_idxmin_idxmax_extremes(how, any_real_numpy_dtype):
    # GH#57040
    if any_real_numpy_dtype is int or any_real_numpy_dtype is float:
        # No need to test
        return
    info = np.iinfo if "int" in any_real_numpy_dtype else np.finfo
    min_value = info(any_real_numpy_dtype).min
    max_value = info(any_real_numpy_dtype).max
    df = DataFrame(
        {"a": [2, 1, 1, 2], "b": [min_value, max_value, max_value, min_value]},
        dtype=any_real_numpy_dtype,
    )
    gb = df.groupby("a")
    result = getattr(gb, how)()
    expected = DataFrame(
        {"b": [1, 0]}, index=pd.Index([1, 2], name="a", dtype=any_real_numpy_dtype)
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("how", ["idxmin", "idxmax"])
def test_idxmin_idxmax_extremes_skipna(skipna, how, float_numpy_dtype):
    # GH#57040
    min_value = np.finfo(float_numpy_dtype).min
    max_value = np.finfo(float_numpy_dtype).max
    df = DataFrame(
        {
            "a": Series(np.repeat(range(1, 5), repeats=2), dtype="intp"),
            "b": Series(
                [
                    np.nan,
                    min_value,
                    np.nan,
                    max_value,
                    min_value,
                    np.nan,
                    max_value,
                    np.nan,
                ],
                dtype=float_numpy_dtype,
            ),
        },
    )
    gb = df.groupby("a")

    if not skipna:
        msg = f"DataFrameGroupBy.{how} with skipna=False"
        with pytest.raises(ValueError, match=msg):
            getattr(gb, how)(skipna=skipna)
        return
    result = getattr(gb, how)(skipna=skipna)
    expected = DataFrame(
        {"b": [1, 3, 4, 6]}, index=pd.Index(range(1, 5), name="a", dtype="intp")
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "func, values",
    [
        ("idxmin", {"c_int": [0, 2], "c_float": [1, 3], "c_date": [1, 2]}),
        ("idxmax", {"c_int": [1, 3], "c_float": [0, 2], "c_date": [0, 3]}),
    ],
)
@pytest.mark.parametrize("numeric_only", [True, False])
def test_idxmin_idxmax_returns_int_types(func, values, numeric_only):
    # GH 25444
    df = DataFrame(
        {
            "name": ["A", "A", "B", "B"],
            "c_int": [1, 2, 3, 4],
            "c_float": [4.02, 3.03, 2.04, 1.05],
            "c_date": ["2019", "2018", "2016", "2017"],
        }
    )
    df["c_date"] = pd.to_datetime(df["c_date"])
    df["c_date_tz"] = df["c_date"].dt.tz_localize("US/Pacific")
    df["c_timedelta"] = df["c_date"] - df["c_date"].iloc[0]
    df["c_period"] = df["c_date"].dt.to_period("W")
    df["c_Integer"] = df["c_int"].astype("Int64")
    df["c_Floating"] = df["c_float"].astype("Float64")

    result = getattr(df.groupby("name"), func)(numeric_only=numeric_only)

    expected = DataFrame(values, index=pd.Index(["A", "B"], name="name"))
    if numeric_only:
        expected = expected.drop(columns=["c_date"])
    else:
        expected["c_date_tz"] = expected["c_date"]
        expected["c_timedelta"] = expected["c_date"]
        expected["c_period"] = expected["c_date"]
    expected["c_Integer"] = expected["c_int"]
    expected["c_Floating"] = expected["c_float"]

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "data",
    [
        (
            Timestamp("2011-01-15 12:50:28.502376"),
            Timestamp("2011-01-20 12:50:28.593448"),
        ),
        (24650000000000001, 24650000000000002),
    ],
)
@pytest.mark.parametrize("method", ["count", "min", "max", "first", "last"])
def test_groupby_non_arithmetic_agg_int_like_precision(method, data):
    # GH#6620, GH#9311
    df = DataFrame({"a": [1, 1], "b": data})

    grouped = df.groupby("a")
    result = getattr(grouped, method)()
    if method == "count":
        expected_value = 2
    elif method == "first":
        expected_value = data[0]
    elif method == "last":
        expected_value = data[1]
    else:
        expected_value = getattr(df["b"], method)()
    expected = DataFrame({"b": [expected_value]}, index=pd.Index([1], name="a"))

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("how", ["first", "last"])
def test_first_last_skipna(any_real_nullable_dtype, sort, skipna, how):
    # GH#57019
    na_value = na_value_for_dtype(pandas_dtype(any_real_nullable_dtype))
    df = DataFrame(
        {
            "a": [2, 1, 1, 2, 3, 3],
            "b": [na_value, 3.0, na_value, 4.0, np.nan, np.nan],
            "c": [na_value, 3.0, na_value, 4.0, np.nan, np.nan],
        },
        dtype=any_real_nullable_dtype,
    )
    gb = df.groupby("a", sort=sort)
    method = getattr(gb, how)
    result = method(skipna=skipna)

    ilocs = {
        ("first", True): [3, 1, 4],
        ("first", False): [0, 1, 4],
        ("last", True): [3, 1, 5],
        ("last", False): [3, 2, 5],
    }[how, skipna]
    expected = df.iloc[ilocs].set_index("a")
    if sort:
        expected = expected.sort_index()
    tm.assert_frame_equal(result, expected)


def test_groupby_mean_no_overflow():
    # Regression test for (#22487)
    df = DataFrame(
        {
            "user": ["A", "A", "A", "A", "A"],
            "connections": [4970, 4749, 4719, 4704, 18446744073699999744],
        }
    )
    assert df.groupby("user")["connections"].mean()["A"] == 3689348814740003840


def test_mean_on_timedelta():
    # GH 17382
    df = DataFrame({"time": pd.to_timedelta(range(10)), "cat": ["A", "B"] * 5})
    result = df.groupby("cat")["time"].mean()
    expected = Series(
        pd.to_timedelta([4, 5]), name="time", index=pd.Index(["A", "B"], name="cat")
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "values, dtype, result_dtype",
    [
        ([0, 1, np.nan, 3, 4, 5, 6, 7, 8, 9], "float64", "float64"),
        ([0, 1, pd.NA, 3, 4, 5, 6, 7, 8, 9], "Float64", "Float64"),
        ([0, 1, pd.NA, 3, 4, 5, 6, 7, 8, 9], "Int64", "Float64"),
        ([0, 1, np.nan, 3, 4, 5, 6, 7, 8, 9], "timedelta64[ns]", "timedelta64[ns]"),
        (
            pd.to_datetime(
                [
                    "2019-05-09",
                    pd.NaT,
                    "2019-05-11",
                    "2019-05-12",
                    "2019-05-13",
                    "2019-05-14",
                    "2019-05-15",
                    "2019-05-16",
                    "2019-05-17",
                    "2019-05-18",
                ]
            ),
            "datetime64[ns]",
            "datetime64[ns]",
        ),
    ],
)
def test_mean_skipna(values, dtype, result_dtype, skipna):
    # GH#15675
    df = DataFrame(
        {
            "val": values,
            "cat": ["A", "B"] * 5,
        }
    ).astype({"val": dtype})
    # We need to recast the expected values to the result_dtype because
    # Series.mean() changes the dtype to float64/object depending on the input dtype
    expected = (
        df.groupby("cat")["val"]
        .apply(lambda x: x.mean(skipna=skipna))
        .astype(result_dtype)
    )
    result = df.groupby("cat")["val"].mean(skipna=skipna)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "values, dtype",
    [
        ([0, 1, np.nan, 3, 4, 5, 6, 7, 8, 9], "float64"),
        ([0, 1, pd.NA, 3, 4, 5, 6, 7, 8, 9], "Float64"),
        ([0, 1, pd.NA, 3, 4, 5, 6, 7, 8, 9], "Int64"),
        ([0, 1, np.nan, 3, 4, 5, 6, 7, 8, 9], "timedelta64[ns]"),
    ],
)
def test_sum_skipna(values, dtype, skipna):
    # GH#15675
    df = DataFrame(
        {
            "val": values,
            "cat": ["A", "B"] * 5,
        }
    ).astype({"val": dtype})
    # We need to recast the expected values to the original dtype because
    # Series.sum() changes the dtype
    expected = (
        df.groupby("cat")["val"].apply(lambda x: x.sum(skipna=skipna)).astype(dtype)
    )
    result = df.groupby("cat")["val"].sum(skipna=skipna)
    tm.assert_series_equal(result, expected)


def test_sum_skipna_object(skipna):
    # GH#15675
    df = DataFrame(
        {
            "val": ["a", "b", np.nan, "d", "e", "f", "g", "h", "i", "j"],
            "cat": ["A", "B"] * 5,
        }
    ).astype({"val": object})
    if skipna:
        expected = Series(
            ["aegi", "bdfhj"], index=pd.Index(["A", "B"], name="cat"), name="val"
        ).astype(object)
    else:
        expected = Series(
            [np.nan, "bdfhj"], index=pd.Index(["A", "B"], name="cat"), name="val"
        ).astype(object)
    result = df.groupby("cat")["val"].sum(skipna=skipna)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "func, values, dtype, result_dtype",
    [
        ("prod", [0, 1, 3, np.nan, 4, 5, 6, 7, -8, 9], "float64", "float64"),
        ("prod", [0, -1, 3, 4, 5, pd.NA, 6, 7, 8, 9], "Float64", "Float64"),
        ("prod", [0, 1, 3, -4, 5, 6, 7, -8, pd.NA, 9], "Int64", "Int64"),
        ("prod", [np.nan] * 10, "float64", "float64"),
        ("prod", [pd.NA] * 10, "Float64", "Float64"),
        ("prod", [pd.NA] * 10, "Int64", "Int64"),
        ("var", [0, -1, 3, 4, np.nan, 5, 6, 7, 8, 9], "float64", "float64"),
        ("var", [0, 1, 3, -4, 5, 6, 7, -8, 9, pd.NA], "Float64", "Float64"),
        ("var", [0, -1, 3, 4, 5, -6, 7, pd.NA, 8, 9], "Int64", "Float64"),
        ("var", [np.nan] * 10, "float64", "float64"),
        ("var", [pd.NA] * 10, "Float64", "Float64"),
        ("var", [pd.NA] * 10, "Int64", "Float64"),
        ("std", [0, 1, 3, -4, 5, 6, 7, -8, np.nan, 9], "float64", "float64"),
        ("std", [0, -1, 3, 4, 5, -6, 7, pd.NA, 8, 9], "Float64", "Float64"),
        ("std", [0, 1, 3, -4, 5, 6, 7, -8, 9, pd.NA], "Int64", "Float64"),
        ("std", [np.nan] * 10, "float64", "float64"),
        ("std", [pd.NA] * 10, "Float64", "Float64"),
        ("std", [pd.NA] * 10, "Int64", "Float64"),
        ("sem", [0, -1, 3, 4, 5, -6, 7, np.nan, 8, 9], "float64", "float64"),
        ("sem", [0, 1, 3, -4, 5, 6, 7, -8, pd.NA, 9], "Float64", "Float64"),
        ("sem", [0, -1, 3, 4, 5, -6, 7, 8, 9, pd.NA], "Int64", "Float64"),
        ("sem", [np.nan] * 10, "float64", "float64"),
        ("sem", [pd.NA] * 10, "Float64", "Float64"),
        ("sem", [pd.NA] * 10, "Int64", "Float64"),
        ("min", [0, -1, 3, 4, 5, -6, 7, np.nan, 8, 9], "float64", "float64"),
        ("min", [0, 1, 3, -4, 5, 6, 7, -8, pd.NA, 9], "Float64", "Float64"),
        ("min", [0, -1, 3, 4, 5, -6, 7, 8, 9, pd.NA], "Int64", "Int64"),
        (
            "min",
            [0, 1, np.nan, 3, 4, 5, 6, 7, 8, 9],
            "timedelta64[ns]",
            "timedelta64[ns]",
        ),
        (
            "min",
            pd.to_datetime(
                [
                    "2019-05-09",
                    pd.NaT,
                    "2019-05-11",
                    "2019-05-12",
                    "2019-05-13",
                    "2019-05-14",
                    "2019-05-15",
                    "2019-05-16",
                    "2019-05-17",
                    "2019-05-18",
                ]
            ),
            "datetime64[ns]",
            "datetime64[ns]",
        ),
        ("min", [np.nan] * 10, "float64", "float64"),
        ("min", [pd.NA] * 10, "Float64", "Float64"),
        ("min", [pd.NA] * 10, "Int64", "Int64"),
        ("max", [0, -1, 3, 4, 5, -6, 7, np.nan, 8, 9], "float64", "float64"),
        ("max", [0, 1, 3, -4, 5, 6, 7, -8, pd.NA, 9], "Float64", "Float64"),
        ("max", [0, -1, 3, 4, 5, -6, 7, 8, 9, pd.NA], "Int64", "Int64"),
        (
            "max",
            [0, 1, np.nan, 3, 4, 5, 6, 7, 8, 9],
            "timedelta64[ns]",
            "timedelta64[ns]",
        ),
        (
            "max",
            pd.to_datetime(
                [
                    "2019-05-09",
                    pd.NaT,
                    "2019-05-11",
                    "2019-05-12",
                    "2019-05-13",
                    "2019-05-14",
                    "2019-05-15",
                    "2019-05-16",
                    "2019-05-17",
                    "2019-05-18",
                ]
            ),
            "datetime64[ns]",
            "datetime64[ns]",
        ),
        ("max", [np.nan] * 10, "float64", "float64"),
        ("max", [pd.NA] * 10, "Float64", "Float64"),
        ("max", [pd.NA] * 10, "Int64", "Int64"),
        ("median", [0, -1, 3, 4, 5, -6, 7, np.nan, 8, 9], "float64", "float64"),
        ("median", [0, 1, 3, -4, 5, 6, 7, -8, pd.NA, 9], "Float64", "Float64"),
        ("median", [0, -1, 3, 4, 5, -6, 7, 8, 9, pd.NA], "Int64", "Float64"),
        (
            "median",
            [0, 1, np.nan, 3, 4, 5, 6, 7, 8, 9],
            "timedelta64[ns]",
            "timedelta64[ns]",
        ),
        (
            "median",
            pd.to_datetime(
                [
                    "2019-05-09",
                    pd.NaT,
                    "2019-05-11",
                    "2019-05-12",
                    "2019-05-13",
                    "2019-05-14",
                    "2019-05-15",
                    "2019-05-16",
                    "2019-05-17",
                    "2019-05-18",
                ]
            ),
            "datetime64[ns]",
            "datetime64[ns]",
        ),
        ("median", [np.nan] * 10, "float64", "float64"),
        ("median", [pd.NA] * 10, "Float64", "Float64"),
        ("median", [pd.NA] * 10, "Int64", "Float64"),
    ],
)
def test_multifunc_skipna(func, values, dtype, result_dtype, skipna):
    # GH#15675
    df = DataFrame(
        {
            "val": values,
            "cat": ["A", "B"] * 5,
        }
    ).astype({"val": dtype})
    # We need to recast the expected values to the result_dtype as some operations
    # change the dtype
    expected = (
        df.groupby("cat")["val"]
        .apply(lambda x: getattr(x, func)(skipna=skipna))
        .astype(result_dtype)
    )
    result = getattr(df.groupby("cat")["val"], func)(skipna=skipna)
    tm.assert_series_equal(result, expected)


def test_cython_median():
    arr = np.random.default_rng(2).standard_normal(1000)
    arr[::2] = np.nan
    df = DataFrame(arr)

    labels = np.random.default_rng(2).integers(0, 50, size=1000).astype(float)
    labels[::17] = np.nan

    result = df.groupby(labels).median()
    exp = df.groupby(labels).agg(np.nanmedian)
    tm.assert_frame_equal(result, exp)

    df = DataFrame(np.random.default_rng(2).standard_normal((1000, 5)))
    rs = df.groupby(labels).agg(np.median)
    xp = df.groupby(labels).median()
    tm.assert_frame_equal(rs, xp)


def test_median_empty_bins(observed):
    df = DataFrame(np.random.default_rng(2).integers(0, 44, 500))

    grps = range(0, 55, 5)
    bins = pd.cut(df[0], grps)

    result = df.groupby(bins, observed=observed).median()
    expected = df.groupby(bins, observed=observed).agg(lambda x: x.median())
    tm.assert_frame_equal(result, expected)


def test_max_min_non_numeric():
    # #2700
    aa = DataFrame({"nn": [11, 11, 22, 22], "ii": [1, 2, 3, 4], "ss": 4 * ["mama"]})

    result = aa.groupby("nn").max()
    assert "ss" in result

    result = aa.groupby("nn").max(numeric_only=False)
    assert "ss" in result

    result = aa.groupby("nn").min()
    assert "ss" in result

    result = aa.groupby("nn").min(numeric_only=False)
    assert "ss" in result


def test_max_min_object_multiple_columns(using_infer_string):
    # GH#41111 case where the aggregation is valid for some columns but not
    # others; we split object blocks column-wise, consistent with
    # DataFrame._reduce

    df = DataFrame(
        {
            "A": [1, 1, 2, 2, 3],
            "B": [1, "foo", 2, "bar", False],
            "C": ["a", "b", "c", "d", "e"],
        }
    )
    df._consolidate_inplace()  # should already be consolidate, but double-check
    assert len(df._mgr.blocks) == 3 if using_infer_string else 2

    gb = df.groupby("A")

    result = gb[["C"]].max()
    # "max" is valid for column "C" but not for "B"
    ei = pd.Index([1, 2, 3], name="A")
    expected = DataFrame({"C": ["b", "d", "e"]}, index=ei)
    tm.assert_frame_equal(result, expected)

    result = gb[["C"]].min()
    # "min" is valid for column "C" but not for "B"
    ei = pd.Index([1, 2, 3], name="A")
    expected = DataFrame({"C": ["a", "c", "e"]}, index=ei)
    tm.assert_frame_equal(result, expected)


def test_min_date_with_nans():
    # GH26321
    dates = pd.to_datetime(
        Series(["2019-05-09", "2019-05-09", "2019-05-09"]), format="%Y-%m-%d"
    ).dt.date
    df = DataFrame({"a": [np.nan, "1", np.nan], "b": [0, 1, 1], "c": dates})

    result = df.groupby("b", as_index=False)["c"].min()["c"]
    expected = pd.to_datetime(
        Series(["2019-05-09", "2019-05-09"], name="c"), format="%Y-%m-%d"
    ).dt.date
    tm.assert_series_equal(result, expected)

    result = df.groupby("b")["c"].min()
    expected.index.name = "b"
    tm.assert_series_equal(result, expected)


def test_max_inat():
    # GH#40767 dont interpret iNaT as NaN
    ser = Series([1, iNaT])
    key = np.array([1, 1], dtype=np.int64)
    gb = ser.groupby(key)

    result = gb.max(min_count=2)
    expected = Series({1: 1}, dtype=np.int64)
    tm.assert_series_equal(result, expected, check_exact=True)

    result = gb.min(min_count=2)
    expected = Series({1: iNaT}, dtype=np.int64)
    tm.assert_series_equal(result, expected, check_exact=True)

    # not enough entries -> gets masked to NaN
    result = gb.min(min_count=3)
    expected = Series({1: np.nan})
    tm.assert_series_equal(result, expected, check_exact=True)


def test_max_inat_not_all_na():
    # GH#40767 dont interpret iNaT as NaN

    # make sure we dont round iNaT+1 to iNaT
    ser = Series([1, iNaT, 2, iNaT + 1])
    gb = ser.groupby([1, 2, 3, 3])
    result = gb.min(min_count=2)

    # Note: in converting to float64, the iNaT + 1 maps to iNaT, i.e. is lossy
    expected = Series({1: np.nan, 2: np.nan, 3: iNaT + 1})
    expected.index = expected.index.astype(int)
    tm.assert_series_equal(result, expected, check_exact=True)


@pytest.mark.parametrize("func", ["min", "max"])
def test_groupby_aggregate_period_column(func):
    # GH 31471
    groups = [1, 2]
    periods = pd.period_range("2020", periods=2, freq="Y")
    df = DataFrame({"a": groups, "b": periods})

    result = getattr(df.groupby("a")["b"], func)()
    idx = pd.Index([1, 2], name="a")
    expected = Series(periods, index=idx, name="b")

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("func", ["min", "max"])
def test_groupby_aggregate_period_frame(func):
    # GH 31471
    groups = [1, 2]
    periods = pd.period_range("2020", periods=2, freq="Y")
    df = DataFrame({"a": groups, "b": periods})

    result = getattr(df.groupby("a"), func)()
    idx = pd.Index([1, 2], name="a")
    expected = DataFrame({"b": periods}, index=idx)

    tm.assert_frame_equal(result, expected)


def test_aggregate_numeric_object_dtype():
    # https://github.com/pandas-dev/pandas/issues/39329
    # simplified case: multiple object columns where one is all-NaN
    # -> gets split as the all-NaN is inferred as float
    df = DataFrame(
        {"key": ["A", "A", "B", "B"], "col1": list("abcd"), "col2": [np.nan] * 4},
    ).astype(object)
    result = df.groupby("key").min()
    expected = (
        DataFrame(
            {"key": ["A", "B"], "col1": ["a", "c"], "col2": [np.nan, np.nan]},
        )
        .set_index("key")
        .astype(object)
    )
    tm.assert_frame_equal(result, expected)

    # same but with numbers
    df = DataFrame(
        {"key": ["A", "A", "B", "B"], "col1": list("abcd"), "col2": range(4)},
    ).astype(object)
    result = df.groupby("key").min()
    expected = (
        DataFrame({"key": ["A", "B"], "col1": ["a", "c"], "col2": [0, 2]})
        .set_index("key")
        .astype(object)
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("func", ["min", "max"])
def test_aggregate_categorical_lost_index(func: str):
    # GH: 28641 groupby drops index, when grouping over categorical column with min/max
    ds = Series(["b"], dtype="category").cat.as_ordered()
    df = DataFrame({"A": [1997], "B": ds})
    result = df.groupby("A").agg({"B": func})
    expected = DataFrame({"B": ["b"]}, index=pd.Index([1997], name="A"))

    # ordered categorical dtype should be preserved
    expected["B"] = expected["B"].astype(ds.dtype)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("dtype", ["Int64", "Int32", "Float64", "Float32", "boolean"])
def test_groupby_min_max_nullable(dtype):
    if dtype == "Int64":
        # GH#41743 avoid precision loss
        ts = 1618556707013635762
    elif dtype == "boolean":
        ts = 0
    else:
        ts = 4.0

    df = DataFrame({"id": [2, 2], "ts": [ts, ts + 1]})
    df["ts"] = df["ts"].astype(dtype)

    gb = df.groupby("id")

    result = gb.min()
    expected = df.iloc[:1].set_index("id")
    tm.assert_frame_equal(result, expected)

    res_max = gb.max()
    expected_max = df.iloc[1:].set_index("id")
    tm.assert_frame_equal(res_max, expected_max)

    result2 = gb.min(min_count=3)
    expected2 = DataFrame({"ts": [pd.NA]}, index=expected.index, dtype=dtype)
    tm.assert_frame_equal(result2, expected2)

    res_max2 = gb.max(min_count=3)
    tm.assert_frame_equal(res_max2, expected2)

    # Case with NA values
    df2 = DataFrame({"id": [2, 2, 2], "ts": [ts, pd.NA, ts + 1]})
    df2["ts"] = df2["ts"].astype(dtype)
    gb2 = df2.groupby("id")

    result3 = gb2.min()
    tm.assert_frame_equal(result3, expected)

    res_max3 = gb2.max()
    tm.assert_frame_equal(res_max3, expected_max)

    result4 = gb2.min(min_count=100)
    tm.assert_frame_equal(result4, expected2)

    res_max4 = gb2.max(min_count=100)
    tm.assert_frame_equal(res_max4, expected2)


def test_min_max_nullable_uint64_empty_group():
    # don't raise NotImplementedError from libgroupby
    cat = pd.Categorical([0] * 10, categories=[0, 1])
    df = DataFrame({"A": cat, "B": pd.array(np.arange(10, dtype=np.uint64))})
    gb = df.groupby("A", observed=False)

    res = gb.min()

    idx = pd.CategoricalIndex([0, 1], dtype=cat.dtype, name="A")
    expected = DataFrame({"B": pd.array([0, pd.NA], dtype="UInt64")}, index=idx)
    tm.assert_frame_equal(res, expected)

    res = gb.max()
    expected.iloc[0, 0] = 9
    tm.assert_frame_equal(res, expected)


@pytest.mark.parametrize("func", ["first", "last", "min", "max"])
def test_groupby_min_max_categorical(func):
    # GH: 52151
    df = DataFrame(
        {
            "col1": pd.Categorical(["A"], categories=list("AB"), ordered=True),
            "col2": pd.Categorical([1], categories=[1, 2], ordered=True),
            "value": 0.1,
        }
    )
    result = getattr(df.groupby("col1", observed=False), func)()

    idx = pd.CategoricalIndex(data=["A", "B"], name="col1", ordered=True)
    expected = DataFrame(
        {
            "col2": pd.Categorical([1, None], categories=[1, 2], ordered=True),
            "value": [0.1, None],
        },
        index=idx,
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("func", ["min", "max"])
def test_min_empty_string_dtype(func, string_dtype_no_object):
    # GH#55619
    dtype = string_dtype_no_object
    df = DataFrame({"a": ["a"], "b": "a", "c": "a"}, dtype=dtype).iloc[:0]
    result = getattr(df.groupby("a"), func)()
    expected = DataFrame(
        columns=["b", "c"], dtype=dtype, index=pd.Index([], dtype=dtype, name="a")
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("min_count", [0, 1])
@pytest.mark.parametrize("test_series", [True, False])
def test_string_dtype_all_na(
    string_dtype_no_object, reduction_func, skipna, min_count, test_series
):
    # https://github.com/pandas-dev/pandas/issues/60985
    if reduction_func == "corrwith":
        # corrwith is deprecated.
        return

    dtype = string_dtype_no_object

    if reduction_func in [
        "any",
        "all",
        "idxmin",
        "idxmax",
        "mean",
        "median",
        "std",
        "var",
    ]:
        kwargs = {"skipna": skipna}
    elif reduction_func in ["kurt"]:
        kwargs = {"min_count": min_count}
    elif reduction_func in ["count", "nunique", "quantile", "sem", "size"]:
        kwargs = {}
    else:
        kwargs = {"skipna": skipna, "min_count": min_count}

    expected_dtype, expected_value = dtype, pd.NA
    if reduction_func in ["all", "any"]:
        expected_dtype = "bool"
        # TODO: For skipna=False, bool(pd.NA) raises; should groupby?
        expected_value = not skipna if reduction_func == "any" else True
    elif reduction_func in ["count", "nunique", "size"]:
        # TODO: Should be more consistent - return Int64 when dtype.na_value is pd.NA?
        if (
            test_series
            and reduction_func == "size"
            and dtype.storage == "pyarrow"
            and dtype.na_value is pd.NA
        ):
            expected_dtype = "Int64"
        else:
            expected_dtype = "int64"
        expected_value = 1 if reduction_func == "size" else 0
    elif not skipna or min_count > 0:
        expected_value = pd.NA
    elif reduction_func == "sum":
        # https://github.com/pandas-dev/pandas/pull/60936
        expected_value = ""

    df = DataFrame({"a": ["x"], "b": [pd.NA]}, dtype=dtype)
    obj = df["b"] if test_series else df
    args = get_groupby_method_args(reduction_func, obj)
    gb = obj.groupby(df["a"])
    method = getattr(gb, reduction_func)

    if reduction_func in [
        "mean",
        "median",
        "kurt",
        "prod",
        "quantile",
        "sem",
        "skew",
        "std",
        "var",
    ]:
        msg = f"dtype '{dtype}' does not support operation '{reduction_func}'"
        with pytest.raises(TypeError, match=msg):
            method(*args, **kwargs)
        return
    elif reduction_func in ["idxmin", "idxmax"]:
        if skipna:
            msg = f"{reduction_func} with skipna=True encountered all NA values"
        else:
            msg = f"{reduction_func} with skipna=False encountered an NA value."
        with pytest.raises(ValueError, match=msg):
            method(*args, **kwargs)
        return

    result = method(*args, **kwargs)
    index = pd.Index(["x"], name="a", dtype=dtype)
    if test_series or reduction_func == "size":
        name = None if not test_series and reduction_func == "size" else "b"
        expected = Series(expected_value, index=index, dtype=expected_dtype, name=name)
    else:
        expected = DataFrame({"b": expected_value}, index=index, dtype=expected_dtype)
    tm.assert_equal(result, expected)


def test_max_nan_bug():
    df = DataFrame(
        {
            "Unnamed: 0": ["-04-23", "-05-06", "-05-07"],
            "Date": [
                "2013-04-23 00:00:00",
                "2013-05-06 00:00:00",
                "2013-05-07 00:00:00",
            ],
            "app": Series([np.nan, np.nan, "OE"]),
            "File": ["log080001.log", "log.log", "xlsx"],
        }
    )
    gb = df.groupby("Date")
    r = gb[["File"]].max()
    e = gb["File"].max().to_frame()
    tm.assert_frame_equal(r, e)
    assert not r["File"].isna().any()


@pytest.mark.slow
@pytest.mark.parametrize("with_nan", [True, False])
@pytest.mark.parametrize("keys", [["joe"], ["joe", "jim"]])
def test_series_groupby_nunique(sort, dropna, as_index, with_nan, keys):
    n = 100
    m = 10
    days = date_range("2015-08-23", periods=10)
    df = DataFrame(
        {
            "jim": np.random.default_rng(2).choice(list(ascii_lowercase), n),
            "joe": np.random.default_rng(2).choice(days, n),
            "julie": np.random.default_rng(2).integers(0, m, n),
        }
    )
    if with_nan:
        df = df.astype({"julie": float})  # Explicit cast to avoid implicit cast below
        df.loc[1::17, "jim"] = None
        df.loc[3::37, "joe"] = None
        df.loc[7::19, "julie"] = None
        df.loc[8::19, "julie"] = None
        df.loc[9::19, "julie"] = None
    original_df = df.copy()
    gr = df.groupby(keys, as_index=as_index, sort=sort)
    left = gr["julie"].nunique(dropna=dropna)

    gr = df.groupby(keys, as_index=as_index, sort=sort)
    right = gr["julie"].apply(Series.nunique, dropna=dropna)
    if not as_index:
        right = right.reset_index(drop=True)

    if as_index:
        tm.assert_series_equal(left, right, check_names=False)
    else:
        tm.assert_frame_equal(left, right, check_names=False)
    tm.assert_frame_equal(df, original_df)


def test_nunique():
    df = DataFrame({"A": list("abbacc"), "B": list("abxacc"), "C": list("abbacx")})

    expected = DataFrame({"A": list("abc"), "B": [1, 2, 1], "C": [1, 1, 2]})
    result = df.groupby("A", as_index=False).nunique()
    tm.assert_frame_equal(result, expected)

    # as_index
    expected.index = list("abc")
    expected.index.name = "A"
    expected = expected.drop(columns="A")
    result = df.groupby("A").nunique()
    tm.assert_frame_equal(result, expected)

    # with na
    result = df.replace({"x": None}).groupby("A").nunique(dropna=False)
    tm.assert_frame_equal(result, expected)

    # dropna
    expected = DataFrame({"B": [1] * 3, "C": [1] * 3}, index=list("abc"))
    expected.index.name = "A"
    result = df.replace({"x": None}).groupby("A").nunique()
    tm.assert_frame_equal(result, expected)


def test_nunique_with_object():
    # GH 11077
    data = DataFrame(
        [
            [100, 1, "Alice"],
            [200, 2, "Bob"],
            [300, 3, "Charlie"],
            [-400, 4, "Dan"],
            [500, 5, "Edith"],
        ],
        columns=["amount", "id", "name"],
    )

    result = data.groupby(["id", "amount"])["name"].nunique()
    index = MultiIndex.from_arrays([data.id, data.amount])
    expected = Series([1] * 5, name="name", index=index)
    tm.assert_series_equal(result, expected)


def test_nunique_with_empty_series():
    # GH 12553
    data = Series(name="name", dtype=object)
    result = data.groupby(level=0).nunique()
    expected = Series(name="name", dtype="int64")
    tm.assert_series_equal(result, expected)


def test_nunique_with_timegrouper():
    # GH 13453
    test = DataFrame(
        {
            "time": [
                Timestamp("2016-06-28 09:35:35"),
                Timestamp("2016-06-28 16:09:30"),
                Timestamp("2016-06-28 16:46:28"),
            ],
            "data": ["1", "2", "3"],
        }
    ).set_index("time")
    result = test.groupby(pd.Grouper(freq="h"))["data"].nunique()
    expected = test.groupby(pd.Grouper(freq="h"))["data"].apply(Series.nunique)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "key, data, dropna, expected",
    [
        (
            ["x", "x", "x"],
            [Timestamp("2019-01-01"), pd.NaT, Timestamp("2019-01-01")],
            True,
            Series([1], index=pd.Index(["x"], name="key"), name="data"),
        ),
        (
            ["x", "x", "x"],
            [dt.date(2019, 1, 1), pd.NaT, dt.date(2019, 1, 1)],
            True,
            Series([1], index=pd.Index(["x"], name="key"), name="data"),
        ),
        (
            ["x", "x", "x", "y", "y"],
            [
                dt.date(2019, 1, 1),
                pd.NaT,
                dt.date(2019, 1, 1),
                pd.NaT,
                dt.date(2019, 1, 1),
            ],
            False,
            Series([2, 2], index=pd.Index(["x", "y"], name="key"), name="data"),
        ),
        (
            ["x", "x", "x", "x", "y"],
            [
                dt.date(2019, 1, 1),
                pd.NaT,
                dt.date(2019, 1, 1),
                pd.NaT,
                dt.date(2019, 1, 1),
            ],
            False,
            Series([2, 1], index=pd.Index(["x", "y"], name="key"), name="data"),
        ),
    ],
)
def test_nunique_with_NaT(key, data, dropna, expected):
    # GH 27951
    df = DataFrame({"key": key, "data": data})
    result = df.groupby(["key"])["data"].nunique(dropna=dropna)
    tm.assert_series_equal(result, expected)


def test_nunique_preserves_column_level_names():
    # GH 23222
    test = DataFrame([1, 2, 2], columns=pd.Index(["A"], name="level_0"))
    result = test.groupby([0, 0, 0]).nunique()
    expected = DataFrame([2], index=np.array([0]), columns=test.columns)
    tm.assert_frame_equal(result, expected)


def test_nunique_transform_with_datetime():
    # GH 35109 - transform with nunique on datetimes results in integers
    df = DataFrame(date_range("2008-12-31", "2009-01-02"), columns=["date"])
    result = df.groupby([0, 0, 1])["date"].transform("nunique")
    expected = Series([2, 2, 1], name="date")
    tm.assert_series_equal(result, expected)


def test_empty_categorical(observed):
    # GH#21334
    cat = Series([1]).astype("category")
    ser = cat[:0]
    gb = ser.groupby(ser, observed=observed)
    result = gb.nunique()
    if observed:
        expected = Series([], index=cat[:0], dtype="int64")
    else:
        expected = Series([0], index=cat, dtype="int64")
    tm.assert_series_equal(result, expected)


def test_intercept_builtin_sum():
    s = Series([1.0, 2.0, np.nan, 3.0])
    grouped = s.groupby([0, 1, 2, 2])

    # GH#53425
    result = grouped.agg(builtins.sum)
    # GH#53425
    result2 = grouped.apply(builtins.sum)
    expected = Series([1.0, 2.0, np.nan], index=np.array([0, 1, 2]))
    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(result2, expected)


@pytest.mark.parametrize("min_count", [0, 10])
def test_groupby_sum_mincount_boolean(min_count):
    b = True
    a = False
    na = np.nan
    dfg = pd.array([b, b, na, na, a, a, b], dtype="boolean")

    df = DataFrame({"A": [1, 1, 2, 2, 3, 3, 1], "B": dfg})
    result = df.groupby("A").sum(min_count=min_count)
    if min_count == 0:
        expected = DataFrame(
            {"B": pd.array([3, 0, 0], dtype="Int64")},
            index=pd.Index([1, 2, 3], name="A"),
        )
        tm.assert_frame_equal(result, expected)
    else:
        expected = DataFrame(
            {"B": pd.array([pd.NA] * 3, dtype="Int64")},
            index=pd.Index([1, 2, 3], name="A"),
        )
        tm.assert_frame_equal(result, expected)


def test_groupby_sum_below_mincount_nullable_integer():
    # https://github.com/pandas-dev/pandas/issues/32861
    df = DataFrame({"a": [0, 1, 2], "b": [0, 1, 2], "c": [0, 1, 2]}, dtype="Int64")
    grouped = df.groupby("a")
    idx = pd.Index([0, 1, 2], name="a", dtype="Int64")

    result = grouped["b"].sum(min_count=2)
    expected = Series([pd.NA] * 3, dtype="Int64", index=idx, name="b")
    tm.assert_series_equal(result, expected)

    result = grouped.sum(min_count=2)
    expected = DataFrame({"b": [pd.NA] * 3, "c": [pd.NA] * 3}, dtype="Int64", index=idx)
    tm.assert_frame_equal(result, expected)


def test_groupby_sum_timedelta_with_nat():
    # GH#42659
    df = DataFrame(
        {
            "a": [1, 1, 2, 2],
            "b": [pd.Timedelta("1D"), pd.Timedelta("2D"), pd.Timedelta("3D"), pd.NaT],
        }
    )
    td3 = pd.Timedelta(days=3)

    gb = df.groupby("a")

    res = gb.sum()
    expected = DataFrame({"b": [td3, td3]}, index=pd.Index([1, 2], name="a"))
    tm.assert_frame_equal(res, expected)

    res = gb["b"].sum()
    tm.assert_series_equal(res, expected["b"])

    res = gb["b"].sum(min_count=2)
    expected = Series([td3, pd.NaT], dtype="m8[ns]", name="b", index=expected.index)
    tm.assert_series_equal(res, expected)


@pytest.mark.parametrize(
    "dtype", ["int8", "int16", "int32", "int64", "float32", "float64", "uint64"]
)
@pytest.mark.parametrize(
    "method,data",
    [
        ("first", {"df": [{"a": 1, "b": 1}, {"a": 2, "b": 3}]}),
        ("last", {"df": [{"a": 1, "b": 2}, {"a": 2, "b": 4}]}),
        ("min", {"df": [{"a": 1, "b": 1}, {"a": 2, "b": 3}]}),
        ("max", {"df": [{"a": 1, "b": 2}, {"a": 2, "b": 4}]}),
        ("count", {"df": [{"a": 1, "b": 2}, {"a": 2, "b": 2}], "out_type": "int64"}),
    ],
)
def test_groupby_non_arithmetic_agg_types(dtype, method, data):
    # GH9311, GH6620
    df = DataFrame(
        [{"a": 1, "b": 1}, {"a": 1, "b": 2}, {"a": 2, "b": 3}, {"a": 2, "b": 4}]
    )

    df["b"] = df.b.astype(dtype)

    if "args" not in data:
        data["args"] = []

    if "out_type" in data:
        out_type = data["out_type"]
    else:
        out_type = dtype

    exp = data["df"]
    df_out = DataFrame(exp)

    df_out["b"] = df_out.b.astype(out_type)
    df_out.set_index("a", inplace=True)

    grpd = df.groupby("a")
    t = getattr(grpd, method)(*data["args"])
    tm.assert_frame_equal(t, df_out)


def scipy_sem(*args, **kwargs):
    from scipy.stats import sem

    return sem(*args, ddof=1, **kwargs)


@pytest.mark.parametrize(
    "op,targop",
    [
        ("mean", np.mean),
        ("median", np.median),
        ("std", np.std),
        ("var", np.var),
        ("sum", np.sum),
        ("prod", np.prod),
        ("min", np.min),
        ("max", np.max),
        ("first", lambda x: x.iloc[0]),
        ("last", lambda x: x.iloc[-1]),
        ("count", np.size),
        pytest.param("sem", scipy_sem, marks=td.skip_if_no("scipy")),
    ],
)
def test_ops_general(op, targop):
    df = DataFrame(np.random.default_rng(2).standard_normal(1000))
    labels = np.random.default_rng(2).integers(0, 50, size=1000).astype(float)

    result = getattr(df.groupby(labels), op)()
    kwargs = {"ddof": 1, "axis": 0} if op in ["std", "var"] else {}
    expected = df.groupby(labels).agg(targop, **kwargs)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "values",
    [
        {
            "a": [1, 1, 1, 2, 2, 2, 3, 3, 3],
            "b": [1, pd.NA, 2, 1, pd.NA, 2, 1, pd.NA, 2],
        },
        {"a": [1, 1, 2, 2, 3, 3], "b": [1, 2, 1, 2, 1, 2]},
    ],
)
@pytest.mark.parametrize("function", ["mean", "median", "var"])
def test_apply_to_nullable_integer_returns_float(values, function):
    # https://github.com/pandas-dev/pandas/issues/32219
    output = 0.5 if function == "var" else 1.5
    arr = np.array([output] * 3, dtype=float)
    idx = pd.Index([1, 2, 3], name="a", dtype="Int64")
    expected = DataFrame({"b": arr}, index=idx).astype("Float64")

    groups = DataFrame(values, dtype="Int64").groupby("a")

    result = getattr(groups, function)()
    tm.assert_frame_equal(result, expected)

    result = groups.agg(function)
    tm.assert_frame_equal(result, expected)

    result = groups.agg([function])
    expected.columns = MultiIndex.from_tuples([("b", function)])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "op",
    [
        "sum",
        "prod",
        "min",
        "max",
        "median",
        "mean",
        "skew",
        "kurt",
        "std",
        "var",
        "sem",
    ],
)
def test_regression_allowlist_methods(op, skipna, sort):
    # GH6944
    # GH 17537
    # explicitly test the allowlist methods
    frame = DataFrame([0])

    grouped = frame.groupby(level=0, sort=sort)

    if op in ["skew", "kurt", "sum", "mean"]:
        # skew, kurt, sum, mean have skipna
        result = getattr(grouped, op)(skipna=skipna)
        expected = frame.groupby(level=0).apply(lambda h: getattr(h, op)(skipna=skipna))
        if sort:
            expected = expected.sort_index()
        tm.assert_frame_equal(result, expected)
    else:
        result = getattr(grouped, op)()
        expected = frame.groupby(level=0).apply(lambda h: getattr(h, op)())
        if sort:
            expected = expected.sort_index()
        tm.assert_frame_equal(result, expected)


def test_groupby_prod_with_int64_dtype():
    # GH#46573
    data = [
        [1, 11],
        [1, 41],
        [1, 17],
        [1, 37],
        [1, 7],
        [1, 29],
        [1, 31],
        [1, 2],
        [1, 3],
        [1, 43],
        [1, 5],
        [1, 47],
        [1, 19],
        [1, 88],
    ]
    df = DataFrame(data, columns=["A", "B"], dtype="int64")
    result = df.groupby(["A"]).prod().reset_index()
    expected = DataFrame({"A": [1], "B": [180970905912331920]}, dtype="int64")
    tm.assert_frame_equal(result, expected)


def test_groupby_std_datetimelike():
    # GH#48481
    tdi = pd.timedelta_range("1 Day", periods=10000)
    ser = Series(tdi)
    ser[::5] *= 2  # get different std for different groups

    df = ser.to_frame("A").copy()

    df["B"] = ser + Timestamp(0)
    df["C"] = ser + Timestamp(0, tz="UTC")
    df.iloc[-1] = pd.NaT  # last group includes NaTs

    gb = df.groupby(list(range(5)) * 2000)

    result = gb.std()

    # Note: this does not _exactly_ match what we would get if we did
    # [gb.get_group(i).std() for i in gb.groups]
    #  but it _does_ match the floating point error we get doing the
    #  same operation on int64 data xref GH#51332
    td1 = pd.Timedelta("2887 days 11:21:02.326710176")
    td4 = pd.Timedelta("2886 days 00:42:34.664668096")
    exp_ser = Series([td1 * 2, td1, td1, td1, td4], index=np.arange(5))
    expected = DataFrame({"A": exp_ser, "B": exp_ser, "C": exp_ser})
    tm.assert_frame_equal(result, expected)
