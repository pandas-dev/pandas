from datetime import (
    date,
    datetime,
)

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    bdate_range,
)
import pandas._testing as tm
from pandas.tests.groupby import get_groupby_method_args


def test_apply_func_that_appends_group_to_list_without_copy():
    # GH: 17718

    df = DataFrame(1, index=list(range(10)) * 10, columns=[0]).reset_index()
    groups = []

    def store(group):
        groups.append(group)

    df.groupby("index").apply(store)
    expected_value = DataFrame({0: [1] * 10}, index=pd.RangeIndex(0, 100, 10))
    expected_value.columns = expected_value.columns.astype(object)

    tm.assert_frame_equal(groups[0], expected_value)


def test_apply_index_date(using_infer_string):
    # GH 5788
    ts = [
        "2011-05-16 00:00",
        "2011-05-16 01:00",
        "2011-05-16 02:00",
        "2011-05-16 03:00",
        "2011-05-17 02:00",
        "2011-05-17 03:00",
        "2011-05-17 04:00",
        "2011-05-17 05:00",
        "2011-05-18 02:00",
        "2011-05-18 03:00",
        "2011-05-18 04:00",
        "2011-05-18 05:00",
    ]
    df = DataFrame(
        {
            "value": [
                1.40893,
                1.40760,
                1.40750,
                1.40649,
                1.40893,
                1.40760,
                1.40750,
                1.40649,
                1.40893,
                1.40760,
                1.40750,
                1.40649,
            ],
        },
        index=Index(pd.to_datetime(ts), name="date_time"),
    )
    expected = df.groupby(df.index.date).idxmax()
    result = df.groupby(df.index.date).apply(lambda x: x.idxmax())
    tm.assert_frame_equal(result, expected)


def test_apply_index_date_object():
    # GH 5789
    # don't auto coerce dates
    ts = [
        "2011-05-16 00:00",
        "2011-05-16 01:00",
        "2011-05-16 02:00",
        "2011-05-16 03:00",
        "2011-05-17 02:00",
        "2011-05-17 03:00",
        "2011-05-17 04:00",
        "2011-05-17 05:00",
        "2011-05-18 02:00",
        "2011-05-18 03:00",
        "2011-05-18 04:00",
        "2011-05-18 05:00",
    ]
    df = DataFrame([row.split() for row in ts], columns=["date", "time"])
    df["value"] = [
        1.40893,
        1.40760,
        1.40750,
        1.40649,
        1.40893,
        1.40760,
        1.40750,
        1.40649,
        1.40893,
        1.40760,
        1.40750,
        1.40649,
    ]
    exp_idx = Index(["2011-05-16", "2011-05-17", "2011-05-18"], name="date")
    expected = Series(["00:00", "02:00", "02:00"], index=exp_idx)
    result = df.groupby("date").apply(lambda x: x["time"][x["value"].idxmax()])
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "df, group_names",
    [
        (DataFrame({"a": [1, 1, 1, 2, 3], "b": ["a", "a", "a", "b", "c"]}), [1, 2, 3]),
        (DataFrame({"a": [0, 0, 1, 1], "b": [0, 1, 0, 1]}), [0, 1]),
        (DataFrame({"a": [1]}), [1]),
        (DataFrame({"a": [1, 1, 1, 2, 2, 1, 1, 2], "b": range(8)}), [1, 2]),
        (DataFrame({"a": [1, 2, 3, 1, 2, 3], "two": [4, 5, 6, 7, 8, 9]}), [1, 2, 3]),
        (
            DataFrame(
                {
                    "a": list("aaabbbcccc"),
                    "B": [3, 4, 3, 6, 5, 2, 1, 9, 5, 4],
                    "C": [4, 0, 2, 2, 2, 7, 8, 6, 2, 8],
                }
            ),
            ["a", "b", "c"],
        ),
        (DataFrame([[1, 2, 3], [2, 2, 3]], columns=["a", "b", "c"]), [1, 2]),
    ],
    ids=[
        "GH2936",
        "GH7739 & GH10519",
        "GH10519",
        "GH2656",
        "GH12155",
        "GH20084",
        "GH21417",
    ],
)
def test_group_apply_once_per_group(df, group_names):
    # GH2936, GH7739, GH10519, GH2656, GH12155, GH20084, GH21417

    # This test should ensure that a function is only evaluated
    # once per group. Previously the function has been evaluated twice
    # on the first group to check if the Cython index slider is safe to use
    # This test ensures that the side effect (append to list) is only triggered
    # once per group

    names = []
    # cannot parameterize over the functions since they need external
    # `names` to detect side effects

    def f_copy(group):
        # this takes the fast apply path
        names.append(group.name)
        return group.copy()

    def f_nocopy(group):
        # this takes the slow apply path
        names.append(group.name)
        return group

    def f_scalar(group):
        # GH7739, GH2656
        names.append(group.name)
        return 0

    def f_none(group):
        # GH10519, GH12155, GH21417
        names.append(group.name)

    def f_constant_df(group):
        # GH2936, GH20084
        names.append(group.name)
        return DataFrame({"a": [1], "b": [1]})

    for func in [f_copy, f_nocopy, f_scalar, f_none, f_constant_df]:
        del names[:]

        df.groupby("a").apply(func)
        assert names == group_names


def test_group_apply_once_per_group2(capsys):
    # GH: 31111
    # groupby-apply need to execute len(set(group_by_columns)) times

    expected = 2  # Number of times `apply` should call a function for the current test

    df = DataFrame(
        {
            "group_by_column": [0, 0, 0, 0, 1, 1, 1, 1],
            "test_column": ["0", "2", "4", "6", "8", "10", "12", "14"],
        },
        index=["0", "2", "4", "6", "8", "10", "12", "14"],
    )

    df.groupby("group_by_column", group_keys=False).apply(
        lambda df: print("function_called")
    )

    result = capsys.readouterr().out.count("function_called")
    # If `groupby` behaves unexpectedly, this test will break
    assert result == expected


def test_apply_fast_slow_identical():
    # GH 31613

    df = DataFrame({"A": [0, 0, 1], "b": range(3)})

    # For simple index structures we check for fast/slow apply using
    # an identity check on in/output
    def slow(group):
        return group

    def fast(group):
        return group.copy()

    fast_df = df.groupby("A", group_keys=False).apply(fast)
    slow_df = df.groupby("A", group_keys=False).apply(slow)
    tm.assert_frame_equal(fast_df, slow_df)


def test_apply_fast_slow_identical_index():
    # GH#44803
    df = DataFrame(
        {
            "name": ["Alice", "Bob", "Carl"],
            "age": [20, 21, 20],
        }
    ).set_index("name")

    grp_by_same_value = df.groupby(["age"], group_keys=False).apply(lambda group: group)
    grp_by_copy = df.groupby(["age"], group_keys=False).apply(
        lambda group: group.copy()
    )
    tm.assert_frame_equal(grp_by_same_value, grp_by_copy)


@pytest.mark.parametrize(
    "func",
    [
        lambda x: x,
        lambda x: x[:],
        lambda x: x.copy(deep=False),
        lambda x: x.copy(deep=True),
    ],
)
def test_groupby_apply_identity_maybecopy_index_identical(func):
    # GH 14927
    # Whether the function returns a copy of the input data or not should not
    # have an impact on the index structure of the result since this is not
    # transparent to the user

    df = DataFrame({"g": [1, 2, 2, 2], "a": [1, 2, 3, 4], "b": [5, 6, 7, 8]})
    result = df.groupby("g", group_keys=False).apply(func)
    tm.assert_frame_equal(result, df[["a", "b"]])


def test_apply_with_mixed_dtype():
    # GH3480, apply with mixed dtype on axis=1 breaks in 0.11
    df = DataFrame(
        {
            "foo1": np.random.default_rng(2).standard_normal(6),
            "foo2": ["one", "two", "two", "three", "one", "two"],
        }
    )
    result = df.apply(lambda x: x, axis=1)
    expected = df
    tm.assert_frame_equal(result, expected)

    # GH 3610 incorrect dtype conversion with as_index=False
    df = DataFrame({"c1": [1, 2, 6, 6, 8]})
    df["c2"] = df.c1 / 2.0
    result1 = df.groupby("c2").mean().reset_index()
    result2 = df.groupby("c2", as_index=False).mean()
    tm.assert_frame_equal(result1, result2)


def test_groupby_as_index_apply(as_index):
    # GH #4648 and #3417
    df = DataFrame(
        {
            "item_id": ["b", "b", "a", "c", "a", "b"],
            "user_id": [1, 2, 1, 1, 3, 1],
            "time": range(6),
        }
    )
    gb = df.groupby("user_id", as_index=as_index)

    expected = DataFrame(
        {
            "item_id": ["b", "b", "a", "a"],
            "user_id": [1, 2, 1, 3],
            "time": [0, 1, 2, 4],
        },
        index=[0, 1, 2, 4],
    )
    result = gb.head(2)
    tm.assert_frame_equal(result, expected)

    # apply doesn't maintain the original ordering
    # changed in GH5610 as the as_index=False returns a MI here
    if as_index:
        tp = [(1, 0), (1, 2), (2, 1), (3, 4)]
        index = MultiIndex.from_tuples(tp, names=["user_id", None])
    else:
        index = Index([0, 2, 1, 4])
    expected = DataFrame(
        {
            "item_id": list("baba"),
            "time": [0, 2, 1, 4],
        },
        index=index,
    )
    result = gb.apply(lambda x: x.head(2))
    tm.assert_frame_equal(result, expected)


def test_groupby_as_index_apply_str():
    ind = Index(list("abcde"))
    df = DataFrame([[1, 2], [2, 3], [1, 4], [1, 5], [2, 6]], index=ind)
    res = df.groupby(0, as_index=False, group_keys=False).apply(lambda x: x).index
    tm.assert_index_equal(res, ind)


def test_apply_concat_preserve_names(three_group):
    grouped = three_group.groupby(["A", "B"])

    def desc(group):
        result = group.describe()
        result.index.name = "stat"
        return result

    def desc2(group):
        result = group.describe()
        result.index.name = "stat"
        result = result[: len(group)]
        # weirdo
        return result

    def desc3(group):
        result = group.describe()

        # names are different
        result.index.name = f"stat_{len(group):d}"

        result = result[: len(group)]
        # weirdo
        return result

    result = grouped.apply(desc)
    assert result.index.names == ("A", "B", "stat")

    result2 = grouped.apply(desc2)
    assert result2.index.names == ("A", "B", "stat")

    result3 = grouped.apply(desc3)
    assert result3.index.names == ("A", "B", None)


def test_apply_series_to_frame():
    def f(piece):
        with np.errstate(invalid="ignore"):
            logged = np.log(piece)
        return DataFrame(
            {"value": piece, "demeaned": piece - piece.mean(), "logged": logged}
        )

    dr = bdate_range("1/1/2000", periods=10)
    ts = Series(np.random.default_rng(2).standard_normal(10), index=dr)

    grouped = ts.groupby(lambda x: x.month, group_keys=False)
    result = grouped.apply(f)

    assert isinstance(result, DataFrame)
    assert not hasattr(result, "name")  # GH49907
    tm.assert_index_equal(result.index, ts.index)


def test_apply_series_yield_constant(df):
    result = df.groupby(["A", "B"])["C"].apply(len)
    assert result.index.names[:2] == ("A", "B")


def test_apply_frame_yield_constant(df):
    # GH13568
    result = df.groupby(["A", "B"]).apply(len)
    assert isinstance(result, Series)
    assert result.name is None

    result = df.groupby(["A", "B"])[["C", "D"]].apply(len)
    assert isinstance(result, Series)
    assert result.name is None


def test_apply_frame_to_series(df):
    grouped = df.groupby(["A", "B"])
    result = grouped.apply(len)
    expected = grouped.count()["C"]
    tm.assert_index_equal(result.index, expected.index)
    tm.assert_numpy_array_equal(result.values, expected.values)


def test_apply_frame_not_as_index_column_name(df):
    # GH 35964 - path within _wrap_applied_output not hit by a test
    grouped = df.groupby(["A", "B"], as_index=False)
    result = grouped.apply(len)
    expected = grouped.count().rename(columns={"C": np.nan}).drop(columns="D")
    # TODO(GH#34306): Use assert_frame_equal when column name is not np.nan
    tm.assert_index_equal(result.index, expected.index)
    tm.assert_numpy_array_equal(result.values, expected.values)


def test_apply_frame_concat_series():
    def trans(group):
        return group.groupby("B")["C"].sum().sort_values().iloc[:2]

    def trans2(group):
        grouped = group.groupby(df.reindex(group.index)["B"])
        return grouped.sum().sort_values().iloc[:2]

    df = DataFrame(
        {
            "A": np.random.default_rng(2).integers(0, 5, 1000),
            "B": np.random.default_rng(2).integers(0, 5, 1000),
            "C": np.random.default_rng(2).standard_normal(1000),
        }
    )

    result = df.groupby("A").apply(trans)
    exp = df.groupby("A")["C"].apply(trans2)
    tm.assert_series_equal(result, exp, check_names=False)
    assert result.name == "C"


def test_apply_transform(ts):
    grouped = ts.groupby(lambda x: x.month, group_keys=False)
    result = grouped.apply(lambda x: x * 2)
    expected = grouped.transform(lambda x: x * 2)
    tm.assert_series_equal(result, expected)


def test_apply_multikey_corner(tsframe):
    grouped = tsframe.groupby([lambda x: x.year, lambda x: x.month])

    def f(group):
        return group.sort_values("A")[-5:]

    result = grouped.apply(f)
    for key, group in grouped:
        tm.assert_frame_equal(result.loc[key], f(group))


@pytest.mark.parametrize("group_keys", [True, False])
def test_apply_chunk_view(group_keys):
    # Low level tinkering could be unsafe, make sure not
    df = DataFrame({"key": [1, 1, 1, 2, 2, 2, 3, 3, 3], "value": range(9)})

    result = df.groupby("key", group_keys=group_keys).apply(lambda x: x.iloc[:2])
    expected = df[["value"]].take([0, 1, 3, 4, 6, 7])
    if group_keys:
        expected.index = MultiIndex.from_arrays(
            [[1, 1, 2, 2, 3, 3], expected.index], names=["key", None]
        )

    tm.assert_frame_equal(result, expected)


def test_apply_no_name_column_conflict():
    df = DataFrame(
        {
            "name": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
            "name2": [0, 0, 0, 1, 1, 1, 0, 0, 1, 1],
            "value": range(9, -1, -1),
        }
    )

    # it works! #2605
    grouped = df.groupby(["name", "name2"])
    grouped.apply(lambda x: x.sort_values("value", inplace=True))


def test_apply_typecast_fail():
    df = DataFrame(
        {
            "d": [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
            "c": np.tile(["a", "b", "c"], 2),
            "v": np.arange(1.0, 7.0),
        }
    )

    def f(group):
        v = group["v"]
        group["v2"] = (v - v.min()) / (v.max() - v.min())
        return group

    result = df.groupby("d", group_keys=False).apply(f)

    expected = df[["c", "v"]]
    expected["v2"] = np.tile([0.0, 0.5, 1], 2)

    tm.assert_frame_equal(result, expected)


def test_apply_multiindex_fail():
    index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1], [1, 2, 3, 1, 2, 3]])
    df = DataFrame(
        {
            "d": [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
            "c": np.tile(["a", "b", "c"], 2),
            "v": np.arange(1.0, 7.0),
        },
        index=index,
    )

    def f(group):
        v = group["v"]
        group["v2"] = (v - v.min()) / (v.max() - v.min())
        return group

    result = df.groupby("d", group_keys=False).apply(f)

    expected = df[["c", "v"]]
    expected["v2"] = np.tile([0.0, 0.5, 1], 2)
    tm.assert_frame_equal(result, expected)


def test_apply_corner(tsframe):
    result = tsframe.groupby(lambda x: x.year, group_keys=False).apply(lambda x: x * 2)
    expected = tsframe * 2
    tm.assert_frame_equal(result, expected)


def test_apply_without_copy():
    # GH 5545
    # returning a non-copy in an applied function fails

    data = DataFrame(
        {
            "id_field": [100, 100, 200, 300],
            "category": ["a", "b", "c", "c"],
            "value": [1, 2, 3, 4],
        }
    )

    def filt1(x):
        if x.shape[0] == 1:
            return x.copy()
        else:
            return x[x.category == "c"]

    def filt2(x):
        if x.shape[0] == 1:
            return x
        else:
            return x[x.category == "c"]

    expected = data.groupby("id_field").apply(filt1)
    result = data.groupby("id_field").apply(filt2)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("test_series", [True, False])
def test_apply_with_duplicated_non_sorted_axis(test_series):
    # GH 30667
    df = DataFrame(
        [["x", "p"], ["x", "p"], ["x", "o"]], columns=["X", "Y"], index=[1, 2, 2]
    )
    if test_series:
        ser = df.set_index("Y")["X"]
        result = ser.groupby(level=0, group_keys=False).apply(lambda x: x)
        expected = ser
        tm.assert_series_equal(result, expected)
    else:
        result = df.groupby("Y", group_keys=False).apply(lambda x: x)
        expected = df[["X"]]
        tm.assert_frame_equal(result, expected)


def test_apply_reindex_values():
    # GH: 26209
    # reindexing from a single column of a groupby object with duplicate indices caused
    # a ValueError (cannot reindex from duplicate axis) in 0.24.2, the problem was
    # solved in #30679
    values = [1, 2, 3, 4]
    indices = [1, 1, 2, 2]
    df = DataFrame({"group": ["Group1", "Group2"] * 2, "value": values}, index=indices)
    expected = Series(values, index=indices, name="value")

    def reindex_helper(x):
        return x.reindex(np.arange(x.index.min(), x.index.max() + 1))

    # the following group by raised a ValueError
    result = df.groupby("group", group_keys=False).value.apply(reindex_helper)
    tm.assert_series_equal(expected, result)


def test_apply_corner_cases():
    # #535, can't use sliding iterator

    N = 10
    labels = np.random.default_rng(2).integers(0, 100, size=N)
    df = DataFrame(
        {
            "key": labels,
            "value1": np.random.default_rng(2).standard_normal(N),
            "value2": ["foo", "bar", "baz", "qux", "a"] * (N // 5),
        }
    )

    grouped = df.groupby("key", group_keys=False)

    def f(g):
        g["value3"] = g["value1"] * 2
        return g

    result = grouped.apply(f)
    assert "value3" in result


def test_apply_numeric_coercion_when_datetime():
    # In the past, group-by/apply operations have been over-eager
    # in converting dtypes to numeric, in the presence of datetime
    # columns.  Various GH issues were filed, the reproductions
    # for which are here.

    # GH 15670
    df = DataFrame(
        {"Number": [1, 2], "Date": ["2017-03-02"] * 2, "Str": ["foo", "inf"]}
    )
    expected = df.groupby(["Number"]).apply(lambda x: x.iloc[0])
    df.Date = pd.to_datetime(df.Date)
    result = df.groupby(["Number"]).apply(lambda x: x.iloc[0])
    tm.assert_series_equal(result["Str"], expected["Str"])


def test_apply_numeric_coercion_when_datetime_getitem():
    # GH 15421
    df = DataFrame(
        {"A": [10, 20, 30], "B": ["foo", "3", "4"], "T": [pd.Timestamp("12:31:22")] * 3}
    )

    def get_B(g):
        return g.iloc[0][["B"]]

    result = df.groupby("A").apply(get_B)["B"]
    expected = df.B
    expected.index = df.A
    tm.assert_series_equal(result, expected)


def test_apply_numeric_coercion_when_datetime_with_nat():
    # GH 14423
    def predictions(tool):
        out = Series(index=["p1", "p2", "useTime"], dtype=object)
        if "step1" in list(tool.State):
            out["p1"] = str(tool[tool.State == "step1"].Machine.values[0])
        if "step2" in list(tool.State):
            out["p2"] = str(tool[tool.State == "step2"].Machine.values[0])
            out["useTime"] = str(tool[tool.State == "step2"].oTime.values[0])
        return out

    df1 = DataFrame(
        {
            "Key": ["B", "B", "A", "A"],
            "State": ["step1", "step2", "step1", "step2"],
            "oTime": ["", "2016-09-19 05:24:33", "", "2016-09-19 23:59:04"],
            "Machine": ["23", "36L", "36R", "36R"],
        }
    )
    df2 = df1.copy()
    df2.oTime = pd.to_datetime(df2.oTime)
    expected = df1.groupby("Key").apply(predictions).p1
    result = df2.groupby("Key").apply(predictions).p1
    tm.assert_series_equal(expected, result)


def test_apply_aggregating_timedelta_and_datetime():
    # Regression test for GH 15562
    # The following groupby caused ValueErrors and IndexErrors pre 0.20.0

    df = DataFrame(
        {
            "clientid": ["A", "B", "C"],
            "datetime": [np.datetime64("2017-02-01 00:00:00")] * 3,
        }
    )
    df["time_delta_zero"] = df.datetime - df.datetime
    result = df.groupby("clientid").apply(
        lambda ddf: Series(
            {"clientid_age": ddf.time_delta_zero.min(), "date": ddf.datetime.min()}
        )
    )
    expected = DataFrame(
        {
            "clientid": ["A", "B", "C"],
            "clientid_age": [np.timedelta64(0, "D")] * 3,
            "date": [np.datetime64("2017-02-01 00:00:00")] * 3,
        }
    ).set_index("clientid")

    tm.assert_frame_equal(result, expected)


def test_apply_groupby_datetimeindex():
    # GH 26182
    # groupby apply failed on dataframe with DatetimeIndex

    data = [["A", 10], ["B", 20], ["B", 30], ["C", 40], ["C", 50]]
    df = DataFrame(
        data, columns=["Name", "Value"], index=pd.date_range("2020-09-01", "2020-09-05")
    )

    result = df.groupby("Name").sum()

    expected = DataFrame({"Name": ["A", "B", "C"], "Value": [10, 50, 90]})
    expected.set_index("Name", inplace=True)

    tm.assert_frame_equal(result, expected)


def test_time_field_bug():
    # Test a fix for the following error related to GH issue 11324 When
    # non-key fields in a group-by dataframe contained time-based fields
    # that were not returned by the apply function, an exception would be
    # raised.

    df = DataFrame({"a": 1, "b": [datetime.now() for nn in range(10)]})

    def func_with_no_date(batch):
        return Series({"c": 2})

    def func_with_date(batch):
        return Series({"b": datetime(2015, 1, 1), "c": 2})

    dfg_no_conversion = df.groupby(by=["a"]).apply(func_with_no_date)
    dfg_no_conversion_expected = DataFrame({"c": 2}, index=[1])
    dfg_no_conversion_expected.index.name = "a"

    dfg_conversion = df.groupby(by=["a"]).apply(func_with_date)
    dfg_conversion_expected = DataFrame(
        {"b": pd.Timestamp(2015, 1, 1), "c": 2}, index=[1]
    )
    dfg_conversion_expected.index.name = "a"

    tm.assert_frame_equal(dfg_no_conversion, dfg_no_conversion_expected)
    tm.assert_frame_equal(dfg_conversion, dfg_conversion_expected)


def test_gb_apply_list_of_unequal_len_arrays():
    # GH1738
    df = DataFrame(
        {
            "group1": ["a", "a", "a", "b", "b", "b", "a", "a", "a", "b", "b", "b"],
            "group2": ["c", "c", "d", "d", "d", "e", "c", "c", "d", "d", "d", "e"],
            "weight": [1.1, 2, 3, 4, 5, 6, 2, 4, 6, 8, 1, 2],
            "value": [7.1, 8, 9, 10, 11, 12, 8, 7, 6, 5, 4, 3],
        }
    )
    df = df.set_index(["group1", "group2"])
    df_grouped = df.groupby(level=["group1", "group2"], sort=True)

    def noddy(value, weight):
        out = np.array(value * weight).repeat(3)
        return out

    # the kernel function returns arrays of unequal length
    # pandas sniffs the first one, sees it's an array and not
    # a list, and assumed the rest are of equal length
    # and so tries a vstack

    # don't die
    df_grouped.apply(lambda x: noddy(x.value, x.weight))


def test_groupby_apply_all_none():
    # Tests to make sure no errors if apply function returns all None
    # values. Issue 9684.
    test_df = DataFrame({"groups": [0, 0, 1, 1], "random_vars": [8, 7, 4, 5]})

    def test_func(x):
        pass

    result = test_df.groupby("groups").apply(test_func)
    expected = DataFrame(columns=["random_vars"], dtype="int64")
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "in_data, out_idx, out_data",
    [
        [
            {"groups": [1, 1, 1, 2], "vars": [0, 1, 2, 3]},
            [[1, 1], [0, 2]],
            {"vars": [0, 2]},
        ],
        [
            {"groups": [1, 2, 2, 2], "vars": [0, 1, 2, 3]},
            [[2, 2], [1, 3]],
            {"vars": [1, 3]},
        ],
    ],
)
def test_groupby_apply_none_first(in_data, out_idx, out_data):
    # GH 12824. Tests if apply returns None first.
    test_df1 = DataFrame(in_data)

    def test_func(x):
        if x.shape[0] < 2:
            return None
        return x.iloc[[0, -1]]

    result1 = test_df1.groupby("groups").apply(test_func)
    index1 = MultiIndex.from_arrays(out_idx, names=["groups", None])
    expected1 = DataFrame(out_data, index=index1)
    tm.assert_frame_equal(result1, expected1)


def test_groupby_apply_return_empty_chunk():
    # GH 22221: apply filter which returns some empty groups
    df = DataFrame({"value": [0, 1], "group": ["filled", "empty"]})
    groups = df.groupby("group")
    result = groups.apply(lambda group: group[group.value != 1]["value"])
    expected = Series(
        [0],
        name="value",
        index=MultiIndex.from_product(
            [["empty", "filled"], [0]], names=["group", None]
        ).drop("empty"),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("meth", ["apply", "transform"])
def test_apply_with_mixed_types(meth):
    # gh-20949
    df = DataFrame({"A": "a a b".split(), "B": [1, 2, 3], "C": [4, 6, 5]})
    g = df.groupby("A", group_keys=False)

    result = getattr(g, meth)(lambda x: x / x.sum())
    expected = DataFrame({"B": [1 / 3.0, 2 / 3.0, 1], "C": [0.4, 0.6, 1.0]})
    tm.assert_frame_equal(result, expected)


def test_func_returns_object():
    # GH 28652
    df = DataFrame({"a": [1, 2]}, index=Index([1, 2]))
    result = df.groupby("a").apply(lambda g: g.index)
    expected = Series([Index([1]), Index([2])], index=Index([1, 2], name="a"))

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "group_column_dtlike",
    [datetime.today(), datetime.today().date(), datetime.today().time()],
)
def test_apply_datetime_issue(group_column_dtlike):
    # GH-28247
    # groupby-apply throws an error if one of the columns in the DataFrame
    #   is a datetime object and the column labels are different from
    #   standard int values in range(len(num_columns))

    df = DataFrame({"a": ["foo"], "b": [group_column_dtlike]})
    result = df.groupby("a").apply(lambda x: Series(["spam"], index=[42]))

    expected = DataFrame(["spam"], Index(["foo"], dtype="str", name="a"), columns=[42])
    tm.assert_frame_equal(result, expected)


def test_apply_series_return_dataframe_groups():
    # GH 10078
    tdf = DataFrame(
        {
            "day": {
                0: pd.Timestamp("2015-02-24 00:00:00"),
                1: pd.Timestamp("2015-02-24 00:00:00"),
                2: pd.Timestamp("2015-02-24 00:00:00"),
                3: pd.Timestamp("2015-02-24 00:00:00"),
                4: pd.Timestamp("2015-02-24 00:00:00"),
            },
            "userAgent": {
                0: "some UA string",
                1: "some UA string",
                2: "some UA string",
                3: "another UA string",
                4: "some UA string",
            },
            "userId": {
                0: "17661101",
                1: "17661101",
                2: "17661101",
                3: "17661101",
                4: "17661101",
            },
        }
    )

    def most_common_values(df):
        return Series({c: s.value_counts().index[0] for c, s in df.items()})

    result = tdf.groupby("day").apply(most_common_values)["userId"]
    expected = Series(
        ["17661101"], index=pd.DatetimeIndex(["2015-02-24"], name="day"), name="userId"
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("category", [False, True])
def test_apply_multi_level_name(category):
    # https://github.com/pandas-dev/pandas/issues/31068
    b = [1, 2] * 5
    if category:
        b = pd.Categorical(b, categories=[1, 2, 3])
        expected_index = pd.CategoricalIndex([1, 2, 3], categories=[1, 2, 3], name="B")
        expected_values = [20, 25, 0]
    else:
        expected_index = Index([1, 2], name="B")
        expected_values = [20, 25]
    expected = DataFrame(
        {"C": expected_values, "D": expected_values}, index=expected_index
    )

    df = DataFrame(
        {"A": np.arange(10), "B": b, "C": list(range(10)), "D": list(range(10))}
    ).set_index(["A", "B"])
    result = df.groupby("B", observed=False).apply(lambda x: x.sum())
    tm.assert_frame_equal(result, expected)
    assert df.index.names == ["A", "B"]


def test_groupby_apply_datetime_result_dtypes(using_infer_string):
    # GH 14849
    data = DataFrame.from_records(
        [
            (pd.Timestamp(2016, 1, 1), "red", "dark", 1, "8"),
            (pd.Timestamp(2015, 1, 1), "green", "stormy", 2, "9"),
            (pd.Timestamp(2014, 1, 1), "blue", "bright", 3, "10"),
            (pd.Timestamp(2013, 1, 1), "blue", "calm", 4, "potato"),
        ],
        columns=["observation", "color", "mood", "intensity", "score"],
    )
    result = data.groupby("color").apply(lambda g: g.iloc[0]).dtypes
    dtype = pd.StringDtype(na_value=np.nan) if using_infer_string else object
    expected = Series(
        [np.dtype("datetime64[us]"), dtype, np.int64, dtype],
        index=["observation", "mood", "intensity", "score"],
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "index",
    [
        pd.CategoricalIndex(list("abc")),
        pd.interval_range(0, 3),
        pd.period_range("2020", periods=3, freq="D"),
        MultiIndex.from_tuples([("a", 0), ("a", 1), ("b", 0)]),
    ],
)
def test_apply_index_has_complex_internals(index):
    # GH 31248
    df = DataFrame({"group": [1, 1, 2], "value": [0, 1, 0]}, index=index)
    result = df.groupby("group", group_keys=False).apply(lambda x: x)
    tm.assert_frame_equal(result, df[["value"]])


@pytest.mark.parametrize(
    "function, expected_values",
    [
        (lambda x: x.index.to_list(), [[0, 1], [2, 3]]),
        (lambda x: set(x.index.to_list()), [{0, 1}, {2, 3}]),
        (lambda x: tuple(x.index.to_list()), [(0, 1), (2, 3)]),
        (
            lambda x: dict(enumerate(x.index.to_list())),
            [{0: 0, 1: 1}, {0: 2, 1: 3}],
        ),
        (
            lambda x: [{n: i} for (n, i) in enumerate(x.index.to_list())],
            [[{0: 0}, {1: 1}], [{0: 2}, {1: 3}]],
        ),
    ],
)
def test_apply_function_returns_non_pandas_non_scalar(function, expected_values):
    # GH 31441
    df = DataFrame(["A", "A", "B", "B"], columns=["groups"])
    result = df.groupby("groups").apply(function)
    expected = Series(expected_values, index=Index(["A", "B"], name="groups"))
    tm.assert_series_equal(result, expected)


def test_apply_function_returns_numpy_array():
    # GH 31605
    def fct(group):
        return group["B"].values.flatten()

    df = DataFrame({"A": ["a", "a", "b", "none"], "B": [1, 2, 3, np.nan]})

    result = df.groupby("A").apply(fct)
    expected = Series(
        [[1.0, 2.0], [3.0], [np.nan]], index=Index(["a", "b", "none"], name="A")
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("function", [lambda gr: gr.index, lambda gr: gr.index + 1 - 1])
def test_apply_function_index_return(function):
    # GH: 22541
    df = DataFrame([1, 2, 2, 2, 1, 2, 3, 1, 3, 1], columns=["id"])
    result = df.groupby("id").apply(function)
    expected = Series(
        [Index([0, 4, 7, 9]), Index([1, 2, 3, 5]), Index([6, 8])],
        index=Index([1, 2, 3], name="id"),
    )
    tm.assert_series_equal(result, expected)


def test_apply_function_with_indexing_return_column():
    # GH#7002, GH#41480, GH#49256
    df = DataFrame(
        {
            "foo1": ["one", "two", "two", "three", "one", "two"],
            "foo2": [1, 2, 4, 4, 5, 6],
        }
    )
    result = df.groupby("foo1", as_index=False).apply(lambda x: x.mean())
    expected = DataFrame(
        {
            "foo1": ["one", "three", "two"],
            "foo2": [3.0, 4.0, 4.0],
        }
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "udf",
    [lambda x: x.copy(), lambda x: x.copy().rename(lambda y: y + 1)],
)
@pytest.mark.parametrize("group_keys", [True, False])
def test_apply_result_type(group_keys, udf):
    # https://github.com/pandas-dev/pandas/issues/34809
    # We'd like to control whether the group keys end up in the index
    # regardless of whether the UDF happens to be a transform.
    df = DataFrame({"A": ["a", "b"], "B": [1, 2]})
    df_result = df.groupby("A", group_keys=group_keys).apply(udf)
    series_result = df.B.groupby(df.A, group_keys=group_keys).apply(udf)

    if group_keys:
        assert df_result.index.nlevels == 2
        assert series_result.index.nlevels == 2
    else:
        assert df_result.index.nlevels == 1
        assert series_result.index.nlevels == 1


def test_result_order_group_keys_false():
    # GH 34998
    # apply result order should not depend on whether index is the same or just equal
    df = DataFrame({"A": [2, 1, 2], "B": [1, 2, 3]})
    result = df.groupby("A", group_keys=False).apply(lambda x: x)
    expected = df.groupby("A", group_keys=False).apply(lambda x: x.copy())
    tm.assert_frame_equal(result, expected)


def test_apply_with_timezones_aware():
    # GH: 27212
    dates = ["2001-01-01"] * 2 + ["2001-01-02"] * 2 + ["2001-01-03"] * 2
    index_no_tz = pd.DatetimeIndex(dates)
    index_tz = pd.DatetimeIndex(dates, tz="UTC")
    df1 = DataFrame({"x": list(range(2)) * 3, "y": range(6), "t": index_no_tz})
    df2 = DataFrame({"x": list(range(2)) * 3, "y": range(6), "t": index_tz})

    result1 = df1.groupby("x", group_keys=False).apply(lambda df: df[["y"]].copy())
    result2 = df2.groupby("x", group_keys=False).apply(lambda df: df[["y"]].copy())

    tm.assert_frame_equal(result1, result2)


def test_apply_is_unchanged_when_other_methods_are_called_first(reduction_func):
    # GH #34656
    # GH #34271
    df = DataFrame(
        {
            "a": [99, 99, 99, 88, 88, 88],
            "b": [1, 2, 3, 4, 5, 6],
            "c": [10, 20, 30, 40, 50, 60],
        }
    )

    expected = DataFrame(
        {"b": [15, 6], "c": [150, 60]},
        index=Index([88, 99], name="a"),
    )

    # Check output when no other methods are called before .apply()
    grp = df.groupby(by="a")
    result = grp.apply(np.sum, axis=0)
    tm.assert_frame_equal(result, expected)

    # Check output when another method is called before .apply()
    grp = df.groupby(by="a")
    args = get_groupby_method_args(reduction_func, df)
    if reduction_func == "corrwith":
        warn = Pandas4Warning
        msg = "DataFrameGroupBy.corrwith is deprecated"
    else:
        warn = None
        msg = ""
    with tm.assert_produces_warning(warn, match=msg):
        _ = getattr(grp, reduction_func)(*args)
    result = grp.apply(np.sum, axis=0)
    tm.assert_frame_equal(result, expected)


def test_apply_with_date_in_multiindex_does_not_convert_to_timestamp():
    # GH 29617

    df = DataFrame(
        {
            "A": ["a", "a", "a", "b"],
            "B": [
                date(2020, 1, 10),
                date(2020, 1, 10),
                date(2020, 2, 10),
                date(2020, 2, 10),
            ],
            "C": [1, 2, 3, 4],
        },
        index=Index([100, 101, 102, 103], name="idx"),
    )

    grp = df.groupby(["A", "B"])
    result = grp.apply(lambda x: x.head(1))

    expected = df.iloc[[0, 2, 3]]
    expected = expected.reset_index()
    expected.index = MultiIndex.from_frame(expected[["A", "B", "idx"]])
    expected = expected.drop(columns=["A", "B", "idx"])

    tm.assert_frame_equal(result, expected)
    for val in result.index.levels[1]:
        assert type(val) is date


def test_apply_dropna_with_indexed_same(dropna):
    # GH 38227
    # GH#43205
    df = DataFrame(
        {
            "col": [1, 2, 3, 4, 5],
            "group": ["a", np.nan, np.nan, "b", "b"],
        },
        index=list("xxyxz"),
    )
    result = df.groupby("group", dropna=dropna, group_keys=False).apply(lambda x: x)
    expected = df.dropna()[["col"]] if dropna else df[["col"]].iloc[[0, 3, 1, 2, 4]]
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "as_index, expected",
    [
        [
            False,
            DataFrame(
                [[1, 1, 1], [2, 2, 1]], columns=Index(["a", "b", None], dtype=object)
            ),
        ],
        [
            True,
            Series(
                [1, 1], index=MultiIndex.from_tuples([(1, 1), (2, 2)], names=["a", "b"])
            ),
        ],
    ],
)
def test_apply_as_index_constant_lambda(as_index, expected):
    # GH 13217
    df = DataFrame({"a": [1, 1, 2, 2], "b": [1, 1, 2, 2], "c": [1, 1, 1, 1]})
    result = df.groupby(["a", "b"], as_index=as_index).apply(lambda x: 1)
    tm.assert_equal(result, expected)


def test_sort_index_groups():
    # GH 20420
    df = DataFrame(
        {"A": [1, 2, 3, 4, 5], "B": [6, 7, 8, 9, 0], "C": [1, 1, 1, 2, 2]},
        index=range(5),
    )
    result = df.groupby("C").apply(lambda x: x.A.sort_index())
    expected = Series(
        range(1, 6),
        index=MultiIndex.from_tuples(
            [(1, 0), (1, 1), (1, 2), (2, 3), (2, 4)], names=["C", None]
        ),
        name="A",
    )
    tm.assert_series_equal(result, expected)


def test_positional_slice_groups_datetimelike():
    # GH 21651
    expected = DataFrame(
        {
            "date": pd.date_range("2010-01-01", freq="12h", periods=5),
            "vals": range(5),
            "let": list("abcde"),
        }
    )
    result = expected.groupby(
        [expected.let, expected.date.dt.date], group_keys=False
    ).apply(lambda x: x.iloc[0:])
    tm.assert_frame_equal(result, expected[["date", "vals"]])


def test_groupby_apply_shape_cache_safety():
    # GH#42702 this fails if we cache_readonly Block.shape
    df = DataFrame({"A": ["a", "a", "b"], "B": [1, 2, 3], "C": [4, 6, 5]})
    gb = df.groupby("A")
    result = gb[["B", "C"]].apply(lambda x: x.astype(float).max() - x.min())

    expected = DataFrame(
        {"B": [1.0, 0.0], "C": [2.0, 0.0]}, index=Index(["a", "b"], name="A")
    )
    tm.assert_frame_equal(result, expected)


def test_groupby_apply_to_series_name():
    # GH52444
    df = DataFrame.from_dict(
        {
            "a": ["a", "b", "a", "b"],
            "b1": ["aa", "ac", "ac", "ad"],
            "b2": ["aa", "aa", "aa", "ac"],
        }
    )
    grp = df.groupby("a")[["b1", "b2"]]
    result = grp.apply(lambda x: x.unstack().value_counts())

    expected_idx = MultiIndex.from_arrays(
        arrays=[["a", "a", "b", "b", "b"], ["aa", "ac", "ac", "ad", "aa"]],
        names=["a", None],
    )
    expected = Series([3, 1, 2, 1, 1], index=expected_idx, name="count")
    tm.assert_series_equal(result, expected)


def test_apply_na(dropna):
    # GH#28984
    df = DataFrame(
        {"grp": [1, 1, 2, 2], "y": [1, 0, 2, 5], "z": [1, 2, np.nan, np.nan]}
    )
    dfgrp = df.groupby("grp", dropna=dropna)
    result = dfgrp.apply(lambda grp_df: grp_df.nlargest(1, "z"))
    expected = dfgrp.apply(lambda x: x.sort_values("z", ascending=False).head(1))
    tm.assert_frame_equal(result, expected)


def test_apply_empty_string_nan_coerce_bug():
    # GH#24903
    result = (
        DataFrame(
            {
                "a": [1, 1, 2, 2],
                "b": ["", "", "", ""],
                "c": pd.to_datetime([1, 2, 3, 4], unit="s"),
            }
        )
        .groupby(["a", "b"])
        .apply(lambda df: df.iloc[-1])
    )
    expected = DataFrame(
        [[pd.to_datetime(2, unit="s")], [pd.to_datetime(4, unit="s")]],
        columns=["c"],
        index=MultiIndex.from_tuples([(1, ""), (2, "")], names=["a", "b"]),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("index_values", [[1, 2, 3], [1.0, 2.0, 3.0]])
def test_apply_index_key_error_bug(index_values):
    # GH 44310
    result = DataFrame(
        {
            "a": ["aa", "a2", "a3"],
            "b": [1, 2, 3],
        },
        index=Index(index_values),
    )
    expected = DataFrame(
        {
            "b_mean": [2.0, 3.0, 1.0],
        },
        index=Index(["a2", "a3", "aa"], name="a"),
    )
    result = result.groupby("a").apply(
        lambda df: Series([df["b"].mean()], index=["b_mean"])
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "arg,idx",
    [
        [
            [
                1,
                2,
                3,
            ],
            [
                0.1,
                0.3,
                0.2,
            ],
        ],
        [
            [
                1,
                2,
                3,
            ],
            [
                0.1,
                0.2,
                0.3,
            ],
        ],
        [
            [
                1,
                4,
                3,
            ],
            [
                0.1,
                0.4,
                0.2,
            ],
        ],
    ],
)
def test_apply_nonmonotonic_float_index(arg, idx):
    # GH 34455
    df = DataFrame({"grp": arg, "col": arg}, index=idx)
    result = df.groupby("grp", group_keys=False).apply(lambda x: x)
    expected = df[["col"]]
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("args, kwargs", [([True], {}), ([], {"numeric_only": True})])
def test_apply_str_with_args(df, args, kwargs):
    # GH#46479
    gb = df.groupby("A")
    result = gb.apply("sum", *args, **kwargs)
    expected = gb.sum(numeric_only=True)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("name", ["some_name", None])
def test_result_name_when_one_group(name):
    # GH 46369
    ser = Series([1, 2], name=name)
    result = ser.groupby(["a", "a"], group_keys=False).apply(lambda x: x)
    expected = Series([1, 2], name=name)

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, op",
    [
        ("apply", lambda gb: gb.values[-1]),
        ("apply", lambda gb: gb["b"].iloc[0]),
        ("agg", "skew"),
        ("agg", "kurt"),
        ("agg", "prod"),
        ("agg", "sum"),
    ],
)
def test_empty_df(method, op):
    # GH 47985
    empty_df = DataFrame({"a": [], "b": []})
    gb = empty_df.groupby("a", group_keys=True)
    group = gb.b

    result = getattr(group, method)(op)
    expected = Series(
        [], name="b", dtype="float64", index=Index([], dtype="float64", name="a")
    )

    tm.assert_series_equal(result, expected)


def test_include_groups():
    # GH#7155
    df = DataFrame({"a": [1, 1, 2], "b": [3, 4, 5]})
    gb = df.groupby("a")
    with pytest.raises(ValueError, match="include_groups=True is no longer allowed"):
        gb.apply(lambda x: x.sum(), include_groups=True)


@pytest.mark.parametrize("func, value", [(max, 2), (min, 1), (sum, 3)])
def test_builtins_apply(func, value):
    # GH#8155, GH#53974
    # Builtins act as e.g. sum(group), which sums the column labels of group
    df = DataFrame({0: [1, 1, 2], 1: [3, 4, 5], 2: [3, 4, 5]})
    gb = df.groupby(0)
    result = gb.apply(func)

    expected = Series([value, value], index=Index([1, 2], name=0))
    tm.assert_series_equal(result, expected)


def test_inconsistent_return_type():
    # GH5592
    # inconsistent return type
    df = DataFrame(
        {
            "A": ["Tiger", "Tiger", "Tiger", "Lamb", "Lamb", "Pony", "Pony"],
            "B": Series(np.arange(7), dtype="int64"),
            "C": pd.date_range("20130101", periods=7),
        }
    )

    def f_0(grp):
        return grp.iloc[0]

    expected = df.groupby("A").first()[["B"]]
    result = df.groupby("A").apply(f_0)[["B"]]
    tm.assert_frame_equal(result, expected)

    def f_1(grp):
        if grp.name == "Tiger":
            return None
        return grp.iloc[0]

    result = df.groupby("A").apply(f_1)[["B"]]
    e = expected.copy()
    e.loc["Tiger"] = np.nan
    tm.assert_frame_equal(result, e)

    def f_2(grp):
        if grp.name == "Pony":
            return None
        return grp.iloc[0]

    result = df.groupby("A").apply(f_2)[["B"]]
    e = expected.copy()
    e.loc["Pony"] = np.nan
    tm.assert_frame_equal(result, e)

    # 5592 revisited, with datetimes
    def f_3(grp):
        if grp.name == "Pony":
            return None
        return grp.iloc[0]

    result = df.groupby("A").apply(f_3)[["C"]]
    e = df.groupby("A").first()[["C"]]
    e.loc["Pony"] = pd.NaT
    tm.assert_frame_equal(result, e)

    # scalar outputs
    def f_4(grp):
        if grp.name == "Pony":
            return None
        return grp.iloc[0].loc["C"]

    result = df.groupby("A").apply(f_4)
    e = df.groupby("A").first()["C"].copy()
    e.loc["Pony"] = np.nan
    e.name = None
    tm.assert_series_equal(result, e)


def test_nonreducer_nonstransform():
    # GH3380, GH60619
    # Was originally testing mutating in a UDF; now kept as an example
    # of using apply with a nonreducer and nontransformer.
    df = DataFrame(
        {
            "cat1": ["a"] * 8 + ["b"] * 6,
            "cat2": ["c"] * 2
            + ["d"] * 2
            + ["e"] * 2
            + ["f"] * 2
            + ["c"] * 2
            + ["d"] * 2
            + ["e"] * 2,
            "val": np.random.default_rng(2).integers(100, size=14),
        }
    )

    def f(x):
        x = x.copy()
        x["rank"] = x.val.rank(method="min")
        return x.groupby("cat2")["rank"].min()

    expected = DataFrame(
        {
            "cat1": list("aaaabbb"),
            "cat2": list("cdefcde"),
            "rank": [3.0, 2.0, 5.0, 1.0, 2.0, 4.0, 1.0],
        }
    ).set_index(["cat1", "cat2"])["rank"]
    result = df.groupby("cat1").apply(f)
    tm.assert_series_equal(result, expected)
