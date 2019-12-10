""" test with the .transform """
from io import StringIO

import numpy as np
import pytest

from pandas._libs import groupby

from pandas.core.dtypes.common import ensure_platform_int, is_timedelta64_dtype

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    MultiIndex,
    Series,
    Timestamp,
    concat,
    date_range,
)
from pandas.core.groupby.groupby import DataError
import pandas.util.testing as tm


def assert_fp_equal(a, b):
    assert (np.abs(a - b) < 1e-12).all()


def test_transform():
    data = Series(np.arange(9) // 3, index=np.arange(9))

    index = np.arange(9)
    np.random.shuffle(index)
    data = data.reindex(index)

    grouped = data.groupby(lambda x: x // 3)

    transformed = grouped.transform(lambda x: x * x.sum())
    assert transformed[7] == 12

    # GH 8046
    # make sure that we preserve the input order

    df = DataFrame(
        np.arange(6, dtype="int64").reshape(3, 2), columns=["a", "b"], index=[0, 2, 1]
    )
    key = [0, 0, 1]
    expected = (
        df.sort_index()
        .groupby(key)
        .transform(lambda x: x - x.mean())
        .groupby(key)
        .mean()
    )
    result = df.groupby(key).transform(lambda x: x - x.mean()).groupby(key).mean()
    tm.assert_frame_equal(result, expected)

    def demean(arr):
        return arr - arr.mean()

    people = DataFrame(
        np.random.randn(5, 5),
        columns=["a", "b", "c", "d", "e"],
        index=["Joe", "Steve", "Wes", "Jim", "Travis"],
    )
    key = ["one", "two", "one", "two", "one"]
    result = people.groupby(key).transform(demean).groupby(key).mean()
    expected = people.groupby(key).apply(demean).groupby(key).mean()
    tm.assert_frame_equal(result, expected)

    # GH 8430
    df = tm.makeTimeDataFrame()
    g = df.groupby(pd.Grouper(freq="M"))
    g.transform(lambda x: x - 1)

    # GH 9700
    df = DataFrame({"a": range(5, 10), "b": range(5)})
    result = df.groupby("a").transform(max)
    expected = DataFrame({"b": range(5)})
    tm.assert_frame_equal(result, expected)


def test_transform_fast():

    df = DataFrame({"id": np.arange(100000) / 3, "val": np.random.randn(100000)})

    grp = df.groupby("id")["val"]

    values = np.repeat(grp.mean().values, ensure_platform_int(grp.count().values))
    expected = pd.Series(values, index=df.index, name="val")

    result = grp.transform(np.mean)
    tm.assert_series_equal(result, expected)

    result = grp.transform("mean")
    tm.assert_series_equal(result, expected)

    # GH 12737
    df = pd.DataFrame(
        {
            "grouping": [0, 1, 1, 3],
            "f": [1.1, 2.1, 3.1, 4.5],
            "d": pd.date_range("2014-1-1", "2014-1-4"),
            "i": [1, 2, 3, 4],
        },
        columns=["grouping", "f", "i", "d"],
    )
    result = df.groupby("grouping").transform("first")

    dates = [
        pd.Timestamp("2014-1-1"),
        pd.Timestamp("2014-1-2"),
        pd.Timestamp("2014-1-2"),
        pd.Timestamp("2014-1-4"),
    ]
    expected = pd.DataFrame(
        {"f": [1.1, 2.1, 2.1, 4.5], "d": dates, "i": [1, 2, 2, 4]},
        columns=["f", "i", "d"],
    )
    tm.assert_frame_equal(result, expected)

    # selection
    result = df.groupby("grouping")[["f", "i"]].transform("first")
    expected = expected[["f", "i"]]
    tm.assert_frame_equal(result, expected)

    # dup columns
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=["g", "a", "a"])
    result = df.groupby("g").transform("first")
    expected = df.drop("g", axis=1)
    tm.assert_frame_equal(result, expected)


def test_transform_broadcast(tsframe, ts):
    grouped = ts.groupby(lambda x: x.month)
    result = grouped.transform(np.mean)

    tm.assert_index_equal(result.index, ts.index)
    for _, gp in grouped:
        assert_fp_equal(result.reindex(gp.index), gp.mean())

    grouped = tsframe.groupby(lambda x: x.month)
    result = grouped.transform(np.mean)
    tm.assert_index_equal(result.index, tsframe.index)
    for _, gp in grouped:
        agged = gp.mean()
        res = result.reindex(gp.index)
        for col in tsframe:
            assert_fp_equal(res[col], agged[col])

    # group columns
    grouped = tsframe.groupby({"A": 0, "B": 0, "C": 1, "D": 1}, axis=1)
    result = grouped.transform(np.mean)
    tm.assert_index_equal(result.index, tsframe.index)
    tm.assert_index_equal(result.columns, tsframe.columns)
    for _, gp in grouped:
        agged = gp.mean(1)
        res = result.reindex(columns=gp.columns)
        for idx in gp.index:
            assert_fp_equal(res.xs(idx), agged[idx])


def test_transform_axis(tsframe):

    # make sure that we are setting the axes
    # correctly when on axis=0 or 1
    # in the presence of a non-monotonic indexer
    # GH12713

    base = tsframe.iloc[0:5]
    r = len(base.index)
    c = len(base.columns)
    tso = DataFrame(
        np.random.randn(r, c), index=base.index, columns=base.columns, dtype="float64"
    )
    # monotonic
    ts = tso
    grouped = ts.groupby(lambda x: x.weekday())
    result = ts - grouped.transform("mean")
    expected = grouped.apply(lambda x: x - x.mean())
    tm.assert_frame_equal(result, expected)

    ts = ts.T
    grouped = ts.groupby(lambda x: x.weekday(), axis=1)
    result = ts - grouped.transform("mean")
    expected = grouped.apply(lambda x: (x.T - x.mean(1)).T)
    tm.assert_frame_equal(result, expected)

    # non-monotonic
    ts = tso.iloc[[1, 0] + list(range(2, len(base)))]
    grouped = ts.groupby(lambda x: x.weekday())
    result = ts - grouped.transform("mean")
    expected = grouped.apply(lambda x: x - x.mean())
    tm.assert_frame_equal(result, expected)

    ts = ts.T
    grouped = ts.groupby(lambda x: x.weekday(), axis=1)
    result = ts - grouped.transform("mean")
    expected = grouped.apply(lambda x: (x.T - x.mean(1)).T)
    tm.assert_frame_equal(result, expected)


def test_transform_dtype():
    # GH 9807
    # Check transform dtype output is preserved
    df = DataFrame([[1, 3], [2, 3]])
    result = df.groupby(1).transform("mean")
    expected = DataFrame([[1.5], [1.5]])
    tm.assert_frame_equal(result, expected)


def test_transform_bug():
    # GH 5712
    # transforming on a datetime column
    df = DataFrame(dict(A=Timestamp("20130101"), B=np.arange(5)))
    result = df.groupby("A")["B"].transform(lambda x: x.rank(ascending=False))
    expected = Series(np.arange(5, 0, step=-1), name="B")
    tm.assert_series_equal(result, expected)


def test_transform_numeric_to_boolean():
    # GH 16875
    # inconsistency in transforming boolean values
    expected = pd.Series([True, True], name="A")

    df = pd.DataFrame({"A": [1.1, 2.2], "B": [1, 2]})
    result = df.groupby("B").A.transform(lambda x: True)
    tm.assert_series_equal(result, expected)

    df = pd.DataFrame({"A": [1, 2], "B": [1, 2]})
    result = df.groupby("B").A.transform(lambda x: True)
    tm.assert_series_equal(result, expected)


def test_transform_datetime_to_timedelta():
    # GH 15429
    # transforming a datetime to timedelta
    df = DataFrame(dict(A=Timestamp("20130101"), B=np.arange(5)))
    expected = pd.Series([Timestamp("20130101") - Timestamp("20130101")] * 5, name="A")

    # this does date math without changing result type in transform
    base_time = df["A"][0]
    result = (
        df.groupby("A")["A"].transform(lambda x: x.max() - x.min() + base_time)
        - base_time
    )
    tm.assert_series_equal(result, expected)

    # this does date math and causes the transform to return timedelta
    result = df.groupby("A")["A"].transform(lambda x: x.max() - x.min())
    tm.assert_series_equal(result, expected)


def test_transform_datetime_to_numeric():
    # GH 10972
    # convert dt to float
    df = DataFrame({"a": 1, "b": date_range("2015-01-01", periods=2, freq="D")})
    result = df.groupby("a").b.transform(
        lambda x: x.dt.dayofweek - x.dt.dayofweek.mean()
    )

    expected = Series([-0.5, 0.5], name="b")
    tm.assert_series_equal(result, expected)

    # convert dt to int
    df = DataFrame({"a": 1, "b": date_range("2015-01-01", periods=2, freq="D")})
    result = df.groupby("a").b.transform(
        lambda x: x.dt.dayofweek - x.dt.dayofweek.min()
    )

    expected = Series([0, 1], name="b")
    tm.assert_series_equal(result, expected)


def test_transform_casting():
    # 13046
    data = """
    idx     A         ID3              DATETIME
    0   B-028  b76cd912ff "2014-10-08 13:43:27"
    1   B-054  4a57ed0b02 "2014-10-08 14:26:19"
    2   B-076  1a682034f8 "2014-10-08 14:29:01"
    3   B-023  b76cd912ff "2014-10-08 18:39:34"
    4   B-023  f88g8d7sds "2014-10-08 18:40:18"
    5   B-033  b76cd912ff "2014-10-08 18:44:30"
    6   B-032  b76cd912ff "2014-10-08 18:46:00"
    7   B-037  b76cd912ff "2014-10-08 18:52:15"
    8   B-046  db959faf02 "2014-10-08 18:59:59"
    9   B-053  b76cd912ff "2014-10-08 19:17:48"
    10  B-065  b76cd912ff "2014-10-08 19:21:38"
    """
    df = pd.read_csv(
        StringIO(data), sep=r"\s+", index_col=[0], parse_dates=["DATETIME"]
    )

    result = df.groupby("ID3")["DATETIME"].transform(lambda x: x.diff())
    assert is_timedelta64_dtype(result.dtype)

    result = df[["ID3", "DATETIME"]].groupby("ID3").transform(lambda x: x.diff())
    assert is_timedelta64_dtype(result.DATETIME.dtype)


def test_transform_multiple(ts):
    grouped = ts.groupby([lambda x: x.year, lambda x: x.month])

    grouped.transform(lambda x: x * 2)
    grouped.transform(np.mean)


def test_dispatch_transform(tsframe):
    df = tsframe[::5].reindex(tsframe.index)

    grouped = df.groupby(lambda x: x.month)

    filled = grouped.fillna(method="pad")
    fillit = lambda x: x.fillna(method="pad")
    expected = df.groupby(lambda x: x.month).transform(fillit)
    tm.assert_frame_equal(filled, expected)


def test_transform_select_columns(df):
    f = lambda x: x.mean()
    result = df.groupby("A")["C", "D"].transform(f)

    selection = df[["C", "D"]]
    expected = selection.groupby(df["A"]).transform(f)

    tm.assert_frame_equal(result, expected)


def test_transform_exclude_nuisance(df):

    # this also tests orderings in transform between
    # series/frame to make sure it's consistent
    expected = {}
    grouped = df.groupby("A")
    expected["C"] = grouped["C"].transform(np.mean)
    expected["D"] = grouped["D"].transform(np.mean)
    expected = DataFrame(expected)
    result = df.groupby("A").transform(np.mean)

    tm.assert_frame_equal(result, expected)


def test_transform_function_aliases(df):
    result = df.groupby("A").transform("mean")
    expected = df.groupby("A").transform(np.mean)
    tm.assert_frame_equal(result, expected)

    result = df.groupby("A")["C"].transform("mean")
    expected = df.groupby("A")["C"].transform(np.mean)
    tm.assert_series_equal(result, expected)


def test_series_fast_transform_date():
    # GH 13191
    df = pd.DataFrame(
        {"grouping": [np.nan, 1, 1, 3], "d": pd.date_range("2014-1-1", "2014-1-4")}
    )
    result = df.groupby("grouping")["d"].transform("first")
    dates = [
        pd.NaT,
        pd.Timestamp("2014-1-2"),
        pd.Timestamp("2014-1-2"),
        pd.Timestamp("2014-1-4"),
    ]
    expected = pd.Series(dates, name="d")
    tm.assert_series_equal(result, expected)


def test_transform_length():
    # GH 9697
    df = pd.DataFrame({"col1": [1, 1, 2, 2], "col2": [1, 2, 3, np.nan]})
    expected = pd.Series([3.0] * 4)

    def nsum(x):
        return np.nansum(x)

    results = [
        df.groupby("col1").transform(sum)["col2"],
        df.groupby("col1")["col2"].transform(sum),
        df.groupby("col1").transform(nsum)["col2"],
        df.groupby("col1")["col2"].transform(nsum),
    ]
    for result in results:
        tm.assert_series_equal(result, expected, check_names=False)


def test_transform_coercion():

    # 14457
    # when we are transforming be sure to not coerce
    # via assignment
    df = pd.DataFrame(dict(A=["a", "a"], B=[0, 1]))
    g = df.groupby("A")

    expected = g.transform(np.mean)
    result = g.transform(lambda x: np.mean(x))
    tm.assert_frame_equal(result, expected)


def test_groupby_transform_with_int():

    # GH 3740, make sure that we might upcast on item-by-item transform

    # floats
    df = DataFrame(
        dict(
            A=[1, 1, 1, 2, 2, 2],
            B=Series(1, dtype="float64"),
            C=Series([1, 2, 3, 1, 2, 3], dtype="float64"),
            D="foo",
        )
    )
    with np.errstate(all="ignore"):
        result = df.groupby("A").transform(lambda x: (x - x.mean()) / x.std())
    expected = DataFrame(
        dict(B=np.nan, C=Series([-1, 0, 1, -1, 0, 1], dtype="float64"))
    )
    tm.assert_frame_equal(result, expected)

    # int case
    df = DataFrame(dict(A=[1, 1, 1, 2, 2, 2], B=1, C=[1, 2, 3, 1, 2, 3], D="foo"))
    with np.errstate(all="ignore"):
        result = df.groupby("A").transform(lambda x: (x - x.mean()) / x.std())
    expected = DataFrame(dict(B=np.nan, C=[-1, 0, 1, -1, 0, 1]))
    tm.assert_frame_equal(result, expected)

    # int that needs float conversion
    s = Series([2, 3, 4, 10, 5, -1])
    df = DataFrame(dict(A=[1, 1, 1, 2, 2, 2], B=1, C=s, D="foo"))
    with np.errstate(all="ignore"):
        result = df.groupby("A").transform(lambda x: (x - x.mean()) / x.std())

    s1 = s.iloc[0:3]
    s1 = (s1 - s1.mean()) / s1.std()
    s2 = s.iloc[3:6]
    s2 = (s2 - s2.mean()) / s2.std()
    expected = DataFrame(dict(B=np.nan, C=concat([s1, s2])))
    tm.assert_frame_equal(result, expected)

    # int downcasting
    result = df.groupby("A").transform(lambda x: x * 2 / 2)
    expected = DataFrame(dict(B=1, C=[2, 3, 4, 10, 5, -1]))
    tm.assert_frame_equal(result, expected)


def test_groupby_transform_with_nan_group():
    # GH 9941
    df = pd.DataFrame({"a": range(10), "b": [1, 1, 2, 3, np.nan, 4, 4, 5, 5, 5]})
    result = df.groupby(df.b)["a"].transform(max)
    expected = pd.Series(
        [1.0, 1.0, 2.0, 3.0, np.nan, 6.0, 6.0, 9.0, 9.0, 9.0], name="a"
    )
    tm.assert_series_equal(result, expected)


def test_transform_mixed_type():
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
        group["g"] = group["d"] * 2
        return group[:1]

    grouped = df.groupby("c")
    result = grouped.apply(f)

    assert result["d"].dtype == np.float64

    # this is by definition a mutating operation!
    with pd.option_context("mode.chained_assignment", None):
        for key, group in grouped:
            res = f(group)
            tm.assert_frame_equal(res, result.loc[key])


def _check_cython_group_transform_cumulative(pd_op, np_op, dtype):
    """
    Check a group transform that executes a cumulative function.

    Parameters
    ----------
    pd_op : callable
        The pandas cumulative function.
    np_op : callable
        The analogous one in NumPy.
    dtype : type
        The specified dtype of the data.
    """

    is_datetimelike = False

    data = np.array([[1], [2], [3], [4]], dtype=dtype)
    ans = np.zeros_like(data)

    labels = np.array([0, 0, 0, 0], dtype=np.int64)
    ngroups = 1
    pd_op(ans, data, labels, ngroups, is_datetimelike)

    tm.assert_numpy_array_equal(np_op(data), ans[:, 0], check_dtype=False)


def test_cython_group_transform_cumsum(any_real_dtype):
    # see gh-4095
    dtype = np.dtype(any_real_dtype).type
    pd_op, np_op = groupby.group_cumsum, np.cumsum
    _check_cython_group_transform_cumulative(pd_op, np_op, dtype)


def test_cython_group_transform_cumprod():
    # see gh-4095
    dtype = np.float64
    pd_op, np_op = groupby.group_cumprod_float64, np.cumproduct
    _check_cython_group_transform_cumulative(pd_op, np_op, dtype)


def test_cython_group_transform_algos():
    # see gh-4095
    is_datetimelike = False

    # with nans
    labels = np.array([0, 0, 0, 0, 0], dtype=np.int64)
    ngroups = 1

    data = np.array([[1], [2], [3], [np.nan], [4]], dtype="float64")
    actual = np.zeros_like(data)
    actual.fill(np.nan)
    groupby.group_cumprod_float64(actual, data, labels, ngroups, is_datetimelike)
    expected = np.array([1, 2, 6, np.nan, 24], dtype="float64")
    tm.assert_numpy_array_equal(actual[:, 0], expected)

    actual = np.zeros_like(data)
    actual.fill(np.nan)
    groupby.group_cumsum(actual, data, labels, ngroups, is_datetimelike)
    expected = np.array([1, 3, 6, np.nan, 10], dtype="float64")
    tm.assert_numpy_array_equal(actual[:, 0], expected)

    # timedelta
    is_datetimelike = True
    data = np.array([np.timedelta64(1, "ns")] * 5, dtype="m8[ns]")[:, None]
    actual = np.zeros_like(data, dtype="int64")
    groupby.group_cumsum(actual, data.view("int64"), labels, ngroups, is_datetimelike)
    expected = np.array(
        [
            np.timedelta64(1, "ns"),
            np.timedelta64(2, "ns"),
            np.timedelta64(3, "ns"),
            np.timedelta64(4, "ns"),
            np.timedelta64(5, "ns"),
        ]
    )
    tm.assert_numpy_array_equal(actual[:, 0].view("m8[ns]"), expected)


@pytest.mark.parametrize(
    "op, args, targop",
    [
        ("cumprod", (), lambda x: x.cumprod()),
        ("cumsum", (), lambda x: x.cumsum()),
        ("shift", (-1,), lambda x: x.shift(-1)),
        ("shift", (1,), lambda x: x.shift()),
    ],
)
def test_cython_transform_series(op, args, targop):
    # GH 4095
    s = Series(np.random.randn(1000))
    s_missing = s.copy()
    s_missing.iloc[2:10] = np.nan
    labels = np.random.randint(0, 50, size=1000).astype(float)

    # series
    for data in [s, s_missing]:
        # print(data.head())
        expected = data.groupby(labels).transform(targop)

        tm.assert_series_equal(expected, data.groupby(labels).transform(op, *args))
        tm.assert_series_equal(expected, getattr(data.groupby(labels), op)(*args))


@pytest.mark.parametrize("op", ["cumprod", "cumsum"])
@pytest.mark.parametrize("skipna", [False, True])
@pytest.mark.parametrize(
    "input, exp",
    [
        # When everything is NaN
        ({"key": ["b"] * 10, "value": np.nan}, pd.Series([np.nan] * 10, name="value")),
        # When there is a single NaN
        (
            {"key": ["b"] * 10 + ["a"] * 2, "value": [3] * 3 + [np.nan] + [3] * 8},
            {
                ("cumprod", False): [3.0, 9.0, 27.0] + [np.nan] * 7 + [3.0, 9.0],
                ("cumprod", True): [
                    3.0,
                    9.0,
                    27.0,
                    np.nan,
                    81.0,
                    243.0,
                    729.0,
                    2187.0,
                    6561.0,
                    19683.0,
                    3.0,
                    9.0,
                ],
                ("cumsum", False): [3.0, 6.0, 9.0] + [np.nan] * 7 + [3.0, 6.0],
                ("cumsum", True): [
                    3.0,
                    6.0,
                    9.0,
                    np.nan,
                    12.0,
                    15.0,
                    18.0,
                    21.0,
                    24.0,
                    27.0,
                    3.0,
                    6.0,
                ],
            },
        ),
    ],
)
def test_groupby_cum_skipna(op, skipna, input, exp):
    df = pd.DataFrame(input)
    result = df.groupby("key")["value"].transform(op, skipna=skipna)
    if isinstance(exp, dict):
        expected = exp[(op, skipna)]
    else:
        expected = exp
    expected = pd.Series(expected, name="value")
    tm.assert_series_equal(expected, result)


@pytest.mark.parametrize(
    "op, args, targop",
    [
        ("cumprod", (), lambda x: x.cumprod()),
        ("cumsum", (), lambda x: x.cumsum()),
        ("shift", (-1,), lambda x: x.shift(-1)),
        ("shift", (1,), lambda x: x.shift()),
    ],
)
def test_cython_transform_frame(op, args, targop):
    s = Series(np.random.randn(1000))
    s_missing = s.copy()
    s_missing.iloc[2:10] = np.nan
    labels = np.random.randint(0, 50, size=1000).astype(float)
    strings = list("qwertyuiopasdfghjklz")
    strings_missing = strings[:]
    strings_missing[5] = np.nan
    df = DataFrame(
        {
            "float": s,
            "float_missing": s_missing,
            "int": [1, 1, 1, 1, 2] * 200,
            "datetime": pd.date_range("1990-1-1", periods=1000),
            "timedelta": pd.timedelta_range(1, freq="s", periods=1000),
            "string": strings * 50,
            "string_missing": strings_missing * 50,
        },
        columns=[
            "float",
            "float_missing",
            "int",
            "datetime",
            "timedelta",
            "string",
            "string_missing",
        ],
    )
    df["cat"] = df["string"].astype("category")

    df2 = df.copy()
    df2.index = pd.MultiIndex.from_product([range(100), range(10)])

    # DataFrame - Single and MultiIndex,
    # group by values, index level, columns
    for df in [df, df2]:
        for gb_target in [
            dict(by=labels),
            dict(level=0),
            dict(by="string"),
        ]:  # dict(by='string_missing')]:
            # dict(by=['int','string'])]:

            gb = df.groupby(**gb_target)
            # whitelisted methods set the selection before applying
            # bit a of hack to make sure the cythonized shift
            # is equivalent to pre 0.17.1 behavior
            if op == "shift":
                gb._set_group_selection()

            if op != "shift" and "int" not in gb_target:
                # numeric apply fastpath promotes dtype so have
                # to apply separately and concat
                i = gb[["int"]].apply(targop)
                f = gb[["float", "float_missing"]].apply(targop)
                expected = pd.concat([f, i], axis=1)
            else:
                expected = gb.apply(targop)

            expected = expected.sort_index(axis=1)
            tm.assert_frame_equal(expected, gb.transform(op, *args).sort_index(axis=1))
            tm.assert_frame_equal(expected, getattr(gb, op)(*args).sort_index(axis=1))
            # individual columns
            for c in df:
                if c not in ["float", "int", "float_missing"] and op != "shift":
                    msg = "No numeric types to aggregate"
                    with pytest.raises(DataError, match=msg):
                        gb[c].transform(op)
                    with pytest.raises(DataError, match=msg):
                        getattr(gb[c], op)()
                else:
                    expected = gb[c].apply(targop)
                    expected.name = c
                    tm.assert_series_equal(expected, gb[c].transform(op, *args))
                    tm.assert_series_equal(expected, getattr(gb[c], op)(*args))


def test_transform_with_non_scalar_group():
    # GH 10165
    cols = pd.MultiIndex.from_tuples(
        [
            ("syn", "A"),
            ("mis", "A"),
            ("non", "A"),
            ("syn", "C"),
            ("mis", "C"),
            ("non", "C"),
            ("syn", "T"),
            ("mis", "T"),
            ("non", "T"),
            ("syn", "G"),
            ("mis", "G"),
            ("non", "G"),
        ]
    )
    df = pd.DataFrame(
        np.random.randint(1, 10, (4, 12)), columns=cols, index=["A", "C", "G", "T"]
    )

    msg = "transform must return a scalar value for each group.*"
    with pytest.raises(ValueError, match=msg):
        df.groupby(axis=1, level=1).transform(lambda z: z.div(z.sum(axis=1), axis=0))


@pytest.mark.parametrize(
    "cols,exp,comp_func",
    [
        ("a", pd.Series([1, 1, 1], name="a"), tm.assert_series_equal),
        (
            ["a", "c"],
            pd.DataFrame({"a": [1, 1, 1], "c": [1, 1, 1]}),
            tm.assert_frame_equal,
        ),
    ],
)
@pytest.mark.parametrize("agg_func", ["count", "rank", "size"])
def test_transform_numeric_ret(cols, exp, comp_func, agg_func):
    if agg_func == "size" and isinstance(cols, list):
        pytest.xfail("'size' transformation not supported with NDFrameGroupy")

    # GH 19200
    df = pd.DataFrame(
        {"a": pd.date_range("2018-01-01", periods=3), "b": range(3), "c": range(7, 10)}
    )

    result = df.groupby("b")[cols].transform(agg_func)

    if agg_func == "rank":
        exp = exp.astype("float")

    comp_func(result, exp)


@pytest.mark.parametrize("mix_groupings", [True, False])
@pytest.mark.parametrize("as_series", [True, False])
@pytest.mark.parametrize("val1,val2", [("foo", "bar"), (1, 2), (1.0, 2.0)])
@pytest.mark.parametrize(
    "fill_method,limit,exp_vals",
    [
        (
            "ffill",
            None,
            [np.nan, np.nan, "val1", "val1", "val1", "val2", "val2", "val2"],
        ),
        ("ffill", 1, [np.nan, np.nan, "val1", "val1", np.nan, "val2", "val2", np.nan]),
        (
            "bfill",
            None,
            ["val1", "val1", "val1", "val2", "val2", "val2", np.nan, np.nan],
        ),
        ("bfill", 1, [np.nan, "val1", "val1", np.nan, "val2", "val2", np.nan, np.nan]),
    ],
)
def test_group_fill_methods(
    mix_groupings, as_series, val1, val2, fill_method, limit, exp_vals
):
    vals = [np.nan, np.nan, val1, np.nan, np.nan, val2, np.nan, np.nan]
    _exp_vals = list(exp_vals)
    # Overwrite placeholder values
    for index, exp_val in enumerate(_exp_vals):
        if exp_val == "val1":
            _exp_vals[index] = val1
        elif exp_val == "val2":
            _exp_vals[index] = val2

    # Need to modify values and expectations depending on the
    # Series / DataFrame that we ultimately want to generate
    if mix_groupings:  # ['a', 'b', 'a, 'b', ...]
        keys = ["a", "b"] * len(vals)

        def interweave(list_obj):
            temp = list()
            for x in list_obj:
                temp.extend([x, x])

            return temp

        _exp_vals = interweave(_exp_vals)
        vals = interweave(vals)
    else:  # ['a', 'a', 'a', ... 'b', 'b', 'b']
        keys = ["a"] * len(vals) + ["b"] * len(vals)
        _exp_vals = _exp_vals * 2
        vals = vals * 2

    df = DataFrame({"key": keys, "val": vals})
    if as_series:
        result = getattr(df.groupby("key")["val"], fill_method)(limit=limit)
        exp = Series(_exp_vals, name="val")
        tm.assert_series_equal(result, exp)
    else:
        result = getattr(df.groupby("key"), fill_method)(limit=limit)
        exp = DataFrame({"val": _exp_vals})
        tm.assert_frame_equal(result, exp)


@pytest.mark.parametrize("fill_method", ["ffill", "bfill"])
def test_pad_stable_sorting(fill_method):
    # GH 21207
    x = [0] * 20
    y = [np.nan] * 10 + [1] * 10

    if fill_method == "bfill":
        y = y[::-1]

    df = pd.DataFrame({"x": x, "y": y})
    expected = df.drop("x", 1)

    result = getattr(df.groupby("x"), fill_method)()

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("test_series", [True, False])
@pytest.mark.parametrize(
    "freq",
    [
        None,
        pytest.param(
            "D",
            marks=pytest.mark.xfail(
                reason="GH#23918 before method uses freq in vectorized approach"
            ),
        ),
    ],
)
@pytest.mark.parametrize(
    "periods,fill_method,limit",
    [
        (1, "ffill", None),
        (1, "ffill", 1),
        (1, "bfill", None),
        (1, "bfill", 1),
        (-1, "ffill", None),
        (-1, "ffill", 1),
        (-1, "bfill", None),
        (-1, "bfill", 1),
    ],
)
def test_pct_change(test_series, freq, periods, fill_method, limit):
    # GH  21200, 21621
    vals = [3, np.nan, np.nan, np.nan, 1, 2, 4, 10, np.nan, 4]
    keys = ["a", "b"]
    key_v = np.repeat(keys, len(vals))
    df = DataFrame({"key": key_v, "vals": vals * 2})

    df_g = getattr(df.groupby("key"), fill_method)(limit=limit)
    grp = df_g.groupby(df.key)

    expected = grp["vals"].obj / grp["vals"].shift(periods) - 1

    if test_series:
        result = df.groupby("key")["vals"].pct_change(
            periods=periods, fill_method=fill_method, limit=limit, freq=freq
        )
        tm.assert_series_equal(result, expected)
    else:
        result = df.groupby("key").pct_change(
            periods=periods, fill_method=fill_method, limit=limit, freq=freq
        )
        tm.assert_frame_equal(result, expected.to_frame("vals"))


@pytest.mark.parametrize(
    "func, expected_status",
    [
        ("ffill", ["shrt", "shrt", "lng", np.nan, "shrt", "ntrl", "ntrl"]),
        ("bfill", ["shrt", "lng", "lng", "shrt", "shrt", "ntrl", np.nan]),
    ],
)
def test_ffill_bfill_non_unique_multilevel(func, expected_status):
    # GH 19437
    date = pd.to_datetime(
        [
            "2018-01-01",
            "2018-01-01",
            "2018-01-01",
            "2018-01-01",
            "2018-01-02",
            "2018-01-01",
            "2018-01-02",
        ]
    )
    symbol = ["MSFT", "MSFT", "MSFT", "AAPL", "AAPL", "TSLA", "TSLA"]
    status = ["shrt", np.nan, "lng", np.nan, "shrt", "ntrl", np.nan]

    df = DataFrame({"date": date, "symbol": symbol, "status": status})
    df = df.set_index(["date", "symbol"])
    result = getattr(df.groupby("symbol")["status"], func)()

    index = MultiIndex.from_tuples(
        tuples=list(zip(*[date, symbol])), names=["date", "symbol"]
    )
    expected = Series(expected_status, index=index, name="status")

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("func", [np.any, np.all])
def test_any_all_np_func(func):
    # GH 20653
    df = pd.DataFrame(
        [["foo", True], [np.nan, True], ["foo", True]], columns=["key", "val"]
    )

    exp = pd.Series([True, np.nan, True], name="val")

    res = df.groupby("key")["val"].transform(func)
    tm.assert_series_equal(res, exp)


def test_groupby_transform_rename():
    # https://github.com/pandas-dev/pandas/issues/23461
    def demean_rename(x):
        result = x - x.mean()

        if isinstance(x, pd.Series):
            return result

        result = result.rename(
            columns={c: "{}_demeaned".format(c) for c in result.columns}
        )

        return result

    df = pd.DataFrame({"group": list("ababa"), "value": [1, 1, 1, 2, 2]})
    expected = pd.DataFrame({"value": [-1.0 / 3, -0.5, -1.0 / 3, 0.5, 2.0 / 3]})

    result = df.groupby("group").transform(demean_rename)
    tm.assert_frame_equal(result, expected)
    result_single = df.groupby("group").value.transform(demean_rename)
    tm.assert_series_equal(result_single, expected["value"])


@pytest.mark.parametrize("func", [min, max, np.min, np.max, "first", "last"])
def test_groupby_transform_timezone_column(func):
    # GH 24198
    ts = pd.to_datetime("now", utc=True).tz_convert("Asia/Singapore")
    result = pd.DataFrame({"end_time": [ts], "id": [1]})
    result["max_end_time"] = result.groupby("id").end_time.transform(func)
    expected = pd.DataFrame([[ts, 1, ts]], columns=["end_time", "id", "max_end_time"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "func, values",
    [
        ("idxmin", ["1/1/2011"] * 2 + ["1/3/2011"] * 7 + ["1/10/2011"]),
        ("idxmax", ["1/2/2011"] * 2 + ["1/9/2011"] * 7 + ["1/10/2011"]),
    ],
)
def test_groupby_transform_with_datetimes(func, values):
    # GH 15306
    dates = pd.date_range("1/1/2011", periods=10, freq="D")

    stocks = pd.DataFrame({"price": np.arange(10.0)}, index=dates)
    stocks["week_id"] = pd.to_datetime(stocks.index).week

    result = stocks.groupby(stocks["week_id"])["price"].transform(func)

    expected = pd.Series(data=pd.to_datetime(values), index=dates, name="price")

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("func", ["cumsum", "cumprod", "cummin", "cummax"])
def test_transform_absent_categories(func):
    # GH 16771
    # cython transforms with more groups than rows
    x_vals = [1]
    x_cats = range(2)
    y = [1]
    df = DataFrame(dict(x=Categorical(x_vals, x_cats), y=y))
    result = getattr(df.y.groupby(df.x), func)()
    expected = df.y
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("func", ["ffill", "bfill", "shift"])
@pytest.mark.parametrize("key, val", [("level", 0), ("by", Series([0]))])
def test_ffill_not_in_axis(func, key, val):
    # GH 21521
    df = pd.DataFrame([[np.nan]])
    result = getattr(df.groupby(**{key: val}), func)()
    expected = df

    tm.assert_frame_equal(result, expected)


def test_transform_invalid_name_raises():
    # GH#27486
    df = DataFrame(dict(a=[0, 1, 1, 2]))
    g = df.groupby(["a", "b", "b", "c"])
    with pytest.raises(ValueError, match="not a valid function name"):
        g.transform("some_arbitrary_name")

    # method exists on the object, but is not a valid transformation/agg
    assert hasattr(g, "aggregate")  # make sure the method exists
    with pytest.raises(ValueError, match="not a valid function name"):
        g.transform("aggregate")

    # Test SeriesGroupBy
    g = df["a"].groupby(["a", "b", "b", "c"])
    with pytest.raises(ValueError, match="not a valid function name"):
        g.transform("some_arbitrary_name")


@pytest.mark.parametrize(
    "obj",
    [
        DataFrame(
            dict(a=[0, 0, 0, 1, 1, 1], b=range(6)), index=["A", "B", "C", "D", "E", "F"]
        ),
        Series([0, 0, 0, 1, 1, 1], index=["A", "B", "C", "D", "E", "F"]),
    ],
)
def test_transform_agg_by_name(reduction_func, obj):
    func = reduction_func
    g = obj.groupby(np.repeat([0, 1], 3))

    if func == "ngroup":  # GH#27468
        pytest.xfail("TODO: g.transform('ngroup') doesn't work")
    if func == "size":  # GH#27469
        pytest.xfail("TODO: g.transform('size') doesn't work")

    args = {"nth": [0], "quantile": [0.5]}.get(func, [])

    result = g.transform(func, *args)

    # this is the *definition* of a transformation
    tm.assert_index_equal(result.index, obj.index)
    if hasattr(obj, "columns"):
        tm.assert_index_equal(result.columns, obj.columns)

    # verify that values were broadcasted across each group
    assert len(set(DataFrame(result).iloc[-3:, -1])) == 1


def test_transform_lambda_with_datetimetz():
    # GH 27496
    df = DataFrame(
        {
            "time": [
                Timestamp("2010-07-15 03:14:45"),
                Timestamp("2010-11-19 18:47:06"),
            ],
            "timezone": ["Etc/GMT+4", "US/Eastern"],
        }
    )
    result = df.groupby(["timezone"])["time"].transform(
        lambda x: x.dt.tz_localize(x.name)
    )
    expected = Series(
        [
            Timestamp("2010-07-15 03:14:45", tz="Etc/GMT+4"),
            Timestamp("2010-11-19 18:47:06", tz="US/Eastern"),
        ],
        name="time",
    )
    tm.assert_series_equal(result, expected)


def test_transform_fastpath_raises():
    # GH#29631 case where fastpath defined in groupby.generic _choose_path
    #  raises, but slow_path does not

    df = pd.DataFrame({"A": [1, 1, 2, 2], "B": [1, -1, 1, 2]})
    gb = df.groupby("A")

    def func(grp):
        # we want a function such that func(frame) fails but func.apply(frame)
        #  works
        if grp.ndim == 2:
            # Ensure that fast_path fails
            raise NotImplementedError("Don't cross the streams")
        return grp * 2

    # Check that the fastpath raises, see _transform_general
    obj = gb._obj_with_exclusions
    gen = gb.grouper.get_iterator(obj, axis=gb.axis)
    fast_path, slow_path = gb._define_paths(func)
    _, group = next(gen)

    with pytest.raises(NotImplementedError, match="Don't cross the streams"):
        fast_path(group)

    result = gb.transform(func)

    expected = pd.DataFrame([2, -2, 2, 4], columns=["B"])
    tm.assert_frame_equal(result, expected)
