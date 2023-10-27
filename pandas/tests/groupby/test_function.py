import builtins
import re

import numpy as np
import pytest

from pandas._libs import lib

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Series,
    Timestamp,
    date_range,
)
import pandas._testing as tm
from pandas.tests.groupby import get_groupby_method_args
from pandas.util import _test_decorators as td


def test_intercept_builtin_sum():
    s = Series([1.0, 2.0, np.nan, 3.0])
    grouped = s.groupby([0, 1, 2, 2])

    msg = "using SeriesGroupBy.sum"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        # GH#53425
        result = grouped.agg(builtins.sum)
    msg = "using np.sum"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        # GH#53425
        result2 = grouped.apply(builtins.sum)
    expected = grouped.sum()
    tm.assert_series_equal(result, expected)
    tm.assert_series_equal(result2, expected)


@pytest.mark.parametrize("f", [max, min, sum])
@pytest.mark.parametrize("keys", ["jim", ["jim", "joe"]])  # Single key  # Multi-key
def test_builtins_apply(keys, f):
    # see gh-8155
    rs = np.random.default_rng(2)
    df = DataFrame(rs.integers(1, 7, (10, 2)), columns=["jim", "joe"])
    df["jolie"] = rs.standard_normal(10)

    gb = df.groupby(keys)

    fname = f.__name__

    warn = None if f is not sum else FutureWarning
    msg = "The behavior of DataFrame.sum with axis=None is deprecated"
    with tm.assert_produces_warning(
        warn, match=msg, check_stacklevel=False, raise_on_extra_warnings=False
    ):
        # Also warns on deprecation GH#53425
        result = gb.apply(f)
    ngroups = len(df.drop_duplicates(subset=keys))

    assert_msg = f"invalid frame shape: {result.shape} (expected ({ngroups}, 3))"
    assert result.shape == (ngroups, 3), assert_msg

    npfunc = lambda x: getattr(np, fname)(x, axis=0)  # numpy's equivalent function
    msg = "DataFrameGroupBy.apply operated on the grouping columns"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        expected = gb.apply(npfunc)
    tm.assert_frame_equal(result, expected)

    with tm.assert_produces_warning(FutureWarning, match=msg):
        expected2 = gb.apply(lambda x: npfunc(x))
    tm.assert_frame_equal(result, expected2)

    if f != sum:
        expected = gb.agg(fname).reset_index()
        expected.set_index(keys, inplace=True, drop=False)
        tm.assert_frame_equal(result, expected, check_dtype=False)

    tm.assert_series_equal(getattr(result, fname)(axis=0), getattr(df, fname)(axis=0))


class TestNumericOnly:
    # make sure that we are passing thru kwargs to our agg functions

    @pytest.fixture
    def df(self):
        # GH3668
        # GH5724
        df = DataFrame(
            {
                "group": [1, 1, 2],
                "int": [1, 2, 3],
                "float": [4.0, 5.0, 6.0],
                "string": list("abc"),
                "category_string": Series(list("abc")).astype("category"),
                "category_int": [7, 8, 9],
                "datetime": date_range("20130101", periods=3),
                "datetimetz": date_range("20130101", periods=3, tz="US/Eastern"),
                "timedelta": pd.timedelta_range("1 s", periods=3, freq="s"),
            },
            columns=[
                "group",
                "int",
                "float",
                "string",
                "category_string",
                "category_int",
                "datetime",
                "datetimetz",
                "timedelta",
            ],
        )
        return df

    @pytest.mark.parametrize("method", ["mean", "median"])
    def test_averages(self, df, method):
        # mean / median
        expected_columns_numeric = Index(["int", "float", "category_int"])

        gb = df.groupby("group")
        expected = DataFrame(
            {
                "category_int": [7.5, 9],
                "float": [4.5, 6.0],
                "timedelta": [pd.Timedelta("1.5s"), pd.Timedelta("3s")],
                "int": [1.5, 3],
                "datetime": [
                    Timestamp("2013-01-01 12:00:00"),
                    Timestamp("2013-01-03 00:00:00"),
                ],
                "datetimetz": [
                    Timestamp("2013-01-01 12:00:00", tz="US/Eastern"),
                    Timestamp("2013-01-03 00:00:00", tz="US/Eastern"),
                ],
            },
            index=Index([1, 2], name="group"),
            columns=[
                "int",
                "float",
                "category_int",
            ],
        )

        result = getattr(gb, method)(numeric_only=True)
        tm.assert_frame_equal(result.reindex_like(expected), expected)

        expected_columns = expected.columns

        self._check(df, method, expected_columns, expected_columns_numeric)

    @pytest.mark.parametrize("method", ["min", "max"])
    def test_extrema(self, df, method):
        # TODO: min, max *should* handle
        # categorical (ordered) dtype

        expected_columns = Index(
            [
                "int",
                "float",
                "string",
                "category_int",
                "datetime",
                "datetimetz",
                "timedelta",
            ]
        )
        expected_columns_numeric = expected_columns

        self._check(df, method, expected_columns, expected_columns_numeric)

    @pytest.mark.parametrize("method", ["first", "last"])
    def test_first_last(self, df, method):
        expected_columns = Index(
            [
                "int",
                "float",
                "string",
                "category_string",
                "category_int",
                "datetime",
                "datetimetz",
                "timedelta",
            ]
        )
        expected_columns_numeric = expected_columns

        self._check(df, method, expected_columns, expected_columns_numeric)

    @pytest.mark.parametrize("method", ["sum", "cumsum"])
    def test_sum_cumsum(self, df, method):
        expected_columns_numeric = Index(["int", "float", "category_int"])
        expected_columns = Index(
            ["int", "float", "string", "category_int", "timedelta"]
        )
        if method == "cumsum":
            # cumsum loses string
            expected_columns = Index(["int", "float", "category_int", "timedelta"])

        self._check(df, method, expected_columns, expected_columns_numeric)

    @pytest.mark.parametrize("method", ["prod", "cumprod"])
    def test_prod_cumprod(self, df, method):
        expected_columns = Index(["int", "float", "category_int"])
        expected_columns_numeric = expected_columns

        self._check(df, method, expected_columns, expected_columns_numeric)

    @pytest.mark.parametrize("method", ["cummin", "cummax"])
    def test_cummin_cummax(self, df, method):
        # like min, max, but don't include strings
        expected_columns = Index(
            ["int", "float", "category_int", "datetime", "datetimetz", "timedelta"]
        )

        # GH#15561: numeric_only=False set by default like min/max
        expected_columns_numeric = expected_columns

        self._check(df, method, expected_columns, expected_columns_numeric)

    def _check(self, df, method, expected_columns, expected_columns_numeric):
        gb = df.groupby("group")

        # object dtypes for transformations are not implemented in Cython and
        # have no Python fallback
        exception = NotImplementedError if method.startswith("cum") else TypeError

        if method in ("min", "max", "cummin", "cummax", "cumsum", "cumprod"):
            # The methods default to numeric_only=False and raise TypeError
            msg = "|".join(
                [
                    "Categorical is not ordered",
                    f"Cannot perform {method} with non-ordered Categorical",
                    re.escape(f"agg function failed [how->{method},dtype->object]"),
                    # cumsum/cummin/cummax/cumprod
                    "function is not implemented for this dtype",
                ]
            )
            with pytest.raises(exception, match=msg):
                getattr(gb, method)()
        elif method in ("sum", "mean", "median", "prod"):
            msg = "|".join(
                [
                    "category type does not support sum operations",
                    re.escape(f"agg function failed [how->{method},dtype->object]"),
                ]
            )
            with pytest.raises(exception, match=msg):
                getattr(gb, method)()
        else:
            result = getattr(gb, method)()
            tm.assert_index_equal(result.columns, expected_columns_numeric)

        if method not in ("first", "last"):
            msg = "|".join(
                [
                    "Categorical is not ordered",
                    "category type does not support",
                    "function is not implemented for this dtype",
                    f"Cannot perform {method} with non-ordered Categorical",
                    re.escape(f"agg function failed [how->{method},dtype->object]"),
                ]
            )
            with pytest.raises(exception, match=msg):
                getattr(gb, method)(numeric_only=False)
        else:
            result = getattr(gb, method)(numeric_only=False)
            tm.assert_index_equal(result.columns, expected_columns)


class TestGroupByNonCythonPaths:
    # GH#5610 non-cython calls should not include the grouper
    # Tests for code not expected to go through cython paths.

    @pytest.fixture
    def df(self):
        df = DataFrame(
            [[1, 2, "foo"], [1, np.nan, "bar"], [3, np.nan, "baz"]],
            columns=["A", "B", "C"],
        )
        return df

    @pytest.fixture
    def gb(self, df):
        gb = df.groupby("A")
        return gb

    @pytest.fixture
    def gni(self, df):
        gni = df.groupby("A", as_index=False)
        return gni

    def test_describe(self, df, gb, gni):
        # describe
        expected_index = Index([1, 3], name="A")
        expected_col = MultiIndex(
            levels=[["B"], ["count", "mean", "std", "min", "25%", "50%", "75%", "max"]],
            codes=[[0] * 8, list(range(8))],
        )
        expected = DataFrame(
            [
                [1.0, 2.0, np.nan, 2.0, 2.0, 2.0, 2.0, 2.0],
                [0.0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            ],
            index=expected_index,
            columns=expected_col,
        )
        result = gb.describe()
        tm.assert_frame_equal(result, expected)

        expected = expected.reset_index()
        result = gni.describe()
        tm.assert_frame_equal(result, expected)


def test_cython_api2():
    # this takes the fast apply path

    # cumsum (GH5614)
    df = DataFrame([[1, 2, np.nan], [1, np.nan, 9], [3, 4, 9]], columns=["A", "B", "C"])
    expected = DataFrame([[2, np.nan], [np.nan, 9], [4, 9]], columns=["B", "C"])
    result = df.groupby("A").cumsum()
    tm.assert_frame_equal(result, expected)

    # GH 5755 - cumsum is a transformer and should ignore as_index
    result = df.groupby("A", as_index=False).cumsum()
    tm.assert_frame_equal(result, expected)

    # GH 13994
    msg = "DataFrameGroupBy.cumsum with axis=1 is deprecated"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        result = df.groupby("A").cumsum(axis=1)
    expected = df.cumsum(axis=1)
    tm.assert_frame_equal(result, expected)

    msg = "DataFrameGroupBy.cumprod with axis=1 is deprecated"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        result = df.groupby("A").cumprod(axis=1)
    expected = df.cumprod(axis=1)
    tm.assert_frame_equal(result, expected)


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


@pytest.mark.parametrize(
    "i",
    [
        (
            Timestamp("2011-01-15 12:50:28.502376"),
            Timestamp("2011-01-20 12:50:28.593448"),
        ),
        (24650000000000001, 24650000000000002),
    ],
)
def test_groupby_non_arithmetic_agg_int_like_precision(i):
    # see gh-6620, gh-9311
    df = DataFrame([{"a": 1, "b": i[0]}, {"a": 1, "b": i[1]}])

    grp_exp = {
        "first": {"expected": i[0]},
        "last": {"expected": i[1]},
        "min": {"expected": i[0]},
        "max": {"expected": i[1]},
        "nth": {"expected": i[1], "args": [1]},
        "count": {"expected": 2},
    }

    for method, data in grp_exp.items():
        if "args" not in data:
            data["args"] = []

        grouped = df.groupby("a")
        res = getattr(grouped, method)(*data["args"])

        assert res.iloc[0].b == data["expected"]


@pytest.mark.parametrize("numeric_only", [True, False, None])
def test_axis1_numeric_only(request, groupby_func, numeric_only):
    if groupby_func in ("idxmax", "idxmin"):
        pytest.skip("idxmax and idx_min tested in test_idxmin_idxmax_axis1")
    if groupby_func in ("corrwith", "skew"):
        msg = "GH#47723 groupby.corrwith and skew do not correctly implement axis=1"
        request.applymarker(pytest.mark.xfail(reason=msg))

    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)), columns=["A", "B", "C", "D"]
    )
    df["E"] = "x"
    groups = [1, 2, 3, 1, 2, 3, 1, 2, 3, 4]
    gb = df.groupby(groups)
    method = getattr(gb, groupby_func)
    args = get_groupby_method_args(groupby_func, df)
    kwargs = {"axis": 1}
    if numeric_only is not None:
        # when numeric_only is None we don't pass any argument
        kwargs["numeric_only"] = numeric_only

    # Functions without numeric_only and axis args
    no_args = ("cumprod", "cumsum", "diff", "fillna", "pct_change", "rank", "shift")
    # Functions with axis args
    has_axis = (
        "cumprod",
        "cumsum",
        "diff",
        "pct_change",
        "rank",
        "shift",
        "cummax",
        "cummin",
        "idxmin",
        "idxmax",
        "fillna",
    )
    warn_msg = f"DataFrameGroupBy.{groupby_func} with axis=1 is deprecated"
    if numeric_only is not None and groupby_func in no_args:
        msg = "got an unexpected keyword argument 'numeric_only'"
        if groupby_func in ["cumprod", "cumsum"]:
            with pytest.raises(TypeError, match=msg):
                with tm.assert_produces_warning(FutureWarning, match=warn_msg):
                    method(*args, **kwargs)
        else:
            with pytest.raises(TypeError, match=msg):
                method(*args, **kwargs)
    elif groupby_func not in has_axis:
        msg = "got an unexpected keyword argument 'axis'"
        with pytest.raises(TypeError, match=msg):
            method(*args, **kwargs)
    # fillna and shift are successful even on object dtypes
    elif (numeric_only is None or not numeric_only) and groupby_func not in (
        "fillna",
        "shift",
    ):
        msgs = (
            # cummax, cummin, rank
            "not supported between instances of",
            # cumprod
            "can't multiply sequence by non-int of type 'float'",
            # cumsum, diff, pct_change
            "unsupported operand type",
        )
        with pytest.raises(TypeError, match=f"({'|'.join(msgs)})"):
            with tm.assert_produces_warning(FutureWarning, match=warn_msg):
                method(*args, **kwargs)
    else:
        with tm.assert_produces_warning(FutureWarning, match=warn_msg):
            result = method(*args, **kwargs)

        df_expected = df.drop(columns="E").T if numeric_only else df.T
        expected = getattr(df_expected, groupby_func)(*args).T
        if groupby_func == "shift" and not numeric_only:
            # shift with axis=1 leaves the leftmost column as numeric
            # but transposing for expected gives us object dtype
            expected = expected.astype(float)

        tm.assert_equal(result, expected)


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
        pytest.param("sem", scipy_sem, marks=td.skip_if_no_scipy),
    ],
)
def test_ops_general(op, targop):
    df = DataFrame(np.random.default_rng(2).standard_normal(1000))
    labels = np.random.default_rng(2).integers(0, 50, size=1000).astype(float)

    result = getattr(df.groupby(labels), op)()
    warn = None if op in ("first", "last", "count", "sem") else FutureWarning
    msg = f"using DataFrameGroupBy.{op}"
    with tm.assert_produces_warning(warn, match=msg):
        expected = df.groupby(labels).agg(targop)
    tm.assert_frame_equal(result, expected)


def test_series_index_name(df):
    grouped = df.loc[:, ["C"]].groupby(df["A"])
    result = grouped.agg(lambda x: x.mean())
    assert result.index.name == "A"


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
    idx = Index([1, 2, 3], name="a", dtype="Int64")
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
    "kernel, has_arg",
    [
        ("all", False),
        ("any", False),
        ("bfill", False),
        ("corr", True),
        ("corrwith", True),
        ("cov", True),
        ("cummax", True),
        ("cummin", True),
        ("cumprod", True),
        ("cumsum", True),
        ("diff", False),
        ("ffill", False),
        ("fillna", False),
        ("first", True),
        ("idxmax", True),
        ("idxmin", True),
        ("last", True),
        ("max", True),
        ("mean", True),
        ("median", True),
        ("min", True),
        ("nth", False),
        ("nunique", False),
        ("pct_change", False),
        ("prod", True),
        ("quantile", True),
        ("sem", True),
        ("skew", True),
        ("std", True),
        ("sum", True),
        ("var", True),
    ],
)
@pytest.mark.parametrize("numeric_only", [True, False, lib.no_default])
@pytest.mark.parametrize("keys", [["a1"], ["a1", "a2"]])
def test_numeric_only(kernel, has_arg, numeric_only, keys):
    # GH#46072
    # drops_nuisance: Whether the op drops nuisance columns even when numeric_only=False
    # has_arg: Whether the op has a numeric_only arg
    df = DataFrame({"a1": [1, 1], "a2": [2, 2], "a3": [5, 6], "b": 2 * [object]})

    args = get_groupby_method_args(kernel, df)
    kwargs = {} if numeric_only is lib.no_default else {"numeric_only": numeric_only}

    gb = df.groupby(keys)
    method = getattr(gb, kernel)
    if has_arg and numeric_only is True:
        # Cases where b does not appear in the result
        result = method(*args, **kwargs)
        assert "b" not in result.columns
    elif (
        # kernels that work on any dtype and have numeric_only arg
        kernel in ("first", "last")
        or (
            # kernels that work on any dtype and don't have numeric_only arg
            kernel in ("any", "all", "bfill", "ffill", "fillna", "nth", "nunique")
            and numeric_only is lib.no_default
        )
    ):
        result = method(*args, **kwargs)
        assert "b" in result.columns
    elif has_arg:
        assert numeric_only is not True
        # kernels that are successful on any dtype were above; this will fail

        # object dtypes for transformations are not implemented in Cython and
        # have no Python fallback
        exception = NotImplementedError if kernel.startswith("cum") else TypeError

        msg = "|".join(
            [
                "not allowed for this dtype",
                "cannot be performed against 'object' dtypes",
                # On PY39 message is "a number"; on PY310 and after is "a real number"
                "must be a string or a.* number",
                "unsupported operand type",
                "function is not implemented for this dtype",
                re.escape(f"agg function failed [how->{kernel},dtype->object]"),
            ]
        )
        if kernel == "idxmin":
            msg = "'<' not supported between instances of 'type' and 'type'"
        elif kernel == "idxmax":
            msg = "'>' not supported between instances of 'type' and 'type'"
        with pytest.raises(exception, match=msg):
            method(*args, **kwargs)
    elif not has_arg and numeric_only is not lib.no_default:
        with pytest.raises(
            TypeError, match="got an unexpected keyword argument 'numeric_only'"
        ):
            method(*args, **kwargs)
    else:
        assert kernel in ("diff", "pct_change")
        assert numeric_only is lib.no_default
        # Doesn't have numeric_only argument and fails on nuisance columns
        with pytest.raises(TypeError, match=r"unsupported operand type"):
            method(*args, **kwargs)


@pytest.mark.filterwarnings("ignore:Downcasting object dtype arrays:FutureWarning")
@pytest.mark.parametrize("dtype", [bool, int, float, object])
def test_deprecate_numeric_only_series(dtype, groupby_func, request):
    # GH#46560
    grouper = [0, 0, 1]

    ser = Series([1, 0, 0], dtype=dtype)
    gb = ser.groupby(grouper)

    if groupby_func == "corrwith":
        # corrwith is not implemented on SeriesGroupBy
        assert not hasattr(gb, groupby_func)
        return

    method = getattr(gb, groupby_func)

    expected_ser = Series([1, 0, 0])
    expected_gb = expected_ser.groupby(grouper)
    expected_method = getattr(expected_gb, groupby_func)

    args = get_groupby_method_args(groupby_func, ser)

    fails_on_numeric_object = (
        "corr",
        "cov",
        "cummax",
        "cummin",
        "cumprod",
        "cumsum",
        "quantile",
    )
    # ops that give an object result on object input
    obj_result = (
        "first",
        "last",
        "nth",
        "bfill",
        "ffill",
        "shift",
        "sum",
        "diff",
        "pct_change",
        "var",
        "mean",
        "median",
        "min",
        "max",
        "prod",
        "skew",
    )

    # Test default behavior; kernels that fail may be enabled in the future but kernels
    # that succeed should not be allowed to fail (without deprecation, at least)
    if groupby_func in fails_on_numeric_object and dtype is object:
        if groupby_func == "quantile":
            msg = "cannot be performed against 'object' dtypes"
        else:
            msg = "is not supported for object dtype"
        with pytest.raises(TypeError, match=msg):
            method(*args)
    elif dtype is object:
        result = method(*args)
        expected = expected_method(*args)
        if groupby_func in obj_result:
            expected = expected.astype(object)
        tm.assert_series_equal(result, expected)

    has_numeric_only = (
        "first",
        "last",
        "max",
        "mean",
        "median",
        "min",
        "prod",
        "quantile",
        "sem",
        "skew",
        "std",
        "sum",
        "var",
        "cummax",
        "cummin",
        "cumprod",
        "cumsum",
    )
    if groupby_func not in has_numeric_only:
        msg = "got an unexpected keyword argument 'numeric_only'"
        with pytest.raises(TypeError, match=msg):
            method(*args, numeric_only=True)
    elif dtype is object:
        msg = "|".join(
            [
                "SeriesGroupBy.sem called with numeric_only=True and dtype object",
                "Series.skew does not allow numeric_only=True with non-numeric",
                "cum(sum|prod|min|max) is not supported for object dtype",
                r"Cannot use numeric_only=True with SeriesGroupBy\..* and non-numeric",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            method(*args, numeric_only=True)
    elif dtype == bool and groupby_func == "quantile":
        msg = "Allowing bool dtype in SeriesGroupBy.quantile"
        with tm.assert_produces_warning(FutureWarning, match=msg):
            # GH#51424
            result = method(*args, numeric_only=True)
            expected = method(*args, numeric_only=False)
        tm.assert_series_equal(result, expected)
    else:
        result = method(*args, numeric_only=True)
        expected = method(*args, numeric_only=False)
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("dtype", [int, float, object])
@pytest.mark.parametrize(
    "kwargs",
    [
        {"percentiles": [0.10, 0.20, 0.30], "include": "all", "exclude": None},
        {"percentiles": [0.10, 0.20, 0.30], "include": None, "exclude": ["int"]},
        {"percentiles": [0.10, 0.20, 0.30], "include": ["int"], "exclude": None},
    ],
)
def test_groupby_empty_dataset(dtype, kwargs):
    # GH#41575
    df = DataFrame([[1, 2, 3]], columns=["A", "B", "C"], dtype=dtype)
    df["B"] = df["B"].astype(int)
    df["C"] = df["C"].astype(float)

    result = df.iloc[:0].groupby("A").describe(**kwargs)
    expected = df.groupby("A").describe(**kwargs).reset_index(drop=True).iloc[:0]
    tm.assert_frame_equal(result, expected)

    result = df.iloc[:0].groupby("A").B.describe(**kwargs)
    expected = df.groupby("A").B.describe(**kwargs).reset_index(drop=True).iloc[:0]
    expected.index = Index([])
    tm.assert_frame_equal(result, expected)


def test_multiindex_group_all_columns_when_empty(groupby_func):
    # GH 32464
    df = DataFrame({"a": [], "b": [], "c": []}).set_index(["a", "b", "c"])
    gb = df.groupby(["a", "b", "c"], group_keys=False)
    method = getattr(gb, groupby_func)
    args = get_groupby_method_args(groupby_func, df)

    result = method(*args).index
    expected = df.index
    tm.assert_index_equal(result, expected)


def test_duplicate_columns(request, groupby_func, as_index):
    # GH#50806
    if groupby_func == "corrwith":
        msg = "GH#50845 - corrwith fails when there are duplicate columns"
        request.applymarker(pytest.mark.xfail(reason=msg))
    df = DataFrame([[1, 3, 6], [1, 4, 7], [2, 5, 8]], columns=list("abb"))
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby("a", as_index=as_index)
    result = getattr(gb, groupby_func)(*args)

    expected_df = df.set_axis(["a", "b", "c"], axis=1)
    expected_args = get_groupby_method_args(groupby_func, expected_df)
    expected_gb = expected_df.groupby("a", as_index=as_index)
    expected = getattr(expected_gb, groupby_func)(*expected_args)
    if groupby_func not in ("size", "ngroup", "cumcount"):
        expected = expected.rename(columns={"c": "b"})
    tm.assert_equal(result, expected)


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
        "std",
        "var",
        "sem",
    ],
)
@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize("sort", [True, False])
def test_regression_allowlist_methods(op, axis, skipna, sort):
    # GH6944
    # GH 17537
    # explicitly test the allowlist methods
    raw_frame = DataFrame([0])
    if axis == 0:
        frame = raw_frame
        msg = "The 'axis' keyword in DataFrame.groupby is deprecated and will be"
    else:
        frame = raw_frame.T
        msg = "DataFrame.groupby with axis=1 is deprecated"

    with tm.assert_produces_warning(FutureWarning, match=msg):
        grouped = frame.groupby(level=0, axis=axis, sort=sort)

    if op == "skew":
        # skew has skipna
        result = getattr(grouped, op)(skipna=skipna)
        expected = frame.groupby(level=0).apply(
            lambda h: getattr(h, op)(axis=axis, skipna=skipna)
        )
        if sort:
            expected = expected.sort_index(axis=axis)
        tm.assert_frame_equal(result, expected)
    else:
        result = getattr(grouped, op)()
        expected = frame.groupby(level=0).apply(lambda h: getattr(h, op)(axis=axis))
        if sort:
            expected = expected.sort_index(axis=axis)
        tm.assert_frame_equal(result, expected)
