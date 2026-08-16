import re

import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def df():
    return pd.DataFrame(
        {
            "g": ["x", "x", "y"],
            "a": [1, 2, 3],
            "b": [4.0, 5.0, 6.0],
        }
    )


def test_select(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    gb = df.groupby("g")
    result = gb.select("a", "b").sum()
    expected = gb[["a", "b"]].sum()
    tm.assert_frame_equal(result, expected)


def test_select_single_label_returns_frame_groupby(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    gb = df.groupby("g")
    result = gb.select("a")
    assert isinstance(result, pd.api.typing.DataFrameGroupBy)
    tm.assert_frame_equal(result.sum(), gb[["a"]].sum())


def test_select_list(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    gb = df.groupby("g")
    result = gb.select(["b", "a"]).sum()
    expected = gb[["b", "a"]].sum()
    tm.assert_frame_equal(result, expected)


def test_select_list_mixed_with_args_raises(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    msg = "`select` supports individual columns .* or a single list .*, but not both"
    with pytest.raises(TypeError, match=msg):
        df.groupby("g").select(["a"], "b")


def test_select_missing_raises(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    with pytest.raises(KeyError, match=re.escape("['nope'] not in index")):
        df.groupby("g").select("a", "nope")


def test_select_after_getitem_selection_raises(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    msg = re.escape("Column(s) ['a'] already selected")
    with pytest.raises(IndexError, match=msg):
        df.groupby("g")[["a"]].select("b")


def test_select_chained(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    gb = df.groupby("g")
    result = gb.select("a", "b").select("a").sum()
    expected = gb[["a"]].sum()
    tm.assert_frame_equal(result, expected)


def test_select_chained_missing_raises(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    # dtype depends on PANDAS_FUTURE_INFER_STRING
    msg = r"None of \[Index\(\['b'\], dtype=.*\)\] are in the \[columns\]"
    with pytest.raises(KeyError, match=msg):
        df.groupby("g").select("a").select("b")


def test_select_then_getitem(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    gb = df.groupby("g")
    result = gb.select("a", "b")["a"].sum()
    expected = gb["a"].sum()
    tm.assert_series_equal(result, expected)


def test_select_multiindex_first_level():
    # https://github.com/pandas-dev/pandas/issues/40322
    df = pd.DataFrame(
        {
            ("g", ""): ["x", "x", "y"],
            ("one", "a"): [1, 2, 3],
            ("one", "b"): [4.0, 5.0, 6.0],
        }
    )
    gb = df.groupby(("g", ""))
    result = gb.select("one").sum()
    expected = gb[[("one", "a"), ("one", "b")]].sum()
    tm.assert_frame_equal(result, expected)


def test_select_expression(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    result = df.groupby("g").select(pd.col("a") * 2).sum()
    expected = df.assign(a=df["a"] * 2).groupby("g")[["a"]].sum()
    tm.assert_frame_equal(result, expected)


def test_select_mixed_labels_and_expressions(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    result = df.groupby("g").select("b", pd.col("a") * 2).sum()
    expected = df.assign(a=df["a"] * 2).groupby("g")[["b", "a"]].sum()
    tm.assert_frame_equal(result, expected)


def test_select_kwargs(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    result = df.groupby("g").select("a", double=pd.col("a") * 2).sum()
    expected = df.assign(double=df["a"] * 2).groupby("g")[["a", "double"]].sum()
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_callable(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    result = df.groupby("g").select(double=lambda x: x["a"] * 2).sum()
    expected = df.assign(double=df["a"] * 2).groupby("g")[["double"]].sum()
    tm.assert_frame_equal(result, expected)


def test_select_expression_uses_grouping_key(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    result = df.groupby("g").select(key_upper=pd.col("g").str.upper()).first()
    expected = (
        df.assign(key_upper=df["g"].str.upper()).groupby("g")[["key_upper"]].first()
    )
    tm.assert_frame_equal(result, expected)


def test_select_expression_unnamed_raises(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    msg = re.escape("expression col('a') + col('b') evaluated to an unnamed Series")
    with pytest.raises(TypeError, match=msg):
        df.groupby("g").select(pd.col("a") + pd.col("b"))


def test_select_expression_scalar_raises(df, using_python_scalars):
    # https://github.com/pandas-dev/pandas/issues/40322
    if using_python_scalars:
        result_type = "float"
    else:
        result_type = "float64"
    msg = re.escape(
        f"expression col('b').sum() evaluated to an object of type {result_type}"
    )
    with pytest.raises(TypeError, match=msg):
        df.groupby("g").select(pd.col("b").sum())


def test_select_expression_transform(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    result = df.groupby("g").select(double=pd.col("a") * 2).cumsum()
    expected = df.assign(double=df["a"] * 2).groupby("g")[["double"]].cumsum()
    tm.assert_frame_equal(result, expected)


def test_select_duplicate_output_names(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    result = df.groupby("g").select("a", a=pd.col("b")).sum()
    expected = df[["g", "a", "b"]].rename(columns={"b": "a"}).groupby("g").sum()
    tm.assert_frame_equal(result, expected)


def test_select_expression_agg(df):
    # https://github.com/pandas-dev/pandas/issues/40322
    result = df.groupby("g").select(double=pd.col("a") * 2).agg("mean")
    expected = df.assign(double=df["a"] * 2).groupby("g")[["double"]].agg("mean")
    tm.assert_frame_equal(result, expected)
