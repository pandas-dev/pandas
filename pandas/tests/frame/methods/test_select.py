import re

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def df():
    return pd.DataFrame(
        {
            "a": [1, 2, 3],
            "b": [4.0, 5.0, 6.0],
            "c": ["x", "y", "z"],
        }
    )


@pytest.fixture
def mi_df():
    return pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]],
        columns=pd.MultiIndex.from_tuples([("one", "a"), ("one", "b"), ("two", "a")]),
    )


def test_select(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", "c")
    expected = df[["a", "c"]]
    tm.assert_frame_equal(result, expected)


def test_select_single_label_returns_frame(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a")
    expected = df[["a"]]
    tm.assert_frame_equal(result, expected)


def test_select_list(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(["b", "a"])
    expected = df[["b", "a"]]
    tm.assert_frame_equal(result, expected)


def test_select_index(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(df.columns[:2])
    expected = df[df.columns[:2]]
    tm.assert_frame_equal(result, expected)


def test_select_ndarray(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(np.array(["b", "a"]))
    expected = df[["b", "a"]]
    tm.assert_frame_equal(result, expected)


def test_select_generator(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(c for c in ["a", "c"])
    expected = df[["a", "c"]]
    tm.assert_frame_equal(result, expected)


def test_select_set_raises(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    msg = "`select` does not support sets"
    with pytest.raises(TypeError, match=msg):
        df.select({"a", "b"})


def test_select_list_mixed_with_args_raises(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    msg = "`select` supports individual columns .* or a single list .*, but not both"
    with pytest.raises(TypeError, match=msg):
        df.select(["a", "b"], "c")
    with pytest.raises(TypeError, match=msg):
        df.select("c", ["a", "b"])
    with pytest.raises(TypeError, match=msg):
        df.select(df.columns[:2], "c")


def test_select_reorder(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("c", "a", "b")
    expected = df[["c", "a", "b"]]
    tm.assert_frame_equal(result, expected)


def test_select_label_twice(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", "a")
    expected = df[["a", "a"]]
    tm.assert_frame_equal(result, expected)


def test_select_duplicate_columns():
    # https://github.com/pandas-dev/pandas/issues/61522
    df = pd.DataFrame([[1, 2, 3]], columns=["a", "b", "a"])
    result = df.select("a")
    expected = df[["a"]]
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("columns", [["a", "a", "b"], ["a", "b", "a"]])
def test_select_kwarg_duplicate_label(columns):
    # https://github.com/pandas-dev/pandas/issues/61522
    # a keyword argument defines a single output column even when its
    # name is duplicated in the DataFrame's columns
    df = pd.DataFrame([[1, 2, 3]], columns=columns)
    result = df.select("b", a=99)
    expected = pd.DataFrame({"b": df["b"], "a": [99]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("columns", [["b", "b", "a"], ["b", "a", "b"]])
def test_select_expression_duplicate_label(columns):
    # https://github.com/pandas-dev/pandas/issues/61522
    # a renamed expression defines a single output column even when its
    # name is duplicated in the DataFrame's columns
    df = pd.DataFrame([[1, 2, 3]], columns=columns)
    result = df.select(pd.col("a").rename("b"))
    expected = pd.DataFrame({"b": df["a"]})
    tm.assert_frame_equal(result, expected)


def test_select_missing_raises(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    with pytest.raises(KeyError, match=r"\['nope'\] not in index"):
        df.select("a", "nope")


def test_select_empty(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select()
    expected = df[[]]
    tm.assert_frame_equal(result, expected)


def test_select_empty_multiindex_columns(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # the empty result keeps the columns' type, levels, and names
    result = mi_df.select()
    expected = mi_df[[]]
    tm.assert_frame_equal(result, expected)


def test_select_empty_preserves_columns_name(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    df.columns.name = "cols"
    result = df.select()
    expected = df[[]]
    tm.assert_frame_equal(result, expected)


def test_select_index_preserved():
    # https://github.com/pandas-dev/pandas/issues/61522
    df = pd.DataFrame({"a": [1, 2]}, index=["x", "y"])
    result = df.select("a")
    expected = pd.DataFrame({"a": [1, 2]}, index=["x", "y"])
    tm.assert_frame_equal(result, expected)


def test_select_returns_copy(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", "b")
    result.iloc[0, 0] = 100
    expected = pd.DataFrame(
        {"a": [1, 2, 3], "b": [4.0, 5.0, 6.0], "c": ["x", "y", "z"]}
    )
    tm.assert_frame_equal(df, expected)


def test_select_tuple_label():
    # https://github.com/pandas-dev/pandas/issues/61522
    df = pd.DataFrame([[1, 2]], columns=pd.Index([("a", "b"), "c"], dtype=object))
    result = df.select(("a", "b"))
    expected = df[[("a", "b")]]
    tm.assert_frame_equal(result, expected)


def test_select_multiindex_first_level(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = mi_df.select("one")
    expected = mi_df[["one"]]
    tm.assert_frame_equal(result, expected)


def test_select_multiindex_tuples(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = mi_df.select(("one", "b"), ("two", "a"))
    expected = mi_df[[("one", "b"), ("two", "a")]]
    tm.assert_frame_equal(result, expected)


def test_select_multiindex_mixed_depth(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # each label resolves independently, so first-level keys and full
    # tuples can be mixed even though df[list] rejects the combination
    result = mi_df.select("one", ("two", "a"))
    expected = mi_df[[("one", "a"), ("one", "b"), ("two", "a")]]
    tm.assert_frame_equal(result, expected)


def test_select_multiindex_mixed_depth_list(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = mi_df.select([("two", "a"), "one"])
    expected = mi_df[[("two", "a"), ("one", "a"), ("one", "b")]]
    tm.assert_frame_equal(result, expected)


def test_select_multiindex_mixed_depth_expression_position(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # an interleaved expression does not change how labels resolve
    expr = (pd.col(("one", "a")) * 2).rename(("one", "X"))
    result = mi_df.select("one", expr, ("two", "a"))
    expected = mi_df.select("one", ("two", "a"), expr)[
        [("one", "a"), ("one", "b"), ("one", "X"), ("two", "a")]
    ]
    tm.assert_frame_equal(result, expected)


def test_select_expression_multiindex(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = mi_df.select(pd.col(("one", "a")) * 2)
    expected = pd.DataFrame(
        [[2], [8]], columns=pd.MultiIndex.from_tuples([("one", "a")])
    )
    tm.assert_frame_equal(result, expected)


def test_select_multiindex_label_and_expression(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = mi_df.select(("one", "b"), pd.col(("two", "a")) + 1)
    expected = pd.DataFrame(
        [[2, 4], [5, 7]],
        columns=pd.MultiIndex.from_tuples([("one", "b"), ("two", "a")]),
    )
    tm.assert_frame_equal(result, expected)


def test_select_expression_multiindex_rename(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    expr = (pd.col(("one", "a")) + pd.col(("one", "b"))).rename(("one", "total"))
    result = mi_df.select(expr)
    expected = pd.DataFrame(
        [[3], [9]], columns=pd.MultiIndex.from_tuples([("one", "total")])
    )
    tm.assert_frame_equal(result, expected)


def test_select_expression_multiindex_str_name_raises(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    msg = "use .rename(...) with a tuple of length 2"
    with pytest.raises(TypeError, match=re.escape(msg)):
        mi_df.select((pd.col(("one", "a")) * 2).rename("x"))


def test_select_expression_multiindex_wrong_length_name_raises(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    msg = "use .rename(...) with a tuple of length 2"
    with pytest.raises(TypeError, match=re.escape(msg)):
        mi_df.select((pd.col(("one", "a")) * 2).rename(("x",)))


def test_select_expression_multiindex_frame_result(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # to_frame on a tuple-named Series produces MultiIndex columns
    result = mi_df.select(pd.col(("one", "a")).to_frame())
    expected = mi_df[[("one", "a")]]
    tm.assert_frame_equal(result, expected)


def test_select_expression_multiindex_frame_nlevels_mismatch_raises(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    msg = "columns have 1 levels, but the DataFrame's columns have 2 levels"
    with pytest.raises(TypeError, match=re.escape(msg)):
        mi_df.select(pd.col(("one", "a")).astype(str).str.get_dummies())


def test_select_kwargs_multiindex_raises(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    msg = "cannot create column 'extra' from a keyword argument"
    with pytest.raises(TypeError, match=re.escape(msg)):
        mi_df.select(("one", "a"), extra=1)


def test_select_kwargs_preserves_columns_name(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    df.columns.name = "cols"
    result = df.select("a", double=pd.col("a") * 2)
    expected = pd.DataFrame({"a": [1, 2, 3], "double": [2, 4, 6]})
    expected.columns.name = "cols"
    tm.assert_frame_equal(result, expected)


def test_select_expression(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(pd.col("a"))
    expected = df[["a"]]
    tm.assert_frame_equal(result, expected)


def test_select_expression_computed(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # binary ops with scalars keep the column name
    result = df.select(pd.col("a") * 2)
    expected = pd.DataFrame({"a": [2, 4, 6]})
    tm.assert_frame_equal(result, expected)


def test_select_expression_method_call(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(pd.col("c").str.upper())
    expected = pd.DataFrame({"c": ["X", "Y", "Z"]})
    tm.assert_frame_equal(result, expected)


def test_select_expression_rename(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select((pd.col("a") + pd.col("b")).rename("total"))
    expected = pd.DataFrame({"total": [5.0, 7.0, 9.0]})
    tm.assert_frame_equal(result, expected)


def test_select_expression_unnamed_raises(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    msg = re.escape("expression col('a') + col('b') evaluated to an unnamed Series")
    with pytest.raises(TypeError, match=msg):
        df.select(pd.col("a") + pd.col("b"))


def test_select_expression_scalar_raises(df, using_python_scalars):
    # https://github.com/pandas-dev/pandas/issues/61522
    if using_python_scalars:
        result_type = "float"
    else:
        result_type = "float64"
    msg = re.escape(
        f"expression col('b').sum() evaluated to an object of type {result_type}"
    )
    with pytest.raises(TypeError, match=msg):
        df.select(pd.col("b").sum())


def test_select_expression_frame_result(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", pd.col("c").str.get_dummies())
    expected = df[["a"]].join(df["c"].str.get_dummies())
    tm.assert_frame_equal(result, expected)


def test_select_interleaved(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("b", pd.col("a") * 2, "c")
    expected = pd.DataFrame(
        {"b": [4.0, 5.0, 6.0], "a": [2, 4, 6], "c": ["x", "y", "z"]}
    )
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_expression(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", total=pd.col("a") + pd.col("b"))
    expected = pd.DataFrame({"a": [1, 2, 3], "total": [5.0, 7.0, 9.0]})
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_callable(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", double=lambda x: x["a"] * 2)
    expected = pd.DataFrame({"a": [1, 2, 3], "double": [2, 4, 6]})
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_series_aligns(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", extra=pd.Series([10.0], index=[1]))
    expected = pd.DataFrame({"a": [1, 2, 3], "extra": [np.nan, 10.0, np.nan]})
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_scalar_broadcasts(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", flag=True)
    expected = pd.DataFrame({"a": [1, 2, 3], "flag": [True, True, True]})
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_dict_aligns(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select("a", extra={0: 1.0, 2: 3.0})
    expected = pd.DataFrame({"a": [1, 2, 3], "extra": [1.0, np.nan, 3.0]})
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_wrong_length_raises(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    msg = re.escape("Length of values (2) does not match length of index (3)")
    with pytest.raises(ValueError, match=msg):
        df.select("a", extra=[1, 2])


def test_select_kwargs_sequential(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # like assign, later kwargs can refer to columns computed earlier
    result = df.select(double=pd.col("a") * 2, quadruple=pd.col("double") * 2)
    expected = pd.DataFrame({"double": [2, 4, 6], "quadruple": [4, 8, 12]})
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_sequential_callable(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(double=lambda x: x["a"] * 2, quadruple=lambda x: x["double"] * 2)
    expected = pd.DataFrame({"double": [2, 4, 6], "quadruple": [4, 8, 12]})
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_see_positional_expression(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(pd.col("a") * 2, b=pd.col("a") + 1)
    expected = pd.DataFrame({"a": [2, 4, 6], "b": [3, 5, 7]})
    tm.assert_frame_equal(result, expected)


def test_select_label_after_computed_column(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # a label refers to the column as it is at that point in the call
    result = df.select("a", pd.col("a") * 2, "a")
    expected = pd.DataFrame([[1, 2, 2], [2, 4, 4], [3, 6, 6]], columns=["a", "a", "a"])
    tm.assert_frame_equal(result, expected)


def test_select_replaced_column_seen_by_later_kwarg(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(a=pd.col("a") * 2, c=pd.col("a") + 1)
    expected = pd.DataFrame({"a": [2, 4, 6], "c": [3, 5, 7]})
    tm.assert_frame_equal(result, expected)


def test_select_swap_uses_replaced_column(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # like assign, the second kwarg sees the already replaced column
    result = df.select(a=pd.col("b"), b=pd.col("a"))
    expected = pd.DataFrame({"a": [4.0, 5.0, 6.0], "b": [4.0, 5.0, 6.0]})
    tm.assert_frame_equal(result, expected)


def test_select_frame_result_seen_by_later_kwarg(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(pd.col("c").str.get_dummies(), total=pd.col("x") + pd.col("y"))
    expected = df["c"].str.get_dummies().assign(total=lambda d: d["x"] + d["y"])
    tm.assert_frame_equal(result, expected)


def test_select_expression_multiindex_sequential(mi_df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = mi_df.select(
        (pd.col(("one", "a")) * 2).rename(("three", "a")),
        (pd.col(("three", "a")) + 1).rename(("three", "b")),
    )
    expected = pd.DataFrame(
        [[2, 3], [8, 9]],
        columns=pd.MultiIndex.from_tuples([("three", "a"), ("three", "b")]),
    )
    tm.assert_frame_equal(result, expected)


def test_select_does_not_modify_original(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    expected = df.copy()
    df.select(a=pd.col("a") * 2, d=pd.col("a"))
    tm.assert_frame_equal(df, expected)


def test_select_duplicate_output_names(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    # a computed column does not replace an already selected column
    # in the result; both are returned
    result = df.select("a", a=pd.col("b"))
    expected = df[["a", "b"]].set_axis(["a", "a"], axis=1)
    tm.assert_frame_equal(result, expected)


def test_select_kwargs_only(df):
    # https://github.com/pandas-dev/pandas/issues/61522
    result = df.select(total=pd.col("a") + pd.col("b"))
    expected = pd.DataFrame({"total": [5.0, 7.0, 9.0]})
    tm.assert_frame_equal(result, expected)
