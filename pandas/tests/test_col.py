from datetime import datetime

import numpy as np
import pytest

from pandas._libs.properties import cache_readonly

import pandas as pd
import pandas._testing as tm
from pandas.api.typing import Expression
from pandas.tests.test_register_accessor import ensure_removed


@pytest.mark.parametrize(
    ("expr", "expected_values", "expected_str"),
    [
        (pd.col("a"), [1, 2], "col('a')"),
        (pd.col("a") * 2, [2, 4], "col('a') * 2"),
        (pd.col("a").sum(), [3, 3], "col('a').sum()"),
        (pd.col("a") + 1, [2, 3], "col('a') + 1"),
        (1 + pd.col("a"), [2, 3], "1 + col('a')"),
        (pd.col("a") - 1, [0, 1], "col('a') - 1"),
        (1 - pd.col("a"), [0, -1], "1 - col('a')"),
        (pd.col("a") * 1, [1, 2], "col('a') * 1"),
        (1 * pd.col("a"), [1, 2], "1 * col('a')"),
        (pd.col("a") / 1, [1.0, 2.0], "col('a') / 1"),
        (1 / pd.col("a"), [1.0, 0.5], "1 / col('a')"),
        (pd.col("a") // 1, [1, 2], "col('a') // 1"),
        (1 // pd.col("a"), [1, 0], "1 // col('a')"),
        (pd.col("a") % 1, [0, 0], "col('a') % 1"),
        (1 % pd.col("a"), [0, 1], "1 % col('a')"),
        (pd.col("a") > 1, [False, True], "col('a') > 1"),
        (pd.col("a") >= 1, [True, True], "col('a') >= 1"),
        (pd.col("a") < 1, [False, False], "col('a') < 1"),
        (pd.col("a") <= 1, [True, False], "col('a') <= 1"),
        (pd.col("a") == 1, [True, False], "col('a') == 1"),
        (np.power(pd.col("a"), 2), [1, 4], "power(col('a'), 2)"),
        (np.divide(pd.col("a"), pd.col("a")), [1.0, 1.0], "divide(col('a'), col('a'))"),
        (
            (pd.col("a") + 1) * (pd.col("b") + 2),
            [10, 18],
            "(col('a') + 1) * (col('b') + 2)",
        ),
        (
            (pd.col("a") - 1).astype("bool"),
            [False, True],
            "(col('a') - 1).astype('bool')",
        ),
    ],
)
def test_col_simple(
    expr: Expression, expected_values: list[object], expected_str: str
) -> None:
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    result = df.assign(c=expr)
    expected = pd.DataFrame({"a": [1, 2], "b": [3, 4], "c": expected_values})
    tm.assert_frame_equal(result, expected)
    assert str(expr) == expected_str


def test_frame_getitem() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    expr = pd.col("a") == 2
    result = df[expr]
    expected = df.iloc[[1]]
    tm.assert_frame_equal(result, expected)


def test_frame_setitem() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    expr = pd.col("a") == 2

    result = df.copy()
    result[expr] = 100
    expected = pd.DataFrame({"a": [1, 100], "b": [3, 100]})
    tm.assert_frame_equal(result, expected)


def test_frame_loc() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    expr = pd.col("a") == 2
    result = df.copy()
    result.loc[expr, "b"] = 100
    expected = pd.DataFrame({"a": [1, 2], "b": [3, 100]})
    tm.assert_frame_equal(result, expected)


def test_frame_iloc() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    expr = pd.col("a") == 2
    result = df.copy()
    result.iloc[expr, 1] = 100
    expected = pd.DataFrame({"a": [1, 2], "b": [3, 100]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    ("expr", "expected_values", "expected_str"),
    [
        (pd.col("a").dt.year, [2020], "col('a').dt.year"),
        (pd.col("a").dt.strftime("%B"), ["January"], "col('a').dt.strftime('%B')"),
        (pd.col("b").str.upper(), ["FOO"], "col('b').str.upper()"),
    ],
)
def test_namespaces(
    expr: Expression, expected_values: list[object], expected_str: str
) -> None:
    df = pd.DataFrame({"a": [datetime(2020, 1, 1)], "b": ["foo"]})
    result = df.assign(c=expr)
    expected = pd.DataFrame(
        {"a": [datetime(2020, 1, 1)], "b": ["foo"], "c": expected_values}
    )
    tm.assert_frame_equal(result, expected, check_dtype=False)
    assert str(expr) == expected_str


def test_invalid() -> None:
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    with pytest.raises(ValueError, match=r"did you mean one of \['a', 'b'\] instead"):
        df.assign(c=pd.col("c").mean())
    df = pd.DataFrame({f"col_{i}": [0] for i in range(11)})
    msg = (
        "did you mean one of "
        r"\['col_0', 'col_1', 'col_2', 'col_3', "
        "'col_4', 'col_5', 'col_6', 'col_7', "
        r"'col_8', 'col_9',\.\.\.\] instead"
    )
    ""
    with pytest.raises(ValueError, match=msg):
        df.assign(c=pd.col("c").mean())


def test_custom_accessor() -> None:
    df = pd.DataFrame({"a": [1, 2, 3]})

    class XYZAccessor:
        def __init__(self, pandas_obj):
            self._obj = pandas_obj

        def mean(self):
            return self._obj.mean()

    with ensure_removed(pd.Series, "xyz"):
        pd.api.extensions.register_series_accessor("xyz")(XYZAccessor)
        result = df.assign(b=pd.col("a").xyz.mean())
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [2.0, 2.0, 2.0]})
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    ("expr", "expected_values", "expected_str"),
    [
        (
            pd.col("a") & pd.col("b"),
            [False, False, True, False],
            "col('a') & col('b')",
        ),
        (
            pd.col("a") & True,
            [True, False, True, False],
            "col('a') & True",
        ),
        (
            pd.col("a") | pd.col("b"),
            [True, True, True, True],
            "col('a') | col('b')",
        ),
        (
            pd.col("a") | False,
            [True, False, True, False],
            "col('a') | False",
        ),
        (
            pd.col("a") ^ pd.col("b"),
            [True, True, False, True],
            "col('a') ^ col('b')",
        ),
        (
            pd.col("a") ^ True,
            [False, True, False, True],
            "col('a') ^ True",
        ),
        (
            ~pd.col("a"),
            [False, True, False, True],
            "~col('a')",
        ),
    ],
)
def test_col_logical_ops(
    expr: Expression, expected_values: list[bool], expected_str: str
) -> None:
    # https://github.com/pandas-dev/pandas/issues/63322
    df = pd.DataFrame({"a": [True, False, True, False], "b": [False, True, True, True]})
    result = df.assign(c=expr)
    expected = pd.DataFrame(
        {
            "a": [True, False, True, False],
            "b": [False, True, True, True],
            "c": expected_values,
        }
    )
    tm.assert_frame_equal(result, expected)
    assert str(expr) == expected_str

    # Test that the expression works with .loc
    result = df.loc[expr]
    expected = df[expected_values]
    tm.assert_frame_equal(result, expected)


def test_expression_getitem() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    df = pd.DataFrame({"a": [1, 2, 3]})
    expr = pd.col("a")[1]
    expected_str = "col('a')[1]"

    assert str(expr) == expected_str

    result = df.assign(b=expr)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [2, 2, 2]})
    tm.assert_frame_equal(result, expected)


def test_property() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    df = pd.DataFrame({"a": [1, 2, 3]})
    expr = pd.col("a").index
    expected_str = "col('a').index"

    assert str(expr) == expected_str

    result = df.assign(b=expr)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [0, 1, 2]})
    tm.assert_frame_equal(result, expected)


def test_cached_property() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    # Ensure test is valid
    assert isinstance(pd.Index.dtype, cache_readonly)

    df = pd.DataFrame({"a": [1, 2, 3]})
    expr = pd.col("a").index.dtype
    expected_str = "col('a').index.dtype"
    assert str(expr) == expected_str

    result = df.assign(b=expr)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": np.int64})
    tm.assert_frame_equal(result, expected)


def test_qcut() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    df = pd.DataFrame({"a": [1, 2, 3]})
    expr = pd.qcut(pd.col("a"), 3)
    expected_str = "qcut(x=col('a'), q=3, labels=None, retbins=False, precision=3)"
    assert str(expr) == expected_str, str(expr)

    result = df.assign(b=expr)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": pd.qcut(df["a"], 3)})
    tm.assert_frame_equal(result, expected)


def test_where() -> None:
    # https://github.com/pandas-dev/pandas/pull/63439
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    expr = pd.col("a").where(pd.col("b") == 5, 100)
    expected_str = "col('a').where(col('b') == 5, 100)"
    assert str(expr) == expected_str, str(expr)

    result = df.assign(c=expr)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [100, 2, 100]})
    tm.assert_frame_equal(result, expected)

    expr = pd.col("a").where(pd.col("b") == 5, pd.col("a") + 1)
    expected_str = "col('a').where(col('b') == 5, col('a') + 1)"
    assert str(expr) == expected_str, str(expr)

    result = df.assign(c=expr)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [2, 2, 4]})
    tm.assert_frame_equal(result, expected)
