from datetime import datetime

import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.test_register_accessor import ensure_removed


@pytest.mark.parametrize(
    ("expr", "expected_values", "expected_str"),
    [
        (pd.col("a"), [1, 2], "col('a')"),
        (pd.col("a") * 2, [2, 4], "(col('a') * 2)"),
        (pd.col("a").sum(), [3, 3], "col('a').sum()"),
        (pd.col("a") + 1, [2, 3], "(col('a') + 1)"),
        (1 + pd.col("a"), [2, 3], "(1 + col('a'))"),
        (pd.col("a") - 1, [0, 1], "(col('a') - 1)"),
        (1 - pd.col("a"), [0, -1], "(1 - col('a'))"),
        (pd.col("a") * 1, [1, 2], "(col('a') * 1)"),
        (1 * pd.col("a"), [1, 2], "(1 * col('a'))"),
        (pd.col("a") / 1, [1.0, 2.0], "(col('a') / 1)"),
        (1 / pd.col("a"), [1.0, 0.5], "(1 / col('a'))"),
        (pd.col("a") // 1, [1, 2], "(col('a') // 1)"),
        (1 // pd.col("a"), [1, 0], "(1 // col('a'))"),
        (pd.col("a") % 1, [0, 0], "(col('a') % 1)"),
        (1 % pd.col("a"), [0, 1], "(1 % col('a'))"),
        (pd.col("a") > 1, [False, True], "(col('a') > 1)"),
        (pd.col("a") >= 1, [True, True], "(col('a') >= 1)"),
        (pd.col("a") < 1, [False, False], "(col('a') < 1)"),
        (pd.col("a") <= 1, [True, False], "(col('a') <= 1)"),
        (pd.col("a") == 1, [True, False], "(col('a') == 1)"),
    ],
)
def test_col_simple(
    expr: pd.Expr, expected_values: list[object], expected_str: str
) -> None:
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    result = df.assign(c=expr)
    expected = pd.DataFrame({"a": [1, 2], "b": [3, 4], "c": expected_values})
    tm.assert_frame_equal(result, expected)
    assert str(expr) == expected_str


@pytest.mark.parametrize(
    ("expr", "expected_values", "expected_str"),
    [
        (pd.col("a").dt.year, [2020], "col('a').dt.year"),
        (pd.col("a").dt.strftime("%B"), ["January"], "col('a').dt.strftime('%B')"),
        (pd.col("b").str.upper(), ["FOO"], "col('b').str.upper()"),
    ],
)
def test_namespaces(
    expr: pd.Expr, expected_values: list[object], expected_str: str
) -> None:
    df = pd.DataFrame({"a": [datetime(2020, 1, 1)], "b": ["foo"]})
    result = df.assign(c=expr)
    expected = pd.DataFrame(
        {"a": [datetime(2020, 1, 1)], "b": ["foo"], "c": expected_values}
    )
    tm.assert_frame_equal(result, expected, check_dtype=False)
    assert str(expr) == expected_str


def test_invalid() -> None:
    df = pd.DataFrame({"a": [1, 2]})
    with pytest.raises(ValueError, match="did you mean"):
        df.assign(c=pd.col("b").mean())


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
