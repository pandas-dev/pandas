import operator

import numpy as np
import pytest

from pandas.compat import (
    WASM,
)

from pandas.core.dtypes.common import is_number

import pandas as pd
import pandas._testing as tm
from pandas.tests.apply.common import (
    frame_transform_kernels,
    series_transform_kernels,
)


@pytest.mark.parametrize("func", ["sum", "mean", "min", "max", "std"])
@pytest.mark.parametrize(
    "kwds",
    [
        pytest.param({}, id="no_kwds"),
        pytest.param({"axis": 1}, id="on_axis"),
        pytest.param({"numeric_only": True}, id="func_kwds"),
        pytest.param({"axis": 1, "numeric_only": True}, id="axis_and_func_kwds"),
    ],
)
@pytest.mark.parametrize("how", ["agg", "apply"])
def test_apply_with_string_funcs(float_frame, func, kwds, how):
    result = getattr(float_frame, how)(func, **kwds)
    expected = getattr(float_frame, func)(**kwds)
    tm.assert_series_equal(result, expected)


def test_with_string_args(datetime_series, all_numeric_reductions):
    result = datetime_series.apply(all_numeric_reductions)
    expected = getattr(datetime_series, all_numeric_reductions)()
    assert result == expected


@pytest.mark.parametrize("op", ["mean", "median", "std", "var"])
@pytest.mark.parametrize("how", ["agg", "apply"])
def test_apply_np_reducer(op, how):
    # GH 39116
    float_frame = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    result = getattr(float_frame, how)(op)
    # pandas ddof defaults to 1, numpy to 0
    kwargs = {"ddof": 1} if op in ("std", "var") else {}
    expected = pd.Series(
        getattr(np, op)(float_frame, axis=0, **kwargs), index=float_frame.columns
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.skipif(WASM, reason="No fp exception support in wasm")
@pytest.mark.parametrize(
    "op", ["abs", "ceil", "cos", "cumsum", "exp", "log", "sqrt", "square"]
)
@pytest.mark.parametrize("how", ["transform", "apply"])
def test_apply_np_transformer(float_frame, op, how):
    # GH 39116

    # float_frame will _usually_ have negative values, which will
    #  trigger the warning here, but let's put one in just to be sure
    float_frame.iloc[0, 0] = -1.0
    warn = None
    if op in ["log", "sqrt"]:
        warn = RuntimeWarning

    msg = f"invalid value encountered in {op}"
    with tm.assert_produces_warning(warn, check_stacklevel=False, match=msg):
        # float_frame fixture is defined in conftest.py, so we don't check the
        # stacklevel as otherwise the test would fail.
        result = getattr(float_frame, how)(op)
        expected = getattr(np, op)(float_frame)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "series, func, expected",
    [
        (pd.Series(dtype=np.float64), "sum", 0),
        (pd.Series(dtype=np.float64), "max", np.nan),
        (pd.Series(dtype=np.float64), "min", np.nan),
        (pd.Series(dtype=np.float64), "all", True),
        (pd.Series(dtype=np.float64), "any", False),
        (pd.Series(dtype=np.float64), "mean", np.nan),
        (pd.Series(dtype=np.float64), "prod", 1),
        (pd.Series(dtype=np.float64), "std", np.nan),
        (pd.Series(dtype=np.float64), "var", np.nan),
        (pd.Series(dtype=np.float64), "median", np.nan),
        (pd.Series([np.nan, 1, 2, 3]), "sum", 6),
        (pd.Series([np.nan, 1, 2, 3]), "max", 3),
        (pd.Series([np.nan, 1, 2, 3]), "min", 1),
        (pd.Series([np.nan, 1, 2, 3]), "all", True),
        (pd.Series([np.nan, 1, 2, 3]), "any", True),
        (pd.Series([np.nan, 1, 2, 3]), "mean", 2),
        (pd.Series([np.nan, 1, 2, 3]), "prod", 6),
        (pd.Series([np.nan, 1, 2, 3]), "std", 1),
        (pd.Series([np.nan, 1, 2, 3]), "var", 1),
        (pd.Series([np.nan, 1, 2, 3]), "median", 2),
        (pd.Series("a b c".split()), "sum", "abc"),
        (pd.Series("a b c".split()), "max", "c"),
        (pd.Series("a b c".split()), "min", "a"),
        (pd.Series("a b c".split()), "all", True),
        (pd.Series("a b c".split()), "any", True),
    ],
)
def test_agg_reduction_series(series, func, expected):
    # GH21224
    result = series.agg(func)
    if is_number(expected):
        assert np.isclose(result, expected, equal_nan=True)
    else:
        assert result == expected


@pytest.mark.parametrize(
    "series, func, expected",
    [
        (pd.Series(dtype=np.float64), "cumprod", pd.Series([], dtype=np.float64)),
        (pd.Series(dtype=np.float64), "cumsum", pd.Series([], dtype=np.float64)),
        (pd.Series([np.nan, 1, 2, 3]), "cumprod", pd.Series([np.nan, 1, 2, 6])),
        (pd.Series([np.nan, 1, 2, 3]), "cumsum", pd.Series([np.nan, 1, 3, 6])),
        (pd.Series("a b c".split()), "cumsum", pd.Series(["a", "ab", "abc"])),
    ],
)
def test_agg_transform_series(series, func, expected):
    # GH21224
    result = series.agg(func)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "df, func, expected",
    [
        (pd.DataFrame(), "sum", pd.Series(dtype="float64")),
        (pd.DataFrame(), "max", pd.Series(dtype="float64")),
        (pd.DataFrame(), "min", pd.Series(dtype="float64")),
        (pd.DataFrame(), "all", pd.Series(dtype=bool)),
        (pd.DataFrame(), "any", pd.Series(dtype=bool)),
        (pd.DataFrame(), "mean", pd.Series(dtype="float64")),
        (pd.DataFrame(), "prod", pd.Series(dtype="float64")),
        (pd.DataFrame(), "std", pd.Series(dtype="float64")),
        (pd.DataFrame(), "var", pd.Series(dtype="float64")),
        (pd.DataFrame(), "median", pd.Series(dtype="float64")),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "sum", pd.Series([1.0, 3])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "max", pd.Series([1.0, 2])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "min", pd.Series([1.0, 1])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "all", pd.Series([True, True])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "any", pd.Series([True, True])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "mean", pd.Series([1, 1.5])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "prod", pd.Series([1.0, 2])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "std", pd.Series([np.nan, 0.707107])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "var", pd.Series([np.nan, 0.5])),
        (pd.DataFrame([[np.nan, 1], [1, 2]]), "median", pd.Series([1, 1.5])),
    ],
)
def test_agg_reduction_frame(df, func, expected, axis):
    # GH 21224
    result = df.agg(func, axis=axis)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "df, func, expected",
    [
        (pd.DataFrame(), "cumprod", pd.DataFrame()),
        (pd.DataFrame(), "cumsum", pd.DataFrame()),
        (
            pd.DataFrame([[np.nan, 1], [1, 2]]),
            "cumprod",
            pd.DataFrame([[np.nan, 1], [1, 2]]),
        ),
        (
            pd.DataFrame([[np.nan, 1], [1, 2]]),
            "cumsum",
            pd.DataFrame([[np.nan, 1], [1, 3]]),
        ),
    ],
)
def test_agg_transform_frame(df, func, expected, axis):
    # GH 21224
    if axis in ("columns", 1):
        # operating blockwise doesn't let us preserve dtypes
        expected = expected.astype("float64")

    result = df.agg(func, axis=axis)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("op", series_transform_kernels)
def test_transform_groupby_kernel_series(request, string_series, op):
    # GH 35964
    if op == "ngroup":
        request.applymarker(
            pytest.mark.xfail(raises=ValueError, reason="ngroup not valid for NDFrame")
        )
    ones = np.ones(string_series.shape[0])
    expected = string_series.groupby(ones).transform(op)
    result = string_series.transform(op, 0)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("op", frame_transform_kernels)
def test_transform_groupby_kernel_frame(request, float_frame, op):
    if op == "ngroup":
        request.applymarker(
            pytest.mark.xfail(raises=ValueError, reason="ngroup not valid for NDFrame")
        )

    # GH 35964

    ones = np.ones(float_frame.shape[0])
    gb = float_frame.groupby(ones)
    expected = gb.transform(op)

    result = float_frame.transform(op, 0)
    tm.assert_frame_equal(result, expected)

    # same thing, but ensuring we have multiple blocks
    assert "E" not in float_frame.columns
    float_frame["E"] = float_frame["A"].copy()
    assert len(float_frame._mgr.blocks) > 1

    ones = np.ones(float_frame.shape[0])
    gb2 = float_frame.groupby(ones)
    expected2 = gb2.transform(op)
    result2 = float_frame.transform(op, 0)
    tm.assert_frame_equal(result2, expected2)


@pytest.mark.parametrize("method", ["abs", "shift", "pct_change", "cumsum", "rank"])
def test_transform_method_name(method):
    # GH 19760
    df = pd.DataFrame({"A": [-1, 2]})
    result = df.transform(method)
    expected = operator.methodcaller(method)(df)
    tm.assert_frame_equal(result, expected)
