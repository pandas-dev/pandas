import operator

import numpy as np
import pytest

from pandas.compat import (
    WASM,
)

from pandas.core.dtypes.common import is_number

from pandas import (
    DataFrame,
    Series,
)
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
    float_frame = DataFrame({"a": [1, 2], "b": [3, 4]})
    result = getattr(float_frame, how)(op)
    # pandas ddof defaults to 1, numpy to 0
    kwargs = {"ddof": 1} if op in ("std", "var") else {}
    expected = Series(
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

    with tm.assert_produces_warning(warn, check_stacklevel=False):
        # float_frame fixture is defined in conftest.py, so we don't check the
        # stacklevel as otherwise the test would fail.
        result = getattr(float_frame, how)(op)
        expected = getattr(np, op)(float_frame)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "series, func, expected",
    [
        (Series(dtype=np.float64), "sum", 0),
        (Series(dtype=np.float64), "max", np.nan),
        (Series(dtype=np.float64), "min", np.nan),
        (Series(dtype=np.float64), "all", True),
        (Series(dtype=np.float64), "any", False),
        (Series(dtype=np.float64), "mean", np.nan),
        (Series(dtype=np.float64), "prod", 1),
        (Series(dtype=np.float64), "std", np.nan),
        (Series(dtype=np.float64), "var", np.nan),
        (Series(dtype=np.float64), "median", np.nan),
        (Series([np.nan, 1, 2, 3]), "sum", 6),
        (Series([np.nan, 1, 2, 3]), "max", 3),
        (Series([np.nan, 1, 2, 3]), "min", 1),
        (Series([np.nan, 1, 2, 3]), "all", True),
        (Series([np.nan, 1, 2, 3]), "any", True),
        (Series([np.nan, 1, 2, 3]), "mean", 2),
        (Series([np.nan, 1, 2, 3]), "prod", 6),
        (Series([np.nan, 1, 2, 3]), "std", 1),
        (Series([np.nan, 1, 2, 3]), "var", 1),
        (Series([np.nan, 1, 2, 3]), "median", 2),
        (Series("a b c".split()), "sum", "abc"),
        (Series("a b c".split()), "max", "c"),
        (Series("a b c".split()), "min", "a"),
        (Series("a b c".split()), "all", True),
        (Series("a b c".split()), "any", True),
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
        (Series(dtype=np.float64), "cumprod", Series([], dtype=np.float64)),
        (Series(dtype=np.float64), "cumsum", Series([], dtype=np.float64)),
        (Series([np.nan, 1, 2, 3]), "cumprod", Series([np.nan, 1, 2, 6])),
        (Series([np.nan, 1, 2, 3]), "cumsum", Series([np.nan, 1, 3, 6])),
        (Series("a b c".split()), "cumsum", Series(["a", "ab", "abc"])),
    ],
)
def test_agg_transform_series(series, func, expected):
    # GH21224
    result = series.agg(func)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "df, func, expected",
    [
        (DataFrame(), "sum", Series(dtype="float64")),
        (DataFrame(), "max", Series(dtype="float64")),
        (DataFrame(), "min", Series(dtype="float64")),
        (DataFrame(), "all", Series(dtype=bool)),
        (DataFrame(), "any", Series(dtype=bool)),
        (DataFrame(), "mean", Series(dtype="float64")),
        (DataFrame(), "prod", Series(dtype="float64")),
        (DataFrame(), "std", Series(dtype="float64")),
        (DataFrame(), "var", Series(dtype="float64")),
        (DataFrame(), "median", Series(dtype="float64")),
        (DataFrame([[np.nan, 1], [1, 2]]), "sum", Series([1.0, 3])),
        (DataFrame([[np.nan, 1], [1, 2]]), "max", Series([1.0, 2])),
        (DataFrame([[np.nan, 1], [1, 2]]), "min", Series([1.0, 1])),
        (DataFrame([[np.nan, 1], [1, 2]]), "all", Series([True, True])),
        (DataFrame([[np.nan, 1], [1, 2]]), "any", Series([True, True])),
        (DataFrame([[np.nan, 1], [1, 2]]), "mean", Series([1, 1.5])),
        (DataFrame([[np.nan, 1], [1, 2]]), "prod", Series([1.0, 2])),
        (DataFrame([[np.nan, 1], [1, 2]]), "std", Series([np.nan, 0.707107])),
        (DataFrame([[np.nan, 1], [1, 2]]), "var", Series([np.nan, 0.5])),
        (DataFrame([[np.nan, 1], [1, 2]]), "median", Series([1, 1.5])),
    ],
)
def test_agg_reduction_frame(df, func, expected, axis):
    # GH 21224
    result = df.agg(func, axis=axis)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "df, func, expected",
    [
        (DataFrame(), "cumprod", DataFrame()),
        (DataFrame(), "cumsum", DataFrame()),
        (
            DataFrame([[np.nan, 1], [1, 2]]),
            "cumprod",
            DataFrame([[np.nan, 1], [1, 2]]),
        ),
        (DataFrame([[np.nan, 1], [1, 2]]), "cumsum", DataFrame([[np.nan, 1], [1, 3]])),
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
    df = DataFrame({"A": [-1, 2]})
    result = df.transform(method)
    expected = operator.methodcaller(method)(df)
    tm.assert_frame_equal(result, expected)
