import numpy as np
import pytest

import pandas as pd


@pytest.fixture
def data():
    """Fixture returning boolean array, with valid and missing values."""
    return pd.array(
        [True, False] * 2 + [np.nan] + [True, False] + [np.nan] + [True, False],
        dtype="boolean",
    )


@pytest.mark.parametrize(
    "values, exp_any, exp_all, exp_any_noskip, exp_all_noskip",
    [
        ([True, pd.NA], True, True, True, pd.NA),
        ([False, pd.NA], False, False, pd.NA, False),
        ([pd.NA], False, True, pd.NA, pd.NA),
        ([], False, True, False, True),
        # GH-33253: all True / all False values buggy with skipna=False
        ([True, True], True, True, True, True),
        ([False, False], False, False, False, False),
    ],
)
@pytest.mark.parametrize("con", [pd.array, pd.Series])
def test_any_all(
    values, exp_any, exp_all, exp_any_noskip, exp_all_noskip, using_python_scalars, con
):
    # the methods return numpy scalars
    if not using_python_scalars or con is pd.array:
        exp_any = pd.NA if exp_any is pd.NA else np.bool_(exp_any)
        exp_all = pd.NA if exp_all is pd.NA else np.bool_(exp_all)
        exp_any_noskip = pd.NA if exp_any_noskip is pd.NA else np.bool_(exp_any_noskip)
        exp_all_noskip = pd.NA if exp_all_noskip is pd.NA else np.bool_(exp_all_noskip)

    a = con(values, dtype="boolean")
    assert a.any() is exp_any
    assert a.all() is exp_all
    assert a.any(skipna=False) is exp_any_noskip
    assert a.all(skipna=False) is exp_all_noskip


def test_reductions_return_types(
    dropna, data, all_numeric_reductions, using_python_scalars
):
    op = all_numeric_reductions
    s = pd.Series(data)
    if dropna:
        s = s.dropna()

    if using_python_scalars:
        expected = {
            "sum": int,
            "prod": int,
            "count": int,
            "min": bool,
            "max": bool,
        }.get(op, float)
    else:
        expected = {
            "sum": np.int_,
            "prod": np.int_,
            "count": np.integer,
            "min": np.bool_,
            "max": np.bool_,
        }.get(op, np.float64)
    result = getattr(s, op)()
    assert isinstance(result, expected), f"{type(result)} vs {expected}"
