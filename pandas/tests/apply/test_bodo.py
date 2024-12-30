import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm

pytestmark = [td.skip_if_no("bodo")]


def test_bodo_vs_python_indexing():
    frame = pd.DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": [7.0, 8.0, 9.0]},
    )
    f = lambda x: x["c"]
    result = frame.apply(f, engine="bodo", axis=1)
    expected = frame.apply(f, engine="python", axis=1)

    tm.assert_series_equal(result, expected, check_series_type=False)


@pytest.mark.parametrize(
    "reduction",
    [lambda x: x.mean(), lambda x: x.min(), lambda x: x.max(), lambda x: x.sum()],
)
def test_bodo_vs_python_reductions(reduction):
    df = pd.DataFrame(np.ones((4, 4), dtype=np.float64))
    result = df.apply(reduction, engine="bodo", axis=1)
    expected = df.apply(reduction, engine="python", axis=1)
    tm.assert_series_equal(result, expected, check_series_type=False)


def test_bodo_vs_python_df_output():
    df = pd.DataFrame({"A": np.arange(20), "B": ["hi", "there"] * 10})

    f = lambda a: pd.Series([a["B"], a["A"]])
    result = df.apply(f, engine="bodo", axis=1)
    expected = df.apply(f, engine="python", axis=1)

    tm.assert_frame_equal(result, expected, check_frame_type=False, check_dtype=False)


@pytest.mark.skip(reason="TODO: pass args/kwargs to bodo jitted function")
def test_bodo_vs_python_args_kwargs():
    def f(x, y, z=3):
        return x.A == y + z

    df = pd.DataFrame({"A": np.arange(20)})

    result = df.apply(f, z=2, engine="bodo", axis=1, args=(2,))
    expected = df.apply(f, z=2, axis=1, args=(2,))
    tm.assert_series_equal(result, expected, check_series_type=False)


@pytest.mark.parametrize("axis", [0, 1])
def test_bodo_vs_python_str_apply(axis):
    df = pd.DataFrame({"A": np.arange(20)})

    func = "mean"
    axis = 1
    result = df.apply(func, axis)
    expected = df.apply(func, axis)

    tm.assert_series_equal(result, expected, check_series_type=False)


def test_bodo_unsupported_axis():
    """Tests that a BodoError is raised when trying to apply UDF column-wise"""
    frame = pd.DataFrame(
        {"a": [1, 2, 3]},
    )
    f = lambda x: 1

    with pytest.raises(
        NotImplementedError,
        match=r"the 'bodo' engine only supports axis=1 for user-defined functions",
    ):
        frame.apply(f, engine="bodo", axis=0)


def test_bodo_raw_unsupported():
    """Tests that error gets raised when using raw=True"""
    frame = pd.DataFrame(
        {"a": [1, 2, 3]},
    )
    f = lambda a: 1

    with pytest.raises(
        NotImplementedError, match="the 'bodo' engine does not support raw=True."
    ):
        frame.apply(f, engine="bodo", raw=True, axis=1)


def test_bodo_result_type_unsupported():
    """Tests that error gets raised when passing any value to result_type"""
    frame = pd.DataFrame(
        {"a": [1, 2, 3]},
    )
    f = lambda a: 1

    with pytest.raises(
        NotImplementedError, match="the 'bodo' engine does not support result_type yet."
    ):
        frame.apply(f, engine="bodo", axis=1, result_type="reduce")
