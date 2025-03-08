import numpy as np
import pytest

from pandas.errors import ExecutionError

import pandas as pd
import pandas._testing as tm

pytestmark = [pytest.mark.single_cpu, pytest.mark.bodo_udf_engine]


@pytest.fixture(params=["bodo"])
def engine(request):
    """Test bodo engine by itself to avoid extensions conflicting with numba.

    Note: Using a fixture here to avoid importing at the start of the session.
    """
    if request.param == "bodo":
        pytest.importorskip("bodo")
    return request.param


def test_bodo_vs_python_indexing(engine):
    frame = pd.DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": [7.0, 8.0, 9.0]},
    )

    def f(a):
        return a["c"]

    result = frame.apply(f, engine="bodo", axis=1)
    expected = frame.apply(f, engine="python", axis=1)

    tm.assert_series_equal(result, expected, check_series_type=False)


@pytest.mark.parametrize(
    "reduction",
    [lambda x: x.mean(), lambda x: x.min(), lambda x: x.max(), lambda x: x.sum()],
)
def test_bodo_vs_python_reductions(reduction, engine):
    df = pd.DataFrame(np.ones((4, 4), dtype=np.float64))
    result = df.apply(reduction, engine="bodo", axis=1)
    expected = df.apply(reduction, engine="python", axis=1)
    tm.assert_series_equal(result, expected, check_series_type=False)


def test_bodo_vs_python_df_output(engine):
    df = pd.DataFrame({"A": np.arange(20), "B": ["hi", "there"] * 10})

    def f(a):
        return pd.Series([a["B"], a["A"]])

    result = df.apply(f, engine="bodo", axis=1)
    expected = df.apply(f, engine="python", axis=1)

    tm.assert_frame_equal(result, expected, check_frame_type=False, check_dtype=False)


def test_bodo_vs_python_args(engine):
    msg = (
        "the 'bodo' engine does not support passing additional args/kwargs "
        "to apply function yet."
    )

    def f(x, y):
        return x.A + y

    df = pd.DataFrame({"A": np.arange(20)})

    with pytest.raises(NotImplementedError, match=msg):
        df.apply(f, engine="bodo", axis=1, args=(2,))

    with pytest.raises(NotImplementedError, match=msg):
        df.apply(f, engine="bodo", axis=1, y=2)


@pytest.mark.parametrize("axis", [0, 1])
def test_bodo_vs_python_str_apply(axis, engine):
    df = pd.DataFrame({"A": np.arange(20)})

    func = "mean"
    axis = 1
    result = df.apply(func, axis, engine="bodo")
    expected = df.apply(func, axis)

    tm.assert_series_equal(result, expected, check_series_type=False)


def test_bodo_unsupported_axis(engine):
    """Tests that a BodoError is raised when trying to apply UDF column-wise"""
    frame = pd.DataFrame(
        {"a": [1, 2, 3]},
    )

    def f(a):
        return 1

    with pytest.raises(
        NotImplementedError,
        match=r"the 'bodo' engine only supports axis=1 for user-defined functions",
    ):
        frame.apply(f, engine="bodo", axis=0)


def test_bodo_raw_unsupported(engine):
    """Tests that error gets raised when using raw=True"""
    frame = pd.DataFrame(
        {"a": [1, 2, 3]},
    )

    def f(a):
        return 1

    with pytest.raises(
        NotImplementedError, match="the 'bodo' engine does not support raw=True."
    ):
        frame.apply(f, engine="bodo", raw=True, axis=1)


def test_bodo_result_type_unsupported(engine):
    """Tests that error gets raised when passing any value to result_type"""
    frame = pd.DataFrame(
        {"a": [1, 2, 3]},
    )

    def f(a):
        return 1

    with pytest.raises(
        NotImplementedError, match="the 'bodo' engine does not support result_type yet."
    ):
        frame.apply(f, engine="bodo", axis=1, result_type="reduce")


def test_bodo_engine_execution_error(engine):
    frame = pd.DataFrame(
        {"a": [1, 2, 3], "b": ["1", "2", "3"]},
    )

    def f(x):
        return x.a + x.b

    with pytest.raises(ExecutionError, match="Execution with engine='bodo' failed."):
        frame.apply(f, engine="bodo", axis=1)
