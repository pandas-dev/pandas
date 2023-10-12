import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Index,
)
import pandas._testing as tm

pytestmark = td.skip_if_no("numba")


def test_numba_vs_python_noop(float_frame, apply_axis):
    func = lambda x: x
    result = float_frame.apply(func, engine="numba", axis=apply_axis)
    expected = float_frame.apply(func, engine="python", axis=apply_axis)
    tm.assert_frame_equal(result, expected)


def test_numba_vs_python_indexing(float_frame):
    row_func = lambda x: x["A"]
    result = float_frame.apply(row_func, engine="numba", axis=1)
    expected = float_frame.apply(row_func, engine="python", axis=1)
    tm.assert_series_equal(result, expected)

    row_func = lambda x: x["ZqgszYBfuL"]  # This is a label in the index
    result = float_frame.apply(row_func, engine="numba", axis=0)
    expected = float_frame.apply(row_func, engine="python", axis=0)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "reduction",
    [lambda x: x.mean(), lambda x: x.min(), lambda x: x.max(), lambda x: x.sum()],
)
def test_numba_vs_python_reductions(float_frame, reduction, apply_axis):
    result = float_frame.apply(reduction, engine="numba", axis=apply_axis)
    expected = float_frame.apply(reduction, engine="python", axis=apply_axis)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("colnames", [[1, 2, 3], [1.0, 2.0, 3.0]])
def test_numba_numeric_colnames(colnames):
    # Check that numeric column names lower properly and can be indxed on
    df = DataFrame(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), columns=colnames)
    first_col = colnames[0]
    f = lambda x: x[first_col]  # Get the first column
    result = df.apply(f, engine="numba", axis=1)
    expected = df.apply(f, engine="python", axis=1)
    tm.assert_series_equal(result, expected)


def test_numba_parallel_unsupported(float_frame):
    f = lambda x: x
    with pytest.raises(
        NotImplementedError,
        match="Parallel apply is not supported when raw=False and engine='numba'",
    ):
        float_frame.apply(f, engine="numba", engine_kwargs={"parallel": True})


def test_numba_nonunique_unsupported(apply_axis):
    f = lambda x: x
    df = DataFrame({"a": [1, 2]}, index=Index(["a", "a"]))
    with pytest.raises(
        NotImplementedError,
        match="The index/columns must be unique when raw=False and engine='numba'",
    ):
        df.apply(f, engine="numba", axis=apply_axis)
