import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame
import pandas._testing as tm


@td.skip_if_no("numba", "0.46.0")
def test_correct_function_signature():
    def incorrect_function(x):
        return x + 1

    data = DataFrame(
        {"key": ["a", "a", "b", "b", "a"], "data": [1.0, 2.0, 3.0, 4.0, 5.0]},
        columns=["key", "data"],
    )
    with pytest.raises(ValueError, match=f"The first 3"):
        data.groupby("key").transform(incorrect_function, engine="numba")

    with pytest.raises(ValueError, match=f"The first 2"):
        data.groupby("key")["data"].transform(incorrect_function, engine="numba")


@td.skip_if_no("numba", "0.46.0")
def test_check_nopython_kwargs():
    def incorrect_function(x, **kwargs):
        return x + 1

    data = DataFrame(
        {"key": ["a", "a", "b", "b", "a"], "data": [1.0, 2.0, 3.0, 4.0, 5.0]},
        columns=["key", "data"],
    )
    with pytest.raises(ValueError, match="numba does not support"):
        data.groupby("key").transform(incorrect_function, engine="numba", a=1)

    with pytest.raises(ValueError, match="numba does not support"):
        data.groupby("key")["data"].transform(incorrect_function, engine="numba", a=1)


@td.skip_if_no("numba", "0.46.0")
@pytest.mark.filterwarnings("ignore:\\nThe keyword argument")
# Filter warnings when parallel=True and the function can't be parallelized by Numba
@pytest.mark.parametrize("jit", [True, False])
def test_numba_vs_cython(jit, nogil, parallel, nopython):
    def series_func(values, index):
        return values + 1

    def dataframe_func(values, index, columns):
        return values + 1

    if jit:
        # Test accepted jitted functions
        import numba

        series_func = numba.jit(series_func)
        dataframe_func = numba.jit(dataframe_func)

    data = DataFrame(
        {0: ["a", "a", "b", "b", "a"], 1: [1.0, 2.0, 3.0, 4.0, 5.0]}, columns=[0, 1],
    )
    engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
    grouped = data.groupby(0)

    result = grouped.transform(
        dataframe_func, engine="numba", engine_kwargs=engine_kwargs
    )
    expected = grouped.transform(lambda x: x + 1)

    tm.assert_frame_equal(result, expected)

    result = grouped[1].transform(
        series_func, engine="numba", engine_kwargs=engine_kwargs
    )
    expected = grouped[1].transform(lambda x: x + 1)

    tm.assert_series_equal(result, expected)
