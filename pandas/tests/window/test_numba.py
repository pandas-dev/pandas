import numpy as np
import pytest

from pandas.errors import NumbaUtilError
import pandas.util._test_decorators as td

from pandas import DataFrame, Series, option_context
import pandas._testing as tm
from pandas.core.util.numba_ import NUMBA_FUNC_CACHE


@td.skip_if_no("numba", "0.46.0")
@pytest.mark.filterwarnings("ignore:\\nThe keyword argument")
# Filter warnings when parallel=True and the function can't be parallelized by Numba
class TestRollingApply:
    @pytest.mark.parametrize("jit", [True, False])
    def test_numba_vs_cython(self, jit, nogil, parallel, nopython, center):
        def f(x, *args):
            arg_sum = 0
            for arg in args:
                arg_sum += arg
            return np.mean(x) + arg_sum

        if jit:
            import numba

            f = numba.jit(f)

        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        args = (2,)

        s = Series(range(10))
        result = s.rolling(2, center=center).apply(
            f, args=args, engine="numba", engine_kwargs=engine_kwargs, raw=True
        )
        expected = s.rolling(2, center=center).apply(
            f, engine="cython", args=args, raw=True
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("jit", [True, False])
    def test_cache(self, jit, nogil, parallel, nopython):
        # Test that the functions are cached correctly if we switch functions
        def func_1(x):
            return np.mean(x) + 4

        def func_2(x):
            return np.std(x) * 5

        if jit:
            import numba

            func_1 = numba.jit(func_1)
            func_2 = numba.jit(func_2)

        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}

        roll = Series(range(10)).rolling(2)
        result = roll.apply(
            func_1, engine="numba", engine_kwargs=engine_kwargs, raw=True
        )
        expected = roll.apply(func_1, engine="cython", raw=True)
        tm.assert_series_equal(result, expected)

        # func_1 should be in the cache now
        assert (func_1, "rolling_apply") in NUMBA_FUNC_CACHE

        result = roll.apply(
            func_2, engine="numba", engine_kwargs=engine_kwargs, raw=True
        )
        expected = roll.apply(func_2, engine="cython", raw=True)
        tm.assert_series_equal(result, expected)
        # This run should use the cached func_1
        result = roll.apply(
            func_1, engine="numba", engine_kwargs=engine_kwargs, raw=True
        )
        expected = roll.apply(func_1, engine="cython", raw=True)
        tm.assert_series_equal(result, expected)


@td.skip_if_no("numba", "0.46.0")
class TestGroupbyEWMMean:
    def test_invalid_engine(self):
        df = DataFrame({"A": ["a", "b", "a", "b"], "B": range(4)})
        with pytest.raises(ValueError, match="engine must be either"):
            df.groupby("A").ewm(com=1.0).mean(engine="foo")

    def test_invalid_engine_kwargs(self):
        df = DataFrame({"A": ["a", "b", "a", "b"], "B": range(4)})
        with pytest.raises(ValueError, match="cython engine does not"):
            df.groupby("A").ewm(com=1.0).mean(
                engine="cython", engine_kwargs={"nopython": True}
            )

    def test_cython_vs_numba(self, nogil, parallel, nopython, ignore_na, adjust):
        df = DataFrame({"A": ["a", "b", "a", "b"], "B": range(4)})
        gb_ewm = df.groupby("A").ewm(com=1.0, adjust=adjust, ignore_na=ignore_na)

        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        result = gb_ewm.mean(engine="numba", engine_kwargs=engine_kwargs)
        expected = gb_ewm.mean(engine="cython")

        tm.assert_frame_equal(result, expected)


@td.skip_if_no("numba", "0.46.0")
def test_use_global_config():
    def f(x):
        return np.mean(x) + 2

    s = Series(range(10))
    with option_context("compute.use_numba", True):
        result = s.rolling(2).apply(f, engine=None, raw=True)
    expected = s.rolling(2).apply(f, engine="numba", raw=True)
    tm.assert_series_equal(expected, result)


@td.skip_if_no("numba", "0.46.0")
def test_invalid_kwargs_nopython():
    with pytest.raises(NumbaUtilError, match="numba does not support kwargs with"):
        Series(range(1)).rolling(1).apply(
            lambda x: x, kwargs={"a": 1}, engine="numba", raw=True
        )
