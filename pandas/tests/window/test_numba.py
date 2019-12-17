import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import Series
import pandas.util.testing as tm


@td.skip_if_no("numba", "0.46.0")
class TestApply:
    @pytest.mark.parametrize("nogil", [True, False])
    @pytest.mark.parametrize("parallel", [True, False])
    @pytest.mark.parametrize("nopython", [True, False])
    def test_numba_vs_cython(self, nogil, parallel, nopython):
        def f(x, *args):
            arg_sum = 0
            for arg in args:
                arg_sum += arg
            return np.mean(x) + arg_sum

        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        args = (2,)

        s = Series(range(10))
        result = s.rolling(2).apply(
            f, args=args, engine="numba", engine_kwargs=engine_kwargs, raw=True
        )
        expected = s.rolling(2).apply(f, engine="cython", args=args, raw=True)
        tm.assert_series_equal(result, expected)
