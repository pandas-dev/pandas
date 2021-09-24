import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


@td.skip_if_no("numba")
@pytest.mark.filterwarnings("ignore:\\nThe keyword argument")
# Filter warnings when parallel=True and the function can't be parallelized by Numba
class TestEngine:
    def test_cython_vs_numba_frame(self, nogil, parallel, nopython):
        df = DataFrame({"a": [1, 2, 1, 2], "b": range(4), "c": range(1, 5)})
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        result = df.groupby("a").mean(engine="numba", engine_kwargs=engine_kwargs)
        expected = df.groupby("a").mean()
        tm.assert_frame_equal(result, expected)

    def test_cython_vs_numba_getitem(self, nogil, parallel, nopython):
        df = DataFrame({"a": [1, 2, 1, 2], "b": range(4), "c": range(1, 5)})
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        result = df.groupby("a")["c"].mean(engine="numba", engine_kwargs=engine_kwargs)
        expected = df.groupby("a")["c"].mean()
        tm.assert_series_equal(result, expected)

    def test_cython_vs_numba_series(self, nogil, parallel, nopython):
        ser = Series(range(3), index=[1, 2, 1], name="foo")
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        result = ser.groupby(level=0).mean(engine="numba", engine_kwargs=engine_kwargs)
        expected = ser.groupby(level=0).mean()
        tm.assert_series_equal(result, expected)
