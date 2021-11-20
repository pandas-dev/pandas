import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


@td.skip_if_no("numba")
@pytest.mark.filterwarnings("ignore:\n")
# Filter warnings when parallel=True and the function can't be parallelized by Numba
class TestEngine:
    def test_cython_vs_numba_frame(self, sort, nogil, parallel, nopython):
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        result = df.groupby("a", sort=sort).mean(
            engine="numba", engine_kwargs=engine_kwargs
        )
        expected = df.groupby("a", sort=sort).mean()
        tm.assert_frame_equal(result, expected)

    def test_cython_vs_numba_getitem(self, sort, nogil, parallel, nopython):
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        result = df.groupby("a", sort=sort)["c"].mean(
            engine="numba", engine_kwargs=engine_kwargs
        )
        expected = df.groupby("a", sort=sort)["c"].mean()
        tm.assert_series_equal(result, expected)

    def test_cython_vs_numba_series(self, sort, nogil, parallel, nopython):
        ser = Series(range(3), index=[1, 2, 1], name="foo")
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        result = ser.groupby(level=0, sort=sort).mean(
            engine="numba", engine_kwargs=engine_kwargs
        )
        expected = ser.groupby(level=0, sort=sort).mean()
        tm.assert_series_equal(result, expected)

    def test_as_index_false_unsupported(self):
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        with pytest.raises(NotImplementedError, match="as_index=False"):
            df.groupby("a", as_index=False).mean(engine="numba")

    def test_axis_1_unsupported(self):
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        with pytest.raises(NotImplementedError, match="axis=1"):
            df.groupby("a", axis=1).mean(engine="numba")
