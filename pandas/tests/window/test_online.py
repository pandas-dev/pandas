import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


@td.skip_if_no("numba", "0.46.0")
@pytest.mark.filterwarnings("ignore:\\nThe keyword argument")
class TestEWM:
    def test_online_vs_non_online(self, nogil, parallel, nopython):
        df = DataFrame({"a": range(5), "b": range(5)})
        expected = df.ewm(0.5).mean()
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        online_ewm = df.head(2).ewm(0.5).online(engine_kwargs=engine_kwargs)
        result = online_ewm.mean()
        tm.assert_frame_equal(result, expected.head(2))

        result = online_ewm.update(update=df.tail(3))
        tm.assert_frame_equal(result, expected.tail(3))

    def test_update_times(self, nogil, parallel, nopython):
        times = Series(
            np.array(
                ["2020-01-01", "2020-01-02", "2020-01-04", "2020-01-17", "2020-01-21"],
                dtype="datetime64",
            )
        )
        df = DataFrame({"a": range(5), "b": range(5)})
        expected = df.ewm(0.5, times=times).mean()
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        online_ewm = (
            df.head(2).ewm(0.5, times=times.head(2)).online(engine_kwargs=engine_kwargs)
        )
        result = online_ewm.mean()
        tm.assert_frame_equal(result, expected.head(2))

        result = online_ewm.update(update=df.tail(3), update_times=times.tail(3))
        tm.assert_frame_equal(result, expected.tail(3))
