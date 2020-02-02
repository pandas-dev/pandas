import re

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm

from pandas.tseries.offsets import BDay


class TestXS:
    def test_xs(self, float_frame, datetime_frame):
        idx = float_frame.index[5]
        xs = float_frame.xs(idx)
        for item, value in xs.items():
            if np.isnan(value):
                assert np.isnan(float_frame[item][idx])
            else:
                assert value == float_frame[item][idx]

        # mixed-type xs
        test_data = {"A": {"1": 1, "2": 2}, "B": {"1": "1", "2": "2", "3": "3"}}
        frame = DataFrame(test_data)
        xs = frame.xs("1")
        assert xs.dtype == np.object_
        assert xs["A"] == 1
        assert xs["B"] == "1"

        with pytest.raises(
            KeyError, match=re.escape("Timestamp('1999-12-31 00:00:00', freq='B')")
        ):
            datetime_frame.xs(datetime_frame.index[0] - BDay())

        # xs get column
        series = float_frame.xs("A", axis=1)
        expected = float_frame["A"]
        tm.assert_series_equal(series, expected)

        # view is returned if possible
        series = float_frame.xs("A", axis=1)
        series[:] = 5
        assert (expected == 5).all()

    def test_xs_corner(self):
        # pathological mixed-type reordering case
        df = DataFrame(index=[0])
        df["A"] = 1.0
        df["B"] = "foo"
        df["C"] = 2.0
        df["D"] = "bar"
        df["E"] = 3.0

        xs = df.xs(0)
        exp = pd.Series([1.0, "foo", 2.0, "bar", 3.0], index=list("ABCDE"), name=0)
        tm.assert_series_equal(xs, exp)

        # no columns but Index(dtype=object)
        df = DataFrame(index=["a", "b", "c"])
        result = df.xs("a")
        expected = Series([], name="a", index=pd.Index([]), dtype=np.float64)
        tm.assert_series_equal(result, expected)

    def test_xs_duplicates(self):
        df = DataFrame(np.random.randn(5, 2), index=["b", "b", "c", "b", "a"])

        cross = df.xs("c")
        exp = df.iloc[2]
        tm.assert_series_equal(cross, exp)

    def test_xs_keep_level(self):
        df = DataFrame(
            {
                "day": {0: "sat", 1: "sun"},
                "flavour": {0: "strawberry", 1: "strawberry"},
                "sales": {0: 10, 1: 12},
                "year": {0: 2008, 1: 2008},
            }
        ).set_index(["year", "flavour", "day"])
        result = df.xs("sat", level="day", drop_level=False)
        expected = df[:1]
        tm.assert_frame_equal(result, expected)

        result = df.xs([2008, "sat"], level=["year", "day"], drop_level=False)
        tm.assert_frame_equal(result, expected)

    def test_xs_view(self):
        # in 0.14 this will return a view if possible a copy otherwise, but
        # this is numpy dependent

        dm = DataFrame(np.arange(20.0).reshape(4, 5), index=range(4), columns=range(5))

        dm.xs(2)[:] = 10
        assert (dm.xs(2) == 10).all()
