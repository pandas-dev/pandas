from datetime import datetime

import numpy as np
import pytest

from pandas import Index, Series
import pandas._testing as tm


class TestSeriesAlterAxes:
    def test_setindex(self, string_series):
        # wrong type
        msg = (
            r"Index\(\.\.\.\) must be called with a collection of some "
            r"kind, None was passed"
        )
        with pytest.raises(TypeError, match=msg):
            string_series.index = None

        # wrong length
        msg = (
            "Length mismatch: Expected axis has 30 elements, "
            "new values have 29 elements"
        )
        with pytest.raises(ValueError, match=msg):
            string_series.index = np.arange(len(string_series) - 1)

        # works
        string_series.index = np.arange(len(string_series))
        assert isinstance(string_series.index, Index)

    # Renaming

    def test_set_name_attribute(self):
        s = Series([1, 2, 3])
        s2 = Series([1, 2, 3], name="bar")
        for name in [7, 7.0, "name", datetime(2001, 1, 1), (1,), "\u05D0"]:
            s.name = name
            assert s.name == name
            s2.name = name
            assert s2.name == name

    def test_set_name(self):
        s = Series([1, 2, 3])
        s2 = s._set_name("foo")
        assert s2.name == "foo"
        assert s.name is None
        assert s is not s2

    def test_set_index_makes_timeseries(self):
        idx = tm.makeDateIndex(10)

        s = Series(range(10))
        s.index = idx
        assert s.index.is_all_dates

    def test_set_axis_inplace_axes(self, axis_series):
        # GH14636
        ser = Series(np.arange(4), index=[1, 3, 5, 7], dtype="int64")

        expected = ser.copy()
        expected.index = list("abcd")

        # inplace=True
        # The FutureWarning comes from the fact that we would like to have
        # inplace default to False some day
        result = ser.copy()
        result.set_axis(list("abcd"), axis=axis_series, inplace=True)
        tm.assert_series_equal(result, expected)

    def test_set_axis_inplace(self):
        # GH14636

        s = Series(np.arange(4), index=[1, 3, 5, 7], dtype="int64")

        expected = s.copy()
        expected.index = list("abcd")

        # inplace=False
        result = s.set_axis(list("abcd"), axis=0, inplace=False)
        tm.assert_series_equal(expected, result)

        # omitting the "axis" parameter
        with tm.assert_produces_warning(None):
            result = s.set_axis(list("abcd"), inplace=False)
        tm.assert_series_equal(result, expected)

        # wrong values for the "axis" parameter
        for axis in [2, "foo"]:
            with pytest.raises(ValueError, match="No axis named"):
                s.set_axis(list("abcd"), axis=axis, inplace=False)
