import numpy as np
import pytest

from pandas import DataFrame, Series
import pandas._testing as tm


class SharedSetAxisTests:
    @pytest.fixture
    def obj(self):
        raise NotImplementedError("Implemented by subclasses")

    def test_set_axis(self, obj):
        # GH14636; this tests setting index for both Series and DataFrame
        new_index = list("abcd")[: len(obj)]

        expected = obj.copy()
        expected.index = new_index

        # inplace=False
        result = obj.set_axis(new_index, axis=0, inplace=False)
        tm.assert_equal(expected, result)

    @pytest.mark.parametrize("axis", [0, "index", 1, "columns"])
    def test_set_axis_inplace_axis(self, axis, obj):
        # GH#14636
        if obj.ndim == 1 and axis in [1, "columns"]:
            # Series only has [0, "index"]
            return

        new_index = list("abcd")[: len(obj)]

        expected = obj.copy()
        if axis in [0, "index"]:
            expected.index = new_index
        else:
            expected.columns = new_index

        result = obj.copy()
        result.set_axis(new_index, axis=axis, inplace=True)
        tm.assert_equal(result, expected)

    def test_set_axis_unnamed_kwarg_warns(self, obj):
        # omitting the "axis" parameter
        new_index = list("abcd")[: len(obj)]

        expected = obj.copy()
        expected.index = new_index

        with tm.assert_produces_warning(None):
            result = obj.set_axis(new_index, inplace=False)
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize("axis", [3, "foo"])
    def test_set_axis_invalid_axis_name(self, axis, obj):
        # wrong values for the "axis" parameter
        with pytest.raises(ValueError, match="No axis named"):
            obj.set_axis(list("abc"), axis=axis)


class TestDataFrameSetAxis(SharedSetAxisTests):
    @pytest.fixture
    def obj(self):
        df = DataFrame(
            {"A": [1.1, 2.2, 3.3], "B": [5.0, 6.1, 7.2], "C": [4.4, 5.5, 6.6]},
            index=[2010, 2011, 2012],
        )
        return df


class TestSeriesSetAxis(SharedSetAxisTests):
    @pytest.fixture
    def obj(self):
        ser = Series(np.arange(4), index=[1, 3, 5, 7], dtype="int64")
        return ser
