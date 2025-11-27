import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
)
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
        result = obj.set_axis(new_index, axis=0)
        tm.assert_equal(expected, result)

    def test_set_axis_copy(self, obj):
        # Test copy keyword GH#47932
        new_index = list("abcd")[: len(obj)]

        orig = obj.iloc[:]
        expected = obj.copy()
        expected.index = new_index

        result = obj.set_axis(new_index, axis=0)
        tm.assert_equal(expected, result)
        assert result is not obj
        # check we did NOT make a copy
        if obj.ndim == 1:
            assert tm.shares_memory(result, obj)
        else:
            assert all(
                tm.shares_memory(result.iloc[:, i], obj.iloc[:, i])
                for i in range(obj.shape[1])
            )

        result = obj.set_axis(new_index, axis=0)
        tm.assert_equal(expected, result)
        assert result is not obj
        # check we DID NOT make a copy
        if obj.ndim == 1:
            assert tm.shares_memory(result, obj)
        else:
            assert any(
                tm.shares_memory(result.iloc[:, i], obj.iloc[:, i])
                for i in range(obj.shape[1])
            )

        res = obj.set_axis(new_index)
        tm.assert_equal(expected, res)
        # check we did NOT make a copy
        if res.ndim == 1:
            assert tm.shares_memory(res, orig)
        else:
            assert all(
                tm.shares_memory(res.iloc[:, i], orig.iloc[:, i])
                for i in range(res.shape[1])
            )

    def test_set_axis_unnamed_kwarg_warns(self, obj):
        # omitting the "axis" parameter
        new_index = list("abcd")[: len(obj)]

        expected = obj.copy()
        expected.index = new_index

        result = obj.set_axis(new_index)
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize("axis", [3, "foo"])
    def test_set_axis_invalid_axis_name(self, axis, obj):
        # wrong values for the "axis" parameter
        with pytest.raises(ValueError, match="No axis named"):
            obj.set_axis(list("abc"), axis=axis)

    def test_set_axis_setattr_index_not_collection(self, obj):
        # wrong type
        msg = (
            r"Index\(\.\.\.\) must be called with a collection of some "
            r"kind, None was passed"
        )
        with pytest.raises(TypeError, match=msg):
            obj.index = None

    def test_set_axis_setattr_index_wrong_length(self, obj):
        # wrong length
        msg = (
            f"Length mismatch: Expected axis has {len(obj)} elements, "
            f"new values have {len(obj) - 1} elements"
        )
        with pytest.raises(ValueError, match=msg):
            obj.index = np.arange(len(obj) - 1)

        if obj.ndim == 2:
            with pytest.raises(ValueError, match="Length mismatch"):
                obj.columns = obj.columns[::2]


class TestDataFrameSetAxis(SharedSetAxisTests):
    @pytest.fixture
    def obj(self):
        df = DataFrame(
            {"A": [1.1, 2.2, 3.3], "B": [5.0, 6.1, 7.2], "C": [4.4, 5.5, 6.6]},
            index=[2010, 2011, 2012],
        )
        return df

    def test_set_axis_with_allows_duplicate_labels_false(self):
        # GH#44958
        df = DataFrame([[1, 2], [3, 4]], columns=["a", "b"]).set_flags(
            allows_duplicate_labels=False
        )

        result = df.set_axis(labels=["x", "y"], axis=0)
        expected = DataFrame([[1, 2], [3, 4]], index=["x", "y"], columns=["a", "b"])
        tm.assert_frame_equal(result, expected, check_flags=False)


class TestSeriesSetAxis(SharedSetAxisTests):
    @pytest.fixture
    def obj(self):
        ser = Series(np.arange(4), index=[1, 3, 5, 7], dtype="int64")
        return ser
