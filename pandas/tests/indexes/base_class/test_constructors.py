import numpy as np
import pytest

import pandas as pd
from pandas import (
    Index,
    MultiIndex,
    Series,
)
import pandas._testing as tm


class TestIndexConstructor:
    # Tests for the Index constructor, specifically for cases that do
    #  not return a subclass

    @pytest.mark.parametrize("value", [1, np.int64(1)])
    def test_constructor_corner(self, value):
        # corner case
        msg = (
            r"Index\(\.\.\.\) must be called with a collection of some "
            f"kind, {value} was passed"
        )
        with pytest.raises(TypeError, match=msg):
            Index(value)

    @pytest.mark.parametrize("index_vals", [[("A", 1), "B"], ["B", ("A", 1)]])
    def test_construction_list_mixed_tuples(self, index_vals):
        # see gh-10697: if we are constructing from a mixed list of tuples,
        # make sure that we are independent of the sorting order.
        index = Index(index_vals)
        assert isinstance(index, Index)
        assert not isinstance(index, MultiIndex)

    def test_constructor_cast(self):
        msg = "could not convert string to float"
        with pytest.raises(ValueError, match=msg):
            Index(["a", "b", "c"], dtype=float)

    @pytest.mark.parametrize("tuple_list", [[()], [(), ()]])
    def test_construct_empty_tuples(self, tuple_list):
        # GH #45608
        result = Index(tuple_list)
        expected = MultiIndex.from_tuples(tuple_list)

        tm.assert_index_equal(result, expected)

    def test_index_string_inference(self):
        # GH#54430
        expected = Index(["a", "b"], dtype=pd.StringDtype(na_value=np.nan))
        with pd.option_context("future.infer_string", True):
            ser = Index(["a", "b"])
        tm.assert_index_equal(ser, expected)

        expected = Index(["a", 1], dtype="object")
        with pd.option_context("future.infer_string", True):
            ser = Index(["a", 1])
        tm.assert_index_equal(ser, expected)

    @pytest.mark.parametrize("klass", [Series, Index])
    def test_inference_on_pandas_objects(self, klass):
        # GH#56012
        obj = klass([pd.Timestamp("2019-12-31")], dtype=object)
        result = Index(obj)
        assert result.dtype == np.object_

    def test_constructor_not_read_only(self):
        # GH#57130
        ser = Series([1, 2], dtype=object)
        idx = Index(ser)
        assert idx._values.flags.writeable
