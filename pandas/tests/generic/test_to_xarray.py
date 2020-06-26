import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm


class TestDataFrameToXArray:
    @td.skip_if_no("xarray", "0.10.0")
    def test_to_xarray_index_types(self, index):
        if isinstance(index, pd.MultiIndex):
            pytest.skip("MultiIndex is tested separately")
        if len(index) == 0:
            pytest.skip("Test doesn't make sense for empty index")

        from xarray import Dataset

        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.Categorical(list("abc")),
                "g": pd.date_range("20130101", periods=3),
                "h": pd.date_range("20130101", periods=3, tz="US/Eastern"),
            }
        )

        df.index = index[:3]
        df.index.name = "foo"
        df.columns.name = "bar"
        result = df.to_xarray()
        assert result.dims["foo"] == 3
        assert len(result.coords) == 1
        assert len(result.data_vars) == 8
        tm.assert_almost_equal(list(result.coords.keys()), ["foo"])
        assert isinstance(result, Dataset)

        # idempotency
        # datetimes w/tz are preserved
        # column names are lost
        expected = df.copy()
        expected["f"] = expected["f"].astype(object)
        expected.columns.name = None
        tm.assert_frame_equal(
            result.to_dataframe(), expected,
        )

    @td.skip_if_no("xarray", min_version="0.7.0")
    def test_to_xarray(self):
        from xarray import Dataset

        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.Categorical(list("abc")),
                "g": pd.date_range("20130101", periods=3),
                "h": pd.date_range("20130101", periods=3, tz="US/Eastern"),
            }
        )

        df.index.name = "foo"
        result = df[0:0].to_xarray()
        assert result.dims["foo"] == 0
        assert isinstance(result, Dataset)

        # available in 0.7.1
        # MultiIndex
        df.index = pd.MultiIndex.from_product([["a"], range(3)], names=["one", "two"])
        result = df.to_xarray()
        assert result.dims["one"] == 1
        assert result.dims["two"] == 3
        assert len(result.coords) == 2
        assert len(result.data_vars) == 8
        tm.assert_almost_equal(list(result.coords.keys()), ["one", "two"])
        assert isinstance(result, Dataset)

        result = result.to_dataframe()
        expected = df.copy()
        expected["f"] = expected["f"].astype(object)
        expected.columns.name = None
        tm.assert_frame_equal(result, expected, check_index_type=False)


class TestSeriesToXArray:
    @td.skip_if_no("xarray", "0.10.0")
    def test_to_xarray_index_types(self, index):
        if isinstance(index, pd.MultiIndex):
            pytest.skip("MultiIndex is tested separately")

        from xarray import DataArray

        s = Series(range(len(index)), index=index, dtype="int64")
        s.index.name = "foo"
        result = s.to_xarray()
        repr(result)
        assert len(result) == len(index)
        assert len(result.coords) == 1
        tm.assert_almost_equal(list(result.coords.keys()), ["foo"])
        assert isinstance(result, DataArray)

        # idempotency
        tm.assert_series_equal(result.to_series(), s, check_index_type=False)

    @td.skip_if_no("xarray", min_version="0.7.0")
    def test_to_xarray(self):
        from xarray import DataArray

        s = Series([], dtype=object)
        s.index.name = "foo"
        result = s.to_xarray()
        assert len(result) == 0
        assert len(result.coords) == 1
        tm.assert_almost_equal(list(result.coords.keys()), ["foo"])
        assert isinstance(result, DataArray)

        s = Series(range(6), dtype="int64")
        s.index.name = "foo"
        s.index = pd.MultiIndex.from_product(
            [["a", "b"], range(3)], names=["one", "two"]
        )
        result = s.to_xarray()
        assert len(result) == 2
        tm.assert_almost_equal(list(result.coords.keys()), ["one", "two"])
        assert isinstance(result, DataArray)
        tm.assert_series_equal(result.to_series(), s)
