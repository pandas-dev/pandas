from copy import deepcopy
from distutils.version import LooseVersion
from operator import methodcaller

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas.core.dtypes.generic import ABCMultiIndex

import pandas as pd
from pandas import DataFrame, MultiIndex, Series, date_range
import pandas._testing as tm

from .test_generic import Generic

try:
    import xarray

    _XARRAY_INSTALLED = True
except ImportError:
    _XARRAY_INSTALLED = False


class TestDataFrame(Generic):
    _typ = DataFrame
    _comparator = lambda self, x, y: tm.assert_frame_equal(x, y)

    def test_rename_mi(self):
        df = DataFrame(
            [11, 21, 31],
            index=MultiIndex.from_tuples([("A", x) for x in ["a", "B", "c"]]),
        )
        df.rename(str.lower)

    @pytest.mark.parametrize("func", ["_set_axis_name", "rename_axis"])
    def test_set_axis_name(self, func):
        df = pd.DataFrame([[1, 2], [3, 4]])

        result = methodcaller(func, "foo")(df)
        assert df.index.name is None
        assert result.index.name == "foo"

        result = methodcaller(func, "cols", axis=1)(df)
        assert df.columns.name is None
        assert result.columns.name == "cols"

    @pytest.mark.parametrize("func", ["_set_axis_name", "rename_axis"])
    def test_set_axis_name_mi(self, func):
        df = DataFrame(
            np.empty((3, 3)),
            index=MultiIndex.from_tuples([("A", x) for x in list("aBc")]),
            columns=MultiIndex.from_tuples([("C", x) for x in list("xyz")]),
        )

        level_names = ["L1", "L2"]

        result = methodcaller(func, level_names)(df)
        assert result.index.names == level_names
        assert result.columns.names == [None, None]

        result = methodcaller(func, level_names, axis=1)(df)
        assert result.columns.names == ["L1", "L2"]
        assert result.index.names == [None, None]

    def test_nonzero_single_element(self):

        # allow single item via bool method
        df = DataFrame([[True]])
        assert df.bool()

        df = DataFrame([[False]])
        assert not df.bool()

        df = DataFrame([[False, False]])
        msg = "The truth value of a DataFrame is ambiguous"
        with pytest.raises(ValueError, match=msg):
            df.bool()
        with pytest.raises(ValueError, match=msg):
            bool(df)

    def test_get_numeric_data_preserve_dtype(self):

        # get the numeric data
        o = DataFrame({"A": [1, "2", 3.0]})
        result = o._get_numeric_data()
        expected = DataFrame(index=[0, 1, 2], dtype=object)
        self._compare(result, expected)

    def test_metadata_propagation_indiv(self):

        # groupby
        df = DataFrame(
            {
                "A": ["foo", "bar", "foo", "bar", "foo", "bar", "foo", "foo"],
                "B": ["one", "one", "two", "three", "two", "two", "one", "three"],
                "C": np.random.randn(8),
                "D": np.random.randn(8),
            }
        )
        result = df.groupby("A").sum()
        self.check_metadata(df, result)

        # resample
        df = DataFrame(
            np.random.randn(1000, 2),
            index=date_range("20130101", periods=1000, freq="s"),
        )
        result = df.resample("1T")
        self.check_metadata(df, result)

        # merging with override
        # GH 6923
        _metadata = DataFrame._metadata
        _finalize = DataFrame.__finalize__

        np.random.seed(10)
        df1 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=["a", "b"])
        df2 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=["c", "d"])
        DataFrame._metadata = ["filename"]
        df1.filename = "fname1.csv"
        df2.filename = "fname2.csv"

        def finalize(self, other, method=None, **kwargs):

            for name in self._metadata:
                if method == "merge":
                    left, right = other.left, other.right
                    value = getattr(left, name, "") + "|" + getattr(right, name, "")
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, ""))

            return self

        DataFrame.__finalize__ = finalize
        result = df1.merge(df2, left_on=["a"], right_on=["c"], how="inner")
        assert result.filename == "fname1.csv|fname2.csv"

        # concat
        # GH 6927
        DataFrame._metadata = ["filename"]
        df1 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=list("ab"))
        df1.filename = "foo"

        def finalize(self, other, method=None, **kwargs):
            for name in self._metadata:
                if method == "concat":
                    value = "+".join(
                        [getattr(o, name) for o in other.objs if getattr(o, name, None)]
                    )
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, None))

            return self

        DataFrame.__finalize__ = finalize

        result = pd.concat([df1, df1])
        assert result.filename == "foo+foo"

        # reset
        DataFrame._metadata = _metadata
        DataFrame.__finalize__ = _finalize  # FIXME: use monkeypatch

    def test_set_attribute(self):
        # Test for consistent setattr behavior when an attribute and a column
        # have the same name (Issue #8994)
        df = DataFrame({"x": [1, 2, 3]})

        df.y = 2
        df["y"] = [2, 4, 6]
        df.y = 5

        assert df.y == 5
        tm.assert_series_equal(df["y"], Series([2, 4, 6], name="y"))

    def test_deepcopy_empty(self):
        # This test covers empty frame copying with non-empty column sets
        # as reported in issue GH15370
        empty_frame = DataFrame(data=[], index=[], columns=["A"])
        empty_frame_copy = deepcopy(empty_frame)

        self._compare(empty_frame_copy, empty_frame)


# formerly in Generic but only test DataFrame
class TestDataFrame2:
    @pytest.mark.parametrize("value", [1, "True", [1, 2, 3], 5.0])
    def test_validate_bool_args(self, value):
        df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

        msg = 'For argument "inplace" expected type bool, received type'
        with pytest.raises(ValueError, match=msg):
            super(DataFrame, df).rename_axis(
                mapper={"a": "x", "b": "y"}, axis=1, inplace=value
            )

        with pytest.raises(ValueError, match=msg):
            super(DataFrame, df).drop("a", axis=1, inplace=value)

        with pytest.raises(ValueError, match=msg):
            super(DataFrame, df)._consolidate(inplace=value)

        with pytest.raises(ValueError, match=msg):
            super(DataFrame, df).fillna(value=0, inplace=value)

        with pytest.raises(ValueError, match=msg):
            super(DataFrame, df).replace(to_replace=1, value=7, inplace=value)

        with pytest.raises(ValueError, match=msg):
            super(DataFrame, df).interpolate(inplace=value)

        with pytest.raises(ValueError, match=msg):
            super(DataFrame, df)._where(cond=df.a > 2, inplace=value)

        with pytest.raises(ValueError, match=msg):
            super(DataFrame, df).mask(cond=df.a > 2, inplace=value)

    def test_unexpected_keyword(self):
        # GH8597
        df = DataFrame(np.random.randn(5, 2), columns=["jim", "joe"])
        ca = pd.Categorical([0, 0, 2, 2, 3, np.nan])
        ts = df["joe"].copy()
        ts[2] = np.nan

        msg = "unexpected keyword"
        with pytest.raises(TypeError, match=msg):
            df.drop("joe", axis=1, in_place=True)

        with pytest.raises(TypeError, match=msg):
            df.reindex([1, 0], inplace=True)

        with pytest.raises(TypeError, match=msg):
            ca.fillna(0, inplace=True)

        with pytest.raises(TypeError, match=msg):
            ts.fillna(0, in_place=True)


class TestToXArray:
    @pytest.mark.skipif(
        not _XARRAY_INSTALLED
        or _XARRAY_INSTALLED
        and LooseVersion(xarray.__version__) < LooseVersion("0.10.0"),
        reason="xarray >= 0.10.0 required",
    )
    def test_to_xarray_index_types(self, indices):
        if isinstance(indices, ABCMultiIndex):
            pytest.skip("MultiIndex is tested separately")
        if len(indices) == 0:
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

        df.index = indices[:3]
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
