from __future__ import annotations

import copy
from datetime import datetime
from typing import Any

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray.coding.cftimeindex import CFTimeIndex
from xarray.core.indexes import (
    Hashable,
    Index,
    Indexes,
    PandasIndex,
    PandasMultiIndex,
    _asarray_tuplesafe,
    safe_cast_to_index,
)
from xarray.core.variable import IndexVariable, Variable
from xarray.tests import assert_array_equal, assert_identical, requires_cftime
from xarray.tests.test_coding_times import _all_cftime_date_types


def test_asarray_tuplesafe() -> None:
    res = _asarray_tuplesafe(("a", 1))
    assert isinstance(res, np.ndarray)
    assert res.ndim == 0
    assert res.item() == ("a", 1)

    res = _asarray_tuplesafe([(0,), (1,)])
    assert res.shape == (2,)
    assert res[0] == (0,)
    assert res[1] == (1,)


class CustomIndex(Index):
    def __init__(self, dims) -> None:
        self.dims = dims


class TestIndex:
    @pytest.fixture
    def index(self) -> CustomIndex:
        return CustomIndex({"x": 2})

    def test_from_variables(self) -> None:
        with pytest.raises(NotImplementedError):
            Index.from_variables({}, options={})

    def test_concat(self) -> None:
        with pytest.raises(NotImplementedError):
            Index.concat([], "x")

    def test_stack(self) -> None:
        with pytest.raises(NotImplementedError):
            Index.stack({}, "x")

    def test_unstack(self, index) -> None:
        with pytest.raises(NotImplementedError):
            index.unstack()

    def test_create_variables(self, index) -> None:
        assert index.create_variables() == {}
        assert index.create_variables({"x": "var"}) == {"x": "var"}

    def test_to_pandas_index(self, index) -> None:
        with pytest.raises(TypeError):
            index.to_pandas_index()

    def test_isel(self, index) -> None:
        assert index.isel({}) is None

    def test_sel(self, index) -> None:
        with pytest.raises(NotImplementedError):
            index.sel({})

    def test_join(self, index) -> None:
        with pytest.raises(NotImplementedError):
            index.join(CustomIndex({"y": 2}))

    def test_reindex_like(self, index) -> None:
        with pytest.raises(NotImplementedError):
            index.reindex_like(CustomIndex({"y": 2}))

    def test_equals(self, index) -> None:
        with pytest.raises(NotImplementedError):
            index.equals(CustomIndex({"y": 2}))

    def test_roll(self, index) -> None:
        assert index.roll({}) is None

    def test_rename(self, index) -> None:
        assert index.rename({}, {}) is index

    @pytest.mark.parametrize("deep", [True, False])
    def test_copy(self, index, deep) -> None:
        copied = index.copy(deep=deep)
        assert isinstance(copied, CustomIndex)
        assert copied is not index

        copied.dims["x"] = 3
        if deep:
            assert copied.dims != index.dims
            assert copied.dims != copy.deepcopy(index).dims
        else:
            assert copied.dims is index.dims
            assert copied.dims is copy.copy(index).dims

    def test_getitem(self, index) -> None:
        with pytest.raises(NotImplementedError):
            index[:]


class TestPandasIndex:
    def test_constructor(self) -> None:
        pd_idx = pd.Index([1, 2, 3])
        index = PandasIndex(pd_idx, "x")

        assert index.index.equals(pd_idx)
        # makes a shallow copy
        assert index.index is not pd_idx
        assert index.dim == "x"

        # test no name set for pd.Index
        pd_idx.name = None
        index = PandasIndex(pd_idx, "x")
        assert index.index.name == "x"

    def test_from_variables(self) -> None:
        # pandas has only Float64Index but variable dtype should be preserved
        data = np.array([1.1, 2.2, 3.3], dtype=np.float32)
        var = xr.Variable(
            "x", data, attrs={"unit": "m"}, encoding={"dtype": np.float64}
        )

        index = PandasIndex.from_variables({"x": var}, options={})
        assert index.dim == "x"
        assert index.index.equals(pd.Index(data))
        assert index.coord_dtype == data.dtype

        var2 = xr.Variable(("x", "y"), [[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError, match=r".*only accepts one variable.*"):
            PandasIndex.from_variables({"x": var, "foo": var2}, options={})

        with pytest.raises(
            ValueError, match=r".*cannot set a PandasIndex.*scalar variable.*"
        ):
            PandasIndex.from_variables({"foo": xr.Variable((), 1)}, options={})

        with pytest.raises(
            ValueError, match=r".*only accepts a 1-dimensional variable.*"
        ):
            PandasIndex.from_variables({"foo": var2}, options={})

    def test_from_variables_index_adapter(self) -> None:
        # test index type is preserved when variable wraps a pd.Index
        data = pd.Series(["foo", "bar"], dtype="category")
        pd_idx = pd.Index(data)
        var = xr.Variable("x", pd_idx)

        index = PandasIndex.from_variables({"x": var}, options={})
        assert isinstance(index.index, pd.CategoricalIndex)

    def test_concat_periods(self):
        periods = pd.period_range("2000-01-01", periods=10)
        indexes = [PandasIndex(periods[:5], "t"), PandasIndex(periods[5:], "t")]
        expected = PandasIndex(periods, "t")
        actual = PandasIndex.concat(indexes, dim="t")
        assert actual.equals(expected)
        assert isinstance(actual.index, pd.PeriodIndex)

        positions = [list(range(5)), list(range(5, 10))]
        actual = PandasIndex.concat(indexes, dim="t", positions=positions)
        assert actual.equals(expected)
        assert isinstance(actual.index, pd.PeriodIndex)

    @pytest.mark.parametrize("dtype", [str, bytes])
    def test_concat_str_dtype(self, dtype) -> None:
        a = PandasIndex(np.array(["a"], dtype=dtype), "x", coord_dtype=dtype)
        b = PandasIndex(np.array(["b"], dtype=dtype), "x", coord_dtype=dtype)
        expected = PandasIndex(
            np.array(["a", "b"], dtype=dtype), "x", coord_dtype=dtype
        )

        actual = PandasIndex.concat([a, b], "x")
        assert actual.equals(expected)
        assert np.issubdtype(actual.coord_dtype, dtype)

    def test_concat_empty(self) -> None:
        idx = PandasIndex.concat([], "x")
        assert idx.coord_dtype is np.dtype("O")

    def test_concat_dim_error(self) -> None:
        indexes = [PandasIndex([0, 1], "x"), PandasIndex([2, 3], "y")]

        with pytest.raises(ValueError, match=r"Cannot concatenate.*dimensions.*"):
            PandasIndex.concat(indexes, "x")

    def test_create_variables(self) -> None:
        # pandas has only Float64Index but variable dtype should be preserved
        data = np.array([1.1, 2.2, 3.3], dtype=np.float32)
        pd_idx = pd.Index(data, name="foo")
        index = PandasIndex(pd_idx, "x", coord_dtype=data.dtype)
        index_vars = {
            "foo": IndexVariable(
                "x", data, attrs={"unit": "m"}, encoding={"fill_value": 0.0}
            )
        }

        actual = index.create_variables(index_vars)
        assert_identical(actual["foo"], index_vars["foo"])
        assert actual["foo"].dtype == index_vars["foo"].dtype
        assert actual["foo"].dtype == index.coord_dtype

    def test_to_pandas_index(self) -> None:
        pd_idx = pd.Index([1, 2, 3], name="foo")
        index = PandasIndex(pd_idx, "x")
        assert index.to_pandas_index() is index.index

    def test_sel(self) -> None:
        # TODO: add tests that aren't just for edge cases
        index = PandasIndex(pd.Index([1, 2, 3]), "x")
        with pytest.raises(KeyError, match=r"not all values found"):
            index.sel({"x": [0]})
        with pytest.raises(KeyError):
            index.sel({"x": 0})
        with pytest.raises(ValueError, match=r"does not have a MultiIndex"):
            index.sel({"x": {"one": 0}})

    def test_sel_boolean(self) -> None:
        # index should be ignored and indexer dtype should not be coerced
        # see https://github.com/pydata/xarray/issues/5727
        index = PandasIndex(pd.Index([0.0, 2.0, 1.0, 3.0]), "x")
        actual = index.sel({"x": [False, True, False, True]})
        expected_dim_indexers = {"x": [False, True, False, True]}
        np.testing.assert_array_equal(
            actual.dim_indexers["x"], expected_dim_indexers["x"]
        )

    def test_sel_datetime(self) -> None:
        index = PandasIndex(
            pd.to_datetime(["2000-01-01", "2001-01-01", "2002-01-01"]), "x"
        )
        actual = index.sel({"x": "2001-01-01"})
        expected_dim_indexers = {"x": 1}
        assert actual.dim_indexers == expected_dim_indexers

        actual = index.sel({"x": index.to_pandas_index().to_numpy()[1]})
        assert actual.dim_indexers == expected_dim_indexers

    def test_sel_unsorted_datetime_index_raises(self) -> None:
        index = PandasIndex(pd.to_datetime(["2001", "2000", "2002"]), "x")
        with pytest.raises(KeyError):
            # pandas will try to convert this into an array indexer. We should
            # raise instead, so we can be sure the result of indexing with a
            # slice is always a view.
            index.sel({"x": slice("2001", "2002")})

    def test_equals(self) -> None:
        index1 = PandasIndex([1, 2, 3], "x")
        index2 = PandasIndex([1, 2, 3], "x")
        assert index1.equals(index2) is True

    def test_join(self) -> None:
        index1 = PandasIndex(["a", "aa", "aaa"], "x", coord_dtype="<U3")
        index2 = PandasIndex(["aa", "aaa", "aaaa"], "x", coord_dtype="<U4")

        expected = PandasIndex(["aa", "aaa"], "x")
        actual = index1.join(index2)
        print(actual.index)
        assert actual.equals(expected)
        assert actual.coord_dtype == "=U4"

        expected = PandasIndex(["a", "aa", "aaa", "aaaa"], "x")
        actual = index1.join(index2, how="outer")
        print(actual.index)
        assert actual.equals(expected)
        assert actual.coord_dtype == "=U4"

    def test_reindex_like(self) -> None:
        index1 = PandasIndex([0, 1, 2], "x")
        index2 = PandasIndex([1, 2, 3, 4], "x")

        expected = {"x": [1, 2, -1, -1]}
        actual = index1.reindex_like(index2)
        assert actual.keys() == expected.keys()
        np.testing.assert_array_equal(actual["x"], expected["x"])

        index3 = PandasIndex([1, 1, 2], "x")
        with pytest.raises(ValueError, match=r".*index has duplicate values"):
            index3.reindex_like(index2)

    def test_rename(self) -> None:
        index = PandasIndex(pd.Index([1, 2, 3], name="a"), "x", coord_dtype=np.int32)

        # shortcut
        new_index = index.rename({}, {})
        assert new_index is index

        new_index = index.rename({"a": "b"}, {})
        assert new_index.index.name == "b"
        assert new_index.dim == "x"
        assert new_index.coord_dtype == np.int32

        new_index = index.rename({}, {"x": "y"})
        assert new_index.index.name == "a"
        assert new_index.dim == "y"
        assert new_index.coord_dtype == np.int32

    def test_copy(self) -> None:
        expected = PandasIndex([1, 2, 3], "x", coord_dtype=np.int32)
        actual = expected.copy()

        assert actual.index.equals(expected.index)
        assert actual.index is not expected.index
        assert actual.dim == expected.dim
        assert actual.coord_dtype == expected.coord_dtype

    def test_getitem(self) -> None:
        pd_idx = pd.Index([1, 2, 3])
        expected = PandasIndex(pd_idx, "x", coord_dtype=np.int32)
        actual = expected[1:]

        assert actual.index.equals(pd_idx[1:])
        assert actual.dim == expected.dim
        assert actual.coord_dtype == expected.coord_dtype


class TestPandasMultiIndex:
    def test_constructor(self) -> None:
        foo_data = np.array([0, 0, 1], dtype="int64")
        bar_data = np.array([1.1, 1.2, 1.3], dtype="float64")
        pd_idx = pd.MultiIndex.from_arrays([foo_data, bar_data], names=("foo", "bar"))

        index = PandasMultiIndex(pd_idx, "x")

        assert index.dim == "x"
        assert index.index.equals(pd_idx)
        assert index.index.names == ("foo", "bar")
        assert index.index.name == "x"
        assert index.level_coords_dtype == {
            "foo": foo_data.dtype,
            "bar": bar_data.dtype,
        }

        with pytest.raises(ValueError, match=r".*conflicting multi-index level name.*"):
            PandasMultiIndex(pd_idx, "foo")

        # default level names
        pd_idx = pd.MultiIndex.from_arrays([foo_data, bar_data])
        index = PandasMultiIndex(pd_idx, "x")
        assert list(index.index.names) == ["x_level_0", "x_level_1"]

    def test_from_variables(self) -> None:
        v_level1 = xr.Variable(
            "x", [1, 2, 3], attrs={"unit": "m"}, encoding={"dtype": np.int32}
        )
        v_level2 = xr.Variable(
            "x", ["a", "b", "c"], attrs={"unit": "m"}, encoding={"dtype": "U"}
        )

        index = PandasMultiIndex.from_variables(
            {"level1": v_level1, "level2": v_level2}, options={}
        )

        expected_idx = pd.MultiIndex.from_arrays([v_level1.data, v_level2.data])
        assert index.dim == "x"
        assert index.index.equals(expected_idx)
        assert index.index.name == "x"
        assert list(index.index.names) == ["level1", "level2"]

        var = xr.Variable(("x", "y"), [[1, 2, 3], [4, 5, 6]])
        with pytest.raises(
            ValueError, match=r".*only accepts 1-dimensional variables.*"
        ):
            PandasMultiIndex.from_variables({"var": var}, options={})

        v_level3 = xr.Variable("y", [4, 5, 6])
        with pytest.raises(
            ValueError, match=r"unmatched dimensions for multi-index variables.*"
        ):
            PandasMultiIndex.from_variables(
                {"level1": v_level1, "level3": v_level3}, options={}
            )

    def test_concat(self) -> None:
        pd_midx = pd.MultiIndex.from_product(
            [[0, 1, 2], ["a", "b"]], names=("foo", "bar")
        )
        level_coords_dtype = {"foo": np.int32, "bar": "=U1"}

        midx1 = PandasMultiIndex(
            pd_midx[:2], "x", level_coords_dtype=level_coords_dtype
        )
        midx2 = PandasMultiIndex(
            pd_midx[2:], "x", level_coords_dtype=level_coords_dtype
        )
        expected = PandasMultiIndex(pd_midx, "x", level_coords_dtype=level_coords_dtype)

        actual = PandasMultiIndex.concat([midx1, midx2], "x")
        assert actual.equals(expected)
        assert actual.level_coords_dtype == expected.level_coords_dtype

    def test_stack(self) -> None:
        prod_vars = {
            "x": xr.Variable("x", pd.Index(["b", "a"]), attrs={"foo": "bar"}),
            "y": xr.Variable("y", pd.Index([1, 3, 2])),
        }

        index_xr = PandasMultiIndex.stack(prod_vars, "z")

        assert index_xr.dim == "z"
        index_pd = index_xr.index
        assert isinstance(index_pd, pd.MultiIndex)
        # TODO: change to tuple when pandas 3 is minimum
        assert list(index_pd.names) == ["x", "y"]
        np.testing.assert_array_equal(
            index_pd.codes, [[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]]
        )

        with pytest.raises(
            ValueError, match=r"conflicting dimensions for multi-index product.*"
        ):
            PandasMultiIndex.stack(
                {"x": xr.Variable("x", ["a", "b"]), "x2": xr.Variable("x", [1, 2])},
                "z",
            )

    def test_stack_non_unique(self) -> None:
        prod_vars = {
            "x": xr.Variable("x", pd.Index(["b", "a"]), attrs={"foo": "bar"}),
            "y": xr.Variable("y", pd.Index([1, 1, 2])),
        }

        index_xr = PandasMultiIndex.stack(prod_vars, "z")
        index_pd = index_xr.index
        assert isinstance(index_pd, pd.MultiIndex)

        np.testing.assert_array_equal(
            index_pd.codes, [[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 0, 1]]
        )
        np.testing.assert_array_equal(index_pd.levels[0], ["b", "a"])
        np.testing.assert_array_equal(index_pd.levels[1], [1, 2])

    def test_unstack(self) -> None:
        pd_midx = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2, 3]], names=["one", "two"]
        )
        index = PandasMultiIndex(pd_midx, "x")

        new_indexes, new_pd_idx = index.unstack()
        assert list(new_indexes) == ["one", "two"]
        assert new_indexes["one"].equals(PandasIndex(["a", "b"], "one"))
        assert new_indexes["two"].equals(PandasIndex([1, 2, 3], "two"))
        assert new_pd_idx.equals(pd_midx)

    def test_unstack_requires_unique(self) -> None:
        pd_midx = pd.MultiIndex.from_product([["a", "a"], [1, 2]], names=["one", "two"])
        index = PandasMultiIndex(pd_midx, "x")

        with pytest.raises(
            ValueError, match="Cannot unstack MultiIndex containing duplicates"
        ):
            index.unstack()

    def test_create_variables(self) -> None:
        foo_data = np.array([0, 0, 1], dtype="int64")
        bar_data = np.array([1.1, 1.2, 1.3], dtype="float64")
        pd_idx = pd.MultiIndex.from_arrays([foo_data, bar_data], names=("foo", "bar"))
        index_vars = {
            "x": IndexVariable("x", pd_idx),
            "foo": IndexVariable("x", foo_data, attrs={"unit": "m"}),
            "bar": IndexVariable("x", bar_data, encoding={"fill_value": 0}),
        }

        index = PandasMultiIndex(pd_idx, "x")
        actual = index.create_variables(index_vars)

        for k, expected in index_vars.items():
            assert_identical(actual[k], expected)
            assert actual[k].dtype == expected.dtype
            if k != "x":
                assert actual[k].dtype == index.level_coords_dtype[k]

    def test_sel(self) -> None:
        index = PandasMultiIndex(
            pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=("one", "two")), "x"
        )

        # test tuples inside slice are considered as scalar indexer values
        actual = index.sel({"x": slice(("a", 1), ("b", 2))})
        expected_dim_indexers = {"x": slice(0, 4)}
        assert actual.dim_indexers == expected_dim_indexers

        with pytest.raises(KeyError, match=r"not all values found"):
            index.sel({"x": [0]})
        with pytest.raises(KeyError):
            index.sel({"x": 0})
        with pytest.raises(ValueError, match=r"cannot provide labels for both.*"):
            index.sel({"one": 0, "x": "a"})
        with pytest.raises(
            ValueError,
            match=r"multi-index level names \('three',\) not found in indexes",
        ):
            index.sel({"x": {"three": 0}})
        with pytest.raises(IndexError):
            index.sel({"x": (slice(None), 1, "no_level")})

    def test_join(self):
        midx = pd.MultiIndex.from_product([["a", "aa"], [1, 2]], names=("one", "two"))
        level_coords_dtype = {"one": "=U2", "two": "i"}
        index1 = PandasMultiIndex(midx, "x", level_coords_dtype=level_coords_dtype)
        index2 = PandasMultiIndex(midx[0:2], "x", level_coords_dtype=level_coords_dtype)

        actual = index1.join(index2)
        assert actual.equals(index2)
        assert actual.level_coords_dtype == level_coords_dtype

        actual = index1.join(index2, how="outer")
        assert actual.equals(index1)
        assert actual.level_coords_dtype == level_coords_dtype

    def test_rename(self) -> None:
        level_coords_dtype = {"one": "<U1", "two": np.int32}
        index = PandasMultiIndex(
            pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=("one", "two")),
            "x",
            level_coords_dtype=level_coords_dtype,
        )

        # shortcut
        new_index = index.rename({}, {})
        assert new_index is index

        new_index = index.rename({"two": "three"}, {})
        assert list(new_index.index.names) == ["one", "three"]
        assert new_index.dim == "x"
        assert new_index.level_coords_dtype == {"one": "<U1", "three": np.int32}

        new_index = index.rename({}, {"x": "y"})
        assert list(new_index.index.names) == ["one", "two"]
        assert new_index.dim == "y"
        assert new_index.level_coords_dtype == level_coords_dtype

    def test_copy(self) -> None:
        level_coords_dtype = {"one": "U<1", "two": np.int32}
        expected = PandasMultiIndex(
            pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=("one", "two")),
            "x",
            level_coords_dtype=level_coords_dtype,
        )
        actual = expected.copy()

        assert actual.index.equals(expected.index)
        assert actual.index is not expected.index
        assert actual.dim == expected.dim
        assert actual.level_coords_dtype == expected.level_coords_dtype


class TestIndexes:
    @pytest.fixture
    def indexes_and_vars(self) -> tuple[list[PandasIndex], dict[Hashable, Variable]]:
        x_idx = PandasIndex(pd.Index([1, 2, 3], name="x"), "x")
        y_idx = PandasIndex(pd.Index([4, 5, 6], name="y"), "y")
        z_pd_midx = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2]], names=["one", "two"]
        )
        z_midx = PandasMultiIndex(z_pd_midx, "z")

        indexes = [x_idx, y_idx, z_midx]

        variables = {}
        for idx in indexes:
            variables.update(idx.create_variables())

        return indexes, variables

    @pytest.fixture(params=["pd_index", "xr_index"])
    def unique_indexes(
        self, request, indexes_and_vars
    ) -> list[PandasIndex] | list[pd.Index]:
        xr_indexes, _ = indexes_and_vars

        if request.param == "pd_index":
            return [idx.index for idx in xr_indexes]
        else:
            return xr_indexes

    @pytest.fixture
    def indexes(
        self, unique_indexes, indexes_and_vars
    ) -> Indexes[Index] | Indexes[pd.Index]:
        x_idx, y_idx, z_midx = unique_indexes
        indexes: dict[Any, Index] = {
            "x": x_idx,
            "y": y_idx,
            "z": z_midx,
            "one": z_midx,
            "two": z_midx,
        }

        _, variables = indexes_and_vars

        index_type = Index if isinstance(x_idx, Index) else pd.Index

        return Indexes(indexes, variables, index_type=index_type)

    def test_interface(self, unique_indexes, indexes) -> None:
        x_idx = unique_indexes[0]
        assert list(indexes) == ["x", "y", "z", "one", "two"]
        assert len(indexes) == 5
        assert "x" in indexes
        assert indexes["x"] is x_idx

    def test_variables(self, indexes) -> None:
        assert tuple(indexes.variables) == ("x", "y", "z", "one", "two")

    def test_dims(self, indexes) -> None:
        assert indexes.dims == {"x": 3, "y": 3, "z": 4}

    def test_get_unique(self, unique_indexes, indexes) -> None:
        assert indexes.get_unique() == unique_indexes

    def test_is_multi(self, indexes) -> None:
        assert indexes.is_multi("one") is True
        assert indexes.is_multi("x") is False

    def test_get_all_coords(self, indexes) -> None:
        expected = {
            "z": indexes.variables["z"],
            "one": indexes.variables["one"],
            "two": indexes.variables["two"],
        }
        assert indexes.get_all_coords("one") == expected

        with pytest.raises(ValueError, match=r"errors must be.*"):
            indexes.get_all_coords("x", errors="invalid")

        with pytest.raises(ValueError, match=r"no index found.*"):
            indexes.get_all_coords("no_coord")

        assert indexes.get_all_coords("no_coord", errors="ignore") == {}

    def test_get_all_dims(self, indexes) -> None:
        expected = {"z": 4}
        assert indexes.get_all_dims("one") == expected

    def test_group_by_index(self, unique_indexes, indexes):
        expected = [
            (unique_indexes[0], {"x": indexes.variables["x"]}),
            (unique_indexes[1], {"y": indexes.variables["y"]}),
            (
                unique_indexes[2],
                {
                    "z": indexes.variables["z"],
                    "one": indexes.variables["one"],
                    "two": indexes.variables["two"],
                },
            ),
        ]

        assert indexes.group_by_index() == expected

    def test_to_pandas_indexes(self, indexes) -> None:
        pd_indexes = indexes.to_pandas_indexes()
        assert isinstance(pd_indexes, Indexes)
        assert all(isinstance(idx, pd.Index) for idx in pd_indexes.values())
        assert indexes.variables == pd_indexes.variables

    def test_copy_indexes(self, indexes) -> None:
        copied, index_vars = indexes.copy_indexes()

        assert copied.keys() == indexes.keys()
        for new, original in zip(copied.values(), indexes.values(), strict=True):
            assert new.equals(original)
        # check unique index objects preserved
        assert copied["z"] is copied["one"] is copied["two"]

        assert index_vars.keys() == indexes.variables.keys()
        for new, original in zip(
            index_vars.values(), indexes.variables.values(), strict=True
        ):
            assert_identical(new, original)


def test_safe_cast_to_index():
    dates = pd.date_range("2000-01-01", periods=10)
    x = np.arange(5)
    td = x * np.timedelta64(1, "D")
    for expected, array in [
        (dates, dates.values),
        (pd.Index(x, dtype=object), x.astype(object)),
        (pd.Index(td), td),
        (pd.Index(td, dtype=object), td.astype(object)),
    ]:
        actual = safe_cast_to_index(array)
        assert_array_equal(expected, actual)
        assert expected.dtype == actual.dtype


@requires_cftime
def test_safe_cast_to_index_cftimeindex():
    date_types = _all_cftime_date_types()
    for date_type in date_types.values():
        dates = [date_type(1, 1, day) for day in range(1, 20)]
        expected = CFTimeIndex(dates)
        actual = safe_cast_to_index(np.array(dates))
        assert_array_equal(expected, actual)
        assert expected.dtype == actual.dtype
        assert isinstance(actual, type(expected))


# Test that datetime.datetime objects are never used in a CFTimeIndex
@requires_cftime
def test_safe_cast_to_index_datetime_datetime():
    dates = [datetime(1, 1, day) for day in range(1, 20)]

    expected = pd.Index(dates)
    actual = safe_cast_to_index(np.array(dates))
    assert_array_equal(expected, actual)
    assert isinstance(actual, pd.Index)


@pytest.mark.parametrize("dtype", ["int32", "float32"])
def test_restore_dtype_on_multiindexes(dtype: str) -> None:
    foo = xr.Dataset(coords={"bar": ("bar", np.array([0, 1], dtype=dtype))})
    foo = foo.stack(baz=("bar",))
    assert str(foo["bar"].values.dtype) == dtype


class IndexWithExtraVariables(Index):
    @classmethod
    def from_variables(cls, variables, *, options=None):
        return cls()

    def create_variables(self, variables=None):
        if variables is None:
            # For Coordinates.from_xindex(), return all variables the index can create
            return {
                "time": Variable(dims=("time",), data=[1, 2, 3]),
                "valid_time": Variable(
                    dims=("time",),
                    data=[2, 3, 4],  # time + 1
                    attrs={"description": "time + 1"},
                ),
            }

        result = dict(variables)
        if "time" in variables:
            result["valid_time"] = Variable(
                dims=("time",),
                data=variables["time"].data + 1,
                attrs={"description": "time + 1"},
            )
        return result


def test_set_xindex_with_extra_variables() -> None:
    """Test that set_xindex raises an error when custom index creates extra variables."""

    ds = xr.Dataset(coords={"time": [1, 2, 3]}).reset_index("time")

    # Test that set_xindex raises error for extra variables
    with pytest.raises(ValueError, match="extra variables 'valid_time'"):
        ds.set_xindex("time", IndexWithExtraVariables)


def test_set_xindex_factory_method_pattern() -> None:
    ds = xr.Dataset(coords={"time": [1, 2, 3]}).reset_index("time")

    # Test the recommended factory method pattern
    coord_vars = {"time": ds._variables["time"]}
    index = IndexWithExtraVariables.from_variables(coord_vars)
    coords = xr.Coordinates.from_xindex(index)
    result = ds.assign_coords(coords)

    assert "time" in result.variables
    assert "valid_time" in result.variables
    assert_array_equal(result.valid_time.data, result.time.data + 1)
