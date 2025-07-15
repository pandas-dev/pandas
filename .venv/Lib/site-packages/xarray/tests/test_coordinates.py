from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
import pytest

from xarray.core.coordinates import Coordinates
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.indexes import Index, PandasIndex, PandasMultiIndex
from xarray.core.variable import IndexVariable, Variable
from xarray.structure.alignment import align
from xarray.tests import assert_identical, source_ndarray


class TestCoordinates:
    def test_init_noindex(self) -> None:
        coords = Coordinates(coords={"foo": ("x", [0, 1, 2])})
        expected = Dataset(coords={"foo": ("x", [0, 1, 2])})
        assert_identical(coords.to_dataset(), expected)

    def test_init_default_index(self) -> None:
        coords = Coordinates(coords={"x": [1, 2]})
        expected = Dataset(coords={"x": [1, 2]})
        assert_identical(coords.to_dataset(), expected)
        assert "x" in coords.xindexes

    @pytest.mark.filterwarnings("error:IndexVariable")
    def test_init_no_default_index(self) -> None:
        # dimension coordinate with no default index (explicit)
        coords = Coordinates(coords={"x": [1, 2]}, indexes={})
        assert "x" not in coords.xindexes
        assert not isinstance(coords["x"], IndexVariable)

    def test_init_from_coords(self) -> None:
        expected = Dataset(coords={"foo": ("x", [0, 1, 2])})
        coords = Coordinates(coords=expected.coords)
        assert_identical(coords.to_dataset(), expected)

        # test variables copied
        assert coords.variables["foo"] is not expected.variables["foo"]

        # test indexes are extracted
        expected = Dataset(coords={"x": [0, 1, 2]})
        coords = Coordinates(coords=expected.coords)
        assert_identical(coords.to_dataset(), expected)
        assert expected.xindexes == coords.xindexes

        # coords + indexes not supported
        with pytest.raises(
            ValueError, match="passing both.*Coordinates.*indexes.*not allowed"
        ):
            coords = Coordinates(
                coords=expected.coords, indexes={"x": PandasIndex([0, 1, 2], "x")}
            )

    def test_init_empty(self) -> None:
        coords = Coordinates()
        assert len(coords) == 0

    def test_init_index_error(self) -> None:
        idx = PandasIndex([1, 2, 3], "x")
        with pytest.raises(ValueError, match="no coordinate variables found"):
            Coordinates(indexes={"x": idx})

        with pytest.raises(TypeError, match=".* is not an `xarray.indexes.Index`"):
            Coordinates(
                coords={"x": ("x", [1, 2, 3])},
                indexes={"x": "not_an_xarray_index"},  # type: ignore[dict-item]
            )

    def test_init_dim_sizes_conflict(self) -> None:
        with pytest.raises(ValueError):
            Coordinates(coords={"foo": ("x", [1, 2]), "bar": ("x", [1, 2, 3, 4])})

    def test_from_xindex(self) -> None:
        idx = PandasIndex([1, 2, 3], "x")
        coords = Coordinates.from_xindex(idx)

        assert isinstance(coords.xindexes["x"], PandasIndex)
        assert coords.xindexes["x"].equals(idx)

        expected = PandasIndex(idx, "x").create_variables()
        assert list(coords.variables) == list(expected)
        assert_identical(expected["x"], coords.variables["x"])

    def test_from_xindex_error(self) -> None:
        class CustomIndexNoCoordsGenerated(Index):
            def create_variables(self, variables: Mapping | None = None):
                return {}

        idx = CustomIndexNoCoordsGenerated()

        with pytest.raises(ValueError, match=".*index.*did not create any coordinate"):
            Coordinates.from_xindex(idx)

    def test_from_pandas_multiindex(self) -> None:
        midx = pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=("one", "two"))
        coords = Coordinates.from_pandas_multiindex(midx, "x")

        assert isinstance(coords.xindexes["x"], PandasMultiIndex)
        assert coords.xindexes["x"].index.equals(midx)
        assert coords.xindexes["x"].dim == "x"

        expected = PandasMultiIndex(midx, "x").create_variables()
        assert list(coords.variables) == list(expected)
        for name in ("x", "one", "two"):
            assert_identical(expected[name], coords.variables[name])

    @pytest.mark.filterwarnings("ignore:return type")
    def test_dims(self) -> None:
        coords = Coordinates(coords={"x": [0, 1, 2]})
        assert set(coords.dims) == {"x"}

    def test_sizes(self) -> None:
        coords = Coordinates(coords={"x": [0, 1, 2]})
        assert coords.sizes == {"x": 3}

    def test_dtypes(self) -> None:
        coords = Coordinates(coords={"x": [0, 1, 2]})
        assert coords.dtypes == {"x": int}

    def test_getitem(self) -> None:
        coords = Coordinates(coords={"x": [0, 1, 2]})
        assert_identical(
            coords["x"],
            DataArray([0, 1, 2], coords={"x": [0, 1, 2]}, name="x"),
        )

    def test_delitem(self) -> None:
        coords = Coordinates(coords={"x": [0, 1, 2]})
        del coords["x"]
        assert "x" not in coords

        with pytest.raises(
            KeyError, match="'nonexistent' is not in coordinate variables"
        ):
            del coords["nonexistent"]

    def test_update(self) -> None:
        coords = Coordinates(coords={"x": [0, 1, 2]})

        coords.update({"y": ("y", [4, 5, 6])})
        assert "y" in coords
        assert "y" in coords.xindexes
        expected = DataArray([4, 5, 6], coords={"y": [4, 5, 6]}, name="y")
        assert_identical(coords["y"], expected)

    def test_equals(self):
        coords = Coordinates(coords={"x": [0, 1, 2]})

        assert coords.equals(coords)
        assert not coords.equals("not_a_coords")

    def test_identical(self):
        coords = Coordinates(coords={"x": [0, 1, 2]})

        assert coords.identical(coords)
        assert not coords.identical("not_a_coords")

    def test_assign(self) -> None:
        coords = Coordinates(coords={"x": [0, 1, 2]})
        expected = Coordinates(coords={"x": [0, 1, 2], "y": [3, 4]})

        actual = coords.assign(y=[3, 4])
        assert_identical(actual, expected)

        actual = coords.assign({"y": [3, 4]})
        assert_identical(actual, expected)

    def test_copy(self) -> None:
        no_index_coords = Coordinates({"foo": ("x", [1, 2, 3])})
        copied = no_index_coords.copy()
        assert_identical(no_index_coords, copied)
        v0 = no_index_coords.variables["foo"]
        v1 = copied.variables["foo"]
        assert v0 is not v1
        assert source_ndarray(v0.data) is source_ndarray(v1.data)

        deep_copied = no_index_coords.copy(deep=True)
        assert_identical(no_index_coords.to_dataset(), deep_copied.to_dataset())
        v0 = no_index_coords.variables["foo"]
        v1 = deep_copied.variables["foo"]
        assert v0 is not v1
        assert source_ndarray(v0.data) is not source_ndarray(v1.data)

    def test_align(self) -> None:
        coords = Coordinates(coords={"x": [0, 1, 2]})

        left = coords

        # test Coordinates._reindex_callback
        right = coords.to_dataset().isel(x=[0, 1]).coords
        left2, right2 = align(left, right, join="inner")
        assert_identical(left2, right2)

        # test Coordinates._overwrite_indexes
        right.update({"x": ("x", [4, 5, 6])})
        left2, right2 = align(left, right, join="override")
        assert_identical(left2, left)
        assert_identical(left2, right2)

    def test_dataset_from_coords_with_multidim_var_same_name(self):
        # regression test for GH #8883
        var = Variable(data=np.arange(6).reshape(2, 3), dims=["x", "y"])
        coords = Coordinates(coords={"x": var}, indexes={})
        ds = Dataset(coords=coords)
        assert ds.coords["x"].dims == ("x", "y")
