from collections.abc import Hashable
from typing import Any

import numpy as np
import pytest

import xarray as xr
from xarray.core.coordinate_transform import CoordinateTransform
from xarray.core.indexes import CoordinateTransformIndex
from xarray.tests import assert_equal, assert_identical


class SimpleCoordinateTransform(CoordinateTransform):
    """Simple uniform scale transform in a 2D space (x/y coordinates)."""

    def __init__(self, shape: tuple[int, int], scale: float, dtype: Any = None):
        super().__init__(("x", "y"), {"x": shape[1], "y": shape[0]}, dtype=dtype)

        self.scale = scale

        # array dimensions in reverse order (y = rows, x = cols)
        self.xy_dims = tuple(self.dims)
        self.dims = (self.dims[1], self.dims[0])

    def forward(self, dim_positions: dict[str, Any]) -> dict[Hashable, Any]:
        assert set(dim_positions) == set(self.dims)
        return {
            name: dim_positions[dim] * self.scale
            for name, dim in zip(self.coord_names, self.xy_dims, strict=False)
        }

    def reverse(self, coord_labels: dict[Hashable, Any]) -> dict[str, Any]:
        return {dim: coord_labels[dim] / self.scale for dim in self.xy_dims}

    def equals(
        self, other: CoordinateTransform, exclude: frozenset[Hashable] | None = None
    ) -> bool:
        if not isinstance(other, SimpleCoordinateTransform):
            return False
        return self.scale == other.scale

    def __repr__(self) -> str:
        return f"Scale({self.scale})"


def test_abstract_coordinate_transform() -> None:
    tr = CoordinateTransform(["x"], {"x": 5})

    with pytest.raises(NotImplementedError):
        tr.forward({"x": [1, 2]})

    with pytest.raises(NotImplementedError):
        tr.reverse({"x": [3.0, 4.0]})

    with pytest.raises(NotImplementedError):
        tr.equals(CoordinateTransform(["x"], {"x": 5}))


def test_coordinate_transform_init() -> None:
    tr = SimpleCoordinateTransform((4, 4), 2.0)

    assert tr.coord_names == ("x", "y")
    # array dimensions in reverse order (y = rows, x = cols)
    assert tr.dims == ("y", "x")
    assert tr.dim_size == {"x": 4, "y": 4}
    assert tr.dtype == np.dtype(np.float64)

    tr2 = SimpleCoordinateTransform((4, 4), 2.0, dtype=np.int64)
    assert tr2.dtype == np.dtype(np.int64)


@pytest.mark.parametrize("dims", [None, ("y", "x")])
def test_coordinate_transform_generate_coords(dims) -> None:
    tr = SimpleCoordinateTransform((2, 2), 2.0)

    actual = tr.generate_coords(dims)
    expected = {"x": [[0.0, 2.0], [0.0, 2.0]], "y": [[0.0, 0.0], [2.0, 2.0]]}
    assert set(actual) == set(expected)
    np.testing.assert_array_equal(actual["x"], expected["x"])
    np.testing.assert_array_equal(actual["y"], expected["y"])


def create_coords(scale: float, shape: tuple[int, int]) -> xr.Coordinates:
    """Create x/y Xarray coordinate variables from a simple coordinate transform."""
    tr = SimpleCoordinateTransform(shape, scale)
    index = CoordinateTransformIndex(tr)
    return xr.Coordinates.from_xindex(index)


def test_coordinate_transform_variable() -> None:
    coords = create_coords(scale=2.0, shape=(2, 2))

    assert coords["x"].dtype == np.dtype(np.float64)
    assert coords["y"].dtype == np.dtype(np.float64)
    assert coords["x"].shape == (2, 2)
    assert coords["y"].shape == (2, 2)

    np.testing.assert_array_equal(np.array(coords["x"]), [[0.0, 2.0], [0.0, 2.0]])
    np.testing.assert_array_equal(np.array(coords["y"]), [[0.0, 0.0], [2.0, 2.0]])

    def assert_repr(var: xr.Variable):
        assert (
            repr(var._data)
            == "CoordinateTransformIndexingAdapter(transform=Scale(2.0))"
        )

    assert_repr(coords["x"].variable)
    assert_repr(coords["y"].variable)


def test_coordinate_transform_variable_repr_inline() -> None:
    var = create_coords(scale=2.0, shape=(2, 2))["x"].variable

    actual = var._data._repr_inline_(70)  # type: ignore[union-attr]
    assert actual == "0.0 2.0 0.0 2.0"

    # truncated inline repr
    var2 = create_coords(scale=2.0, shape=(10, 10))["x"].variable

    actual2 = var2._data._repr_inline_(70)  # type: ignore[union-attr]
    assert (
        actual2 == "0.0 2.0 4.0 6.0 8.0 10.0 12.0 ... 6.0 8.0 10.0 12.0 14.0 16.0 18.0"
    )


def test_coordinate_transform_variable_repr() -> None:
    var = create_coords(scale=2.0, shape=(2, 2))["x"].variable

    actual = repr(var)
    expected = """
<xarray.Variable (y: 2, x: 2)> Size: 32B
[4 values with dtype=float64]
    """.strip()
    assert actual == expected


def test_coordinate_transform_variable_basic_outer_indexing() -> None:
    var = create_coords(scale=2.0, shape=(4, 4))["x"].variable

    assert var[0, 0] == 0.0
    assert var[0, 1] == 2.0
    assert var[0, -1] == 6.0
    np.testing.assert_array_equal(var[:, 0:2], [[0.0, 2.0]] * 4)

    with pytest.raises(IndexError, match="out of bounds index"):
        var[5]

    with pytest.raises(IndexError, match="out of bounds index"):
        var[-5]


def test_coordinate_transform_variable_vectorized_indexing() -> None:
    var = create_coords(scale=2.0, shape=(4, 4))["x"].variable

    actual = var[{"x": xr.Variable("z", [0]), "y": xr.Variable("z", [0])}]
    expected = xr.Variable("z", [0.0])
    assert_equal(actual, expected)

    with pytest.raises(IndexError, match="out of bounds index"):
        var[{"x": xr.Variable("z", [5]), "y": xr.Variable("z", [5])}]


def test_coordinate_transform_setitem_error() -> None:
    var = create_coords(scale=2.0, shape=(4, 4))["x"].variable

    # basic indexing
    with pytest.raises(TypeError, match="setting values is not supported"):
        var[0, 0] = 1.0

    # outer indexing
    with pytest.raises(TypeError, match="setting values is not supported"):
        var[[0, 2], 0] = [1.0, 2.0]

    # vectorized indexing
    with pytest.raises(TypeError, match="setting values is not supported"):
        var[{"x": xr.Variable("z", [0]), "y": xr.Variable("z", [0])}] = 1.0


def test_coordinate_transform_transpose() -> None:
    coords = create_coords(scale=2.0, shape=(2, 2))

    actual = coords["x"].transpose().values
    expected = [[0.0, 0.0], [2.0, 2.0]]
    np.testing.assert_array_equal(actual, expected)


def test_coordinate_transform_equals() -> None:
    ds1 = create_coords(scale=2.0, shape=(2, 2)).to_dataset()
    ds2 = create_coords(scale=2.0, shape=(2, 2)).to_dataset()
    ds3 = create_coords(scale=4.0, shape=(2, 2)).to_dataset()

    # cannot use `assert_equal()` test utility function here yet
    # (indexes invariant check are still based on IndexVariable, which
    # doesn't work with coordinate transform index coordinate variables)
    assert ds1.equals(ds2)
    assert not ds1.equals(ds3)


def test_coordinate_transform_sel() -> None:
    ds = create_coords(scale=2.0, shape=(4, 4)).to_dataset()

    data = [
        [0.0, 1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0, 7.0],
        [8.0, 9.0, 10.0, 11.0],
        [12.0, 13.0, 14.0, 15.0],
    ]
    ds["data"] = (("y", "x"), data)

    actual = ds.sel(
        x=xr.Variable("z", [0.5, 5.5]), y=xr.Variable("z", [0.0, 0.5]), method="nearest"
    )
    expected = ds.isel(x=xr.Variable("z", [0, 3]), y=xr.Variable("z", [0, 0]))

    # cannot use `assert_equal()` test utility function here yet
    # (indexes invariant check are still based on IndexVariable, which
    # doesn't work with coordinate transform index coordinate variables)
    assert actual.equals(expected)

    with pytest.raises(ValueError, match=".*only supports selection.*nearest"):
        ds.sel(x=xr.Variable("z", [0.5, 5.5]), y=xr.Variable("z", [0.0, 0.5]))

    with pytest.raises(ValueError, match="missing labels for coordinate.*y"):
        ds.sel(x=[0.5, 5.5], method="nearest")

    with pytest.raises(TypeError, match=".*only supports advanced.*indexing"):
        ds.sel(x=[0.5, 5.5], y=[0.0, 0.5], method="nearest")

    with pytest.raises(ValueError, match=".*only supports advanced.*indexing"):
        ds.sel(
            x=xr.Variable("z", [0.5, 5.5]),
            y=xr.Variable("z", [0.0, 0.5, 1.5]),
            method="nearest",
        )


def test_coordinate_transform_rename() -> None:
    ds = xr.Dataset(coords=create_coords(scale=2.0, shape=(2, 2)))
    roundtripped = ds.rename(x="u", y="v").rename(u="x", v="y")
    assert_identical(ds, roundtripped, check_default_indexes=False)
