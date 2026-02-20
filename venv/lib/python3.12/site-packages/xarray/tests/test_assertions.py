from __future__ import annotations

import warnings

import numpy as np
import pytest

import xarray as xr
from xarray.core.coordinates import Coordinates
from xarray.core.indexes import Index
from xarray.tests import has_dask

try:
    from dask.array import from_array as dask_from_array
except ImportError:
    dask_from_array = lambda x: x  # type: ignore[assignment, misc]

try:
    import pint

    unit_registry: pint.UnitRegistry = pint.UnitRegistry(force_ndarray_like=True)

    def quantity(x):
        return unit_registry.Quantity(x, "m")

    has_pint = True
except ImportError:

    def quantity(x):
        return x

    has_pint = False


def test_allclose_regression() -> None:
    x = xr.DataArray(1.01)
    y = xr.DataArray(1.02)
    xr.testing.assert_allclose(x, y, atol=0.01)


@pytest.mark.parametrize(
    "obj1,obj2",
    (
        pytest.param(
            xr.Variable("x", [1e-17, 2]), xr.Variable("x", [0, 3]), id="Variable"
        ),
        pytest.param(
            xr.DataArray([1e-17, 2], dims="x"),
            xr.DataArray([0, 3], dims="x"),
            id="DataArray",
        ),
        pytest.param(
            xr.Dataset({"a": ("x", [1e-17, 2]), "b": ("y", [-2e-18, 2])}),
            xr.Dataset({"a": ("x", [0, 2]), "b": ("y", [0, 1])}),
            id="Dataset",
        ),
        pytest.param(
            xr.DataArray(np.array("a", dtype="|S1")),
            xr.DataArray(np.array("b", dtype="|S1")),
            id="DataArray_with_character_dtype",
        ),
        pytest.param(
            xr.Coordinates({"x": [1e-17, 2]}),
            xr.Coordinates({"x": [0, 3]}),
            id="Coordinates",
        ),
        pytest.param(
            xr.DataTree.from_dict(
                {
                    "/b": xr.Dataset({"a": ("x", [1e-17, 2]), "b": ("y", [-2e-18, 2])}),
                }
            ),
            xr.DataTree.from_dict(
                {
                    "/b": xr.Dataset({"a": ("x", [0, 2]), "b": ("y", [0, 1])}),
                }
            ),
            id="DataTree",
        ),
    ),
)
def test_assert_allclose(obj1, obj2) -> None:
    with pytest.raises(AssertionError):
        xr.testing.assert_allclose(obj1, obj2)
    with pytest.raises(AssertionError):
        xr.testing.assert_allclose(obj1, obj2, check_dim_order=False)


@pytest.mark.parametrize("func", ["assert_equal", "assert_allclose"])
def test_assert_allclose_equal_transpose(func) -> None:
    """Transposed DataArray raises assertion unless check_dim_order=False."""
    obj1 = xr.DataArray([[0, 1, 2], [2, 3, 4]], dims=["a", "b"])
    obj2 = xr.DataArray([[0, 2], [1, 3], [2, 4]], dims=["b", "a"])
    with pytest.raises(AssertionError):
        getattr(xr.testing, func)(obj1, obj2)
    getattr(xr.testing, func)(obj1, obj2, check_dim_order=False)
    ds1 = obj1.to_dataset(name="varname")
    ds1["var2"] = obj1
    ds2 = obj1.to_dataset(name="varname")
    ds2["var2"] = obj1.transpose()
    with pytest.raises(AssertionError):
        getattr(xr.testing, func)(ds1, ds2)
    getattr(xr.testing, func)(ds1, ds2, check_dim_order=False)


def test_assert_equal_transpose_datatree() -> None:
    """Ensure `check_dim_order=False` works for transposed DataTree"""
    ds = xr.Dataset(data_vars={"data": (("x", "y"), [[1, 2]])})

    a = xr.DataTree.from_dict({"node": ds})
    b = xr.DataTree.from_dict({"node": ds.transpose("y", "x")})

    with pytest.raises(AssertionError):
        xr.testing.assert_equal(a, b)

    xr.testing.assert_equal(a, b, check_dim_order=False)


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize(
    "duckarray",
    (
        pytest.param(np.array, id="numpy"),
        pytest.param(
            dask_from_array,
            id="dask",
            marks=pytest.mark.skipif(not has_dask, reason="requires dask"),
        ),
        pytest.param(
            quantity,
            id="pint",
            marks=pytest.mark.skipif(not has_pint, reason="requires pint"),
        ),
    ),
)
@pytest.mark.parametrize(
    ["obj1", "obj2"],
    (
        pytest.param([1e-10, 2], [0.0, 2.0], id="both arrays"),
        pytest.param([1e-17, 2], 0.0, id="second scalar"),
        pytest.param(0.0, [1e-17, 2], id="first scalar"),
    ),
)
def test_assert_duckarray_equal_failing(duckarray, obj1, obj2) -> None:
    # TODO: actually check the repr
    a = duckarray(obj1)
    b = duckarray(obj2)
    with pytest.raises(AssertionError):
        xr.testing.assert_duckarray_equal(a, b)


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize(
    "duckarray",
    (
        pytest.param(
            np.array,
            id="numpy",
        ),
        pytest.param(
            dask_from_array,
            id="dask",
            marks=pytest.mark.skipif(not has_dask, reason="requires dask"),
        ),
        pytest.param(
            quantity,
            id="pint",
            marks=pytest.mark.skipif(not has_pint, reason="requires pint"),
        ),
    ),
)
@pytest.mark.parametrize(
    ["obj1", "obj2"],
    (
        pytest.param([0, 2], [0.0, 2.0], id="both arrays"),
        pytest.param([0, 0], 0.0, id="second scalar"),
        pytest.param(0.0, [0, 0], id="first scalar"),
    ),
)
def test_assert_duckarray_equal(duckarray, obj1, obj2) -> None:
    a = duckarray(obj1)
    b = duckarray(obj2)

    xr.testing.assert_duckarray_equal(a, b)


@pytest.mark.parametrize(
    "func",
    [
        "assert_equal",
        "assert_identical",
        "assert_allclose",
        "assert_duckarray_equal",
        "assert_duckarray_allclose",
    ],
)
def test_ensure_warnings_not_elevated(func) -> None:
    # make sure warnings are not elevated to errors in the assertion functions
    # e.g. by @pytest.mark.filterwarnings("error")
    # see https://github.com/pydata/xarray/pull/4760#issuecomment-774101639

    # define a custom Variable class that raises a warning in assert_*
    class WarningVariable(xr.Variable):
        @property  # type: ignore[misc]
        def dims(self):
            warnings.warn("warning in test", stacklevel=2)
            return super().dims

        def __array__(
            self,
            dtype: np.typing.DTypeLike | None = None,
            /,
            *,
            copy: bool | None = None,
        ) -> np.ndarray:
            warnings.warn("warning in test", stacklevel=2)
            return super().__array__(dtype, copy=copy)

    a = WarningVariable("x", [1])
    b = WarningVariable("x", [2])

    with warnings.catch_warnings(record=True) as w:
        # elevate warnings to errors
        warnings.filterwarnings("error")
        with pytest.raises(AssertionError):
            getattr(xr.testing, func)(a, b)

        assert len(w) > 0

        # ensure warnings still raise outside of assert_*
        with pytest.raises(UserWarning):
            warnings.warn("test", stacklevel=2)

    # ensure warnings stay ignored in assert_*
    with warnings.catch_warnings(record=True) as w:
        # ignore warnings
        warnings.filterwarnings("ignore")
        with pytest.raises(AssertionError):
            getattr(xr.testing, func)(a, b)

        assert len(w) == 0


class CustomIndex(Index):
    """Custom index without equals() implementation for testing."""

    pass


class CustomIndexWithEquals(Index):
    """Custom index with equals() implementation for testing."""

    def __init__(self, name: str = "default"):
        self.name = name

    def equals(self, other: Index, **kwargs) -> bool:
        if not isinstance(other, CustomIndexWithEquals):
            return False
        return self.name == other.name


class TestAssertIdenticalXindexes:
    """Tests for xindex comparison in assert_identical."""

    @pytest.fixture
    def dataset_with_extra_coord(self) -> xr.Dataset:
        """Dataset with a coordinate that can be indexed.

        Returns a Dataset with 'time' indexed and 'time_metadata' not indexed::

            <xarray.Dataset>
            Dimensions:        (time: 4)
            Coordinates:
              * time           (time) float64 0.1 0.2 0.3 0.4
                time_metadata  (time) int64 10 15 20 25
            Data variables:
                data           (time) int64 0 1 2 3

            xindexes: ['time']
        """
        return xr.Dataset(
            {"data": ("time", [0, 1, 2, 3])},
            coords={
                "time": [0.1, 0.2, 0.3, 0.4],
                "time_metadata": ("time", [10, 15, 20, 25]),
            },
        )

    @pytest.fixture
    def dataset_2d(self) -> xr.Dataset:
        """2D dataset for MultiIndex tests.

        Returns a Dataset with both 'x' and 'y' indexed::

            <xarray.Dataset>
            Dimensions:  (x: 2, y: 2)
            Coordinates:
              * x        (x) int64 10 20
              * y        (y) <U1 'a' 'b'
            Data variables:
                data     (x, y) int64 1 2 3 4

            xindexes: ['x', 'y']
        """
        return xr.Dataset(
            {"data": (("x", "y"), [[1, 2], [3, 4]])},
            coords={"x": [10, 20], "y": ["a", "b"]},
        )

    def test_dataset_xindex_difference(
        self, dataset_with_extra_coord: xr.Dataset
    ) -> None:
        """Test that Dataset.identical() and assert_identical detect different xindexes."""
        ds = dataset_with_extra_coord
        ds_extra_index = ds.set_xindex("time_metadata")

        # equals should pass (indexes not compared)
        assert ds.equals(ds_extra_index)
        xr.testing.assert_equal(ds, ds_extra_index)

        # identical should fail (indexes ARE compared)
        assert not ds.identical(ds_extra_index)
        with pytest.raises(AssertionError, match="Indexes only on the right"):
            xr.testing.assert_identical(ds, ds_extra_index)

    def test_assert_identical_same_xindexes(
        self, dataset_with_extra_coord: xr.Dataset
    ) -> None:
        """Test that assert_identical passes when xindexes match."""
        ds = dataset_with_extra_coord

        # Same base datasets - should pass
        xr.testing.assert_identical(ds, ds.copy())

        # Both with extra index - should pass
        ds_extra1 = ds.set_xindex("time_metadata")
        ds_extra2 = ds.set_xindex("time_metadata")
        xr.testing.assert_identical(ds_extra1, ds_extra2)

    def test_dataarray_xindex_difference(
        self, dataset_with_extra_coord: xr.Dataset
    ) -> None:
        """Test that DataArray.identical() and assert_identical detect different xindexes."""
        ds = dataset_with_extra_coord
        ds_extra_index = ds.set_xindex("time_metadata")

        da = ds["data"]
        da_extra_index = ds_extra_index["data"]

        # equals should pass (indexes not compared)
        assert da.equals(da_extra_index)
        xr.testing.assert_equal(da, da_extra_index)

        # identical should fail (indexes ARE compared)
        assert not da.identical(da_extra_index)
        with pytest.raises(AssertionError, match="Indexes only on the right"):
            xr.testing.assert_identical(da, da_extra_index)

    def test_identical_custom_index_without_equals(self) -> None:
        """Test identical() with custom index that doesn't implement equals().

        When equals() is not implemented, falls back to variable comparison.
        Two different CustomIndex objects with the same coordinate values
        should be considered identical via the fallback mechanism.
        """
        coords1 = Coordinates(
            coords={"x": ("x", [1, 2, 3])}, indexes={"x": CustomIndex()}
        )
        coords2 = Coordinates(
            coords={"x": ("x", [1, 2, 3])}, indexes={"x": CustomIndex()}
        )

        ds1 = xr.Dataset({"data": ("x", [10, 20, 30])}, coords=coords1)
        ds2 = xr.Dataset({"data": ("x", [10, 20, 30])}, coords=coords2)

        # Different index objects but same variable values
        # Should be identical via fallback to variable comparison
        assert ds1.xindexes["x"] is not ds2.xindexes["x"]  # Different objects
        assert ds1.identical(ds2)
        xr.testing.assert_identical(ds1, ds2)

    def test_identical_custom_index_with_equals(self) -> None:
        """Test identical() with custom index that implements equals().

        This tests that a custom Index.equals() implementation is actually
        called and its result determines identity.
        """
        coords1 = Coordinates(
            coords={"x": ("x", [1, 2, 3])},
            indexes={"x": CustomIndexWithEquals("index_a")},
        )
        coords2 = Coordinates(
            coords={"x": ("x", [1, 2, 3])},
            indexes={"x": CustomIndexWithEquals("index_a")},
        )
        coords3 = Coordinates(
            coords={"x": ("x", [1, 2, 3])},
            indexes={"x": CustomIndexWithEquals("index_b")},
        )

        ds1 = xr.Dataset({"data": ("x", [10, 20, 30])}, coords=coords1)
        ds2 = xr.Dataset({"data": ("x", [10, 20, 30])}, coords=coords2)
        ds3 = xr.Dataset({"data": ("x", [10, 20, 30])}, coords=coords3)

        # Same index name - should be identical
        assert ds1.identical(ds2)
        xr.testing.assert_identical(ds1, ds2)

        # Different index name (same coord values) - should not be identical
        # This specifically tests the custom equals() is being called
        assert not ds1.identical(ds3)
        with pytest.raises(AssertionError, match="Differing indexes"):
            xr.testing.assert_identical(ds1, ds3)

    def test_identical_mixed_index_types(self) -> None:
        """Test identical() when comparing different index types.

        Different index types should NOT be considered identical, even if
        the underlying coordinate values are the same.
        """
        # Create dataset with PandasIndex (default)
        ds_pandas = xr.Dataset(
            {"data": ("x", [10, 20, 30])},
            coords={"x": [1, 2, 3]},
        )

        # Create dataset with CustomIndex
        coords_custom = Coordinates(
            coords={"x": ("x", [1, 2, 3])}, indexes={"x": CustomIndex()}
        )
        ds_custom = xr.Dataset({"data": ("x", [10, 20, 30])}, coords=coords_custom)

        # Different index types - should NOT be identical
        assert not ds_pandas.identical(ds_custom)

        # But equals should still pass (compares data values only)
        assert ds_pandas.equals(ds_custom)

    def test_identical_pandas_multiindex(self, dataset_2d: xr.Dataset) -> None:
        """Test identical() with PandasMultiIndex."""
        # Stack to create MultiIndex
        ds_stacked1 = dataset_2d.stack(z=("x", "y"))
        ds_stacked2 = dataset_2d.stack(z=("x", "y"))

        # Same MultiIndex - should be identical
        assert ds_stacked1.identical(ds_stacked2)
        xr.testing.assert_identical(ds_stacked1, ds_stacked2)

        # Different stacking order creates different MultiIndex
        ds_stacked_different = dataset_2d.stack(z=("y", "x"))
        assert not ds_stacked1.identical(ds_stacked_different)

    def test_identical_no_indexes(self) -> None:
        """Test identical() when both objects have no indexes.

        Dimensions without coordinates have no indexes.
        """
        ds1 = xr.Dataset({"data": (("x", "y"), [[1, 2], [3, 4]])})
        ds2 = xr.Dataset({"data": (("x", "y"), [[1, 2], [3, 4]])})

        # Dimensions without coordinates = no indexes
        assert list(ds1.xindexes.keys()) == []
        assert list(ds2.xindexes.keys()) == []
        assert ds1.identical(ds2)
        xr.testing.assert_identical(ds1, ds2)

    def test_identical_index_on_different_coords(self) -> None:
        """Test identical() when indexes are on different coordinates."""
        # Index on 'x'
        coords1 = Coordinates(
            coords={"x": ("x", [1, 2, 3]), "y": ("x", [4, 5, 6])},
            indexes={"x": CustomIndex()},
        )
        # Index on 'y' instead
        coords2 = Coordinates(
            coords={"x": ("x", [1, 2, 3]), "y": ("x", [4, 5, 6])},
            indexes={"y": CustomIndex()},
        )

        ds1 = xr.Dataset({"data": ("x", [10, 20, 30])}, coords=coords1)
        ds2 = xr.Dataset({"data": ("x", [10, 20, 30])}, coords=coords2)

        # Different indexed coordinates - should not be identical
        assert not ds1.identical(ds2)
        with pytest.raises(AssertionError, match="Indexes only on the left"):
            xr.testing.assert_identical(ds1, ds2)
