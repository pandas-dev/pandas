from __future__ import annotations

import warnings

import numpy as np
import pytest

import xarray as xr
from xarray.tests import has_dask

try:
    from dask.array import from_array as dask_from_array
except ImportError:
    dask_from_array = lambda x: x  # type: ignore[assignment, misc]

try:
    import pint

    unit_registry = pint.UnitRegistry(force_ndarray_like=True)

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
            self, dtype: np.typing.DTypeLike = None, /, *, copy: bool | None = None
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
