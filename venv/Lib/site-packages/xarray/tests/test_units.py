from __future__ import annotations

import contextlib
import functools
import operator
from typing import Any

import numpy as np
import pytest

import xarray as xr
from xarray.core import dtypes, duck_array_ops
from xarray.tests import (
    assert_allclose,
    assert_duckarray_allclose,
    assert_equal,
    assert_identical,
    requires_dask,
    requires_matplotlib,
    requires_numbagg,
)
from xarray.tests.test_plot import PlotTestCase
from xarray.tests.test_variable import _PAD_XR_NP_ARGS

with contextlib.suppress(ImportError):
    import matplotlib.pyplot as plt


pint = pytest.importorskip("pint")
DimensionalityError = pint.errors.DimensionalityError


def create_nan_array(values, dtype):
    """Create array with NaN values, handling cast warnings for int dtypes."""
    import warnings

    # When casting float arrays with NaN to integer, NumPy raises a warning
    # This is expected behavior when dtype is int
    with warnings.catch_warnings():
        if np.issubdtype(dtype, np.integer):
            warnings.filterwarnings("ignore", "invalid value encountered in cast")
        return np.array(values).astype(dtype)


# make sure scalars are converted to 0d arrays so quantities can
# always be treated like ndarrays
unit_registry = pint.UnitRegistry(force_ndarray_like=True)
Quantity = unit_registry.Quantity
no_unit_values = ("none", None)


pytestmark = [
    pytest.mark.filterwarnings("error::pint.UnitStrippedWarning"),
]


def is_compatible(unit1, unit2):
    def dimensionality(obj):
        if isinstance(obj, unit_registry.Quantity | unit_registry.Unit):
            unit_like = obj
        else:
            unit_like = unit_registry.dimensionless

        return unit_like.dimensionality

    return dimensionality(unit1) == dimensionality(unit2)


def compatible_mappings(first, second):
    return {
        key: is_compatible(unit1, unit2)
        for key, (unit1, unit2) in zip_mappings(first, second)
    }


def merge_mappings(base, *mappings):
    result = base.copy()
    for m in mappings:
        result.update(m)

    return result


def zip_mappings(*mappings):
    for key in set(mappings[0]).intersection(*mappings[1:]):
        yield key, tuple(m[key] for m in mappings)


def array_extract_units(obj):
    if isinstance(obj, xr.Variable | xr.DataArray | xr.Dataset):
        obj = obj.data

    try:
        return obj.units
    except AttributeError:
        return None


def array_strip_units(array):
    try:
        return array.magnitude
    except AttributeError:
        return array


def array_attach_units(data, unit):
    if isinstance(data, Quantity) and data.units != unit:
        raise ValueError(f"cannot attach unit {unit} to quantity {data}")

    if unit in no_unit_values or (isinstance(unit, int) and unit == 1):
        return data

    quantity = unit_registry.Quantity(data, unit)
    return quantity


def extract_units(obj):
    if isinstance(obj, xr.Dataset):
        vars_units = {
            name: array_extract_units(value) for name, value in obj.data_vars.items()
        }
        coords_units = {
            name: array_extract_units(value) for name, value in obj.coords.items()
        }

        units = {**vars_units, **coords_units}
    elif isinstance(obj, xr.DataArray):
        vars_units = {obj.name: array_extract_units(obj)}
        coords_units = {
            name: array_extract_units(value) for name, value in obj.coords.items()
        }

        units = {**vars_units, **coords_units}
    elif isinstance(obj, xr.Variable):
        vars_units = {None: array_extract_units(obj.data)}

        units = {**vars_units}
    elif isinstance(obj, Quantity):
        vars_units = {None: array_extract_units(obj)}

        units = {**vars_units}
    else:
        units = {}

    return units


def strip_units(obj):
    if isinstance(obj, xr.Dataset):
        data_vars = {
            strip_units(name): strip_units(value)
            for name, value in obj.data_vars.items()
        }
        coords = {
            strip_units(name): strip_units(value) for name, value in obj.coords.items()
        }

        new_obj = xr.Dataset(data_vars=data_vars, coords=coords)
    elif isinstance(obj, xr.DataArray):
        data = array_strip_units(obj.variable._data)
        coords = {
            strip_units(name): (
                (value.dims, array_strip_units(value.variable._data))
                if isinstance(value.data, Quantity)
                else value  # to preserve multiindexes
            )
            for name, value in obj.coords.items()
        }

        new_obj = xr.DataArray(  # type: ignore[assignment]
            name=strip_units(obj.name), data=data, coords=coords, dims=obj.dims
        )
    elif isinstance(obj, xr.Variable):
        data = array_strip_units(obj.data)
        new_obj = obj.copy(data=data)  # type: ignore[assignment]
    elif isinstance(obj, unit_registry.Quantity):
        new_obj = obj.magnitude
    elif isinstance(obj, list | tuple):
        return type(obj)(strip_units(elem) for elem in obj)
    else:
        new_obj = obj

    return new_obj


def attach_units(obj, units):
    if not isinstance(obj, xr.DataArray | xr.Dataset | xr.Variable):
        units = units.get("data", None) or units.get(None, None) or 1
        return array_attach_units(obj, units)

    if isinstance(obj, xr.Dataset):
        data_vars = {
            name: attach_units(value, units) for name, value in obj.data_vars.items()
        }

        coords = {
            name: attach_units(value, units) for name, value in obj.coords.items()
        }

        new_obj = xr.Dataset(data_vars=data_vars, coords=coords, attrs=obj.attrs)
    elif isinstance(obj, xr.DataArray):
        # try the array name, "data" and None, then fall back to dimensionless
        data_units = units.get(obj.name, None) or units.get(None, None) or 1

        data = array_attach_units(obj.data, data_units)

        coords = {
            name: (
                (value.dims, array_attach_units(value.data, units.get(name) or 1))
                if name in units
                else (value.dims, value.data)
            )
            for name, value in obj.coords.items()
        }
        dims = obj.dims
        attrs = obj.attrs

        new_obj = xr.DataArray(  # type: ignore[assignment]
            name=obj.name, data=data, coords=coords, attrs=attrs, dims=dims
        )
    else:
        data_units = units.get("data", None) or units.get(None, None) or 1

        data = array_attach_units(obj.data, data_units)
        new_obj = obj.copy(data=data)  # type: ignore[assignment]

    return new_obj


def convert_units(obj, to):
    # preprocess
    to = {
        key: None if not isinstance(value, unit_registry.Unit) else value
        for key, value in to.items()
    }
    if isinstance(obj, xr.Dataset):
        data_vars = {
            name: convert_units(array.variable, {None: to.get(name)})
            for name, array in obj.data_vars.items()
        }
        coords = {
            name: convert_units(array.variable, {None: to.get(name)})
            for name, array in obj.coords.items()
        }

        new_obj = xr.Dataset(data_vars=data_vars, coords=coords, attrs=obj.attrs)
    elif isinstance(obj, xr.DataArray):
        name = obj.name

        new_units = to.get(name) or to.get("data") or to.get(None) or None
        data = convert_units(obj.variable, {None: new_units})

        coords = {
            name: (array.dims, convert_units(array.variable, {None: to.get(name)}))
            for name, array in obj.coords.items()
            if name != obj.name
        }

        new_obj = xr.DataArray(  # type: ignore[assignment]
            name=name, data=data, coords=coords, attrs=obj.attrs, dims=obj.dims
        )
    elif isinstance(obj, xr.Variable):
        new_data = convert_units(obj.data, to)
        new_obj = obj.copy(data=new_data)  # type: ignore[assignment]
    elif isinstance(obj, unit_registry.Quantity):
        units = to.get(None)
        new_obj = obj.to(units) if units is not None else obj
    else:
        new_obj = obj

    return new_obj


def assert_units_equal(a, b):
    __tracebackhide__ = True
    assert extract_units(a) == extract_units(b)


@pytest.fixture(params=[np.dtype(float), np.dtype(int)], ids=str)
def dtype(request):
    return request.param


def merge_args(default_args, new_args):
    from itertools import zip_longest

    fill_value = object()
    return [
        second if second is not fill_value else first
        for first, second in zip_longest(default_args, new_args, fillvalue=fill_value)
    ]


class method:
    """wrapper class to help with passing methods via parametrize

    This is works a bit similar to using `partial(Class.method, arg, kwarg)`
    """

    def __init__(self, name, *args, fallback_func=None, **kwargs):
        self.name = name
        self.fallback = fallback_func
        self.args = args
        self.kwargs = kwargs

    def __call__(self, obj, *args, **kwargs):
        from functools import partial

        all_args = merge_args(self.args, args)
        all_kwargs = {**self.kwargs, **kwargs}

        from xarray.core.groupby import GroupBy

        xarray_classes = (
            xr.Variable,
            xr.DataArray,
            xr.Dataset,
            GroupBy,
        )

        if not isinstance(obj, xarray_classes):
            # remove typical xarray args like "dim"
            exclude_kwargs = ("dim", "dims")
            # TODO: figure out a way to replace dim / dims with axis
            all_kwargs = {
                key: value
                for key, value in all_kwargs.items()
                if key not in exclude_kwargs
            }
            if self.fallback is not None:
                func = partial(self.fallback, obj)
            else:
                func_attr = getattr(obj, self.name, None)

                if func_attr is None or not callable(func_attr):
                    # fall back to module level numpy functions
                    numpy_func = getattr(np, self.name)
                    func = partial(numpy_func, obj)
                else:
                    func = func_attr
        else:
            func = getattr(obj, self.name)

        return func(*all_args, **all_kwargs)

    def __repr__(self):
        return f"method_{self.name}"


class function:
    """wrapper class for numpy functions

    Same as method, but the name is used for referencing numpy functions
    """

    def __init__(self, name_or_function, *args, function_label=None, **kwargs):
        if callable(name_or_function):
            self.name = (
                function_label
                if function_label is not None
                else name_or_function.__name__
            )
            self.func = name_or_function
        else:
            self.name = name_or_function if function_label is None else function_label
            self.func = getattr(np, name_or_function)
            if self.func is None:
                raise AttributeError(
                    f"module 'numpy' has no attribute named '{self.name}'"
                )

        self.args = args
        self.kwargs = kwargs

    def __call__(self, *args, **kwargs):
        all_args = merge_args(self.args, args)
        all_kwargs = {**self.kwargs, **kwargs}

        return self.func(*all_args, **all_kwargs)

    def __repr__(self):
        return f"function_{self.name}"


@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
def test_apply_ufunc_dataarray(variant, dtype):
    variants = {
        "data": (unit_registry.m, 1, 1),
        "dims": (1, unit_registry.m, 1),
        "coords": (1, 1, unit_registry.m),
    }
    data_unit, dim_unit, coord_unit = variants[variant]
    func = functools.partial(
        xr.apply_ufunc, np.mean, input_core_dims=[["x"]], kwargs={"axis": -1}
    )

    array = np.linspace(0, 10, 20).astype(dtype) * data_unit
    x = np.arange(20) * dim_unit
    u = np.linspace(-1, 1, 20) * coord_unit
    data_array = xr.DataArray(data=array, dims="x", coords={"x": x, "u": ("x", u)})

    expected = attach_units(func(strip_units(data_array)), extract_units(data_array))
    actual = func(data_array)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
def test_apply_ufunc_dataset(variant, dtype):
    variants = {
        "data": (unit_registry.m, 1, 1),
        "dims": (1, unit_registry.m, 1),
        "coords": (1, 1, unit_registry.s),
    }
    data_unit, dim_unit, coord_unit = variants[variant]
    func = functools.partial(
        xr.apply_ufunc, np.mean, input_core_dims=[["x"]], kwargs={"axis": -1}
    )

    array1 = np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype) * data_unit
    array2 = np.linspace(0, 10, 5).astype(dtype) * data_unit

    x = np.arange(5) * dim_unit
    y = np.arange(10) * dim_unit

    u = np.linspace(-1, 1, 10) * coord_unit

    ds = xr.Dataset(
        data_vars={"a": (("x", "y"), array1), "b": ("x", array2)},
        coords={"x": x, "y": y, "u": ("y", u)},
    )

    expected = attach_units(func(strip_units(ds)), extract_units(ds))
    actual = func(ds)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
@pytest.mark.parametrize("value", (10, dtypes.NA))
def test_align_dataarray(value, variant, unit, error, dtype):
    if variant == "coords" and (
        value != dtypes.NA or isinstance(unit, unit_registry.Unit)
    ):
        pytest.xfail(
            reason=(
                "fill_value is used for both data variables and coords. "
                "See https://github.com/pydata/xarray/issues/4165"
            )
        )

    fill_value = dtypes.get_fill_value(dtype) if value == dtypes.NA else value

    original_unit = unit_registry.m

    variants = {
        "data": ((original_unit, unit), (1, 1), (1, 1)),
        "dims": ((1, 1), (original_unit, unit), (1, 1)),
        "coords": ((1, 1), (1, 1), (original_unit, unit)),
    }
    (
        (data_unit1, data_unit2),
        (dim_unit1, dim_unit2),
        (coord_unit1, coord_unit2),
    ) = variants[variant]

    array1 = np.linspace(0, 10, 2 * 5).reshape(2, 5).astype(dtype) * data_unit1
    array2 = np.linspace(0, 8, 2 * 5).reshape(2, 5).astype(dtype) * data_unit2

    x = np.arange(2) * dim_unit1
    y1 = np.arange(5) * dim_unit1
    y2 = np.arange(2, 7) * dim_unit2

    u1 = np.array([3, 5, 7, 8, 9]) * coord_unit1
    u2 = np.array([7, 8, 9, 11, 13]) * coord_unit2

    coords1 = {"x": x, "y": y1}
    coords2 = {"x": x, "y": y2}
    if variant == "coords":
        coords1["y_a"] = ("y", u1)
        coords2["y_a"] = ("y", u2)

    data_array1 = xr.DataArray(data=array1, coords=coords1, dims=("x", "y"))
    data_array2 = xr.DataArray(data=array2, coords=coords2, dims=("x", "y"))

    fill_value = fill_value * data_unit2
    func = function(xr.align, join="outer", fill_value=fill_value)
    if error is not None and (value != dtypes.NA or isinstance(fill_value, Quantity)):
        with pytest.raises(error):
            func(data_array1, data_array2)

        return

    stripped_kwargs = {
        key: strip_units(
            convert_units(value, {None: data_unit1 if data_unit2 != 1 else None})
        )
        for key, value in func.kwargs.items()
    }
    units_a = extract_units(data_array1)
    units_b = extract_units(data_array2)
    expected_a, expected_b = func(
        strip_units(data_array1),
        strip_units(convert_units(data_array2, units_a)),
        **stripped_kwargs,
    )
    expected_a = attach_units(expected_a, units_a)
    if isinstance(array2, Quantity):
        expected_b = convert_units(attach_units(expected_b, units_a), units_b)
    else:
        expected_b = attach_units(expected_b, units_b)

    actual_a, actual_b = func(data_array1, data_array2)

    assert_units_equal(expected_a, actual_a)
    assert_allclose(expected_a, actual_a)
    assert_units_equal(expected_b, actual_b)
    assert_allclose(expected_b, actual_b)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
@pytest.mark.parametrize("value", (10, dtypes.NA))
def test_align_dataset(value, unit, variant, error, dtype):
    if variant == "coords" and (
        value != dtypes.NA or isinstance(unit, unit_registry.Unit)
    ):
        pytest.xfail(
            reason=(
                "fill_value is used for both data variables and coords. "
                "See https://github.com/pydata/xarray/issues/4165"
            )
        )

    fill_value = dtypes.get_fill_value(dtype) if value == dtypes.NA else value

    original_unit = unit_registry.m

    variants = {
        "data": ((original_unit, unit), (1, 1), (1, 1)),
        "dims": ((1, 1), (original_unit, unit), (1, 1)),
        "coords": ((1, 1), (1, 1), (original_unit, unit)),
    }
    (
        (data_unit1, data_unit2),
        (dim_unit1, dim_unit2),
        (coord_unit1, coord_unit2),
    ) = variants[variant]

    array1 = np.linspace(0, 10, 2 * 5).reshape(2, 5).astype(dtype) * data_unit1
    array2 = np.linspace(0, 10, 2 * 5).reshape(2, 5).astype(dtype) * data_unit2

    x = np.arange(2) * dim_unit1
    y1 = np.arange(5) * dim_unit1
    y2 = np.arange(2, 7) * dim_unit2

    u1 = np.array([3, 5, 7, 8, 9]) * coord_unit1
    u2 = np.array([7, 8, 9, 11, 13]) * coord_unit2

    coords1 = {"x": x, "y": y1}
    coords2 = {"x": x, "y": y2}
    if variant == "coords":
        coords1["u"] = ("y", u1)
        coords2["u"] = ("y", u2)

    ds1 = xr.Dataset(data_vars={"a": (("x", "y"), array1)}, coords=coords1)
    ds2 = xr.Dataset(data_vars={"a": (("x", "y"), array2)}, coords=coords2)

    fill_value = fill_value * data_unit2
    func = function(xr.align, join="outer", fill_value=fill_value)
    if error is not None and (value != dtypes.NA or isinstance(fill_value, Quantity)):
        with pytest.raises(error):
            func(ds1, ds2)

        return

    stripped_kwargs = {
        key: strip_units(
            convert_units(value, {None: data_unit1 if data_unit2 != 1 else None})
        )
        for key, value in func.kwargs.items()
    }
    units_a = extract_units(ds1)
    units_b = extract_units(ds2)
    expected_a, expected_b = func(
        strip_units(ds1),
        strip_units(convert_units(ds2, units_a)),
        **stripped_kwargs,
    )
    expected_a = attach_units(expected_a, units_a)
    if isinstance(array2, Quantity):
        expected_b = convert_units(attach_units(expected_b, units_a), units_b)
    else:
        expected_b = attach_units(expected_b, units_b)

    actual_a, actual_b = func(ds1, ds2)

    assert_units_equal(expected_a, actual_a)
    assert_allclose(expected_a, actual_a)
    assert_units_equal(expected_b, actual_b)
    assert_allclose(expected_b, actual_b)


def test_broadcast_dataarray(dtype):
    # uses align internally so more thorough tests are not needed
    array1 = np.linspace(0, 10, 2) * unit_registry.Pa
    array2 = np.linspace(0, 10, 3) * unit_registry.Pa

    a = xr.DataArray(data=array1, dims="x")
    b = xr.DataArray(data=array2, dims="y")

    units_a = extract_units(a)
    units_b = extract_units(b)
    expected_a, expected_b = xr.broadcast(strip_units(a), strip_units(b))
    expected_a = attach_units(expected_a, units_a)
    expected_b = convert_units(attach_units(expected_b, units_a), units_b)

    actual_a, actual_b = xr.broadcast(a, b)

    assert_units_equal(expected_a, actual_a)
    assert_identical(expected_a, actual_a)
    assert_units_equal(expected_b, actual_b)
    assert_identical(expected_b, actual_b)


def test_broadcast_dataset(dtype):
    # uses align internally so more thorough tests are not needed
    array1 = np.linspace(0, 10, 2) * unit_registry.Pa
    array2 = np.linspace(0, 10, 3) * unit_registry.Pa

    x1 = np.arange(2)
    y1 = np.arange(3)

    x2 = np.arange(2, 4)
    y2 = np.arange(3, 6)

    ds = xr.Dataset(
        data_vars={"a": ("x", array1), "b": ("y", array2)}, coords={"x": x1, "y": y1}
    )
    other = xr.Dataset(
        data_vars={
            "a": ("x", array1.to(unit_registry.hPa)),
            "b": ("y", array2.to(unit_registry.hPa)),
        },
        coords={"x": x2, "y": y2},
    )

    units_a = extract_units(ds)
    units_b = extract_units(other)
    expected_a, expected_b = xr.broadcast(strip_units(ds), strip_units(other))
    expected_a = attach_units(expected_a, units_a)
    expected_b = attach_units(expected_b, units_b)

    actual_a, actual_b = xr.broadcast(ds, other)

    assert_units_equal(expected_a, actual_a)
    assert_identical(expected_a, actual_a)
    assert_units_equal(expected_b, actual_b)
    assert_identical(expected_b, actual_b)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
@pytest.mark.filterwarnings(
    "ignore:.*the default value for coords will change:FutureWarning"
)
def test_combine_by_coords(variant, unit, error, dtype):
    original_unit = unit_registry.m

    variants = {
        "data": ((original_unit, unit), (1, 1), (1, 1)),
        "dims": ((1, 1), (original_unit, unit), (1, 1)),
        "coords": ((1, 1), (1, 1), (original_unit, unit)),
    }
    (
        (data_unit1, data_unit2),
        (dim_unit1, dim_unit2),
        (coord_unit1, coord_unit2),
    ) = variants[variant]

    array1 = np.zeros(shape=(2, 3), dtype=dtype) * data_unit1
    array2 = np.zeros(shape=(2, 3), dtype=dtype) * data_unit1
    x = np.arange(1, 4) * 10 * dim_unit1
    y = np.arange(2) * dim_unit1
    u = np.arange(3) * coord_unit1

    other_array1 = np.ones_like(array1) * data_unit2
    other_array2 = np.ones_like(array2) * data_unit2
    other_x = np.arange(1, 4) * 10 * dim_unit2
    other_y = np.arange(2, 4) * dim_unit2
    other_u = np.arange(3, 6) * coord_unit2

    ds = xr.Dataset(
        data_vars={"a": (("y", "x"), array1), "b": (("y", "x"), array2)},
        coords={"x": x, "y": y, "u": ("x", u)},
    )
    other = xr.Dataset(
        data_vars={"a": (("y", "x"), other_array1), "b": (("y", "x"), other_array2)},
        coords={"x": other_x, "y": other_y, "u": ("x", other_u)},
    )

    if error is not None:
        with pytest.raises(error):
            xr.combine_by_coords([ds, other], coords="different", compat="no_conflicts")

        return

    units = extract_units(ds)
    expected = attach_units(
        xr.combine_by_coords(
            [strip_units(ds), strip_units(convert_units(other, units))]
        ),
        units,
    )
    actual = xr.combine_by_coords([ds, other])

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
def test_combine_nested(variant, unit, error, dtype):
    original_unit = unit_registry.m

    variants = {
        "data": ((original_unit, unit), (1, 1), (1, 1)),
        "dims": ((1, 1), (original_unit, unit), (1, 1)),
        "coords": ((1, 1), (1, 1), (original_unit, unit)),
    }
    (
        (data_unit1, data_unit2),
        (dim_unit1, dim_unit2),
        (coord_unit1, coord_unit2),
    ) = variants[variant]

    array1 = np.zeros(shape=(2, 3), dtype=dtype) * data_unit1
    array2 = np.zeros(shape=(2, 3), dtype=dtype) * data_unit1

    x = np.arange(1, 4) * 10 * dim_unit1
    y = np.arange(2) * dim_unit1
    z = np.arange(3) * coord_unit1

    ds1 = xr.Dataset(
        data_vars={"a": (("y", "x"), array1), "b": (("y", "x"), array2)},
        coords={"x": x, "y": y, "z": ("x", z)},
    )
    ds2 = xr.Dataset(
        data_vars={
            "a": (("y", "x"), np.ones_like(array1) * data_unit2),
            "b": (("y", "x"), np.ones_like(array2) * data_unit2),
        },
        coords={
            "x": np.arange(3) * dim_unit2,
            "y": np.arange(2, 4) * dim_unit2,
            "z": ("x", np.arange(-3, 0) * coord_unit2),
        },
    )
    ds3 = xr.Dataset(
        data_vars={
            "a": (("y", "x"), np.full_like(array1, fill_value=np.nan) * data_unit2),
            "b": (("y", "x"), np.full_like(array2, fill_value=np.nan) * data_unit2),
        },
        coords={
            "x": np.arange(3, 6) * dim_unit2,
            "y": np.arange(4, 6) * dim_unit2,
            "z": ("x", np.arange(3, 6) * coord_unit2),
        },
    )
    ds4 = xr.Dataset(
        data_vars={
            "a": (("y", "x"), -1 * np.ones_like(array1) * data_unit2),
            "b": (("y", "x"), -1 * np.ones_like(array2) * data_unit2),
        },
        coords={
            "x": np.arange(6, 9) * dim_unit2,
            "y": np.arange(6, 8) * dim_unit2,
            "z": ("x", np.arange(6, 9) * coord_unit2),
        },
    )

    func = function(xr.combine_nested, concat_dim=["x", "y"], join="outer")
    if error is not None:
        with pytest.raises(error):
            func([[ds1, ds2], [ds3, ds4]])

        return

    units = extract_units(ds1)
    convert_and_strip = lambda ds: strip_units(convert_units(ds, units))
    expected = attach_units(
        func(
            [
                [strip_units(ds1), convert_and_strip(ds2)],
                [convert_and_strip(ds3), convert_and_strip(ds4)],
            ]
        ),
        units,
    )
    actual = func([[ds1, ds2], [ds3, ds4]])

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
def test_concat_dataarray(variant, unit, error, dtype):
    original_unit = unit_registry.m

    variants = {
        "data": ((original_unit, unit), (1, 1), (1, 1)),
        "dims": ((1, 1), (original_unit, unit), (1, 1)),
        "coords": ((1, 1), (1, 1), (original_unit, unit)),
    }
    (
        (data_unit1, data_unit2),
        (dim_unit1, dim_unit2),
        (coord_unit1, coord_unit2),
    ) = variants[variant]

    array1 = np.linspace(0, 5, 10).astype(dtype) * data_unit1
    array2 = np.linspace(-5, 0, 5).astype(dtype) * data_unit2

    x1 = np.arange(5, 15) * dim_unit1
    x2 = np.arange(5) * dim_unit2

    u1 = np.linspace(1, 2, 10).astype(dtype) * coord_unit1
    u2 = np.linspace(0, 1, 5).astype(dtype) * coord_unit2

    arr1 = xr.DataArray(data=array1, coords={"x": x1, "u": ("x", u1)}, dims="x")
    arr2 = xr.DataArray(data=array2, coords={"x": x2, "u": ("x", u2)}, dims="x")

    if error is not None:
        with pytest.raises(error):
            xr.concat([arr1, arr2], dim="x")

        return

    units = extract_units(arr1)
    expected = attach_units(
        xr.concat(
            [strip_units(arr1), strip_units(convert_units(arr2, units))], dim="x"
        ),
        units,
    )
    actual = xr.concat([arr1, arr2], dim="x")

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
def test_concat_dataset(variant, unit, error, dtype):
    original_unit = unit_registry.m

    variants = {
        "data": ((original_unit, unit), (1, 1), (1, 1)),
        "dims": ((1, 1), (original_unit, unit), (1, 1)),
        "coords": ((1, 1), (1, 1), (original_unit, unit)),
    }
    (
        (data_unit1, data_unit2),
        (dim_unit1, dim_unit2),
        (coord_unit1, coord_unit2),
    ) = variants[variant]

    array1 = np.linspace(0, 5, 10).astype(dtype) * data_unit1
    array2 = np.linspace(-5, 0, 5).astype(dtype) * data_unit2

    x1 = np.arange(5, 15) * dim_unit1
    x2 = np.arange(5) * dim_unit2

    u1 = np.linspace(1, 2, 10).astype(dtype) * coord_unit1
    u2 = np.linspace(0, 1, 5).astype(dtype) * coord_unit2

    ds1 = xr.Dataset(data_vars={"a": ("x", array1)}, coords={"x": x1, "u": ("x", u1)})
    ds2 = xr.Dataset(data_vars={"a": ("x", array2)}, coords={"x": x2, "u": ("x", u2)})

    if error is not None:
        with pytest.raises(error):
            xr.concat([ds1, ds2], dim="x")

        return

    units = extract_units(ds1)
    expected = attach_units(
        xr.concat([strip_units(ds1), strip_units(convert_units(ds2, units))], dim="x"),
        units,
    )
    actual = xr.concat([ds1, ds2], dim="x")

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
def test_merge_dataarray(variant, unit, error, dtype):
    original_unit = unit_registry.m

    variants = {
        "data": ((original_unit, unit), (1, 1), (1, 1)),
        "dims": ((1, 1), (original_unit, unit), (1, 1)),
        "coords": ((1, 1), (1, 1), (original_unit, unit)),
    }
    (
        (data_unit1, data_unit2),
        (dim_unit1, dim_unit2),
        (coord_unit1, coord_unit2),
    ) = variants[variant]

    array1 = np.linspace(0, 1, 2 * 3).reshape(2, 3).astype(dtype) * data_unit1
    x1 = np.arange(2) * dim_unit1
    y1 = np.arange(3) * dim_unit1
    u1 = np.linspace(10, 20, 2) * coord_unit1
    v1 = np.linspace(10, 20, 3) * coord_unit1

    array2 = np.linspace(1, 2, 2 * 4).reshape(2, 4).astype(dtype) * data_unit2
    x2 = np.arange(2, 4) * dim_unit2
    z2 = np.arange(4) * dim_unit1
    u2 = np.linspace(20, 30, 2) * coord_unit2
    w2 = np.linspace(10, 20, 4) * coord_unit1

    array3 = np.linspace(0, 2, 3 * 4).reshape(3, 4).astype(dtype) * data_unit2
    y3 = np.arange(3, 6) * dim_unit2
    z3 = np.arange(4, 8) * dim_unit2
    v3 = np.linspace(10, 20, 3) * coord_unit2
    w3 = np.linspace(10, 20, 4) * coord_unit2

    arr1 = xr.DataArray(
        name="a",
        data=array1,
        coords={"x": x1, "y": y1, "u": ("x", u1), "v": ("y", v1)},
        dims=("x", "y"),
    )
    arr2 = xr.DataArray(
        name="a",
        data=array2,
        coords={"x": x2, "z": z2, "u": ("x", u2), "w": ("z", w2)},
        dims=("x", "z"),
    )
    arr3 = xr.DataArray(
        name="a",
        data=array3,
        coords={"y": y3, "z": z3, "v": ("y", v3), "w": ("z", w3)},
        dims=("y", "z"),
    )

    func = function(xr.merge, compat="no_conflicts", join="outer")
    if error is not None:
        with pytest.raises(error):
            func([arr1, arr2, arr3])

        return

    units = {
        "a": data_unit1,
        "u": coord_unit1,
        "v": coord_unit1,
        "w": coord_unit1,
        "x": dim_unit1,
        "y": dim_unit1,
        "z": dim_unit1,
    }
    convert_and_strip = lambda arr: strip_units(convert_units(arr, units))

    expected = attach_units(
        func(
            [convert_and_strip(arr1), convert_and_strip(arr2), convert_and_strip(arr3)]
        ),
        units,
    )

    actual = func([arr1, arr2, arr3])

    assert_units_equal(expected, actual)
    assert_allclose(expected, actual)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
def test_merge_dataset(variant, unit, error, dtype):
    original_unit = unit_registry.m

    variants = {
        "data": ((original_unit, unit), (1, 1), (1, 1)),
        "dims": ((1, 1), (original_unit, unit), (1, 1)),
        "coords": ((1, 1), (1, 1), (original_unit, unit)),
    }
    (
        (data_unit1, data_unit2),
        (dim_unit1, dim_unit2),
        (coord_unit1, coord_unit2),
    ) = variants[variant]

    array1 = np.zeros(shape=(2, 3), dtype=dtype) * data_unit1
    array2 = np.zeros(shape=(2, 3), dtype=dtype) * data_unit1

    x = np.arange(11, 14) * dim_unit1
    y = np.arange(2) * dim_unit1
    u = np.arange(3) * coord_unit1

    ds1 = xr.Dataset(
        data_vars={"a": (("y", "x"), array1), "b": (("y", "x"), array2)},
        coords={"x": x, "y": y, "u": ("x", u)},
    )
    ds2 = xr.Dataset(
        data_vars={
            "a": (("y", "x"), np.ones_like(array1) * data_unit2),
            "b": (("y", "x"), np.ones_like(array2) * data_unit2),
        },
        coords={
            "x": np.arange(3) * dim_unit2,
            "y": np.arange(2, 4) * dim_unit2,
            "u": ("x", np.arange(-3, 0) * coord_unit2),
        },
    )
    ds3 = xr.Dataset(
        data_vars={
            "a": (("y", "x"), np.full_like(array1, np.nan) * data_unit2),
            "b": (("y", "x"), np.full_like(array2, np.nan) * data_unit2),
        },
        coords={
            "x": np.arange(3, 6) * dim_unit2,
            "y": np.arange(4, 6) * dim_unit2,
            "u": ("x", np.arange(3, 6) * coord_unit2),
        },
    )

    func = function(xr.merge, compat="no_conflicts", join="outer")
    if error is not None:
        with pytest.raises(error):
            func([ds1, ds2, ds3])

        return

    units = extract_units(ds1)
    convert_and_strip = lambda ds: strip_units(convert_units(ds, units))
    expected = attach_units(
        func([convert_and_strip(ds1), convert_and_strip(ds2), convert_and_strip(ds3)]),
        units,
    )
    actual = func([ds1, ds2, ds3])

    assert_units_equal(expected, actual)
    assert_allclose(expected, actual)


@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
@pytest.mark.parametrize("func", (xr.zeros_like, xr.ones_like))
def test_replication_dataarray(func, variant, dtype):
    unit = unit_registry.m

    variants = {
        "data": (unit, 1, 1),
        "dims": (1, unit, 1),
        "coords": (1, 1, unit),
    }
    data_unit, dim_unit, coord_unit = variants[variant]
    array = np.linspace(0, 10, 20).astype(dtype) * data_unit
    x = np.arange(20) * dim_unit
    u = np.linspace(0, 1, 20) * coord_unit

    data_array = xr.DataArray(data=array, dims="x", coords={"x": x, "u": ("x", u)})
    units = extract_units(data_array)
    units.pop(data_array.name)

    expected = attach_units(func(strip_units(data_array)), units)
    actual = func(data_array)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        "coords",
    ),
)
@pytest.mark.parametrize("func", (xr.zeros_like, xr.ones_like))
def test_replication_dataset(func, variant, dtype):
    unit = unit_registry.m

    variants = {
        "data": ((unit_registry.m, unit_registry.Pa), 1, 1),
        "dims": ((1, 1), unit, 1),
        "coords": ((1, 1), 1, unit),
    }
    (data_unit1, data_unit2), dim_unit, coord_unit = variants[variant]

    array1 = np.linspace(0, 10, 20).astype(dtype) * data_unit1
    array2 = np.linspace(5, 10, 10).astype(dtype) * data_unit2
    x = np.arange(20).astype(dtype) * dim_unit
    y = np.arange(10).astype(dtype) * dim_unit
    u = np.linspace(0, 1, 10) * coord_unit

    ds = xr.Dataset(
        data_vars={"a": ("x", array1), "b": ("y", array2)},
        coords={"x": x, "y": y, "u": ("y", u)},
    )
    units = {
        name: unit
        for name, unit in extract_units(ds).items()
        if name not in ds.data_vars
    }

    expected = attach_units(func(strip_units(ds)), units)

    actual = func(ds)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        pytest.param(
            "coords",
            marks=pytest.mark.xfail(reason="can't copy quantity into non-quantity"),
        ),
    ),
)
def test_replication_full_like_dataarray(variant, dtype):
    # since full_like will strip units and then use the units of the
    # fill value, we don't need to try multiple units
    unit = unit_registry.m

    variants = {
        "data": (unit, 1, 1),
        "dims": (1, unit, 1),
        "coords": (1, 1, unit),
    }
    data_unit, dim_unit, coord_unit = variants[variant]
    array = np.linspace(0, 5, 10) * data_unit
    x = np.arange(10) * dim_unit
    u = np.linspace(0, 1, 10) * coord_unit
    data_array = xr.DataArray(data=array, dims="x", coords={"x": x, "u": ("x", u)})

    fill_value = -1 * unit_registry.degK

    units = extract_units(data_array)
    units[data_array.name] = fill_value.units
    expected = attach_units(
        xr.full_like(strip_units(data_array), fill_value=strip_units(fill_value)), units
    )
    actual = xr.full_like(data_array, fill_value=fill_value)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "variant",
    (
        "data",
        pytest.param(
            "dims", marks=pytest.mark.skip(reason="indexes don't support units")
        ),
        pytest.param(
            "coords",
            marks=pytest.mark.xfail(reason="can't copy quantity into non-quantity"),
        ),
    ),
)
def test_replication_full_like_dataset(variant, dtype):
    unit = unit_registry.m

    variants = {
        "data": ((unit_registry.s, unit_registry.Pa), 1, 1),
        "dims": ((1, 1), unit, 1),
        "coords": ((1, 1), 1, unit),
    }
    (data_unit1, data_unit2), dim_unit, coord_unit = variants[variant]

    array1 = np.linspace(0, 10, 20).astype(dtype) * data_unit1
    array2 = np.linspace(5, 10, 10).astype(dtype) * data_unit2
    x = np.arange(20).astype(dtype) * dim_unit
    y = np.arange(10).astype(dtype) * dim_unit

    u = np.linspace(0, 1, 10) * coord_unit

    ds = xr.Dataset(
        data_vars={"a": ("x", array1), "b": ("y", array2)},
        coords={"x": x, "y": y, "u": ("y", u)},
    )

    fill_value = -1 * unit_registry.degK

    units = {
        **extract_units(ds),
        **dict.fromkeys(ds.data_vars, unit_registry.degK),
    }
    expected = attach_units(
        xr.full_like(strip_units(ds), fill_value=strip_units(fill_value)), units
    )
    actual = xr.full_like(ds, fill_value=fill_value)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize("fill_value", (np.nan, 10.2))
def test_where_dataarray(fill_value, unit, error, dtype):
    array = np.linspace(0, 5, 10).astype(dtype) * unit_registry.m

    x = xr.DataArray(data=array, dims="x")
    cond = x < 5 * unit_registry.m
    fill_value = fill_value * unit

    if error is not None and not (
        np.isnan(fill_value) and not isinstance(fill_value, Quantity)
    ):
        with pytest.raises(error):
            xr.where(cond, x, fill_value)

        return

    expected = attach_units(
        xr.where(
            cond,
            strip_units(x),
            strip_units(convert_units(fill_value, {None: unit_registry.m})),
        ),
        extract_units(x),
    )
    actual = xr.where(cond, x, fill_value)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "unit,error",
    (
        pytest.param(1, DimensionalityError, id="no_unit"),
        pytest.param(
            unit_registry.dimensionless, DimensionalityError, id="dimensionless"
        ),
        pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
        pytest.param(unit_registry.mm, None, id="compatible_unit"),
        pytest.param(unit_registry.m, None, id="identical_unit"),
    ),
    ids=repr,
)
@pytest.mark.parametrize("fill_value", (np.nan, 10.2))
def test_where_dataset(fill_value, unit, error, dtype):
    array1 = np.linspace(0, 5, 10).astype(dtype) * unit_registry.m
    array2 = np.linspace(-5, 0, 10).astype(dtype) * unit_registry.m

    ds = xr.Dataset(data_vars={"a": ("x", array1), "b": ("x", array2)})
    cond = array1 < 2 * unit_registry.m
    fill_value = fill_value * unit

    if error is not None and not (
        np.isnan(fill_value) and not isinstance(fill_value, Quantity)
    ):
        with pytest.raises(error):
            xr.where(cond, ds, fill_value)

        return

    expected = attach_units(
        xr.where(
            cond,
            strip_units(ds),
            strip_units(convert_units(fill_value, {None: unit_registry.m})),
        ),
        extract_units(ds),
    )
    actual = xr.where(cond, ds, fill_value)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


def test_dot_dataarray(dtype):
    array1 = (
        np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype)
        * unit_registry.m
        / unit_registry.s
    )
    array2 = (
        np.linspace(10, 20, 10 * 20).reshape(10, 20).astype(dtype) * unit_registry.s
    )

    data_array = xr.DataArray(data=array1, dims=("x", "y"))
    other = xr.DataArray(data=array2, dims=("y", "z"))

    with xr.set_options(use_opt_einsum=False):
        expected = attach_units(
            xr.dot(strip_units(data_array), strip_units(other)), {None: unit_registry.m}
        )
        actual = xr.dot(data_array, other)

    assert_units_equal(expected, actual)
    assert_identical(expected, actual)


class TestVariable:
    @pytest.mark.parametrize(
        "func",
        (
            method("all"),
            method("any"),
            method("argmax", dim="x"),
            method("argmin", dim="x"),
            method("argsort"),
            method("cumprod"),
            method("cumsum"),
            method("max"),
            method("mean"),
            method("median"),
            method("min"),
            method("prod"),
            method("std"),
            method("sum"),
            method("var"),
        ),
        ids=repr,
    )
    def test_aggregation(self, func, dtype):
        array = np.linspace(0, 1, 10).astype(dtype) * (
            unit_registry.m if func.name != "cumprod" else unit_registry.dimensionless
        )
        variable = xr.Variable("x", array)

        numpy_kwargs = func.kwargs.copy()
        if "dim" in func.kwargs:
            numpy_kwargs["axis"] = variable.get_axis_num(numpy_kwargs.pop("dim"))

        units = extract_units(func(array, **numpy_kwargs))
        expected = attach_units(func(strip_units(variable)), units)
        actual = func(variable)

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    def test_aggregate_complex(self):
        variable = xr.Variable("x", [1, 2j, np.nan] * unit_registry.m)
        expected = xr.Variable((), (0.5 + 1j) * unit_registry.m)
        actual = variable.mean()

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            method("astype", np.float32),
            method("conj"),
            method("conjugate"),
            method("clip", min=2, max=7),
        ),
        ids=repr,
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_numpy_methods(self, func, unit, error, dtype):
        array = np.linspace(0, 1, 10).astype(dtype) * unit_registry.m
        variable = xr.Variable("x", array)

        args = [
            item * unit if isinstance(item, int | float | list) else item
            for item in func.args
        ]
        kwargs = {
            key: value * unit if isinstance(value, int | float | list) else value
            for key, value in func.kwargs.items()
        }

        if error is not None and func.name in ("searchsorted", "clip"):
            with pytest.raises(error):
                func(variable, *args, **kwargs)

            return

        converted_args = [
            strip_units(convert_units(item, {None: unit_registry.m})) for item in args
        ]
        converted_kwargs = {
            key: strip_units(convert_units(value, {None: unit_registry.m}))
            for key, value in kwargs.items()
        }

        units = extract_units(func(array, *args, **kwargs))
        expected = attach_units(
            func(strip_units(variable), *converted_args, **converted_kwargs), units
        )
        actual = func(variable, *args, **kwargs)

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    @pytest.mark.parametrize(
        "func", (method("item", 5), method("searchsorted", 5)), ids=repr
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_raw_numpy_methods(self, func, unit, error, dtype):
        array = np.linspace(0, 1, 10).astype(dtype) * unit_registry.m
        variable = xr.Variable("x", array)

        args = [
            (
                item * unit
                if isinstance(item, int | float | list) and func.name != "item"
                else item
            )
            for item in func.args
        ]
        kwargs = {
            key: (
                value * unit
                if isinstance(value, int | float | list) and func.name != "item"
                else value
            )
            for key, value in func.kwargs.items()
        }

        if error is not None and func.name != "item":
            with pytest.raises(error):
                func(variable, *args, **kwargs)

            return

        converted_args = [
            (
                strip_units(convert_units(item, {None: unit_registry.m}))
                if func.name != "item"
                else item
            )
            for item in args
        ]
        converted_kwargs = {
            key: (
                strip_units(convert_units(value, {None: unit_registry.m}))
                if func.name != "item"
                else value
            )
            for key, value in kwargs.items()
        }

        units = extract_units(func(array, *args, **kwargs))
        expected = attach_units(
            func(strip_units(variable), *converted_args, **converted_kwargs), units
        )
        actual = func(variable, *args, **kwargs)

        assert_units_equal(expected, actual)
        assert_duckarray_allclose(expected, actual)

    @pytest.mark.parametrize(
        "func", (method("isnull"), method("notnull"), method("count")), ids=repr
    )
    def test_missing_value_detection(self, func):
        array = (
            np.array(
                [
                    [1.4, 2.3, np.nan, 7.2],
                    [np.nan, 9.7, np.nan, np.nan],
                    [2.1, np.nan, np.nan, 4.6],
                    [9.9, np.nan, 7.2, 9.1],
                ]
            )
            * unit_registry.degK
        )
        variable = xr.Variable(("x", "y"), array)

        expected = func(strip_units(variable))
        actual = func(variable)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_missing_value_fillna(self, unit, error):
        value = 10
        array = (
            np.array(
                [
                    [1.4, 2.3, np.nan, 7.2],
                    [np.nan, 9.7, np.nan, np.nan],
                    [2.1, np.nan, np.nan, 4.6],
                    [9.9, np.nan, 7.2, 9.1],
                ]
            )
            * unit_registry.m
        )
        variable = xr.Variable(("x", "y"), array)

        fill_value = value * unit

        if error is not None:
            with pytest.raises(error):
                variable.fillna(value=fill_value)

            return

        expected = attach_units(
            strip_units(variable).fillna(
                value=fill_value.to(unit_registry.m).magnitude
            ),
            extract_units(variable),
        )
        actual = variable.fillna(value=fill_value)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(
                unit_registry.cm,
                id="compatible_unit",
            ),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "convert_data",
        (
            pytest.param(False, id="no_conversion"),
            pytest.param(True, id="with_conversion"),
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("equals"),
            pytest.param(
                method("identical"),
                marks=pytest.mark.skip(reason="behavior of identical is undecided"),
            ),
        ),
        ids=repr,
    )
    def test_comparisons(self, func, unit, convert_data, dtype):
        array = np.linspace(0, 1, 9).astype(dtype)
        quantity1 = array * unit_registry.m
        variable = xr.Variable("x", quantity1)

        if convert_data and is_compatible(unit_registry.m, unit):
            quantity2 = convert_units(array * unit_registry.m, {None: unit})
        else:
            quantity2 = array * unit
        other = xr.Variable("x", quantity2)

        expected = func(
            strip_units(variable),
            strip_units(
                convert_units(other, extract_units(variable))
                if is_compatible(unit_registry.m, unit)
                else other
            ),
        )
        if func.name == "identical":
            expected &= extract_units(variable) == extract_units(other)
        else:
            expected &= all(
                compatible_mappings(
                    extract_units(variable), extract_units(other)
                ).values()
            )

        actual = func(variable, other)

        assert expected == actual

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    def test_broadcast_equals(self, unit, dtype):
        base_unit = unit_registry.m
        left_array = np.ones(shape=(2, 2), dtype=dtype) * base_unit
        value = (
            (1 * base_unit).to(unit).magnitude if is_compatible(unit, base_unit) else 1
        )
        right_array = np.full(shape=(2,), fill_value=value, dtype=dtype) * unit

        left = xr.Variable(("x", "y"), left_array)
        right = xr.Variable("x", right_array)

        units = {
            **extract_units(left),
            **({} if is_compatible(unit, base_unit) else {None: None}),
        }
        expected = strip_units(left).broadcast_equals(
            strip_units(convert_units(right, units))
        ) & is_compatible(unit, base_unit)
        actual = left.broadcast_equals(right)

        assert expected == actual

    @pytest.mark.parametrize("dask", [False, pytest.param(True, marks=[requires_dask])])
    @pytest.mark.parametrize(
        ["variable", "indexers"],
        (
            pytest.param(
                xr.Variable("x", np.linspace(0, 5, 10)),
                {"x": 4},
                id="single value-single indexer",
            ),
            pytest.param(
                xr.Variable("x", np.linspace(0, 5, 10)),
                {"x": [5, 2, 9, 1]},
                id="multiple values-single indexer",
            ),
            pytest.param(
                xr.Variable(("x", "y"), np.linspace(0, 5, 20).reshape(4, 5)),
                {"x": 1, "y": 4},
                id="single value-multiple indexers",
            ),
            pytest.param(
                xr.Variable(("x", "y"), np.linspace(0, 5, 20).reshape(4, 5)),
                {"x": [0, 1, 2], "y": [0, 2, 4]},
                id="multiple values-multiple indexers",
            ),
        ),
    )
    def test_isel(self, variable, indexers, dask, dtype):
        if dask:
            variable = variable.chunk(dict.fromkeys(variable.dims, 2))
        quantified = xr.Variable(
            variable.dims, variable.data.astype(dtype) * unit_registry.s
        )

        expected = attach_units(
            strip_units(quantified).isel(indexers), extract_units(quantified)
        )
        actual = quantified.isel(indexers)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            function(lambda x, *_: +x, function_label="unary_plus"),
            function(lambda x, *_: -x, function_label="unary_minus"),
            function(lambda x, *_: abs(x), function_label="absolute"),
            function(lambda x, y: x + y, function_label="sum"),
            function(lambda x, y: y + x, function_label="commutative_sum"),
            function(lambda x, y: x * y, function_label="product"),
            function(lambda x, y: y * x, function_label="commutative_product"),
        ),
        ids=repr,
    )
    def test_1d_math(self, func, unit, error, dtype):
        base_unit = unit_registry.m
        array = np.arange(5).astype(dtype) * base_unit
        variable = xr.Variable("x", array)

        values = np.ones(5)
        y = values * unit

        if error is not None and func.name in ("sum", "commutative_sum"):
            with pytest.raises(error):
                func(variable, y)

            return

        units = extract_units(func(array, y))
        if all(compatible_mappings(units, extract_units(y)).values()):
            converted_y = convert_units(y, units)
        else:
            converted_y = y

        if all(compatible_mappings(units, extract_units(variable)).values()):
            converted_variable = convert_units(variable, units)
        else:
            converted_variable = variable

        expected = attach_units(
            func(strip_units(converted_variable), strip_units(converted_y)), units
        )
        actual = func(variable, y)

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "func", (method("where"), method("_getitem_with_mask")), ids=repr
    )
    def test_masking(self, func, unit, error, dtype):
        base_unit = unit_registry.m
        array = np.linspace(0, 5, 10).astype(dtype) * base_unit
        variable = xr.Variable("x", array)
        cond = np.array([True, False] * 5)

        other = -1 * unit

        if error is not None:
            with pytest.raises(error):
                func(variable, cond, other)

            return

        expected = attach_units(
            func(
                strip_units(variable),
                cond,
                strip_units(
                    convert_units(
                        other,
                        (
                            {None: base_unit}
                            if is_compatible(base_unit, unit)
                            else {None: None}
                        ),
                    )
                ),
            ),
            extract_units(variable),
        )
        actual = func(variable, cond, other)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("dim", ("x", "y", "z", "t", "all"))
    def test_squeeze(self, dim, dtype):
        shape = (2, 1, 3, 1, 1, 2)
        names = list("abcdef")
        dim_lengths = dict(zip(names, shape, strict=True))
        array = np.ones(shape=shape) * unit_registry.m
        variable = xr.Variable(names, array)

        kwargs = {"dim": dim} if dim != "all" and dim_lengths.get(dim, 0) == 1 else {}
        expected = attach_units(
            strip_units(variable).squeeze(**kwargs), extract_units(variable)
        )
        actual = variable.squeeze(**kwargs)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("compute_backend", ["numbagg", None], indirect=True)
    @pytest.mark.parametrize(
        "func",
        (
            method("coarsen", windows={"y": 2}, func=np.mean),
            method("quantile", q=[0.25, 0.75]),
            pytest.param(
                method("rank", dim="x"),
                marks=pytest.mark.skip(reason="rank not implemented for non-ndarray"),
            ),
            method("roll", {"x": 2}),
            pytest.param(
                method("rolling_window", "x", 3, "window"),
                marks=pytest.mark.xfail(reason="converts to ndarray"),
            ),
            method("reduce", np.std, "x"),
            method("round", 2),
            method("shift", {"x": -2}),
            method("transpose", "y", "x"),
        ),
        ids=repr,
    )
    def test_computation(self, func, dtype, compute_backend):
        base_unit = unit_registry.m
        array = np.linspace(0, 5, 5 * 10).reshape(5, 10).astype(dtype) * base_unit
        variable = xr.Variable(("x", "y"), array)

        expected = attach_units(func(strip_units(variable)), extract_units(variable))

        actual = func(variable)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_searchsorted(self, unit, error, dtype):
        base_unit = unit_registry.m
        array = np.linspace(0, 5, 10).astype(dtype) * base_unit
        variable = xr.Variable("x", array)

        value = 0 * unit

        if error is not None:
            with pytest.raises(error):
                variable.searchsorted(value)  # type: ignore[attr-defined]

            return

        expected = strip_units(variable).searchsorted(
            strip_units(convert_units(value, {None: base_unit}))
        )

        actual = variable.searchsorted(value)  # type: ignore[attr-defined]

        assert_units_equal(expected, actual)
        np.testing.assert_allclose(expected, actual)

    def test_stack(self, dtype):
        array = np.linspace(0, 5, 3 * 10).reshape(3, 10).astype(dtype) * unit_registry.m
        variable = xr.Variable(("x", "y"), array)

        expected = attach_units(
            strip_units(variable).stack(z=("x", "y")), extract_units(variable)
        )
        actual = variable.stack(z=("x", "y"))

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    def test_unstack(self, dtype):
        array = np.linspace(0, 5, 3 * 10).astype(dtype) * unit_registry.m
        variable = xr.Variable("z", array)

        expected = attach_units(
            strip_units(variable).unstack(z={"x": 3, "y": 10}), extract_units(variable)
        )
        actual = variable.unstack(z={"x": 3, "y": 10})

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_concat(self, unit, error, dtype):
        array1 = (
            np.linspace(0, 5, 9 * 10).reshape(3, 6, 5).astype(dtype) * unit_registry.m
        )
        array2 = np.linspace(5, 10, 10 * 3).reshape(3, 2, 5).astype(dtype) * unit

        variable = xr.Variable(("x", "y", "z"), array1)
        other = xr.Variable(("x", "y", "z"), array2)

        if error is not None:
            with pytest.raises(error):
                xr.Variable.concat([variable, other], dim="y")

            return

        units = extract_units(variable)
        expected = attach_units(
            xr.Variable.concat(
                [strip_units(variable), strip_units(convert_units(other, units))],
                dim="y",
            ),
            units,
        )
        actual = xr.Variable.concat([variable, other], dim="y")

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    def test_set_dims(self, dtype):
        array = np.linspace(0, 5, 3 * 10).reshape(3, 10).astype(dtype) * unit_registry.m
        variable = xr.Variable(("x", "y"), array)

        dims = {"z": 6, "x": 3, "a": 1, "b": 4, "y": 10}
        expected = attach_units(
            strip_units(variable).set_dims(dims), extract_units(variable)
        )
        actual = variable.set_dims(dims)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    def test_copy(self, dtype):
        array = np.linspace(0, 5, 10).astype(dtype) * unit_registry.m
        other = np.arange(10).astype(dtype) * unit_registry.s

        variable = xr.Variable("x", array)
        expected = attach_units(
            strip_units(variable).copy(data=strip_units(other)), extract_units(other)
        )
        actual = variable.copy(data=other)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    def test_no_conflicts(self, unit, dtype):
        base_unit = unit_registry.m
        array1 = (
            np.array(
                [
                    [6.3, 0.3, 0.45],
                    [np.nan, 0.3, 0.3],
                    [3.7, np.nan, 0.2],
                    [9.43, 0.3, 0.7],
                ]
            )
            * base_unit
        )
        array2 = np.array([np.nan, 0.3, np.nan]) * unit

        variable = xr.Variable(("x", "y"), array1)
        other = xr.Variable("y", array2)

        expected = strip_units(variable).no_conflicts(
            strip_units(
                convert_units(
                    other, {None: base_unit if is_compatible(base_unit, unit) else None}
                )
            )
        ) & is_compatible(base_unit, unit)
        actual = variable.no_conflicts(other)

        assert expected == actual

    @pytest.mark.parametrize(
        "mode",
        [
            "constant",
            "mean",
            "median",
            "reflect",
            "edge",
            "linear_ramp",
            "maximum",
            "minimum",
            "symmetric",
            "wrap",
        ],
    )
    @pytest.mark.parametrize("xr_arg, np_arg", _PAD_XR_NP_ARGS)
    def test_pad(self, mode, xr_arg, np_arg):
        data = np.arange(4 * 3 * 2).reshape(4, 3, 2) * unit_registry.m
        v = xr.Variable(["x", "y", "z"], data)

        expected = attach_units(
            strip_units(v).pad(mode=mode, **xr_arg),
            extract_units(v),
        )
        actual = v.pad(mode=mode, **xr_arg)

        assert_units_equal(expected, actual)
        assert_equal(actual, expected)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_pad_unit_constant_value(self, unit, error, dtype):
        array = np.linspace(0, 5, 3 * 10).reshape(3, 10).astype(dtype) * unit_registry.m
        variable = xr.Variable(("x", "y"), array)

        fill_value = -100 * unit

        func = method("pad", mode="constant", x=(2, 3), y=(1, 4))
        if error is not None:
            with pytest.raises(error):
                func(variable, constant_values=fill_value)

            return

        units = extract_units(variable)
        expected = attach_units(
            func(
                strip_units(variable),
                constant_values=strip_units(convert_units(fill_value, units)),
            ),
            units,
        )
        actual = func(variable, constant_values=fill_value)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)


class TestDataArray:
    @pytest.mark.parametrize(
        "variant",
        (
            pytest.param(
                "with_dims",
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
            "with_coords",
            "without_coords",
        ),
    )
    def test_init(self, variant, dtype):
        array = np.linspace(1, 2, 10, dtype=dtype) * unit_registry.m

        x = np.arange(len(array)) * unit_registry.s
        y = x.to(unit_registry.ms)

        variants = {
            "with_dims": {"x": x},
            "with_coords": {"y": ("x", y)},
            "without_coords": {},
        }

        kwargs = {"data": array, "dims": "x", "coords": variants[variant]}
        data_array = xr.DataArray(**kwargs)

        assert isinstance(data_array.data, Quantity)
        assert all(
            {
                name: isinstance(coord.data, Quantity)
                for name, coord in data_array.coords.items()
            }.values()
        )

    @pytest.mark.parametrize(
        "func", (pytest.param(str, id="str"), pytest.param(repr, id="repr"))
    )
    @pytest.mark.parametrize(
        "variant",
        (
            pytest.param(
                "with_dims",
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
            pytest.param("with_coords"),
            pytest.param("without_coords"),
        ),
    )
    def test_repr(self, func, variant, dtype):
        array = np.linspace(1, 2, 10, dtype=dtype) * unit_registry.m
        x = np.arange(len(array)) * unit_registry.s
        y = x.to(unit_registry.ms)

        variants = {
            "with_dims": {"x": x},
            "with_coords": {"y": ("x", y)},
            "without_coords": {},
        }

        kwargs = {"data": array, "dims": "x", "coords": variants[variant]}
        data_array = xr.DataArray(**kwargs)

        # FIXME: this just checks that the repr does not raise
        # warnings or errors, but does not check the result
        func(data_array)

    @pytest.mark.parametrize(
        "func",
        (
            function("all"),
            function("any"),
            pytest.param(
                function("argmax"),
                marks=pytest.mark.skip(
                    reason="calling np.argmax as a function on xarray objects is not "
                    "supported"
                ),
            ),
            pytest.param(
                function("argmin"),
                marks=pytest.mark.skip(
                    reason="calling np.argmin as a function on xarray objects is not "
                    "supported"
                ),
            ),
            function("max"),
            function("mean"),
            pytest.param(
                function("median"),
                marks=pytest.mark.skip(
                    reason="median does not work with dataarrays yet"
                ),
            ),
            function("min"),
            function("prod"),
            function("sum"),
            function("std"),
            function("var"),
            function("cumsum"),
            function("cumprod"),
            method("all"),
            method("any"),
            method("argmax", dim="x"),
            method("argmin", dim="x"),
            method("max"),
            method("mean"),
            method("median"),
            method("min"),
            method("prod"),
            method("sum"),
            method("std"),
            method("var"),
            method("cumsum"),
            method("cumprod"),
        ),
        ids=repr,
    )
    def test_aggregation(self, func, dtype):
        array = np.arange(10).astype(dtype) * (
            unit_registry.m if func.name != "cumprod" else unit_registry.dimensionless
        )
        data_array = xr.DataArray(data=array, dims="x")

        numpy_kwargs = func.kwargs.copy()
        if "dim" in numpy_kwargs:
            numpy_kwargs["axis"] = data_array.get_axis_num(numpy_kwargs.pop("dim"))

        # units differ based on the applied function, so we need to
        # first compute the units
        units = extract_units(func(array))
        expected = attach_units(func(strip_units(data_array)), units)
        actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            pytest.param(operator.neg, id="negate"),
            pytest.param(abs, id="absolute"),
            pytest.param(np.round, id="round"),
        ),
    )
    def test_unary_operations(self, func, dtype):
        array = np.arange(10).astype(dtype) * unit_registry.m
        data_array = xr.DataArray(data=array)

        units = extract_units(func(array))
        expected = attach_units(func(strip_units(data_array)), units)
        actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            pytest.param(lambda x: 2 * x, id="multiply"),
            pytest.param(lambda x: x + x, id="add"),
            pytest.param(lambda x: x[0] + x, id="add scalar"),
            pytest.param(lambda x: x.T @ x, id="matrix multiply"),
        ),
    )
    def test_binary_operations(self, func, dtype):
        array = np.arange(10).astype(dtype) * unit_registry.m
        data_array = xr.DataArray(data=array)

        units = extract_units(func(array))
        with xr.set_options(use_opt_einsum=False):
            expected = attach_units(func(strip_units(data_array)), units)
            actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "comparison",
        (
            pytest.param(operator.lt, id="less_than"),
            pytest.param(operator.ge, id="greater_equal"),
            pytest.param(operator.eq, id="equal"),
        ),
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, ValueError, id="without_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.mm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_comparison_operations(self, comparison, unit, error, dtype):
        array = (
            np.array([10.1, 5.2, 6.5, 8.0, 21.3, 7.1, 1.3]).astype(dtype)
            * unit_registry.m
        )
        data_array = xr.DataArray(data=array)

        value = 8
        to_compare_with = value * unit

        # incompatible units are all not equal
        if error is not None and comparison is not operator.eq:
            with pytest.raises(error):
                comparison(array, to_compare_with)

            with pytest.raises(error):
                comparison(data_array, to_compare_with)

            return

        actual = comparison(data_array, to_compare_with)

        expected_units = {None: unit_registry.m if array.check(unit) else None}
        expected = array.check(unit) & comparison(
            strip_units(data_array),
            strip_units(convert_units(to_compare_with, expected_units)),
        )

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "units,error",
        (
            pytest.param(unit_registry.dimensionless, None, id="dimensionless"),
            pytest.param(unit_registry.m, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.degree, None, id="compatible_unit"),
        ),
    )
    def test_univariate_ufunc(self, units, error, dtype):
        array = np.arange(10).astype(dtype) * units
        data_array = xr.DataArray(data=array)

        func = function("sin")

        if error is not None:
            with pytest.raises(error):
                np.sin(data_array)

            return

        expected = attach_units(
            func(strip_units(convert_units(data_array, {None: unit_registry.radians}))),
            {None: unit_registry.dimensionless},
        )
        actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="without_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(
                unit_registry.mm,
                None,
                id="compatible_unit",
                marks=pytest.mark.xfail(reason="pint converts to the wrong units"),
            ),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_bivariate_ufunc(self, unit, error, dtype):
        original_unit = unit_registry.m
        array = np.arange(10).astype(dtype) * original_unit
        data_array = xr.DataArray(data=array)

        if error is not None:
            with pytest.raises(error):
                np.maximum(data_array, 1 * unit)

            return

        expected_units = {None: original_unit}
        expected = attach_units(
            np.maximum(
                strip_units(data_array),
                strip_units(convert_units(1 * unit, expected_units)),
            ),
            expected_units,
        )

        actual = np.maximum(data_array, 1 * unit)
        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

        actual = np.maximum(1 * unit, data_array)
        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("property", ("T", "imag", "real"))
    def test_numpy_properties(self, property, dtype):
        array = (
            np.arange(5 * 10).astype(dtype)
            + 1j * np.linspace(-1, 0, 5 * 10).astype(dtype)
        ).reshape(5, 10) * unit_registry.s

        data_array = xr.DataArray(data=array, dims=("x", "y"))

        expected = attach_units(
            getattr(strip_units(data_array), property), extract_units(data_array)
        )
        actual = getattr(data_array, property)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (method("conj"), method("argsort"), method("conjugate"), method("round")),
        ids=repr,
    )
    def test_numpy_methods(self, func, dtype):
        array = np.arange(10).astype(dtype) * unit_registry.m
        data_array = xr.DataArray(data=array, dims="x")

        units = extract_units(func(array))
        expected = attach_units(strip_units(data_array), units)
        actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    def test_item(self, dtype):
        array = np.arange(10).astype(dtype) * unit_registry.m
        data_array = xr.DataArray(data=array)

        func = method("item", 2)

        expected = func(strip_units(data_array)) * unit_registry.m
        actual = func(data_array)

        assert_duckarray_allclose(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("searchsorted", 5),
            pytest.param(
                function("searchsorted", 5),
                marks=pytest.mark.xfail(
                    reason="xarray does not implement __array_function__"
                ),
            ),
        ),
        ids=repr,
    )
    def test_searchsorted(self, func, unit, error, dtype):
        array = np.arange(10).astype(dtype) * unit_registry.m
        data_array = xr.DataArray(data=array)

        scalar_types = (int, float)
        args = [value * unit for value in func.args]
        kwargs = {
            key: (value * unit if isinstance(value, scalar_types) else value)
            for key, value in func.kwargs.items()
        }

        if error is not None:
            with pytest.raises(error):
                func(data_array, *args, **kwargs)

            return

        units = extract_units(data_array)
        expected_units = extract_units(func(array, *args, **kwargs))
        stripped_args = [strip_units(convert_units(value, units)) for value in args]
        stripped_kwargs = {
            key: strip_units(convert_units(value, units))
            for key, value in kwargs.items()
        }
        expected = attach_units(
            func(strip_units(data_array), *stripped_args, **stripped_kwargs),
            expected_units,
        )
        actual = func(data_array, *args, **kwargs)

        assert_units_equal(expected, actual)
        np.testing.assert_allclose(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            method("clip", min=3, max=8),
            pytest.param(
                function("clip", a_min=3, a_max=8),
                marks=pytest.mark.xfail(
                    reason="xarray does not implement __array_function__"
                ),
            ),
        ),
        ids=repr,
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_numpy_methods_with_args(self, func, unit, error, dtype):
        array = np.arange(10).astype(dtype) * unit_registry.m
        data_array = xr.DataArray(data=array)

        scalar_types = (int, float)
        args = [value * unit for value in func.args]
        kwargs = {
            key: (value * unit if isinstance(value, scalar_types) else value)
            for key, value in func.kwargs.items()
        }
        if error is not None:
            with pytest.raises(error):
                func(data_array, *args, **kwargs)

            return

        units = extract_units(data_array)
        expected_units = extract_units(func(array, *args, **kwargs))
        stripped_args = [strip_units(convert_units(value, units)) for value in args]
        stripped_kwargs = {
            key: strip_units(convert_units(value, units))
            for key, value in kwargs.items()
        }
        expected = attach_units(
            func(strip_units(data_array), *stripped_args, **stripped_kwargs),
            expected_units,
        )
        actual = func(data_array, *args, **kwargs)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "func", (method("isnull"), method("notnull"), method("count")), ids=repr
    )
    def test_missing_value_detection(self, func, dtype):
        array = (
            np.array(
                [
                    [1.4, 2.3, np.nan, 7.2],
                    [np.nan, 9.7, np.nan, np.nan],
                    [2.1, np.nan, np.nan, 4.6],
                    [9.9, np.nan, 7.2, 9.1],
                ]
            )
            * unit_registry.degK
        )
        data_array = xr.DataArray(data=array)

        expected = func(strip_units(data_array))
        actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.xfail(reason="ffill and bfill lose units in data")
    @pytest.mark.parametrize("func", (method("ffill"), method("bfill")), ids=repr)
    def test_missing_value_filling(self, func, dtype):
        array = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * unit_registry.degK
        )
        x = np.arange(len(array))
        data_array = xr.DataArray(data=array, coords={"x": x}, dims="x")

        expected = attach_units(
            func(strip_units(data_array), dim="x"), extract_units(data_array)
        )
        actual = func(data_array, dim="x")

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "fill_value",
        (
            pytest.param(-1, id="python_scalar"),
            pytest.param(np.array(-1), id="numpy_scalar"),
            pytest.param(np.array([-1]), id="numpy_array"),
        ),
    )
    def test_fillna(self, fill_value, unit, error, dtype):
        original_unit = unit_registry.m
        array = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * original_unit
        )
        data_array = xr.DataArray(data=array)

        func = method("fillna")

        value = fill_value * unit
        if error is not None:
            with pytest.raises(error):
                func(data_array, value=value)

            return

        units = extract_units(data_array)
        expected = attach_units(
            func(
                strip_units(data_array), value=strip_units(convert_units(value, units))
            ),
            units,
        )
        actual = func(data_array, value=value)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    def test_dropna(self, dtype):
        array = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * unit_registry.m
        )
        x = np.arange(len(array))
        data_array = xr.DataArray(data=array, coords={"x": x}, dims=["x"])

        units = extract_units(data_array)
        expected = attach_units(strip_units(data_array).dropna(dim="x"), units)
        actual = data_array.dropna(dim="x")

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    def test_isin(self, unit, dtype):
        array = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * unit_registry.m
        )
        data_array = xr.DataArray(data=array, dims="x")

        raw_values = create_nan_array([1.4, np.nan, 2.3], dtype)
        values = raw_values * unit

        units = {None: unit_registry.m if array.check(unit) else None}
        expected = strip_units(data_array).isin(
            strip_units(convert_units(values, units))
        ) & array.check(unit)
        actual = data_array.isin(values)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "variant", ("masking", "replacing_scalar", "replacing_array", "dropping")
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_where(self, variant, unit, error, dtype):
        original_unit = unit_registry.m
        array = np.linspace(0, 1, 10).astype(dtype) * original_unit

        data_array = xr.DataArray(data=array)

        condition = data_array < 0.5 * original_unit
        other = np.linspace(-2, -1, 10).astype(dtype) * unit
        variant_kwargs = {
            "masking": {"cond": condition},
            "replacing_scalar": {"cond": condition, "other": -1 * unit},
            "replacing_array": {"cond": condition, "other": other},
            "dropping": {"cond": condition, "drop": True},
        }
        kwargs = variant_kwargs[variant]
        kwargs_without_units = {
            key: strip_units(
                convert_units(
                    value, {None: original_unit if array.check(unit) else None}
                )
            )
            for key, value in kwargs.items()
        }

        if variant not in ("masking", "dropping") and error is not None:
            with pytest.raises(error):
                data_array.where(**kwargs)

            return

        expected = attach_units(
            strip_units(data_array).where(**kwargs_without_units),
            extract_units(data_array),
        )
        actual = data_array.where(**kwargs)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.xfail(reason="uses numpy.vectorize")
    def test_interpolate_na(self):
        array = (
            np.array([-1.03, 0.1, 1.4, np.nan, 2.3, np.nan, np.nan, 9.1])
            * unit_registry.m
        )
        x = np.arange(len(array))
        data_array = xr.DataArray(data=array, coords={"x": x}, dims="x")

        units = extract_units(data_array)
        expected = attach_units(strip_units(data_array).interpolate_na(dim="x"), units)
        actual = data_array.interpolate_na(dim="x")

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(
                unit_registry.cm,
                None,
                id="compatible_unit",
            ),
            pytest.param(
                unit_registry.m,
                None,
                id="identical_unit",
            ),
        ),
    )
    def test_combine_first(self, unit, error, dtype):
        array = np.zeros(shape=(2, 2), dtype=dtype) * unit_registry.m
        other_array = np.ones_like(array) * unit

        data_array = xr.DataArray(
            data=array, coords={"x": ["a", "b"], "y": [-1, 0]}, dims=["x", "y"]
        )
        other = xr.DataArray(
            data=other_array, coords={"x": ["b", "c"], "y": [0, 1]}, dims=["x", "y"]
        )

        if error is not None:
            with pytest.raises(error):
                data_array.combine_first(other)

            return

        units = extract_units(data_array)
        expected = attach_units(
            strip_units(data_array).combine_first(
                strip_units(convert_units(other, units))
            ),
            units,
        )
        actual = data_array.combine_first(other)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "variation",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("equals"),
            pytest.param(
                method("identical"),
                marks=pytest.mark.skip(reason="the behavior of identical is undecided"),
            ),
        ),
        ids=repr,
    )
    def test_comparisons(self, func, variation, unit, dtype):
        def is_compatible(a, b):
            a = a if a is not None else 1
            b = b if b is not None else 1
            quantity = np.arange(5) * a

            return a == b or quantity.check(b)

        data = np.linspace(0, 5, 10).astype(dtype)
        coord = np.arange(len(data)).astype(dtype)

        base_unit = unit_registry.m
        array = data * (base_unit if variation == "data" else 1)
        x = coord * (base_unit if variation == "dims" else 1)
        y = coord * (base_unit if variation == "coords" else 1)

        variations = {
            "data": (unit, 1, 1),
            "dims": (1, unit, 1),
            "coords": (1, 1, unit),
        }
        data_unit, dim_unit, coord_unit = variations[variation]

        data_array = xr.DataArray(data=array, coords={"x": x, "y": ("x", y)}, dims="x")

        other = attach_units(
            strip_units(data_array), {None: data_unit, "x": dim_unit, "y": coord_unit}
        )

        units = extract_units(data_array)
        other_units = extract_units(other)

        equal_arrays = all(
            is_compatible(units[name], other_units[name]) for name in units.keys()
        ) and (
            strip_units(data_array).equals(
                strip_units(convert_units(other, extract_units(data_array)))
            )
        )
        equal_units = units == other_units
        expected = equal_arrays and (func.name != "identical" or equal_units)

        actual = func(data_array, other)

        assert expected == actual

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    def test_broadcast_like(self, variant, unit, dtype):
        original_unit = unit_registry.m

        variants = {
            "data": ((original_unit, unit), (1, 1), (1, 1)),
            "dims": ((1, 1), (original_unit, unit), (1, 1)),
            "coords": ((1, 1), (1, 1), (original_unit, unit)),
        }
        (
            (data_unit1, data_unit2),
            (dim_unit1, dim_unit2),
            (coord_unit1, coord_unit2),
        ) = variants[variant]

        array1 = np.linspace(1, 2, 2 * 1).reshape(2, 1).astype(dtype) * data_unit1
        array2 = np.linspace(0, 1, 2 * 3).reshape(2, 3).astype(dtype) * data_unit2

        x1 = np.arange(2) * dim_unit1
        x2 = np.arange(2) * dim_unit2
        y1 = np.array([0]) * dim_unit1
        y2 = np.arange(3) * dim_unit2

        u1 = np.linspace(0, 1, 2) * coord_unit1
        u2 = np.linspace(0, 1, 2) * coord_unit2

        arr1 = xr.DataArray(
            data=array1, coords={"x": x1, "y": y1, "u": ("x", u1)}, dims=("x", "y")
        )
        arr2 = xr.DataArray(
            data=array2, coords={"x": x2, "y": y2, "u": ("x", u2)}, dims=("x", "y")
        )

        expected = attach_units(
            strip_units(arr1).broadcast_like(strip_units(arr2)), extract_units(arr1)
        )
        actual = arr1.broadcast_like(arr2)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    def test_broadcast_equals(self, unit, dtype):
        left_array = np.ones(shape=(2, 2), dtype=dtype) * unit_registry.m
        right_array = np.ones(shape=(2,), dtype=dtype) * unit

        left = xr.DataArray(data=left_array, dims=("x", "y"))
        right = xr.DataArray(data=right_array, dims="x")

        units = {
            **extract_units(left),
            **({} if left_array.check(unit) else {None: None}),
        }
        expected = strip_units(left).broadcast_equals(
            strip_units(convert_units(right, units))
        ) & left_array.check(unit)
        actual = left.broadcast_equals(right)

        assert expected == actual

    def test_pad(self, dtype):
        array = np.linspace(0, 5, 10).astype(dtype) * unit_registry.m

        data_array = xr.DataArray(data=array, dims="x")
        units = extract_units(data_array)

        expected = attach_units(strip_units(data_array).pad(x=(2, 3)), units)
        actual = data_array.pad(x=(2, 3))

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("pipe", lambda da: da * 10),
            method("assign_coords", w=("y", np.arange(10) * unit_registry.mm)),
            method("assign_attrs", attr1="value"),
            method("rename", u="v"),
            pytest.param(
                method("swap_dims", {"x": "u"}),
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
            pytest.param(
                method(
                    "expand_dims",
                    dim={"z": np.linspace(10, 20, 12) * unit_registry.s},
                    axis=1,
                ),
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
            method("drop_vars", "x"),
            method("reset_coords", names="u"),
            method("copy"),
            method("astype", np.float32),
        ),
        ids=repr,
    )
    def test_content_manipulation(self, func, variant, dtype):
        unit = unit_registry.m

        variants = {
            "data": (unit, 1, 1),
            "dims": (1, unit, 1),
            "coords": (1, 1, unit),
        }
        data_unit, dim_unit, coord_unit = variants[variant]
        quantity = np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype) * data_unit
        x = np.arange(quantity.shape[0]) * dim_unit
        y = np.arange(quantity.shape[1]) * dim_unit
        u = np.linspace(0, 1, quantity.shape[0]) * coord_unit

        data_array = xr.DataArray(
            name="a",
            data=quantity,
            coords={"x": x, "u": ("x", u), "y": y},
            dims=("x", "y"),
        )

        stripped_kwargs = {
            key: array_strip_units(value) for key, value in func.kwargs.items()
        }
        units = extract_units(data_array)
        units["u"] = getattr(u, "units", None)
        units["v"] = getattr(u, "units", None)

        expected = attach_units(func(strip_units(data_array), **stripped_kwargs), units)
        actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.degK, id="with_unit"),
        ),
    )
    def test_copy(self, unit, dtype):
        quantity = np.linspace(0, 10, 20, dtype=dtype) * unit_registry.pascal
        new_data = np.arange(20)

        data_array = xr.DataArray(data=quantity, dims="x")

        expected = attach_units(
            strip_units(data_array).copy(data=new_data), {None: unit}
        )

        actual = data_array.copy(data=new_data * unit)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "indices",
        (
            pytest.param(4, id="single index"),
            pytest.param([5, 2, 9, 1], id="multiple indices"),
        ),
    )
    def test_isel(self, indices, dtype):
        # TODO: maybe test for units in indexes?
        array = np.arange(10).astype(dtype) * unit_registry.s

        data_array = xr.DataArray(data=array, dims="x")

        expected = attach_units(
            strip_units(data_array).isel(x=indices), extract_units(data_array)
        )
        actual = data_array.isel(x=indices)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "raw_values",
        (
            pytest.param(10, id="single_value"),
            pytest.param([10, 5, 13], id="list_of_values"),
            pytest.param(np.array([9, 3, 7, 12]), id="array_of_values"),
        ),
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, KeyError, id="no_units"),
            pytest.param(unit_registry.dimensionless, KeyError, id="dimensionless"),
            pytest.param(unit_registry.degree, KeyError, id="incompatible_unit"),
            pytest.param(unit_registry.dm, KeyError, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_sel(self, raw_values, unit, error, dtype):
        array = np.linspace(5, 10, 20).astype(dtype) * unit_registry.m
        x = np.arange(len(array)) * unit_registry.m
        data_array = xr.DataArray(data=array, coords={"x": x}, dims="x")

        values = raw_values * unit

        if error is not None and not (
            isinstance(raw_values, int | float) and x.check(unit)
        ):
            with pytest.raises(error):
                data_array.sel(x=values)

            return

        expected = attach_units(
            strip_units(data_array).sel(
                x=strip_units(convert_units(values, {None: array.units}))
            ),
            extract_units(data_array),
        )
        actual = data_array.sel(x=values)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "raw_values",
        (
            pytest.param(10, id="single_value"),
            pytest.param([10, 5, 13], id="list_of_values"),
            pytest.param(np.array([9, 3, 7, 12]), id="array_of_values"),
        ),
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, KeyError, id="no_units"),
            pytest.param(unit_registry.dimensionless, KeyError, id="dimensionless"),
            pytest.param(unit_registry.degree, KeyError, id="incompatible_unit"),
            pytest.param(unit_registry.dm, KeyError, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_loc(self, raw_values, unit, error, dtype):
        array = np.linspace(5, 10, 20).astype(dtype) * unit_registry.m
        x = np.arange(len(array)) * unit_registry.m
        data_array = xr.DataArray(data=array, coords={"x": x}, dims="x")

        values = raw_values * unit

        if error is not None and not (
            isinstance(raw_values, int | float) and x.check(unit)
        ):
            with pytest.raises(error):
                data_array.loc[{"x": values}]

            return

        expected = attach_units(
            strip_units(data_array).loc[
                {"x": strip_units(convert_units(values, {None: array.units}))}
            ],
            extract_units(data_array),
        )
        actual = data_array.loc[{"x": values}]

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "raw_values",
        (
            pytest.param(10, id="single_value"),
            pytest.param([10, 5, 13], id="list_of_values"),
            pytest.param(np.array([9, 3, 7, 12]), id="array_of_values"),
        ),
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, KeyError, id="no_units"),
            pytest.param(unit_registry.dimensionless, KeyError, id="dimensionless"),
            pytest.param(unit_registry.degree, KeyError, id="incompatible_unit"),
            pytest.param(unit_registry.dm, KeyError, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_drop_sel(self, raw_values, unit, error, dtype):
        array = np.linspace(5, 10, 20).astype(dtype) * unit_registry.m
        x = np.arange(len(array)) * unit_registry.m
        data_array = xr.DataArray(data=array, coords={"x": x}, dims="x")

        values = raw_values * unit

        if error is not None and not (
            isinstance(raw_values, int | float) and x.check(unit)
        ):
            with pytest.raises(error):
                data_array.drop_sel(x=values)

            return

        expected = attach_units(
            strip_units(data_array).drop_sel(
                x=strip_units(convert_units(values, {None: x.units}))
            ),
            extract_units(data_array),
        )
        actual = data_array.drop_sel(x=values)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("dim", ("x", "y", "z", "t", "all"))
    @pytest.mark.parametrize(
        "shape",
        (
            pytest.param((10, 20), id="nothing_squeezable"),
            pytest.param((10, 20, 1), id="last_dimension_squeezable"),
            pytest.param((10, 1, 20), id="middle_dimension_squeezable"),
            pytest.param((1, 10, 20), id="first_dimension_squeezable"),
            pytest.param((1, 10, 1, 20), id="first_and_last_dimension_squeezable"),
        ),
    )
    def test_squeeze(self, shape, dim, dtype):
        names = "xyzt"
        dim_lengths = dict(zip(names, shape, strict=False))
        names = "xyzt"
        array = np.arange(10 * 20).astype(dtype).reshape(shape) * unit_registry.J
        data_array = xr.DataArray(data=array, dims=tuple(names[: len(shape)]))

        kwargs = {"dim": dim} if dim != "all" and dim_lengths.get(dim, 0) == 1 else {}

        expected = attach_units(
            strip_units(data_array).squeeze(**kwargs), extract_units(data_array)
        )
        actual = data_array.squeeze(**kwargs)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (method("head", x=7, y=3), method("tail", x=7, y=3), method("thin", x=7, y=3)),
        ids=repr,
    )
    def test_head_tail_thin(self, func, dtype):
        # TODO: works like isel. Maybe also test units in indexes?
        array = np.linspace(1, 2, 10 * 5).reshape(10, 5) * unit_registry.degK

        data_array = xr.DataArray(data=array, dims=("x", "y"))

        expected = attach_units(
            func(strip_units(data_array)), extract_units(data_array)
        )
        actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("variant", ("data", "coords"))
    @pytest.mark.parametrize(
        "func",
        (
            pytest.param(
                method("interp"), marks=pytest.mark.xfail(reason="uses scipy")
            ),
            method("reindex"),
        ),
        ids=repr,
    )
    def test_interp_reindex(self, variant, func, dtype):
        variants = {
            "data": (unit_registry.m, 1),
            "coords": (1, unit_registry.m),
        }
        data_unit, coord_unit = variants[variant]

        array = np.linspace(1, 2, 10).astype(dtype) * data_unit
        y = np.arange(10) * coord_unit

        x = np.arange(10)
        new_x = np.arange(10) + 0.5
        data_array = xr.DataArray(array, coords={"x": x, "y": ("x", y)}, dims="x")

        units = extract_units(data_array)
        expected = attach_units(func(strip_units(data_array), x=new_x), units)
        actual = func(data_array, x=new_x)

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (method("interp"), method("reindex")),
        ids=repr,
    )
    def test_interp_reindex_indexing(self, func, unit, error, dtype):
        array = np.linspace(1, 2, 10).astype(dtype)
        x = np.arange(10) * unit_registry.m
        new_x = (np.arange(10) + 0.5) * unit
        data_array = xr.DataArray(array, coords={"x": x}, dims="x")

        if error is not None:
            with pytest.raises(error):
                func(data_array, x=new_x)

            return

        units = extract_units(data_array)
        expected = attach_units(
            func(
                strip_units(data_array),
                x=strip_units(convert_units(new_x, {None: unit_registry.m})),
            ),
            units,
        )
        actual = func(data_array, x=new_x)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("variant", ("data", "coords"))
    @pytest.mark.parametrize(
        "func",
        (
            pytest.param(
                method("interp_like"), marks=pytest.mark.xfail(reason="uses scipy")
            ),
            method("reindex_like"),
        ),
        ids=repr,
    )
    def test_interp_reindex_like(self, variant, func, dtype):
        variants = {
            "data": (unit_registry.m, 1),
            "coords": (1, unit_registry.m),
        }
        data_unit, coord_unit = variants[variant]

        array = np.linspace(1, 2, 10).astype(dtype) * data_unit
        coord = np.arange(10) * coord_unit

        x = np.arange(10)
        new_x = np.arange(-2, 2) + 0.5
        data_array = xr.DataArray(array, coords={"x": x, "y": ("x", coord)}, dims="x")
        other = xr.DataArray(np.empty_like(new_x), coords={"x": new_x}, dims="x")

        units = extract_units(data_array)
        expected = attach_units(func(strip_units(data_array), other), units)
        actual = func(data_array, other)

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (method("interp_like"), method("reindex_like")),
        ids=repr,
    )
    def test_interp_reindex_like_indexing(self, func, unit, error, dtype):
        array = np.linspace(1, 2, 10).astype(dtype)
        x = np.arange(10) * unit_registry.m
        new_x = (np.arange(-2, 2) + 0.5) * unit

        data_array = xr.DataArray(array, coords={"x": x}, dims="x")
        other = xr.DataArray(np.empty_like(new_x), {"x": new_x}, dims="x")

        if error is not None:
            with pytest.raises(error):
                func(data_array, other)

            return

        units = extract_units(data_array)
        expected = attach_units(
            func(
                strip_units(data_array),
                strip_units(convert_units(other, {None: unit_registry.m})),
            ),
            units,
        )
        actual = func(data_array, other)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (method("unstack"), method("reset_index", "z"), method("reorder_levels")),
        ids=repr,
    )
    def test_stacking_stacked(self, func, dtype):
        array = (
            np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype) * unit_registry.m
        )
        x = np.arange(array.shape[0])
        y = np.arange(array.shape[1])

        data_array = xr.DataArray(
            name="data", data=array, coords={"x": x, "y": y}, dims=("x", "y")
        )
        stacked = data_array.stack(z=("x", "y"))

        expected = attach_units(func(strip_units(stacked)), {"data": unit_registry.m})
        actual = func(stacked)

        assert_units_equal(expected, actual)
        # TODO: strip_units/attach_units reconstruct DataArrays from scratch,
        # losing index structure (e.g., MultiIndex from stack becomes regular Index).
        # Fix these utilities to preserve indexes, then remove check_indexes=False.
        if func.name == "reset_index":
            assert_identical(
                expected, actual, check_default_indexes=False, check_indexes=False
            )
        else:
            assert_identical(expected, actual, check_indexes=False)

    @pytest.mark.skip(reason="indexes don't support units")
    def test_to_unstacked_dataset(self, dtype):
        array = (
            np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype)
            * unit_registry.pascal
        )
        x = np.arange(array.shape[0]) * unit_registry.m
        y = np.arange(array.shape[1]) * unit_registry.s

        data_array = xr.DataArray(
            data=array, coords={"x": x, "y": y}, dims=("x", "y")
        ).stack(z=("x", "y"))

        func = method("to_unstacked_dataset", dim="z")

        expected = attach_units(
            func(strip_units(data_array)),
            {
                "y": y.units,
                **dict(zip(x.magnitude, [array.units] * len(y), strict=True)),
            },
        ).rename({elem.magnitude: elem for elem in x})
        actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            method("transpose", "y", "x", "z"),
            method("stack", a=("x", "y")),
            method("set_index", x="x2"),
            method("shift", x=2),
            pytest.param(
                method("rank", dim="x"),
                marks=pytest.mark.skip(reason="rank not implemented for non-ndarray"),
            ),
            method("roll", x=2, roll_coords=False),
            method("sortby", "x2"),
        ),
        ids=repr,
    )
    def test_stacking_reordering(self, func, dtype):
        array = (
            np.linspace(0, 10, 2 * 5 * 10).reshape(2, 5, 10).astype(dtype)
            * unit_registry.m
        )
        x = np.arange(array.shape[0])
        y = np.arange(array.shape[1])
        z = np.arange(array.shape[2])
        x2 = np.linspace(0, 1, array.shape[0])[::-1]

        data_array = xr.DataArray(
            name="data",
            data=array,
            coords={"x": x, "y": y, "z": z, "x2": ("x", x2)},
            dims=("x", "y", "z"),
        )

        expected = attach_units(func(strip_units(data_array)), {None: unit_registry.m})
        actual = func(data_array)

        assert_units_equal(expected, actual)
        # TODO: strip_units/attach_units reconstruct DataArrays from scratch,
        # losing index structure (e.g., MultiIndex from stack becomes regular Index).
        # Fix these utilities to preserve indexes, then remove check_indexes=False.
        assert_identical(expected, actual, check_indexes=False)

    @pytest.mark.parametrize(
        "variant",
        (
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("differentiate", fallback_func=np.gradient),
            method("integrate", fallback_func=duck_array_ops.cumulative_trapezoid),
            method("cumulative_integrate", fallback_func=duck_array_ops.trapz),
        ),
        ids=repr,
    )
    def test_differentiate_integrate(self, func, variant, dtype):
        data_unit = unit_registry.m
        unit = unit_registry.s

        variants = {
            "dims": ("x", unit, 1),
            "coords": ("u", 1, unit),
        }
        coord, dim_unit, coord_unit = variants[variant]

        array = np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype) * data_unit

        x = np.arange(array.shape[0]) * dim_unit
        y = np.arange(array.shape[1]) * dim_unit

        u = np.linspace(0, 1, array.shape[0]) * coord_unit

        data_array = xr.DataArray(
            data=array, coords={"x": x, "y": y, "u": ("x", u)}, dims=("x", "y")
        )
        # we want to make sure the output unit is correct
        units = extract_units(data_array)
        units.update(
            extract_units(
                func(
                    data_array.data,
                    getattr(data_array, coord).data,
                    axis=0,
                )
            )
        )

        expected = attach_units(
            func(strip_units(data_array), coord=strip_units(coord)),
            units,
        )
        actual = func(data_array, coord=coord)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("compute_backend", ["numbagg", None], indirect=True)
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("diff", dim="x"),
            method("quantile", q=[0.25, 0.75]),
            method("reduce", func=np.sum, dim="x"),
            pytest.param(lambda x: x.dot(x), id="method_dot"),
        ),
        ids=repr,
    )
    def test_computation(self, func, variant, dtype, compute_backend):
        unit = unit_registry.m

        variants = {
            "data": (unit, 1, 1),
            "dims": (1, unit, 1),
            "coords": (1, 1, unit),
        }
        data_unit, dim_unit, coord_unit = variants[variant]
        array = np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype) * data_unit

        x = np.arange(array.shape[0]) * dim_unit
        y = np.arange(array.shape[1]) * dim_unit

        u = np.linspace(0, 1, array.shape[0]) * coord_unit

        data_array = xr.DataArray(
            data=array, coords={"x": x, "y": y, "u": ("x", u)}, dims=("x", "y")
        )

        # we want to make sure the output unit is correct
        units = extract_units(data_array)
        if not isinstance(func, function | method):
            units.update(extract_units(func(array.reshape(-1))))

        with xr.set_options(use_opt_einsum=False):
            expected = attach_units(func(strip_units(data_array)), units)
            actual = func(data_array)

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("groupby", "x"),
            method("groupby_bins", "y", bins=4),
            method("coarsen", y=2),
            method("rolling", y=3),
            pytest.param(method("rolling_exp", y=3), marks=requires_numbagg),
            method("weighted", xr.DataArray(data=np.linspace(0, 1, 10), dims="y")),
        ),
        ids=repr,
    )
    def test_computation_objects(self, func, variant, dtype):
        if variant == "data":
            if func.name == "rolling_exp":
                pytest.xfail(reason="numbagg functions are not supported by pint")
            elif func.name == "rolling":
                pytest.xfail(
                    reason="numpy.lib.stride_tricks.as_strided converts to ndarray"
                )

        unit = unit_registry.m

        variants = {
            "data": (unit, 1, 1),
            "dims": (1, unit, 1),
            "coords": (1, 1, unit),
        }
        data_unit, dim_unit, coord_unit = variants[variant]
        array = np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype) * data_unit

        x = np.array([0, 0, 1, 2, 2]) * dim_unit
        y = np.arange(array.shape[1]) * 3 * dim_unit

        u = np.linspace(0, 1, 5) * coord_unit

        data_array = xr.DataArray(
            data=array, coords={"x": x, "y": y, "u": ("x", u)}, dims=("x", "y")
        )
        units = extract_units(data_array)

        expected = attach_units(func(strip_units(data_array)).mean(), units)
        actual = func(data_array).mean()

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    def test_resample(self, dtype):
        array = np.linspace(0, 5, 10).astype(dtype) * unit_registry.m

        time = xr.date_range("10-09-2010", periods=len(array), freq="YE")
        data_array = xr.DataArray(data=array, coords={"time": time}, dims="time")
        units = extract_units(data_array)

        func = method("resample", time="6ME")

        expected = attach_units(func(strip_units(data_array)).mean(), units)
        actual = func(data_array).mean()

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("compute_backend", ["numbagg", None], indirect=True)
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("assign_coords", z=("x", np.arange(5) * unit_registry.s)),
            method("first"),
            method("last"),
            method("quantile", q=[0.25, 0.5, 0.75], dim="x"),
        ),
        ids=repr,
    )
    def test_grouped_operations(self, func, variant, dtype, compute_backend):
        unit = unit_registry.m

        variants = {
            "data": (unit, 1, 1),
            "dims": (1, unit, 1),
            "coords": (1, 1, unit),
        }
        data_unit, dim_unit, coord_unit = variants[variant]
        array = np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype) * data_unit

        x = np.arange(array.shape[0]) * dim_unit
        y = np.arange(array.shape[1]) * 3 * dim_unit

        u = np.linspace(0, 1, array.shape[0]) * coord_unit

        data_array = xr.DataArray(
            data=array, coords={"x": x, "y": y, "u": ("x", u)}, dims=("x", "y")
        )
        units = {**extract_units(data_array), "z": unit_registry.s, "q": None}

        stripped_kwargs = {
            key: (
                strip_units(value)
                if not isinstance(value, tuple)
                else tuple(strip_units(elem) for elem in value)
            )
            for key, value in func.kwargs.items()
        }
        expected = attach_units(
            func(
                strip_units(data_array).groupby("y", squeeze=False), **stripped_kwargs
            ),
            units,
        )
        actual = func(data_array.groupby("y", squeeze=False))

        assert_units_equal(expected, actual)
        assert_identical(expected, actual)


class TestDataset:
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, xr.MergeError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, xr.MergeError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, xr.MergeError, id="incompatible_unit"),
            pytest.param(unit_registry.mm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="same_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "shared",
        (
            "nothing",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    def test_init(self, shared, unit, error, dtype):
        original_unit = unit_registry.m
        scaled_unit = unit_registry.mm

        a = np.linspace(0, 1, 10).astype(dtype) * unit_registry.Pa
        b = np.linspace(-1, 0, 10).astype(dtype) * unit_registry.degK

        values_a = np.arange(a.shape[0])
        dim_a = values_a * original_unit
        coord_a = dim_a.to(scaled_unit)

        values_b = np.arange(b.shape[0])
        dim_b = values_b * unit
        coord_b = (
            dim_b.to(scaled_unit)
            if unit_registry.is_compatible_with(dim_b, scaled_unit)
            and unit != scaled_unit
            else dim_b * 1000
        )

        variants = {
            "nothing": ({}, {}),
            "dims": ({"x": dim_a}, {"x": dim_b}),
            "coords": (
                {"x": values_a, "y": ("x", coord_a)},
                {"x": values_b, "y": ("x", coord_b)},
            ),
        }
        coords_a, coords_b = variants[shared]

        dims_a, dims_b = ("x", "y") if shared == "nothing" else ("x", "x")

        a = xr.DataArray(data=a, coords=coords_a, dims=dims_a)
        b = xr.DataArray(data=b, coords=coords_b, dims=dims_b)

        if error is not None and shared != "nothing":
            with pytest.raises(error):
                xr.Dataset(data_vars={"a": a, "b": b})

            return

        actual = xr.Dataset(data_vars={"a": a, "b": b})

        units = merge_mappings(
            extract_units(a.rename("a")), extract_units(b.rename("b"))
        )
        expected = attach_units(
            xr.Dataset(data_vars={"a": strip_units(a), "b": strip_units(b)}), units
        )

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "func", (pytest.param(str, id="str"), pytest.param(repr, id="repr"))
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims",
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
            "coords",
        ),
    )
    def test_repr(self, func, variant, dtype):
        unit1, unit2 = (
            (unit_registry.Pa, unit_registry.degK) if variant == "data" else (1, 1)
        )

        array1 = np.linspace(1, 2, 10, dtype=dtype) * unit1
        array2 = np.linspace(0, 1, 10, dtype=dtype) * unit2

        x = np.arange(len(array1)) * unit_registry.s
        y = x.to(unit_registry.ms)

        variants = {
            "dims": {"x": x},
            "coords": {"y": ("x", y)},
            "data": {},
        }

        ds = xr.Dataset(
            data_vars={"a": ("x", array1), "b": ("x", array2)},
            coords=variants[variant],
        )

        # FIXME: this just checks that the repr does not raise
        # warnings or errors, but does not check the result
        func(ds)

    @pytest.mark.parametrize(
        "func",
        (
            method("all"),
            method("any"),
            method("argmax", dim="x"),
            method("argmin", dim="x"),
            method("max"),
            method("min"),
            method("mean"),
            method("median"),
            method("sum"),
            method("prod"),
            method("std"),
            method("var"),
            method("cumsum"),
            method("cumprod"),
        ),
        ids=repr,
    )
    def test_aggregation(self, func, dtype):
        unit_a, unit_b = (
            (unit_registry.Pa, unit_registry.degK)
            if func.name != "cumprod"
            else (unit_registry.dimensionless, unit_registry.dimensionless)
        )

        a = np.linspace(0, 1, 10).astype(dtype) * unit_a
        b = np.linspace(-1, 0, 10).astype(dtype) * unit_b

        ds = xr.Dataset({"a": ("x", a), "b": ("x", b)})

        if "dim" in func.kwargs:
            numpy_kwargs = func.kwargs.copy()
            dim = numpy_kwargs.pop("dim")

            axis_a = ds.a.get_axis_num(dim)
            axis_b = ds.b.get_axis_num(dim)

            numpy_kwargs_a = numpy_kwargs.copy()
            numpy_kwargs_a["axis"] = axis_a
            numpy_kwargs_b = numpy_kwargs.copy()
            numpy_kwargs_b["axis"] = axis_b
        else:
            numpy_kwargs_a = {}
            numpy_kwargs_b = {}

        units_a = array_extract_units(func(a, **numpy_kwargs_a))
        units_b = array_extract_units(func(b, **numpy_kwargs_b))
        units = {"a": units_a, "b": units_b}

        actual = func(ds)
        expected = attach_units(func(strip_units(ds)), units)

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    @pytest.mark.parametrize("property", ("imag", "real"))
    def test_numpy_properties(self, property, dtype):
        a = np.linspace(0, 1, 10) * unit_registry.Pa
        b = np.linspace(-1, 0, 15) * unit_registry.degK
        ds = xr.Dataset({"a": ("x", a), "b": ("y", b)})
        units = extract_units(ds)

        actual = getattr(ds, property)
        expected = attach_units(getattr(strip_units(ds), property), units)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            method("astype", float),
            method("conj"),
            method("argsort"),
            method("conjugate"),
            method("round"),
        ),
        ids=repr,
    )
    def test_numpy_methods(self, func, dtype):
        a = np.linspace(1, -1, 10) * unit_registry.Pa
        b = np.linspace(-1, 1, 15) * unit_registry.degK
        ds = xr.Dataset({"a": ("x", a), "b": ("y", b)})

        units_a = array_extract_units(func(a))
        units_b = array_extract_units(func(b))
        units = {"a": units_a, "b": units_b}

        actual = func(ds)
        expected = attach_units(func(strip_units(ds)), units)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize("func", (method("clip", min=3, max=8),), ids=repr)
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_numpy_methods_with_args(self, func, unit, error, dtype):
        data_unit = unit_registry.m
        a = np.linspace(0, 10, 15) * unit_registry.m
        b = np.linspace(-2, 12, 20) * unit_registry.m
        ds = xr.Dataset({"a": ("x", a), "b": ("y", b)})
        units = extract_units(ds)

        kwargs = {
            key: array_attach_units(value, unit) for key, value in func.kwargs.items()
        }

        if error is not None:
            with pytest.raises(error):
                func(ds, **kwargs)

            return

        stripped_kwargs = {
            key: strip_units(convert_units(value, {None: data_unit}))
            for key, value in kwargs.items()
        }

        actual = func(ds, **kwargs)
        expected = attach_units(func(strip_units(ds), **stripped_kwargs), units)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "func", (method("isnull"), method("notnull"), method("count")), ids=repr
    )
    def test_missing_value_detection(self, func, dtype):
        array1 = (
            np.array(
                [
                    [1.4, 2.3, np.nan, 7.2],
                    [np.nan, 9.7, np.nan, np.nan],
                    [2.1, np.nan, np.nan, 4.6],
                    [9.9, np.nan, 7.2, 9.1],
                ]
            )
            * unit_registry.degK
        )
        array2 = (
            np.array(
                [
                    [np.nan, 5.7, 12.0, 7.2],
                    [np.nan, 12.4, np.nan, 4.2],
                    [9.8, np.nan, 4.6, 1.4],
                    [7.2, np.nan, 6.3, np.nan],
                    [8.4, 3.9, np.nan, np.nan],
                ]
            )
            * unit_registry.Pa
        )

        ds = xr.Dataset({"a": (("x", "y"), array1), "b": (("z", "x"), array2)})

        expected = func(strip_units(ds))
        actual = func(ds)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.xfail(reason="ffill and bfill lose the unit")
    @pytest.mark.parametrize("func", (method("ffill"), method("bfill")), ids=repr)
    def test_missing_value_filling(self, func, dtype):
        array1 = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * unit_registry.degK
        )
        array2 = (
            create_nan_array([4.3, 9.8, 7.5, np.nan, 8.2, np.nan], dtype)
            * unit_registry.Pa
        )

        ds = xr.Dataset({"a": ("x", array1), "b": ("y", array2)})
        units = extract_units(ds)

        expected = attach_units(func(strip_units(ds), dim="x"), units)
        actual = func(ds, dim="x")

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(
                unit_registry.cm,
                None,
                id="compatible_unit",
            ),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "fill_value",
        (
            pytest.param(-1, id="python_scalar"),
            pytest.param(np.array(-1), id="numpy_scalar"),
            pytest.param(np.array([-1]), id="numpy_array"),
        ),
    )
    def test_fillna(self, fill_value, unit, error, dtype):
        array1 = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * unit_registry.m
        )
        array2 = (
            create_nan_array([4.3, 9.8, 7.5, np.nan, 8.2, np.nan], dtype)
            * unit_registry.m
        )
        ds = xr.Dataset({"a": ("x", array1), "b": ("x", array2)})
        value = fill_value * unit
        units = extract_units(ds)

        if error is not None:
            with pytest.raises(error):
                ds.fillna(value=value)

            return

        actual = ds.fillna(value=value)
        expected = attach_units(
            strip_units(ds).fillna(
                value=strip_units(convert_units(value, {None: unit_registry.m}))
            ),
            units,
        )

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    def test_dropna(self, dtype):
        array1 = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * unit_registry.degK
        )
        array2 = (
            create_nan_array([4.3, 9.8, 7.5, np.nan, 8.2, np.nan], dtype)
            * unit_registry.Pa
        )
        ds = xr.Dataset({"a": ("x", array1), "b": ("x", array2)})
        units = extract_units(ds)

        expected = attach_units(strip_units(ds).dropna(dim="x"), units)
        actual = ds.dropna(dim="x")

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="same_unit"),
        ),
    )
    def test_isin(self, unit, dtype):
        array1 = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * unit_registry.m
        )
        array2 = (
            create_nan_array([4.3, 9.8, 7.5, np.nan, 8.2, np.nan], dtype)
            * unit_registry.m
        )
        ds = xr.Dataset({"a": ("x", array1), "b": ("x", array2)})

        raw_values = create_nan_array([1.4, np.nan, 2.3], dtype)
        values = raw_values * unit

        converted_values = (
            convert_units(values, {None: unit_registry.m})
            if is_compatible(unit, unit_registry.m)
            else values
        )

        expected = strip_units(ds).isin(strip_units(converted_values))
        # TODO: use `unit_registry.is_compatible_with(unit, unit_registry.m)` instead.
        # Needs `pint>=0.12.1`, though, so we probably should wait until that is released.
        if not is_compatible(unit, unit_registry.m):
            expected.a[:] = False
            expected.b[:] = False

        actual = ds.isin(values)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "variant", ("masking", "replacing_scalar", "replacing_array", "dropping")
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="same_unit"),
        ),
    )
    def test_where(self, variant, unit, error, dtype):
        original_unit = unit_registry.m
        array1 = np.linspace(0, 1, 10).astype(dtype) * original_unit
        array2 = np.linspace(-1, 0, 10).astype(dtype) * original_unit

        ds = xr.Dataset({"a": ("x", array1), "b": ("x", array2)})
        units = extract_units(ds)

        condition = ds < 0.5 * original_unit
        other = np.linspace(-2, -1, 10).astype(dtype) * unit
        variant_kwargs = {
            "masking": {"cond": condition},
            "replacing_scalar": {"cond": condition, "other": -1 * unit},
            "replacing_array": {"cond": condition, "other": other},
            "dropping": {"cond": condition, "drop": True},
        }
        kwargs = variant_kwargs[variant]
        if variant not in ("masking", "dropping") and error is not None:
            with pytest.raises(error):
                ds.where(**kwargs)

            return

        kwargs_without_units = {
            key: strip_units(convert_units(value, {None: original_unit}))
            for key, value in kwargs.items()
        }

        expected = attach_units(
            strip_units(ds).where(**kwargs_without_units),
            units,
        )
        actual = ds.where(**kwargs)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.xfail(reason="interpolate_na uses numpy.vectorize")
    def test_interpolate_na(self, dtype):
        array1 = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype)
            * unit_registry.degK
        )
        array2 = (
            create_nan_array([4.3, 9.8, 7.5, np.nan, 8.2, np.nan], dtype)
            * unit_registry.Pa
        )
        ds = xr.Dataset({"a": ("x", array1), "b": ("x", array2)})
        units = extract_units(ds)

        expected = attach_units(
            strip_units(ds).interpolate_na(dim="x"),
            units,
        )
        actual = ds.interpolate_na(dim="x")

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="same_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims",
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
        ),
    )
    def test_combine_first(self, variant, unit, error, dtype):
        variants = {
            "data": (unit_registry.m, unit, 1, 1),
            "dims": (1, 1, unit_registry.m, unit),
        }
        data_unit, other_data_unit, dims_unit, other_dims_unit = variants[variant]

        array1 = (
            create_nan_array([1.4, np.nan, 2.3, np.nan, np.nan, 9.1], dtype) * data_unit
        )
        array2 = (
            create_nan_array([4.3, 9.8, 7.5, np.nan, 8.2, np.nan], dtype) * data_unit
        )
        x = np.arange(len(array1)) * dims_unit
        ds = xr.Dataset(
            data_vars={"a": ("x", array1), "b": ("x", array2)},
            coords={"x": x},
        )
        units = extract_units(ds)

        other_array1 = np.ones_like(array1) * other_data_unit
        other_array2 = np.full_like(array2, fill_value=-1) * other_data_unit
        other_x = (np.arange(array1.shape[0]) + 5) * other_dims_unit
        other = xr.Dataset(
            data_vars={"a": ("x", other_array1), "b": ("x", other_array2)},
            coords={"x": other_x},
        )

        if error is not None:
            with pytest.raises(error):
                ds.combine_first(other)

            return

        expected = attach_units(
            strip_units(ds).combine_first(strip_units(convert_units(other, units))),
            units,
        )
        actual = ds.combine_first(other)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    @pytest.mark.parametrize(
        "func",
        (
            method("equals"),
            pytest.param(
                method("identical"),
                marks=pytest.mark.skip("behaviour of identical is unclear"),
            ),
        ),
        ids=repr,
    )
    def test_comparisons(self, func, variant, unit, dtype):
        array1 = np.linspace(0, 5, 10).astype(dtype)
        array2 = np.linspace(-5, 0, 10).astype(dtype)

        coord = np.arange(len(array1)).astype(dtype)

        variants = {
            "data": (unit_registry.m, 1, 1),
            "dims": (1, unit_registry.m, 1),
            "coords": (1, 1, unit_registry.m),
        }
        data_unit, dim_unit, coord_unit = variants[variant]
        a = array1 * data_unit
        b = array2 * data_unit
        x = coord * dim_unit
        y = coord * coord_unit

        ds = xr.Dataset(
            data_vars={"a": ("x", a), "b": ("x", b)},
            coords={"x": x, "y": ("x", y)},
        )
        units = extract_units(ds)

        other_variants = {
            "data": (unit, 1, 1),
            "dims": (1, unit, 1),
            "coords": (1, 1, unit),
        }
        other_data_unit, other_dim_unit, other_coord_unit = other_variants[variant]

        other_units = {
            "a": other_data_unit,
            "b": other_data_unit,
            "x": other_dim_unit,
            "y": other_coord_unit,
        }

        to_convert = {
            key: unit if is_compatible(unit, reference) else None
            for key, (unit, reference) in zip_mappings(units, other_units)
        }
        # convert units where possible, then attach all units to the converted dataset
        other = attach_units(strip_units(convert_units(ds, to_convert)), other_units)
        other_units = extract_units(other)

        # make sure all units are compatible and only then try to
        # convert and compare values
        equal_ds = all(
            is_compatible(unit, other_unit)
            for _, (unit, other_unit) in zip_mappings(units, other_units)
        ) and (strip_units(ds).equals(strip_units(convert_units(other, units))))
        equal_units = units == other_units
        expected = equal_ds and (func.name != "identical" or equal_units)

        actual = func(ds, other)

        assert expected == actual

    # TODO: eventually use another decorator / wrapper function that
    # applies a filter to the parametrize combinations:
    # we only need a single test for data
    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims",
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
        ),
    )
    def test_broadcast_like(self, variant, unit, dtype):
        variants = {
            "data": ((unit_registry.m, unit), (1, 1)),
            "dims": ((1, 1), (unit_registry.m, unit)),
        }
        (data_unit1, data_unit2), (dim_unit1, dim_unit2) = variants[variant]

        array1 = np.linspace(1, 2, 2 * 1).reshape(2, 1).astype(dtype) * data_unit1
        array2 = np.linspace(0, 1, 2 * 3).reshape(2, 3).astype(dtype) * data_unit2

        x1 = np.arange(2) * dim_unit1
        x2 = np.arange(2) * dim_unit2
        y1 = np.array([0]) * dim_unit1
        y2 = np.arange(3) * dim_unit2

        ds1 = xr.Dataset(
            data_vars={"a": (("x", "y"), array1)}, coords={"x": x1, "y": y1}
        )
        ds2 = xr.Dataset(
            data_vars={"a": (("x", "y"), array2)}, coords={"x": x2, "y": y2}
        )

        expected = attach_units(
            strip_units(ds1).broadcast_like(strip_units(ds2)), extract_units(ds1)
        )
        actual = ds1.broadcast_like(ds2)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "unit",
        (
            pytest.param(1, id="no_unit"),
            pytest.param(unit_registry.dimensionless, id="dimensionless"),
            pytest.param(unit_registry.s, id="incompatible_unit"),
            pytest.param(unit_registry.cm, id="compatible_unit"),
            pytest.param(unit_registry.m, id="identical_unit"),
        ),
    )
    def test_broadcast_equals(self, unit, dtype):
        # TODO: does this use indexes?
        left_array1 = np.ones(shape=(2, 3), dtype=dtype) * unit_registry.m
        left_array2 = np.zeros(shape=(3, 6), dtype=dtype) * unit_registry.m

        right_array1 = np.ones(shape=(2,)) * unit
        right_array2 = np.zeros(shape=(3,)) * unit

        left = xr.Dataset(
            {"a": (("x", "y"), left_array1), "b": (("y", "z"), left_array2)},
        )
        right = xr.Dataset({"a": ("x", right_array1), "b": ("y", right_array2)})

        units = merge_mappings(
            extract_units(left),
            {} if is_compatible(left_array1, unit) else {"a": None, "b": None},
        )
        expected = is_compatible(left_array1, unit) and strip_units(
            left
        ).broadcast_equals(strip_units(convert_units(right, units)))
        actual = left.broadcast_equals(right)

        assert expected == actual

    def test_pad(self, dtype):
        a = np.linspace(0, 5, 10).astype(dtype) * unit_registry.Pa
        b = np.linspace(-5, 0, 10).astype(dtype) * unit_registry.degK

        ds = xr.Dataset({"a": ("x", a), "b": ("x", b)})
        units = extract_units(ds)

        expected = attach_units(strip_units(ds).pad(x=(2, 3)), units)
        actual = ds.pad(x=(2, 3))

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (method("unstack"), method("reset_index", "v"), method("reorder_levels")),
        ids=repr,
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims",
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
        ),
    )
    def test_stacking_stacked(self, variant, func, dtype):
        variants = {
            "data": (unit_registry.m, 1),
            "dims": (1, unit_registry.m),
        }
        data_unit, dim_unit = variants[variant]

        array1 = np.linspace(0, 10, 5 * 10).reshape(5, 10).astype(dtype) * data_unit
        array2 = (
            np.linspace(-10, 0, 5 * 10 * 15).reshape(5, 10, 15).astype(dtype)
            * data_unit
        )

        x = np.arange(array1.shape[0]) * dim_unit
        y = np.arange(array1.shape[1]) * dim_unit
        z = np.arange(array2.shape[2]) * dim_unit

        ds = xr.Dataset(
            data_vars={"a": (("x", "y"), array1), "b": (("x", "y", "z"), array2)},
            coords={"x": x, "y": y, "z": z},
        )
        units = extract_units(ds)

        stacked = ds.stack(v=("x", "y"))

        expected = attach_units(func(strip_units(stacked)), units)
        actual = func(stacked)

        assert_units_equal(expected, actual)
        if func.name == "reset_index":
            assert_equal(expected, actual, check_default_indexes=False)
        else:
            assert_equal(expected, actual)

    @pytest.mark.xfail(
        reason="stacked dimension's labels have to be hashable, but is a numpy.array"
    )
    def test_to_stacked_array(self, dtype):
        labels = range(5) * unit_registry.s
        arrays = {
            name: np.linspace(0, 1, 10).astype(dtype) * unit_registry.m
            for name in labels
        }

        ds = xr.Dataset({name: ("x", array) for name, array in arrays.items()})
        units = {None: unit_registry.m, "y": unit_registry.s}

        func = method("to_stacked_array", "z", variable_dim="y", sample_dims=["x"])

        actual = func(ds).rename(None)
        expected = attach_units(
            func(strip_units(ds)).rename(None),
            units,
        )

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            method("transpose", "y", "x", "z1", "z2"),
            method("stack", u=("x", "y")),
            method("set_index", x="x2"),
            method("shift", x=2),
            pytest.param(
                method("rank", dim="x"),
                marks=pytest.mark.skip(reason="rank not implemented for non-ndarray"),
            ),
            method("roll", x=2, roll_coords=False),
            method("sortby", "x2"),
        ),
        ids=repr,
    )
    def test_stacking_reordering(self, func, dtype):
        array1 = (
            np.linspace(0, 10, 2 * 5 * 10).reshape(2, 5, 10).astype(dtype)
            * unit_registry.Pa
        )
        array2 = (
            np.linspace(0, 10, 2 * 5 * 15).reshape(2, 5, 15).astype(dtype)
            * unit_registry.degK
        )

        x = np.arange(array1.shape[0])
        y = np.arange(array1.shape[1])
        z1 = np.arange(array1.shape[2])
        z2 = np.arange(array2.shape[2])

        x2 = np.linspace(0, 1, array1.shape[0])[::-1]

        ds = xr.Dataset(
            data_vars={
                "a": (("x", "y", "z1"), array1),
                "b": (("x", "y", "z2"), array2),
            },
            coords={"x": x, "y": y, "z1": z1, "z2": z2, "x2": ("x", x2)},
        )
        units = extract_units(ds)

        expected = attach_units(func(strip_units(ds)), units)
        actual = func(ds)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "indices",
        (
            pytest.param(4, id="single index"),
            pytest.param([5, 2, 9, 1], id="multiple indices"),
        ),
    )
    def test_isel(self, indices, dtype):
        array1 = np.arange(10).astype(dtype) * unit_registry.s
        array2 = np.linspace(0, 1, 10).astype(dtype) * unit_registry.Pa

        ds = xr.Dataset(data_vars={"a": ("x", array1), "b": ("x", array2)})
        units = extract_units(ds)

        expected = attach_units(strip_units(ds).isel(x=indices), units)
        actual = ds.isel(x=indices)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "raw_values",
        (
            pytest.param(10, id="single_value"),
            pytest.param([10, 5, 13], id="list_of_values"),
            pytest.param(np.array([9, 3, 7, 12]), id="array_of_values"),
        ),
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, KeyError, id="no_units"),
            pytest.param(unit_registry.dimensionless, KeyError, id="dimensionless"),
            pytest.param(unit_registry.degree, KeyError, id="incompatible_unit"),
            pytest.param(unit_registry.mm, KeyError, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_sel(self, raw_values, unit, error, dtype):
        array1 = np.linspace(5, 10, 20).astype(dtype) * unit_registry.degK
        array2 = np.linspace(0, 5, 20).astype(dtype) * unit_registry.Pa
        x = np.arange(len(array1)) * unit_registry.m

        ds = xr.Dataset(
            data_vars={
                "a": xr.DataArray(data=array1, dims="x"),
                "b": xr.DataArray(data=array2, dims="x"),
            },
            coords={"x": x},
        )

        values = raw_values * unit

        # TODO: if we choose dm as compatible unit, single value keys
        # can be found. Should we check that?
        if error is not None:
            with pytest.raises(error):
                ds.sel(x=values)

            return

        expected = attach_units(
            strip_units(ds).sel(
                x=strip_units(convert_units(values, {None: unit_registry.m}))
            ),
            extract_units(ds),
        )
        actual = ds.sel(x=values)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "raw_values",
        (
            pytest.param(10, id="single_value"),
            pytest.param([10, 5, 13], id="list_of_values"),
            pytest.param(np.array([9, 3, 7, 12]), id="array_of_values"),
        ),
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, KeyError, id="no_units"),
            pytest.param(unit_registry.dimensionless, KeyError, id="dimensionless"),
            pytest.param(unit_registry.degree, KeyError, id="incompatible_unit"),
            pytest.param(unit_registry.mm, KeyError, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_drop_sel(self, raw_values, unit, error, dtype):
        array1 = np.linspace(5, 10, 20).astype(dtype) * unit_registry.degK
        array2 = np.linspace(0, 5, 20).astype(dtype) * unit_registry.Pa
        x = np.arange(len(array1)) * unit_registry.m

        ds = xr.Dataset(
            data_vars={
                "a": xr.DataArray(data=array1, dims="x"),
                "b": xr.DataArray(data=array2, dims="x"),
            },
            coords={"x": x},
        )

        values = raw_values * unit

        # TODO: if we choose dm as compatible unit, single value keys
        # can be found. Should we check that?
        if error is not None:
            with pytest.raises(error):
                ds.drop_sel(x=values)

            return

        expected = attach_units(
            strip_units(ds).drop_sel(
                x=strip_units(convert_units(values, {None: unit_registry.m}))
            ),
            extract_units(ds),
        )
        actual = ds.drop_sel(x=values)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "raw_values",
        (
            pytest.param(10, id="single_value"),
            pytest.param([10, 5, 13], id="list_of_values"),
            pytest.param(np.array([9, 3, 7, 12]), id="array_of_values"),
        ),
    )
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, KeyError, id="no_units"),
            pytest.param(unit_registry.dimensionless, KeyError, id="dimensionless"),
            pytest.param(unit_registry.degree, KeyError, id="incompatible_unit"),
            pytest.param(unit_registry.mm, KeyError, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    def test_loc(self, raw_values, unit, error, dtype):
        array1 = np.linspace(5, 10, 20).astype(dtype) * unit_registry.degK
        array2 = np.linspace(0, 5, 20).astype(dtype) * unit_registry.Pa
        x = np.arange(len(array1)) * unit_registry.m

        ds = xr.Dataset(
            data_vars={
                "a": xr.DataArray(data=array1, dims="x"),
                "b": xr.DataArray(data=array2, dims="x"),
            },
            coords={"x": x},
        )

        values = raw_values * unit

        # TODO: if we choose dm as compatible unit, single value keys
        # can be found. Should we check that?
        if error is not None:
            with pytest.raises(error):
                ds.loc[{"x": values}]

            return

        expected = attach_units(
            strip_units(ds).loc[
                {"x": strip_units(convert_units(values, {None: unit_registry.m}))}
            ],
            extract_units(ds),
        )
        actual = ds.loc[{"x": values}]

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            method("head", x=7, y=3, z=6),
            method("tail", x=7, y=3, z=6),
            method("thin", x=7, y=3, z=6),
        ),
        ids=repr,
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    def test_head_tail_thin(self, func, variant, dtype):
        variants = {
            "data": ((unit_registry.degK, unit_registry.Pa), 1, 1),
            "dims": ((1, 1), unit_registry.m, 1),
            "coords": ((1, 1), 1, unit_registry.m),
        }
        (unit_a, unit_b), dim_unit, coord_unit = variants[variant]

        array1 = np.linspace(1, 2, 10 * 5).reshape(10, 5) * unit_a
        array2 = np.linspace(1, 2, 10 * 8).reshape(10, 8) * unit_b

        coords = {
            "x": np.arange(10) * dim_unit,
            "y": np.arange(5) * dim_unit,
            "z": np.arange(8) * dim_unit,
            "u": ("x", np.linspace(0, 1, 10) * coord_unit),
            "v": ("y", np.linspace(1, 2, 5) * coord_unit),
            "w": ("z", np.linspace(-1, 0, 8) * coord_unit),
        }

        ds = xr.Dataset(
            data_vars={
                "a": xr.DataArray(data=array1, dims=("x", "y")),
                "b": xr.DataArray(data=array2, dims=("x", "z")),
            },
            coords=coords,
        )

        expected = attach_units(func(strip_units(ds)), extract_units(ds))
        actual = func(ds)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize("dim", ("x", "y", "z", "t", "all"))
    @pytest.mark.parametrize(
        "shape",
        (
            pytest.param((10, 20), id="nothing squeezable"),
            pytest.param((10, 20, 1), id="last dimension squeezable"),
            pytest.param((10, 1, 20), id="middle dimension squeezable"),
            pytest.param((1, 10, 20), id="first dimension squeezable"),
            pytest.param((1, 10, 1, 20), id="first and last dimension squeezable"),
        ),
    )
    def test_squeeze(self, shape, dim, dtype):
        names = "xyzt"
        dim_lengths = dict(zip(names, shape, strict=False))
        array1 = (
            np.linspace(0, 1, 10 * 20).astype(dtype).reshape(shape) * unit_registry.degK
        )
        array2 = (
            np.linspace(1, 2, 10 * 20).astype(dtype).reshape(shape) * unit_registry.Pa
        )

        ds = xr.Dataset(
            data_vars={
                "a": (tuple(names[: len(shape)]), array1),
                "b": (tuple(names[: len(shape)]), array2),
            },
        )
        units = extract_units(ds)

        kwargs = {"dim": dim} if dim != "all" and dim_lengths.get(dim, 0) == 1 else {}

        expected = attach_units(strip_units(ds).squeeze(**kwargs), units)

        actual = ds.squeeze(**kwargs)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize("variant", ("data", "coords"))
    @pytest.mark.parametrize(
        "func",
        (
            pytest.param(
                method("interp"), marks=pytest.mark.xfail(reason="uses scipy")
            ),
            method("reindex"),
        ),
        ids=repr,
    )
    def test_interp_reindex(self, func, variant, dtype):
        variants = {
            "data": (unit_registry.m, 1),
            "coords": (1, unit_registry.m),
        }
        data_unit, coord_unit = variants[variant]

        array1 = np.linspace(-1, 0, 10).astype(dtype) * data_unit
        array2 = np.linspace(0, 1, 10).astype(dtype) * data_unit

        y = np.arange(10) * coord_unit

        x = np.arange(10)
        new_x = np.arange(8) + 0.5

        ds = xr.Dataset(
            {"a": ("x", array1), "b": ("x", array2)}, coords={"x": x, "y": ("x", y)}
        )
        units = extract_units(ds)

        expected = attach_units(func(strip_units(ds), x=new_x), units)
        actual = func(ds, x=new_x)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize("func", (method("interp"), method("reindex")), ids=repr)
    def test_interp_reindex_indexing(self, func, unit, error, dtype):
        array1 = np.linspace(-1, 0, 10).astype(dtype)
        array2 = np.linspace(0, 1, 10).astype(dtype)

        x = np.arange(10) * unit_registry.m
        new_x = (np.arange(8) + 0.5) * unit

        ds = xr.Dataset({"a": ("x", array1), "b": ("x", array2)}, coords={"x": x})
        units = extract_units(ds)

        if error is not None:
            with pytest.raises(error):
                func(ds, x=new_x)

            return

        expected = attach_units(func(strip_units(ds), x=new_x), units)
        actual = func(ds, x=new_x)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize("variant", ("data", "coords"))
    @pytest.mark.parametrize(
        "func",
        (
            pytest.param(
                method("interp_like"), marks=pytest.mark.xfail(reason="uses scipy")
            ),
            method("reindex_like"),
        ),
        ids=repr,
    )
    def test_interp_reindex_like(self, func, variant, dtype):
        variants = {
            "data": (unit_registry.m, 1),
            "coords": (1, unit_registry.m),
        }
        data_unit, coord_unit = variants[variant]

        array1 = np.linspace(-1, 0, 10).astype(dtype) * data_unit
        array2 = np.linspace(0, 1, 10).astype(dtype) * data_unit

        y = np.arange(10) * coord_unit

        x = np.arange(10)
        new_x = np.arange(8) + 0.5

        ds = xr.Dataset(
            {"a": ("x", array1), "b": ("x", array2)}, coords={"x": x, "y": ("x", y)}
        )
        units = extract_units(ds)

        other = xr.Dataset({"a": ("x", np.empty_like(new_x))}, coords={"x": new_x})

        expected = attach_units(func(strip_units(ds), other), units)
        actual = func(ds, other)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.skip(reason="indexes don't support units")
    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, DimensionalityError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, DimensionalityError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, DimensionalityError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, None, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "func", (method("interp_like"), method("reindex_like")), ids=repr
    )
    def test_interp_reindex_like_indexing(self, func, unit, error, dtype):
        array1 = np.linspace(-1, 0, 10).astype(dtype)
        array2 = np.linspace(0, 1, 10).astype(dtype)

        x = np.arange(10) * unit_registry.m
        new_x = (np.arange(8) + 0.5) * unit

        ds = xr.Dataset({"a": ("x", array1), "b": ("x", array2)}, coords={"x": x})
        units = extract_units(ds)

        other = xr.Dataset({"a": ("x", np.empty_like(new_x))}, coords={"x": new_x})

        if error is not None:
            with pytest.raises(error):
                func(ds, other)

            return

        expected = attach_units(func(strip_units(ds), other), units)
        actual = func(ds, other)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize("compute_backend", ["numbagg", None], indirect=True)
    @pytest.mark.parametrize(
        "func",
        (
            method("diff", dim="x"),
            method("differentiate", coord="x"),
            method("integrate", coord="x"),
            method("quantile", q=[0.25, 0.75]),
            method("reduce", func=np.sum, dim="x"),
            method("map", np.fabs),
        ),
        ids=repr,
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    def test_computation(self, func, variant, dtype, compute_backend):
        variants = {
            "data": ((unit_registry.degK, unit_registry.Pa), 1, 1),
            "dims": ((1, 1), unit_registry.m, 1),
            "coords": ((1, 1), 1, unit_registry.m),
        }
        (unit1, unit2), dim_unit, coord_unit = variants[variant]

        array1 = np.linspace(-5, 5, 4 * 5).reshape(4, 5).astype(dtype) * unit1
        array2 = np.linspace(10, 20, 4 * 3).reshape(4, 3).astype(dtype) * unit2
        x = np.arange(4) * dim_unit
        y = np.arange(5) * dim_unit
        z = np.arange(3) * dim_unit

        ds = xr.Dataset(
            data_vars={
                "a": xr.DataArray(data=array1, dims=("x", "y")),
                "b": xr.DataArray(data=array2, dims=("x", "z")),
            },
            coords={"x": x, "y": y, "z": z, "y2": ("y", np.arange(5) * coord_unit)},
        )

        units = extract_units(ds)

        expected = attach_units(func(strip_units(ds)), units)
        actual = func(ds)

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            method("groupby", "x"),
            method("groupby_bins", "x", bins=2),
            method("coarsen", x=2),
            pytest.param(
                method("rolling", x=3), marks=pytest.mark.xfail(reason="strips units")
            ),
            pytest.param(
                method("rolling_exp", x=3),
                marks=pytest.mark.xfail(
                    reason="numbagg functions are not supported by pint"
                ),
            ),
            method("weighted", xr.DataArray(data=np.linspace(0, 1, 5), dims="y")),
        ),
        ids=repr,
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    def test_computation_objects(self, func, variant, dtype):
        variants = {
            "data": ((unit_registry.degK, unit_registry.Pa), 1, 1),
            "dims": ((1, 1), unit_registry.m, 1),
            "coords": ((1, 1), 1, unit_registry.m),
        }
        (unit1, unit2), dim_unit, coord_unit = variants[variant]

        array1 = np.linspace(-5, 5, 4 * 5).reshape(4, 5).astype(dtype) * unit1
        array2 = np.linspace(10, 20, 4 * 3).reshape(4, 3).astype(dtype) * unit2
        x = np.arange(4) * dim_unit
        y = np.arange(5) * dim_unit
        z = np.arange(3) * dim_unit

        ds = xr.Dataset(
            data_vars={"a": (("x", "y"), array1), "b": (("x", "z"), array2)},
            coords={"x": x, "y": y, "z": z, "y2": ("y", np.arange(5) * coord_unit)},
        )
        units = extract_units(ds)

        args = [] if func.name != "groupby" else ["y"]
        # Doesn't work with flox because pint doesn't implement
        # ufunc.reduceat or np.bincount
        #  kwargs = {"engine": "numpy"} if "groupby" in func.name else {}
        kwargs: dict[str, Any] = {}
        expected = attach_units(func(strip_units(ds)).mean(*args, **kwargs), units)
        actual = func(ds).mean(*args, **kwargs)

        assert_units_equal(expected, actual)
        assert_allclose(expected, actual)

    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    def test_resample(self, variant, dtype):
        # TODO: move this to test_computation_objects
        variants = {
            "data": ((unit_registry.degK, unit_registry.Pa), 1, 1),
            "dims": ((1, 1), unit_registry.m, 1),
            "coords": ((1, 1), 1, unit_registry.m),
        }
        (unit1, unit2), dim_unit, coord_unit = variants[variant]

        array1 = np.linspace(-5, 5, 10 * 5).reshape(10, 5).astype(dtype) * unit1
        array2 = np.linspace(10, 20, 10 * 8).reshape(10, 8).astype(dtype) * unit2

        t = xr.date_range("10-09-2010", periods=array1.shape[0], freq="YE")
        y = np.arange(5) * dim_unit
        z = np.arange(8) * dim_unit

        u = np.linspace(-1, 0, 5) * coord_unit

        ds = xr.Dataset(
            data_vars={"a": (("time", "y"), array1), "b": (("time", "z"), array2)},
            coords={"time": t, "y": y, "z": z, "u": ("y", u)},
        )
        units = extract_units(ds)

        func = method("resample", time="6ME")

        expected = attach_units(func(strip_units(ds)).mean(), units)
        actual = func(ds).mean()

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize("compute_backend", ["numbagg", None], indirect=True)
    @pytest.mark.parametrize(
        "func",
        (
            method("assign", c=lambda ds: 10 * ds.b),
            method("assign_coords", v=("x", np.arange(5) * unit_registry.s)),
            method("first"),
            method("last"),
            method("quantile", q=[0.25, 0.5, 0.75], dim="x"),
        ),
        ids=repr,
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    def test_grouped_operations(self, func, variant, dtype, compute_backend):
        variants = {
            "data": ((unit_registry.degK, unit_registry.Pa), 1, 1),
            "dims": ((1, 1), unit_registry.m, 1),
            "coords": ((1, 1), 1, unit_registry.m),
        }
        (unit1, unit2), dim_unit, coord_unit = variants[variant]

        array1 = np.linspace(-5, 5, 5 * 4).reshape(5, 4).astype(dtype) * unit1
        array2 = np.linspace(10, 20, 5 * 4 * 3).reshape(5, 4, 3).astype(dtype) * unit2
        x = np.arange(5) * dim_unit
        y = np.arange(4) * dim_unit
        z = np.arange(3) * dim_unit

        u = np.linspace(-1, 0, 4) * coord_unit

        ds = xr.Dataset(
            data_vars={"a": (("x", "y"), array1), "b": (("x", "y", "z"), array2)},
            coords={"x": x, "y": y, "z": z, "u": ("y", u)},
        )

        assigned_units = {"c": unit2, "v": unit_registry.s}
        units = merge_mappings(extract_units(ds), assigned_units)

        stripped_kwargs = {
            name: strip_units(value) for name, value in func.kwargs.items()
        }
        expected = attach_units(
            func(strip_units(ds).groupby("y", squeeze=False), **stripped_kwargs), units
        )
        actual = func(ds.groupby("y", squeeze=False))

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "func",
        (
            method("pipe", lambda ds: ds * 10),
            method("assign", d=lambda ds: ds.b * 10),
            method("assign_coords", y2=("y", np.arange(4) * unit_registry.mm)),
            method("assign_attrs", attr1="value"),
            method("rename", x2="x_mm"),
            method("rename_vars", c="temperature"),
            method("rename_dims", x="offset_x"),
            method("swap_dims", {"x": "u"}),
            pytest.param(
                method(
                    "expand_dims", v=np.linspace(10, 20, 12) * unit_registry.s, axis=1
                ),
                marks=pytest.mark.skip(reason="indexes don't support units"),
            ),
            method("drop_vars", "x"),
            method("drop_dims", "z"),
            method("set_coords", names="c"),
            method("reset_coords", names="x2"),
            method("copy"),
        ),
        ids=repr,
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    def test_content_manipulation(self, func, variant, dtype):
        variants = {
            "data": (
                (unit_registry.m**3, unit_registry.Pa, unit_registry.degK),
                1,
                1,
            ),
            "dims": ((1, 1, 1), unit_registry.m, 1),
            "coords": ((1, 1, 1), 1, unit_registry.m),
        }
        (unit1, unit2, unit3), dim_unit, coord_unit = variants[variant]

        array1 = np.linspace(-5, 5, 5 * 4).reshape(5, 4).astype(dtype) * unit1
        array2 = np.linspace(10, 20, 5 * 4 * 3).reshape(5, 4, 3).astype(dtype) * unit2
        array3 = np.linspace(0, 10, 5).astype(dtype) * unit3

        x = np.arange(5) * dim_unit
        y = np.arange(4) * dim_unit
        z = np.arange(3) * dim_unit

        x2 = np.linspace(-1, 0, 5) * coord_unit

        ds = xr.Dataset(
            data_vars={
                "a": (("x", "y"), array1),
                "b": (("x", "y", "z"), array2),
                "c": ("x", array3),
            },
            coords={"x": x, "y": y, "z": z, "x2": ("x", x2)},
        )

        new_units = {
            "y2": unit_registry.mm,
            "x_mm": coord_unit,
            "offset_x": unit_registry.m,
            "d": unit2,
            "temperature": unit3,
        }
        units = merge_mappings(extract_units(ds), new_units)

        stripped_kwargs = {
            key: strip_units(value) for key, value in func.kwargs.items()
        }
        expected = attach_units(func(strip_units(ds), **stripped_kwargs), units)
        actual = func(ds)

        assert_units_equal(expected, actual)
        if func.name == "rename_dims":
            assert_equal(expected, actual, check_default_indexes=False)
        else:
            assert_equal(expected, actual)

    @pytest.mark.parametrize(
        "unit,error",
        (
            pytest.param(1, xr.MergeError, id="no_unit"),
            pytest.param(
                unit_registry.dimensionless, xr.MergeError, id="dimensionless"
            ),
            pytest.param(unit_registry.s, xr.MergeError, id="incompatible_unit"),
            pytest.param(unit_registry.cm, xr.MergeError, id="compatible_unit"),
            pytest.param(unit_registry.m, None, id="identical_unit"),
        ),
    )
    @pytest.mark.parametrize(
        "variant",
        (
            "data",
            pytest.param(
                "dims", marks=pytest.mark.skip(reason="indexes don't support units")
            ),
            "coords",
        ),
    )
    @pytest.mark.filterwarnings(
        "ignore:.*the default value for compat will change:FutureWarning"
    )
    def test_merge(self, variant, unit, error, dtype):
        left_variants = {
            "data": (unit_registry.m, 1, 1),
            "dims": (1, unit_registry.m, 1),
            "coords": (1, 1, unit_registry.m),
        }

        left_data_unit, left_dim_unit, left_coord_unit = left_variants[variant]

        right_variants = {
            "data": (unit, 1, 1),
            "dims": (1, unit, 1),
            "coords": (1, 1, unit),
        }
        right_data_unit, right_dim_unit, right_coord_unit = right_variants[variant]

        left_array = np.arange(10).astype(dtype) * left_data_unit
        right_array = np.arange(-5, 5).astype(dtype) * right_data_unit

        left_dim = np.arange(10, 20) * left_dim_unit
        right_dim = np.arange(5, 15) * right_dim_unit

        left_coord = np.arange(-10, 0) * left_coord_unit
        right_coord = np.arange(-15, -5) * right_coord_unit

        left = xr.Dataset(
            data_vars={"a": ("x", left_array)},
            coords={"x": left_dim, "y": ("x", left_coord)},
        )
        right = xr.Dataset(
            data_vars={"a": ("x", right_array)},
            coords={"x": right_dim, "y": ("x", right_coord)},
        )

        units = extract_units(left)

        if error is not None:
            with pytest.raises(error):
                left.merge(right, compat="no_conflicts", join="outer")

            return

        converted = convert_units(right, units)
        expected = attach_units(
            strip_units(left).merge(strip_units(converted), join="outer"), units
        )
        actual = left.merge(right, join="outer")

        assert_units_equal(expected, actual)
        assert_equal(expected, actual)


@requires_dask
class TestPintWrappingDask:
    def test_duck_array_ops(self):
        import dask.array

        d = dask.array.array([1, 2, 3])
        q = unit_registry.Quantity(d, units="m")
        da = xr.DataArray(q, dims="x")

        actual = da.mean().compute()
        actual.name = None
        expected = xr.DataArray(unit_registry.Quantity(np.array(2.0), units="m"))

        assert_units_equal(expected, actual)
        # Don't use isinstance b/c we don't want to allow subclasses through
        assert type(expected.data) is type(actual.data)


@requires_matplotlib
class TestPlots(PlotTestCase):
    @pytest.mark.parametrize(
        "coord_unit, coord_attrs",
        [
            (1, {"units": "meter"}),
            pytest.param(
                unit_registry.m,
                {},
                marks=pytest.mark.xfail(reason="indexes don't support units"),
            ),
        ],
    )
    def test_units_in_line_plot_labels(self, coord_unit, coord_attrs):
        arr = np.linspace(1, 10, 3) * unit_registry.Pa
        coord_arr = np.linspace(1, 3, 3) * coord_unit
        x_coord = xr.DataArray(coord_arr, dims="x", attrs=coord_attrs)
        da = xr.DataArray(data=arr, dims="x", coords={"x": x_coord}, name="pressure")

        da.plot.line()

        ax = plt.gca()
        assert ax.get_ylabel() == "pressure [pascal]"
        assert ax.get_xlabel() == "x [meter]"

    @pytest.mark.parametrize(
        "coord_unit, coord_attrs",
        [
            (1, {"units": "meter"}),
            pytest.param(
                unit_registry.m,
                {},
                marks=pytest.mark.xfail(reason="indexes don't support units"),
            ),
        ],
    )
    def test_units_in_slice_line_plot_labels_sel(self, coord_unit, coord_attrs):
        arr = xr.DataArray(
            name="var_a",
            data=np.array([[1, 2], [3, 4]]),
            coords=dict(
                a=("a", np.array([5, 6]) * coord_unit, coord_attrs),
                b=("b", np.array([7, 8]) * coord_unit, coord_attrs),
            ),
            dims=("a", "b"),
        )
        arr.sel(a=5).plot(marker="o")  # type: ignore[call-arg]

        assert plt.gca().get_title() == "a = 5 [meter]"

    @pytest.mark.parametrize(
        "coord_unit, coord_attrs",
        [
            (1, {"units": "meter"}),
            pytest.param(
                unit_registry.m,
                {},
                marks=pytest.mark.xfail(reason="pint.errors.UnitStrippedWarning"),
            ),
        ],
    )
    def test_units_in_slice_line_plot_labels_isel(self, coord_unit, coord_attrs):
        arr = xr.DataArray(
            name="var_a",
            data=np.array([[1, 2], [3, 4]]),
            coords=dict(
                a=("x", np.array([5, 6]) * coord_unit, coord_attrs),
                b=("y", np.array([7, 8])),
            ),
            dims=("x", "y"),
        )
        arr.isel(x=0).plot(marker="o")  # type: ignore[call-arg]
        assert plt.gca().get_title() == "a = 5 [meter]"

    def test_units_in_2d_plot_colorbar_label(self):
        arr = np.ones((2, 3)) * unit_registry.Pa
        da = xr.DataArray(data=arr, dims=["x", "y"], name="pressure")

        _fig, (ax, cax) = plt.subplots(1, 2)
        ax = da.plot.contourf(ax=ax, cbar_ax=cax, add_colorbar=True)

        assert cax.get_ylabel() == "pressure [pascal]"

    def test_units_facetgrid_plot_labels(self):
        arr = np.ones((2, 3)) * unit_registry.Pa
        da = xr.DataArray(data=arr, dims=["x", "y"], name="pressure")

        _fig, (_ax, _cax) = plt.subplots(1, 2)
        fgrid = da.plot.line(x="x", col="y")

        assert fgrid.axs[0, 0].get_ylabel() == "pressure [pascal]"

    def test_units_facetgrid_2d_imshow_plot_colorbar_labels(self):
        arr = np.ones((2, 3, 4, 5)) * unit_registry.Pa
        da = xr.DataArray(data=arr, dims=["x", "y", "z", "w"], name="pressure")

        da.plot.imshow(x="x", y="y", col="w")  # no colorbar to check labels of

    def test_units_facetgrid_2d_contourf_plot_colorbar_labels(self):
        arr = np.ones((2, 3, 4)) * unit_registry.Pa
        da = xr.DataArray(data=arr, dims=["x", "y", "z"], name="pressure")

        _fig, (_ax1, _ax2, _ax3, _cax) = plt.subplots(1, 4)
        fgrid = da.plot.contourf(x="x", y="y", col="z")

        assert fgrid.cbar.ax.get_ylabel() == "pressure [pascal]"  # type: ignore[union-attr]
