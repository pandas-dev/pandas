"""Testing functions exposed to the user API"""

import functools
import warnings
from collections.abc import Hashable
from typing import Any

import numpy as np
import pandas as pd

from xarray.core import duck_array_ops, formatting, utils
from xarray.core.coordinates import Coordinates
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.datatree import DataTree
from xarray.core.datatree_mapping import map_over_datasets
from xarray.core.formatting import diff_datatree_repr
from xarray.core.indexes import Index, PandasIndex, PandasMultiIndex, default_indexes
from xarray.core.variable import IndexVariable, Variable


def ensure_warnings(func):
    # sometimes tests elevate warnings to errors
    # -> make sure that does not happen in the assert_* functions
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        __tracebackhide__ = True

        with warnings.catch_warnings():
            # only remove filters that would "error"
            warnings.filters = [f for f in warnings.filters if f[0] != "error"]

            return func(*args, **kwargs)

    return wrapper


def _decode_string_data(data):
    if data.dtype.kind == "S":
        return np.char.decode(data, "utf-8", "replace")
    return data


def _data_allclose_or_equiv(arr1, arr2, rtol=1e-05, atol=1e-08, decode_bytes=True):
    if any(arr.dtype.kind == "S" for arr in [arr1, arr2]) and decode_bytes:
        arr1 = _decode_string_data(arr1)
        arr2 = _decode_string_data(arr2)
    exact_dtypes = ["M", "m", "O", "S", "U"]
    if any(arr.dtype.kind in exact_dtypes for arr in [arr1, arr2]):
        return duck_array_ops.array_equiv(arr1, arr2)
    else:
        return duck_array_ops.allclose_or_equiv(arr1, arr2, rtol=rtol, atol=atol)


@ensure_warnings
def assert_isomorphic(a: DataTree, b: DataTree):
    """
    Two DataTrees are considered isomorphic if the set of paths to their
    descendent nodes are the same.

    Nothing about the data or attrs in each node is checked.

    Isomorphism is a necessary condition for two trees to be used in a nodewise binary operation,
    such as tree1 + tree2.

    Parameters
    ----------
    a : DataTree
        The first object to compare.
    b : DataTree
        The second object to compare.

    See Also
    --------
    DataTree.isomorphic
    assert_equal
    assert_identical
    """
    __tracebackhide__ = True
    assert isinstance(a, type(b))

    if isinstance(a, DataTree):
        assert a.isomorphic(b), diff_datatree_repr(a, b, "isomorphic")
    else:
        raise TypeError(f"{type(a)} not of type DataTree")


def maybe_transpose_dims(a, b, check_dim_order: bool):
    """Helper for assert_equal/allclose/identical"""

    __tracebackhide__ = True

    def _maybe_transpose_dims(a, b):
        if not isinstance(a, Variable | DataArray | Dataset):
            return b
        if set(a.dims) == set(b.dims):
            # Ensure transpose won't fail if a dimension is missing
            # If this is the case, the difference will be caught by the caller
            return b.transpose(*a.dims)
        return b

    if check_dim_order:
        return b

    if isinstance(a, DataTree):
        return map_over_datasets(_maybe_transpose_dims, a, b)

    return _maybe_transpose_dims(a, b)


@ensure_warnings
def assert_equal(a, b, check_dim_order: bool = True):
    """Like :py:func:`numpy.testing.assert_array_equal`, but for xarray
    objects.

    Raises an AssertionError if two objects are not equal. This will match
    data values, dimensions and coordinates, but not names or attributes
    (except for Dataset objects for which the variable names must match).
    Arrays with NaN in the same location are considered equal.

    For DataTree objects, assert_equal is mapped over all Datasets on each node,
    with the DataTrees being equal if both are isomorphic and the corresponding
    Datasets at each node are themselves equal.

    Parameters
    ----------
    a : xarray.Dataset, xarray.DataArray, xarray.Variable, xarray.Coordinates
        or xarray.core.datatree.DataTree. The first object to compare.
    b : xarray.Dataset, xarray.DataArray, xarray.Variable, xarray.Coordinates
        or xarray.core.datatree.DataTree. The second object to compare.
    check_dim_order : bool, optional, default is True
        Whether dimensions must be in the same order.

    See Also
    --------
    assert_identical, assert_allclose, Dataset.equals, DataArray.equals
    numpy.testing.assert_array_equal
    """
    __tracebackhide__ = True
    assert type(a) is type(b) or (
        isinstance(a, Coordinates) and isinstance(b, Coordinates)
    )
    b = maybe_transpose_dims(a, b, check_dim_order)
    if isinstance(a, Variable | DataArray):
        assert a.equals(b), formatting.diff_array_repr(a, b, "equals")
    elif isinstance(a, Dataset):
        assert a.equals(b), formatting.diff_dataset_repr(a, b, "equals")
    elif isinstance(a, Coordinates):
        assert a.equals(b), formatting.diff_coords_repr(a, b, "equals")
    elif isinstance(a, DataTree):
        assert a.equals(b), diff_datatree_repr(a, b, "equals")
    else:
        raise TypeError(f"{type(a)} not supported by assertion comparison")


@ensure_warnings
def assert_identical(a, b):
    """Like :py:func:`xarray.testing.assert_equal`, but also matches the
    objects' names and attributes.

    Raises an AssertionError if two objects are not identical.

    For DataTree objects, assert_identical is mapped over all Datasets on each
    node, with the DataTrees being identical if both are isomorphic and the
    corresponding Datasets at each node are themselves identical.

    Parameters
    ----------
    a : xarray.Dataset, xarray.DataArray, xarray.Variable or xarray.Coordinates
        The first object to compare.
    b : xarray.Dataset, xarray.DataArray, xarray.Variable or xarray.Coordinates
        The second object to compare.

    See Also
    --------
    assert_equal, assert_allclose, Dataset.equals, DataArray.equals
    """
    __tracebackhide__ = True
    assert type(a) is type(b) or (
        isinstance(a, Coordinates) and isinstance(b, Coordinates)
    )
    if isinstance(a, Variable):
        assert a.identical(b), formatting.diff_array_repr(a, b, "identical")
    elif isinstance(a, DataArray):
        assert a.name == b.name, (
            f"DataArray names are different. L: {a.name}, R: {b.name}"
        )
        assert a.identical(b), formatting.diff_array_repr(a, b, "identical")
    elif isinstance(a, Dataset | Variable):
        assert a.identical(b), formatting.diff_dataset_repr(a, b, "identical")
    elif isinstance(a, Coordinates):
        assert a.identical(b), formatting.diff_coords_repr(a, b, "identical")
    elif isinstance(a, DataTree):
        assert a.identical(b), diff_datatree_repr(a, b, "identical")
    else:
        raise TypeError(f"{type(a)} not supported by assertion comparison")


@ensure_warnings
def assert_allclose(
    a, b, rtol=1e-05, atol=1e-08, decode_bytes=True, check_dim_order: bool = True
):
    """Like :py:func:`numpy.testing.assert_allclose`, but for xarray objects.

    Raises an AssertionError if two objects are not equal up to desired
    tolerance.

    Parameters
    ----------
    a : xarray.Dataset, xarray.DataArray or xarray.Variable
        The first object to compare.
    b : xarray.Dataset, xarray.DataArray or xarray.Variable
        The second object to compare.
    rtol : float, optional
        Relative tolerance.
    atol : float, optional
        Absolute tolerance.
    decode_bytes : bool, optional
        Whether byte dtypes should be decoded to strings as UTF-8 or not.
        This is useful for testing serialization methods on Python 3 that
        return saved strings as bytes.
    check_dim_order : bool, optional, default is True
        Whether dimensions must be in the same order.

    See Also
    --------
    assert_identical, assert_equal, numpy.testing.assert_allclose
    """
    __tracebackhide__ = True
    assert type(a) is type(b)
    b = maybe_transpose_dims(a, b, check_dim_order)

    equiv = functools.partial(
        _data_allclose_or_equiv, rtol=rtol, atol=atol, decode_bytes=decode_bytes
    )
    equiv.__name__ = "allclose"  # type: ignore[attr-defined]

    def compat_variable(a, b):
        a = getattr(a, "variable", a)
        b = getattr(b, "variable", b)
        return a.dims == b.dims and (a._data is b._data or equiv(a.data, b.data))

    def compat_node(a, b):
        return a.ds._coord_names == b.ds._coord_names and utils.dict_equiv(
            a.variables, b.variables, compat=compat_variable
        )

    if isinstance(a, Variable):
        allclose = compat_variable(a, b)
        assert allclose, formatting.diff_array_repr(a, b, compat=equiv)
    elif isinstance(a, DataArray):
        allclose = utils.dict_equiv(
            a.coords, b.coords, compat=compat_variable
        ) and compat_variable(a.variable, b.variable)
        assert allclose, formatting.diff_array_repr(a, b, compat=equiv)
    elif isinstance(a, Dataset):
        allclose = a._coord_names == b._coord_names and utils.dict_equiv(
            a.variables, b.variables, compat=compat_variable
        )
        assert allclose, formatting.diff_dataset_repr(a, b, compat=equiv)
    elif isinstance(a, Coordinates):
        allclose = utils.dict_equiv(a.variables, b.variables, compat=compat_variable)
        assert allclose, formatting.diff_coords_repr(a, b, compat=equiv)
    elif isinstance(a, DataTree):
        allclose = utils.dict_equiv(
            dict(a.subtree_with_keys), dict(b.subtree_with_keys), compat=compat_node
        )
        assert allclose, formatting.diff_datatree_repr(a, b, compat=equiv)
    else:
        raise TypeError(f"{type(a)} not supported by assertion comparison")


def _format_message(x, y, err_msg, verbose):
    diff = x - y
    abs_diff = max(abs(diff))
    rel_diff = "not implemented"

    n_diff = np.count_nonzero(diff)
    n_total = diff.size

    fraction = f"{n_diff} / {n_total}"
    percentage = float(n_diff / n_total * 100)

    parts = [
        "Arrays are not equal",
        err_msg,
        f"Mismatched elements: {fraction} ({percentage:.0f}%)",
        f"Max absolute difference: {abs_diff}",
        f"Max relative difference: {rel_diff}",
    ]
    if verbose:
        parts += [
            f" x: {x!r}",
            f" y: {y!r}",
        ]

    return "\n".join(parts)


@ensure_warnings
def assert_duckarray_allclose(
    actual, desired, rtol=1e-07, atol=0, err_msg="", verbose=True
):
    """Like `np.testing.assert_allclose`, but for duckarrays."""
    __tracebackhide__ = True

    allclose = duck_array_ops.allclose_or_equiv(actual, desired, rtol=rtol, atol=atol)
    assert allclose, _format_message(actual, desired, err_msg=err_msg, verbose=verbose)


@ensure_warnings
def assert_duckarray_equal(x, y, err_msg="", verbose=True):
    """Like `np.testing.assert_array_equal`, but for duckarrays"""
    __tracebackhide__ = True

    if not utils.is_duck_array(x) and not utils.is_scalar(x):
        x = np.asarray(x)

    if not utils.is_duck_array(y) and not utils.is_scalar(y):
        y = np.asarray(y)

    if (utils.is_duck_array(x) and utils.is_scalar(y)) or (
        utils.is_scalar(x) and utils.is_duck_array(y)
    ):
        equiv = duck_array_ops.array_all(x == y)
    else:
        equiv = duck_array_ops.array_equiv(x, y)
    assert equiv, _format_message(x, y, err_msg=err_msg, verbose=verbose)


def assert_chunks_equal(a, b):
    """
    Assert that chunksizes along chunked dimensions are equal.

    Parameters
    ----------
    a : xarray.Dataset or xarray.DataArray
        The first object to compare.
    b : xarray.Dataset or xarray.DataArray
        The second object to compare.
    """

    if isinstance(a, DataArray) != isinstance(b, DataArray):
        raise TypeError("a and b have mismatched types")

    left = a.unify_chunks()
    right = b.unify_chunks()
    assert left.chunks == right.chunks


def _assert_indexes_invariants_checks(
    indexes, possible_coord_variables, dims, check_default=True
):
    assert isinstance(indexes, dict), indexes
    assert all(isinstance(v, Index) for v in indexes.values()), {
        k: type(v) for k, v in indexes.items()
    }

    if check_default:
        index_vars = {
            k
            for k, v in possible_coord_variables.items()
            if isinstance(v, IndexVariable)
        }
        assert indexes.keys() <= index_vars, (set(indexes), index_vars)
        assert all(
            k in index_vars
            for k, v in possible_coord_variables.items()
            if v.dims == (k,)
        ), {k: type(v) for k, v in possible_coord_variables.items()}

    assert not any(
        isinstance(v, IndexVariable)
        for k, v in possible_coord_variables.items()
        if k not in indexes.keys()
    ), {k: type(v) for k, v in possible_coord_variables.items()}

    # check pandas index wrappers vs. coordinate data adapters
    for k, index in indexes.items():
        if isinstance(index, PandasIndex):
            pd_index = index.index
            var = possible_coord_variables[k]
            assert (index.dim,) == var.dims, (pd_index, var)
            if k == index.dim:
                # skip multi-index levels here (checked below)
                assert index.coord_dtype == var.dtype, (index.coord_dtype, var.dtype)
            assert isinstance(var._data.array, pd.Index), var._data.array
            # TODO: check identity instead of equality?
            assert pd_index.equals(var._data.array), (pd_index, var)
        if isinstance(index, PandasMultiIndex):
            pd_index = index.index
            for name in index.index.names:
                assert name in possible_coord_variables, (pd_index, index_vars)
                var = possible_coord_variables[name]
                assert (index.dim,) == var.dims, (pd_index, var)
                assert index.level_coords_dtype[name] == var.dtype, (
                    index.level_coords_dtype[name],
                    var.dtype,
                )
                assert isinstance(var._data.array, pd.MultiIndex), var._data.array
                assert pd_index.equals(var._data.array), (pd_index, var)
                # check all all levels are in `indexes`
                assert name in indexes, (name, set(indexes))
                # index identity is used to find unique indexes in `indexes`
                assert index is indexes[name], (pd_index, indexes[name].index)

    if check_default:
        defaults = default_indexes(possible_coord_variables, dims)
        assert indexes.keys() == defaults.keys(), (set(indexes), set(defaults))
        assert all(v.equals(defaults[k]) for k, v in indexes.items()), (
            indexes,
            defaults,
        )


def _assert_variable_invariants(
    var: Variable | Any,
    name: Hashable = None,
) -> None:
    if name is None:
        name_or_empty: tuple = ()
    else:
        name_or_empty = (name,)

    assert isinstance(var, Variable), {name: type(var)}

    assert isinstance(var._dims, tuple), name_or_empty + (var._dims,)
    assert len(var._dims) == len(var._data.shape), name_or_empty + (
        var._dims,
        var._data.shape,
    )
    assert isinstance(var._encoding, type(None) | dict), name_or_empty + (
        var._encoding,
    )
    assert isinstance(var._attrs, type(None) | dict), name_or_empty + (var._attrs,)


def _assert_dataarray_invariants(da: DataArray, check_default_indexes: bool):
    _assert_variable_invariants(da._variable)

    assert isinstance(da._coords, dict), da._coords

    if check_default_indexes:
        assert all(set(v.dims) <= set(da.dims) for v in da._coords.values()), (
            da.dims,
            {k: v.dims for k, v in da._coords.items()},
        )

    for k, v in da._coords.items():
        _assert_variable_invariants(v, k)

    assert da._indexes is not None
    _assert_indexes_invariants_checks(
        da._indexes, da._coords, da.dims, check_default=check_default_indexes
    )


def _assert_dataset_invariants(ds: Dataset, check_default_indexes: bool):
    assert isinstance(ds._variables, dict), type(ds._variables)

    for k, v in ds._variables.items():
        _assert_variable_invariants(v, k)

    assert isinstance(ds._coord_names, set), ds._coord_names
    assert ds._coord_names <= ds._variables.keys(), (
        ds._coord_names,
        set(ds._variables),
    )

    assert type(ds._dims) is dict, ds._dims
    assert all(isinstance(v, int) for v in ds._dims.values()), ds._dims
    var_dims: set[Hashable] = set()
    for v in ds._variables.values():
        var_dims.update(v.dims)
    assert ds._dims.keys() == var_dims, (set(ds._dims), var_dims)
    assert all(
        ds._dims[k] == v.sizes[k] for v in ds._variables.values() for k in v.sizes
    ), (ds._dims, {k: v.sizes for k, v in ds._variables.items()})

    assert ds._indexes is not None
    _assert_indexes_invariants_checks(
        ds._indexes, ds._variables, ds._dims, check_default=check_default_indexes
    )

    assert isinstance(ds._encoding, type(None) | dict)
    assert isinstance(ds._attrs, type(None) | dict)


def _assert_internal_invariants(
    xarray_obj: DataArray | Dataset | Variable, check_default_indexes: bool
):
    """Validate that an xarray object satisfies its own internal invariants.

    This exists for the benefit of xarray's own test suite, but may be useful
    in external projects if they (ill-advisedly) create objects using xarray's
    private APIs.
    """
    if isinstance(xarray_obj, Variable):
        _assert_variable_invariants(xarray_obj)
    elif isinstance(xarray_obj, DataArray):
        _assert_dataarray_invariants(
            xarray_obj, check_default_indexes=check_default_indexes
        )
    elif isinstance(xarray_obj, Dataset):
        _assert_dataset_invariants(
            xarray_obj, check_default_indexes=check_default_indexes
        )
    elif isinstance(xarray_obj, Coordinates):
        _assert_dataset_invariants(
            xarray_obj.to_dataset(), check_default_indexes=check_default_indexes
        )
    else:
        raise TypeError(
            f"{type(xarray_obj)} is not a supported type for xarray invariant checks"
        )
