"""
Functions for applying functions that act on arrays to xarray's labeled data.
"""

from __future__ import annotations

import functools
import itertools
import operator
import warnings
from collections import Counter
from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence, Set
from typing import TYPE_CHECKING, Any, Callable, Literal, TypeVar, Union, cast, overload

import numpy as np

from xarray.core import dtypes, duck_array_ops, utils
from xarray.core.alignment import align, deep_align
from xarray.core.common import zeros_like
from xarray.core.duck_array_ops import datetime_to_numeric
from xarray.core.formatting import limit_lines
from xarray.core.indexes import Index, filter_indexes_from_coords
from xarray.core.merge import merge_attrs, merge_coordinates_without_align
from xarray.core.options import OPTIONS, _get_keep_attrs
from xarray.core.types import Dims, T_DataArray
from xarray.core.utils import is_dict_like, is_duck_dask_array, is_scalar, parse_dims
from xarray.core.variable import Variable
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import is_chunked_array
from xarray.util.deprecation_helpers import deprecate_dims

if TYPE_CHECKING:
    from xarray.core.coordinates import Coordinates
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import CombineAttrsOptions, JoinOptions

    MissingCoreDimOptions = Literal["raise", "copy", "drop"]

_NO_FILL_VALUE = utils.ReprObject("<no-fill-value>")
_DEFAULT_NAME = utils.ReprObject("<default-name>")
_JOINS_WITHOUT_FILL_VALUES = frozenset({"inner", "exact"})


def _first_of_type(args, kind):
    """Return either first object of type 'kind' or raise if not found."""
    for arg in args:
        if isinstance(arg, kind):
            return arg

    raise ValueError("This should be unreachable.")


def _all_of_type(args, kind):
    """Return all objects of type 'kind'"""
    return [arg for arg in args if isinstance(arg, kind)]


class _UFuncSignature:
    """Core dimensions signature for a given function.

    Based on the signature provided by generalized ufuncs in NumPy.

    Attributes
    ----------
    input_core_dims : tuple[tuple]
        Core dimension names on each input variable.
    output_core_dims : tuple[tuple]
        Core dimension names on each output variable.
    """

    __slots__ = (
        "input_core_dims",
        "output_core_dims",
        "_all_input_core_dims",
        "_all_output_core_dims",
        "_all_core_dims",
    )

    def __init__(self, input_core_dims, output_core_dims=((),)):
        self.input_core_dims = tuple(tuple(a) for a in input_core_dims)
        self.output_core_dims = tuple(tuple(a) for a in output_core_dims)
        self._all_input_core_dims = None
        self._all_output_core_dims = None
        self._all_core_dims = None

    @property
    def all_input_core_dims(self):
        if self._all_input_core_dims is None:
            self._all_input_core_dims = frozenset(
                dim for dims in self.input_core_dims for dim in dims
            )
        return self._all_input_core_dims

    @property
    def all_output_core_dims(self):
        if self._all_output_core_dims is None:
            self._all_output_core_dims = frozenset(
                dim for dims in self.output_core_dims for dim in dims
            )
        return self._all_output_core_dims

    @property
    def all_core_dims(self):
        if self._all_core_dims is None:
            self._all_core_dims = self.all_input_core_dims | self.all_output_core_dims
        return self._all_core_dims

    @property
    def dims_map(self):
        return {
            core_dim: f"dim{n}" for n, core_dim in enumerate(sorted(self.all_core_dims))
        }

    @property
    def num_inputs(self):
        return len(self.input_core_dims)

    @property
    def num_outputs(self):
        return len(self.output_core_dims)

    def __eq__(self, other):
        try:
            return (
                self.input_core_dims == other.input_core_dims
                and self.output_core_dims == other.output_core_dims
            )
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return f"{type(self).__name__}({list(self.input_core_dims)!r}, {list(self.output_core_dims)!r})"

    def __str__(self):
        lhs = ",".join("({})".format(",".join(dims)) for dims in self.input_core_dims)
        rhs = ",".join("({})".format(",".join(dims)) for dims in self.output_core_dims)
        return f"{lhs}->{rhs}"

    def to_gufunc_string(self, exclude_dims=frozenset()):
        """Create an equivalent signature string for a NumPy gufunc.

        Unlike __str__, handles dimensions that don't map to Python
        identifiers.

        Also creates unique names for input_core_dims contained in exclude_dims.
        """
        input_core_dims = [
            [self.dims_map[dim] for dim in core_dims]
            for core_dims in self.input_core_dims
        ]
        output_core_dims = [
            [self.dims_map[dim] for dim in core_dims]
            for core_dims in self.output_core_dims
        ]

        # enumerate input_core_dims contained in exclude_dims to make them unique
        if exclude_dims:
            exclude_dims = [self.dims_map[dim] for dim in exclude_dims]

            counter: Counter = Counter()

            def _enumerate(dim):
                if dim in exclude_dims:
                    n = counter[dim]
                    counter.update([dim])
                    dim = f"{dim}_{n}"
                return dim

            input_core_dims = [
                [_enumerate(dim) for dim in arg] for arg in input_core_dims
            ]

        alt_signature = type(self)(input_core_dims, output_core_dims)
        return str(alt_signature)


def result_name(objects: Iterable[Any]) -> Any:
    # use the same naming heuristics as pandas:
    # https://github.com/blaze/blaze/issues/458#issuecomment-51936356
    names = {getattr(obj, "name", _DEFAULT_NAME) for obj in objects}
    names.discard(_DEFAULT_NAME)
    if len(names) == 1:
        (name,) = names
    else:
        name = None
    return name


def _get_coords_list(args: Iterable[Any]) -> list[Coordinates]:
    coords_list = []
    for arg in args:
        try:
            coords = arg.coords
        except AttributeError:
            pass  # skip this argument
        else:
            coords_list.append(coords)
    return coords_list


def build_output_coords_and_indexes(
    args: Iterable[Any],
    signature: _UFuncSignature,
    exclude_dims: Set = frozenset(),
    combine_attrs: CombineAttrsOptions = "override",
) -> tuple[list[dict[Any, Variable]], list[dict[Any, Index]]]:
    """Build output coordinates and indexes for an operation.

    Parameters
    ----------
    args : Iterable
        List of raw operation arguments. Any valid types for xarray operations
        are OK, e.g., scalars, Variable, DataArray, Dataset.
    signature : _UfuncSignature
        Core dimensions signature for the operation.
    exclude_dims : set, optional
        Dimensions excluded from the operation. Coordinates along these
        dimensions are dropped.
    combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts", \
                     "override"} or callable, default: "drop"
        A callable or a string indicating how to combine attrs of the objects being
        merged:

        - "drop": empty attrs on returned Dataset.
        - "identical": all attrs must be the same on every object.
        - "no_conflicts": attrs from all objects are combined, any that have
          the same name must also have the same value.
        - "drop_conflicts": attrs from all objects are combined, any that have
          the same name but different values are dropped.
        - "override": skip comparing and copy attrs from the first dataset to
          the result.

        If a callable, it must expect a sequence of ``attrs`` dicts and a context object
        as its only parameters.

    Returns
    -------
    Dictionaries of Variable and Index objects with merged coordinates.
    """
    coords_list = _get_coords_list(args)

    if len(coords_list) == 1 and not exclude_dims:
        # we can skip the expensive merge
        (unpacked_coords,) = coords_list
        merged_vars = dict(unpacked_coords.variables)
        merged_indexes = dict(unpacked_coords.xindexes)
    else:
        merged_vars, merged_indexes = merge_coordinates_without_align(
            coords_list, exclude_dims=exclude_dims, combine_attrs=combine_attrs
        )

    output_coords = []
    output_indexes = []
    for output_dims in signature.output_core_dims:
        dropped_dims = signature.all_input_core_dims - set(output_dims)
        if dropped_dims:
            filtered_coords = {
                k: v for k, v in merged_vars.items() if dropped_dims.isdisjoint(v.dims)
            }
            filtered_indexes = filter_indexes_from_coords(
                merged_indexes, set(filtered_coords)
            )
        else:
            filtered_coords = merged_vars
            filtered_indexes = merged_indexes
        output_coords.append(filtered_coords)
        output_indexes.append(filtered_indexes)

    return output_coords, output_indexes


def apply_dataarray_vfunc(
    func,
    *args,
    signature: _UFuncSignature,
    join: JoinOptions = "inner",
    exclude_dims=frozenset(),
    keep_attrs="override",
) -> tuple[DataArray, ...] | DataArray:
    """Apply a variable level function over DataArray, Variable and/or ndarray
    objects.
    """
    from xarray.core.dataarray import DataArray

    if len(args) > 1:
        args = tuple(
            deep_align(
                args,
                join=join,
                copy=False,
                exclude=exclude_dims,
                raise_on_invalid=False,
            )
        )

    objs = _all_of_type(args, DataArray)

    if keep_attrs == "drop":
        name = result_name(args)
    else:
        first_obj = _first_of_type(args, DataArray)
        name = first_obj.name
    result_coords, result_indexes = build_output_coords_and_indexes(
        args, signature, exclude_dims, combine_attrs=keep_attrs
    )

    data_vars = [getattr(a, "variable", a) for a in args]
    result_var = func(*data_vars)

    out: tuple[DataArray, ...] | DataArray
    if signature.num_outputs > 1:
        out = tuple(
            DataArray(
                variable, coords=coords, indexes=indexes, name=name, fastpath=True
            )
            for variable, coords, indexes in zip(
                result_var, result_coords, result_indexes
            )
        )
    else:
        (coords,) = result_coords
        (indexes,) = result_indexes
        out = DataArray(
            result_var, coords=coords, indexes=indexes, name=name, fastpath=True
        )

    attrs = merge_attrs([x.attrs for x in objs], combine_attrs=keep_attrs)
    if isinstance(out, tuple):
        for da in out:
            da.attrs = attrs
    else:
        out.attrs = attrs

    return out


def ordered_set_union(all_keys: list[Iterable]) -> Iterable:
    return {key: None for keys in all_keys for key in keys}.keys()


def ordered_set_intersection(all_keys: list[Iterable]) -> Iterable:
    intersection = set(all_keys[0])
    for keys in all_keys[1:]:
        intersection.intersection_update(keys)
    return [key for key in all_keys[0] if key in intersection]


def assert_and_return_exact_match(all_keys):
    first_keys = all_keys[0]
    for keys in all_keys[1:]:
        if keys != first_keys:
            raise ValueError(
                "exact match required for all data variable names, "
                f"but {list(keys)} != {list(first_keys)}: {set(keys) ^ set(first_keys)} are not in both."
            )
    return first_keys


_JOINERS: dict[str, Callable] = {
    "inner": ordered_set_intersection,
    "outer": ordered_set_union,
    "left": operator.itemgetter(0),
    "right": operator.itemgetter(-1),
    "exact": assert_and_return_exact_match,
}


def join_dict_keys(objects: Iterable[Mapping | Any], how: str = "inner") -> Iterable:
    joiner = _JOINERS[how]
    all_keys = [obj.keys() for obj in objects if hasattr(obj, "keys")]
    return joiner(all_keys)


def collect_dict_values(
    objects: Iterable[Mapping | Any], keys: Iterable, fill_value: object = None
) -> list[list]:
    return [
        [obj.get(key, fill_value) if is_dict_like(obj) else obj for obj in objects]
        for key in keys
    ]


def _as_variables_or_variable(arg) -> Variable | tuple[Variable]:
    try:
        return arg.variables
    except AttributeError:
        try:
            return arg.variable
        except AttributeError:
            return arg


def _unpack_dict_tuples(
    result_vars: Mapping[Any, tuple[Variable, ...]], num_outputs: int
) -> tuple[dict[Hashable, Variable], ...]:
    out: tuple[dict[Hashable, Variable], ...] = tuple({} for _ in range(num_outputs))
    for name, values in result_vars.items():
        for value, results_dict in zip(values, out):
            results_dict[name] = value
    return out


def _check_core_dims(signature, variable_args, name):
    """
    Check if an arg has all the core dims required by the signature.

    Slightly awkward design, of returning the error message. But we want to
    give a detailed error message, which requires inspecting the variable in
    the inner loop.
    """
    missing = []
    for i, (core_dims, variable_arg) in enumerate(
        zip(signature.input_core_dims, variable_args)
    ):
        # Check whether all the dims are on the variable. Note that we need the
        # `hasattr` to check for a dims property, to protect against the case where
        # a numpy array is passed in.
        if hasattr(variable_arg, "dims") and set(core_dims) - set(variable_arg.dims):
            missing += [[i, variable_arg, core_dims]]
    if missing:
        message = ""
        for i, variable_arg, core_dims in missing:
            message += f"Missing core dims {set(core_dims) - set(variable_arg.dims)} from arg number {i + 1} on a variable named `{name}`:\n{variable_arg}\n\n"
        message += "Either add the core dimension, or if passing a dataset alternatively pass `on_missing_core_dim` as `copy` or `drop`. "
        return message
    return True


def apply_dict_of_variables_vfunc(
    func,
    *args,
    signature: _UFuncSignature,
    join="inner",
    fill_value=None,
    on_missing_core_dim: MissingCoreDimOptions = "raise",
):
    """Apply a variable level function over dicts of DataArray, DataArray,
    Variable and ndarray objects.
    """
    args = tuple(_as_variables_or_variable(arg) for arg in args)
    names = join_dict_keys(args, how=join)
    grouped_by_name = collect_dict_values(args, names, fill_value)

    result_vars = {}
    for name, variable_args in zip(names, grouped_by_name):
        core_dim_present = _check_core_dims(signature, variable_args, name)
        if core_dim_present is True:
            result_vars[name] = func(*variable_args)
        else:
            if on_missing_core_dim == "raise":
                raise ValueError(core_dim_present)
            elif on_missing_core_dim == "copy":
                result_vars[name] = variable_args[0]
            elif on_missing_core_dim == "drop":
                pass
            else:
                raise ValueError(
                    f"Invalid value for `on_missing_core_dim`: {on_missing_core_dim!r}"
                )

    if signature.num_outputs > 1:
        return _unpack_dict_tuples(result_vars, signature.num_outputs)
    else:
        return result_vars


def _fast_dataset(
    variables: dict[Hashable, Variable],
    coord_variables: Mapping[Hashable, Variable],
    indexes: dict[Hashable, Index],
) -> Dataset:
    """Create a dataset as quickly as possible.

    Beware: the `variables` dict is modified INPLACE.
    """
    from xarray.core.dataset import Dataset

    variables.update(coord_variables)
    coord_names = set(coord_variables)
    return Dataset._construct_direct(variables, coord_names, indexes=indexes)


def apply_dataset_vfunc(
    func,
    *args,
    signature: _UFuncSignature,
    join="inner",
    dataset_join="exact",
    fill_value=_NO_FILL_VALUE,
    exclude_dims=frozenset(),
    keep_attrs="override",
    on_missing_core_dim: MissingCoreDimOptions = "raise",
) -> Dataset | tuple[Dataset, ...]:
    """Apply a variable level function over Dataset, dict of DataArray,
    DataArray, Variable and/or ndarray objects.
    """
    from xarray.core.dataset import Dataset

    if dataset_join not in _JOINS_WITHOUT_FILL_VALUES and fill_value is _NO_FILL_VALUE:
        raise TypeError(
            "to apply an operation to datasets with different "
            "data variables with apply_ufunc, you must supply the "
            "dataset_fill_value argument."
        )

    objs = _all_of_type(args, Dataset)

    if len(args) > 1:
        args = tuple(
            deep_align(
                args,
                join=join,
                copy=False,
                exclude=exclude_dims,
                raise_on_invalid=False,
            )
        )

    list_of_coords, list_of_indexes = build_output_coords_and_indexes(
        args, signature, exclude_dims, combine_attrs=keep_attrs
    )
    args = tuple(getattr(arg, "data_vars", arg) for arg in args)

    result_vars = apply_dict_of_variables_vfunc(
        func,
        *args,
        signature=signature,
        join=dataset_join,
        fill_value=fill_value,
        on_missing_core_dim=on_missing_core_dim,
    )

    out: Dataset | tuple[Dataset, ...]
    if signature.num_outputs > 1:
        out = tuple(
            _fast_dataset(*args)
            for args in zip(result_vars, list_of_coords, list_of_indexes)
        )
    else:
        (coord_vars,) = list_of_coords
        (indexes,) = list_of_indexes
        out = _fast_dataset(result_vars, coord_vars, indexes=indexes)

    attrs = merge_attrs([x.attrs for x in objs], combine_attrs=keep_attrs)
    if isinstance(out, tuple):
        for ds in out:
            ds.attrs = attrs
    else:
        out.attrs = attrs

    return out


def _iter_over_selections(obj, dim, values):
    """Iterate over selections of an xarray object in the provided order."""
    from xarray.core.groupby import _dummy_copy

    dummy = None
    for value in values:
        try:
            obj_sel = obj.sel(**{dim: value})
        except (KeyError, IndexError):
            if dummy is None:
                dummy = _dummy_copy(obj)
            obj_sel = dummy
        yield obj_sel


def apply_groupby_func(func, *args):
    """Apply a dataset or datarray level function over GroupBy, Dataset,
    DataArray, Variable and/or ndarray objects.
    """
    from xarray.core.groupby import GroupBy, peek_at
    from xarray.core.variable import Variable

    groupbys = [arg for arg in args if isinstance(arg, GroupBy)]
    assert groupbys, "must have at least one groupby to iterate over"
    first_groupby = groupbys[0]
    (grouper,) = first_groupby.groupers
    if any(not grouper.group.equals(gb.groupers[0].group) for gb in groupbys[1:]):  # type: ignore[union-attr]
        raise ValueError(
            "apply_ufunc can only perform operations over "
            "multiple GroupBy objects at once if they are all "
            "grouped the same way"
        )

    grouped_dim = grouper.name
    unique_values = grouper.unique_coord.values

    iterators = []
    for arg in args:
        iterator: Iterator[Any]
        if isinstance(arg, GroupBy):
            iterator = (value for _, value in arg)
        elif hasattr(arg, "dims") and grouped_dim in arg.dims:
            if isinstance(arg, Variable):
                raise ValueError(
                    "groupby operations cannot be performed with "
                    "xarray.Variable objects that share a dimension with "
                    "the grouped dimension"
                )
            iterator = _iter_over_selections(arg, grouped_dim, unique_values)
        else:
            iterator = itertools.repeat(arg)
        iterators.append(iterator)

    applied: Iterator = (func(*zipped_args) for zipped_args in zip(*iterators))
    applied_example, applied = peek_at(applied)
    combine = first_groupby._combine  # type: ignore[attr-defined]
    if isinstance(applied_example, tuple):
        combined = tuple(combine(output) for output in zip(*applied))
    else:
        combined = combine(applied)
    return combined


def unified_dim_sizes(
    variables: Iterable[Variable], exclude_dims: Set = frozenset()
) -> dict[Hashable, int]:
    dim_sizes: dict[Hashable, int] = {}

    for var in variables:
        if len(set(var.dims)) < len(var.dims):
            raise ValueError(
                "broadcasting cannot handle duplicate "
                f"dimensions on a variable: {list(var.dims)}"
            )
        for dim, size in zip(var.dims, var.shape):
            if dim not in exclude_dims:
                if dim not in dim_sizes:
                    dim_sizes[dim] = size
                elif dim_sizes[dim] != size:
                    raise ValueError(
                        "operands cannot be broadcast together "
                        "with mismatched lengths for dimension "
                        f"{dim}: {dim_sizes[dim]} vs {size}"
                    )
    return dim_sizes


SLICE_NONE = slice(None)


def broadcast_compat_data(
    variable: Variable,
    broadcast_dims: tuple[Hashable, ...],
    core_dims: tuple[Hashable, ...],
) -> Any:
    data = variable.data

    old_dims = variable.dims
    new_dims = broadcast_dims + core_dims

    if new_dims == old_dims:
        # optimize for the typical case
        return data

    set_old_dims = set(old_dims)
    set_new_dims = set(new_dims)
    unexpected_dims = [d for d in old_dims if d not in set_new_dims]

    if unexpected_dims:
        raise ValueError(
            "operand to apply_ufunc encountered unexpected "
            f"dimensions {unexpected_dims!r} on an input variable: these are core "
            "dimensions on other input or output variables"
        )

    # for consistency with numpy, keep broadcast dimensions to the left
    old_broadcast_dims = tuple(d for d in broadcast_dims if d in set_old_dims)
    reordered_dims = old_broadcast_dims + core_dims
    if reordered_dims != old_dims:
        order = tuple(old_dims.index(d) for d in reordered_dims)
        data = duck_array_ops.transpose(data, order)

    if new_dims != reordered_dims:
        key_parts: list[slice | None] = []
        for dim in new_dims:
            if dim in set_old_dims:
                key_parts.append(SLICE_NONE)
            elif key_parts:
                # no need to insert new axes at the beginning that are already
                # handled by broadcasting
                key_parts.append(np.newaxis)
        data = data[tuple(key_parts)]

    return data


def _vectorize(func, signature, output_dtypes, exclude_dims):
    if signature.all_core_dims:
        func = np.vectorize(
            func,
            otypes=output_dtypes,
            signature=signature.to_gufunc_string(exclude_dims),
        )
    else:
        func = np.vectorize(func, otypes=output_dtypes)

    return func


def apply_variable_ufunc(
    func,
    *args,
    signature: _UFuncSignature,
    exclude_dims=frozenset(),
    dask="forbidden",
    output_dtypes=None,
    vectorize=False,
    keep_attrs="override",
    dask_gufunc_kwargs=None,
) -> Variable | tuple[Variable, ...]:
    """Apply a ndarray level function over Variable and/or ndarray objects."""
    from xarray.core.formatting import short_array_repr
    from xarray.core.variable import Variable, as_compatible_data

    dim_sizes = unified_dim_sizes(
        (a for a in args if hasattr(a, "dims")), exclude_dims=exclude_dims
    )
    broadcast_dims = tuple(
        dim for dim in dim_sizes if dim not in signature.all_core_dims
    )
    output_dims = [broadcast_dims + out for out in signature.output_core_dims]

    input_data = [
        (
            broadcast_compat_data(arg, broadcast_dims, core_dims)
            if isinstance(arg, Variable)
            else arg
        )
        for arg, core_dims in zip(args, signature.input_core_dims)
    ]

    if any(is_chunked_array(array) for array in input_data):
        if dask == "forbidden":
            raise ValueError(
                "apply_ufunc encountered a chunked array on an "
                "argument, but handling for chunked arrays has not "
                "been enabled. Either set the ``dask`` argument "
                "or load your data into memory first with "
                "``.load()`` or ``.compute()``"
            )
        elif dask == "parallelized":
            chunkmanager = get_chunked_array_type(*input_data)

            numpy_func = func

            if dask_gufunc_kwargs is None:
                dask_gufunc_kwargs = {}
            else:
                dask_gufunc_kwargs = dask_gufunc_kwargs.copy()

            allow_rechunk = dask_gufunc_kwargs.get("allow_rechunk", None)
            if allow_rechunk is None:
                for n, (data, core_dims) in enumerate(
                    zip(input_data, signature.input_core_dims)
                ):
                    if is_chunked_array(data):
                        # core dimensions cannot span multiple chunks
                        for axis, dim in enumerate(core_dims, start=-len(core_dims)):
                            if len(data.chunks[axis]) != 1:
                                raise ValueError(
                                    f"dimension {dim} on {n}th function argument to "
                                    "apply_ufunc with dask='parallelized' consists of "
                                    "multiple chunks, but is also a core dimension. To "
                                    "fix, either rechunk into a single array chunk along "
                                    f"this dimension, i.e., ``.chunk(dict({dim}=-1))``, or "
                                    "pass ``allow_rechunk=True`` in ``dask_gufunc_kwargs`` "
                                    "but beware that this may significantly increase memory usage."
                                )
                dask_gufunc_kwargs["allow_rechunk"] = True

            output_sizes = dask_gufunc_kwargs.pop("output_sizes", {})
            if output_sizes:
                output_sizes_renamed = {}
                for key, value in output_sizes.items():
                    if key not in signature.all_output_core_dims:
                        raise ValueError(
                            f"dimension '{key}' in 'output_sizes' must correspond to output_core_dims"
                        )
                    output_sizes_renamed[signature.dims_map[key]] = value
                dask_gufunc_kwargs["output_sizes"] = output_sizes_renamed

            for key in signature.all_output_core_dims:
                if (
                    key not in signature.all_input_core_dims or key in exclude_dims
                ) and key not in output_sizes:
                    raise ValueError(
                        f"dimension '{key}' in 'output_core_dims' needs corresponding (dim, size) in 'output_sizes'"
                    )

            def func(*arrays):
                res = chunkmanager.apply_gufunc(
                    numpy_func,
                    signature.to_gufunc_string(exclude_dims),
                    *arrays,
                    vectorize=vectorize,
                    output_dtypes=output_dtypes,
                    **dask_gufunc_kwargs,
                )

                return res

        elif dask == "allowed":
            pass
        else:
            raise ValueError(
                f"unknown setting for chunked array handling in apply_ufunc: {dask}"
            )
    else:
        if vectorize:
            func = _vectorize(
                func, signature, output_dtypes=output_dtypes, exclude_dims=exclude_dims
            )

    result_data = func(*input_data)

    if signature.num_outputs == 1:
        result_data = (result_data,)
    elif (
        not isinstance(result_data, tuple) or len(result_data) != signature.num_outputs
    ):
        raise ValueError(
            f"applied function does not have the number of "
            f"outputs specified in the ufunc signature. "
            f"Received a {type(result_data)} with {len(result_data)} elements. "
            f"Expected a tuple of {signature.num_outputs} elements:\n\n"
            f"{limit_lines(repr(result_data), limit=10)}"
        )

    objs = _all_of_type(args, Variable)
    attrs = merge_attrs(
        [obj.attrs for obj in objs],
        combine_attrs=keep_attrs,
    )

    output: list[Variable] = []
    for dims, data in zip(output_dims, result_data):
        data = as_compatible_data(data)
        if data.ndim != len(dims):
            raise ValueError(
                "applied function returned data with an unexpected "
                f"number of dimensions. Received {data.ndim} dimension(s) but "
                f"expected {len(dims)} dimensions with names {dims!r}, from:\n\n"
                f"{short_array_repr(data)}"
            )

        var = Variable(dims, data, fastpath=True)
        for dim, new_size in var.sizes.items():
            if dim in dim_sizes and new_size != dim_sizes[dim]:
                raise ValueError(
                    f"size of dimension '{dim}' on inputs was unexpectedly "
                    f"changed by applied function from {dim_sizes[dim]} to {new_size}. Only "
                    "dimensions specified in ``exclude_dims`` with "
                    "xarray.apply_ufunc are allowed to change size. "
                    "The data returned was:\n\n"
                    f"{short_array_repr(data)}"
                )

        var.attrs = attrs
        output.append(var)

    if signature.num_outputs == 1:
        return output[0]
    else:
        return tuple(output)


def apply_array_ufunc(func, *args, dask="forbidden"):
    """Apply a ndarray level function over ndarray objects."""
    if any(is_chunked_array(arg) for arg in args):
        if dask == "forbidden":
            raise ValueError(
                "apply_ufunc encountered a dask array on an "
                "argument, but handling for dask arrays has not "
                "been enabled. Either set the ``dask`` argument "
                "or load your data into memory first with "
                "``.load()`` or ``.compute()``"
            )
        elif dask == "parallelized":
            raise ValueError(
                "cannot use dask='parallelized' for apply_ufunc "
                "unless at least one input is an xarray object"
            )
        elif dask == "allowed":
            pass
        else:
            raise ValueError(f"unknown setting for dask array handling: {dask}")
    return func(*args)


def apply_ufunc(
    func: Callable,
    *args: Any,
    input_core_dims: Sequence[Sequence] | None = None,
    output_core_dims: Sequence[Sequence] | None = ((),),
    exclude_dims: Set = frozenset(),
    vectorize: bool = False,
    join: JoinOptions = "exact",
    dataset_join: str = "exact",
    dataset_fill_value: object = _NO_FILL_VALUE,
    keep_attrs: bool | str | None = None,
    kwargs: Mapping | None = None,
    dask: Literal["forbidden", "allowed", "parallelized"] = "forbidden",
    output_dtypes: Sequence | None = None,
    output_sizes: Mapping[Any, int] | None = None,
    meta: Any = None,
    dask_gufunc_kwargs: dict[str, Any] | None = None,
    on_missing_core_dim: MissingCoreDimOptions = "raise",
) -> Any:
    """Apply a vectorized function for unlabeled arrays on xarray objects.

    The function will be mapped over the data variable(s) of the input
    arguments using xarray's standard rules for labeled computation, including
    alignment, broadcasting, looping over GroupBy/Dataset variables, and
    merging of coordinates.

    Parameters
    ----------
    func : callable
        Function to call like ``func(*args, **kwargs)`` on unlabeled arrays
        (``.data``) that returns an array or tuple of arrays. If multiple
        arguments with non-matching dimensions are supplied, this function is
        expected to vectorize (broadcast) over axes of positional arguments in
        the style of NumPy universal functions [1]_ (if this is not the case,
        set ``vectorize=True``). If this function returns multiple outputs, you
        must set ``output_core_dims`` as well.
    *args : Dataset, DataArray, DataArrayGroupBy, DatasetGroupBy, Variable, \
        numpy.ndarray, dask.array.Array or scalar
        Mix of labeled and/or unlabeled arrays to which to apply the function.
    input_core_dims : sequence of sequence, optional
        List of the same length as ``args`` giving the list of core dimensions
        on each input argument that should not be broadcast. By default, we
        assume there are no core dimensions on any input arguments.

        For example, ``input_core_dims=[[], ['time']]`` indicates that all
        dimensions on the first argument and all dimensions other than 'time'
        on the second argument should be broadcast.

        Core dimensions are automatically moved to the last axes of input
        variables before applying ``func``, which facilitates using NumPy style
        generalized ufuncs [2]_.
    output_core_dims : list of tuple, optional
        List of the same length as the number of output arguments from
        ``func``, giving the list of core dimensions on each output that were
        not broadcast on the inputs. By default, we assume that ``func``
        outputs exactly one array, with axes corresponding to each broadcast
        dimension.

        Core dimensions are assumed to appear as the last dimensions of each
        output in the provided order.
    exclude_dims : set, optional
        Core dimensions on the inputs to exclude from alignment and
        broadcasting entirely. Any input coordinates along these dimensions
        will be dropped. Each excluded dimension must also appear in
        ``input_core_dims`` for at least one argument. Only dimensions listed
        here are allowed to change size between input and output objects.
    vectorize : bool, optional
        If True, then assume ``func`` only takes arrays defined over core
        dimensions as input and vectorize it automatically with
        :py:func:`numpy.vectorize`. This option exists for convenience, but is
        almost always slower than supplying a pre-vectorized function.
    join : {"outer", "inner", "left", "right", "exact"}, default: "exact"
        Method for joining the indexes of the passed objects along each
        dimension, and the variables of Dataset objects with mismatched
        data variables:

        - 'outer': use the union of object indexes
        - 'inner': use the intersection of object indexes
        - 'left': use indexes from the first object with each dimension
        - 'right': use indexes from the last object with each dimension
        - 'exact': raise `ValueError` instead of aligning when indexes to be
          aligned are not equal
    dataset_join : {"outer", "inner", "left", "right", "exact"}, default: "exact"
        Method for joining variables of Dataset objects with mismatched
        data variables.

        - 'outer': take variables from both Dataset objects
        - 'inner': take only overlapped variables
        - 'left': take only variables from the first object
        - 'right': take only variables from the last object
        - 'exact': data variables on all Dataset objects must match exactly
    dataset_fill_value : optional
        Value used in place of missing variables on Dataset inputs when the
        datasets do not share the exact same ``data_vars``. Required if
        ``dataset_join not in {'inner', 'exact'}``, otherwise ignored.
    keep_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts", "override"} or bool, optional
        - 'drop' or False: empty attrs on returned xarray object.
        - 'identical': all attrs must be the same on every object.
        - 'no_conflicts': attrs from all objects are combined, any that have the same name must also have the same value.
        - 'drop_conflicts': attrs from all objects are combined, any that have the same name but different values are dropped.
        - 'override' or True: skip comparing and copy attrs from the first object to the result.
    kwargs : dict, optional
        Optional keyword arguments passed directly on to call ``func``.
    dask : {"forbidden", "allowed", "parallelized"}, default: "forbidden"
        How to handle applying to objects containing lazy data in the form of
        dask arrays:

        - 'forbidden' (default): raise an error if a dask array is encountered.
        - 'allowed': pass dask arrays directly on to ``func``. Prefer this option if
          ``func`` natively supports dask arrays.
        - 'parallelized': automatically parallelize ``func`` if any of the
          inputs are a dask array by using :py:func:`dask.array.apply_gufunc`. Multiple output
          arguments are supported. Only use this option if ``func`` does not natively
          support dask arrays (e.g. converts them to numpy arrays).
    dask_gufunc_kwargs : dict, optional
        Optional keyword arguments passed to :py:func:`dask.array.apply_gufunc` if
        dask='parallelized'. Possible keywords are ``output_sizes``, ``allow_rechunk``
        and ``meta``.
    output_dtypes : list of dtype, optional
        Optional list of output dtypes. Only used if ``dask='parallelized'`` or
        ``vectorize=True``.
    output_sizes : dict, optional
        Optional mapping from dimension names to sizes for outputs. Only used
        if dask='parallelized' and new dimensions (not found on inputs) appear
        on outputs. ``output_sizes`` should be given in the ``dask_gufunc_kwargs``
        parameter. It will be removed as direct parameter in a future version.
    meta : optional
        Size-0 object representing the type of array wrapped by dask array. Passed on to
        :py:func:`dask.array.apply_gufunc`. ``meta`` should be given in the
        ``dask_gufunc_kwargs`` parameter . It will be removed as direct parameter
        a future version.
    on_missing_core_dim : {"raise", "copy", "drop"}, default: "raise"
        How to handle missing core dimensions on input variables.

    Returns
    -------
    Single value or tuple of Dataset, DataArray, Variable, dask.array.Array or
    numpy.ndarray, the first type on that list to appear on an input.

    Notes
    -----
    This function is designed for the more common case where ``func`` can work on numpy
    arrays. If ``func`` needs to manipulate a whole xarray object subset to each block
    it is possible to use :py:func:`xarray.map_blocks`.

    Note that due to the overhead :py:func:`xarray.map_blocks` is considerably slower than ``apply_ufunc``.

    Examples
    --------
    Calculate the vector magnitude of two arguments:

    >>> def magnitude(a, b):
    ...     func = lambda x, y: np.sqrt(x**2 + y**2)
    ...     return xr.apply_ufunc(func, a, b)
    ...

    You can now apply ``magnitude()`` to :py:class:`DataArray` and :py:class:`Dataset`
    objects, with automatically preserved dimensions and coordinates, e.g.,

    >>> array = xr.DataArray([1, 2, 3], coords=[("x", [0.1, 0.2, 0.3])])
    >>> magnitude(array, -array)
    <xarray.DataArray (x: 3)> Size: 24B
    array([1.41421356, 2.82842712, 4.24264069])
    Coordinates:
      * x        (x) float64 24B 0.1 0.2 0.3

    Plain scalars, numpy arrays and a mix of these with xarray objects is also
    supported:

    >>> magnitude(3, 4)
    5.0
    >>> magnitude(3, np.array([0, 4]))
    array([3., 5.])
    >>> magnitude(array, 0)
    <xarray.DataArray (x: 3)> Size: 24B
    array([1., 2., 3.])
    Coordinates:
      * x        (x) float64 24B 0.1 0.2 0.3

    Other examples of how you could use ``apply_ufunc`` to write functions to
    (very nearly) replicate existing xarray functionality:

    Compute the mean (``.mean``) over one dimension:

    >>> def mean(obj, dim):
    ...     # note: apply always moves core dimensions to the end
    ...     return apply_ufunc(
    ...         np.mean, obj, input_core_dims=[[dim]], kwargs={"axis": -1}
    ...     )
    ...

    Inner product over a specific dimension (like :py:func:`dot`):

    >>> def _inner(x, y):
    ...     result = np.matmul(x[..., np.newaxis, :], y[..., :, np.newaxis])
    ...     return result[..., 0, 0]
    ...
    >>> def inner_product(a, b, dim):
    ...     return apply_ufunc(_inner, a, b, input_core_dims=[[dim], [dim]])
    ...

    Stack objects along a new dimension (like :py:func:`concat`):

    >>> def stack(objects, dim, new_coord):
    ...     # note: this version does not stack coordinates
    ...     func = lambda *x: np.stack(x, axis=-1)
    ...     result = apply_ufunc(
    ...         func,
    ...         *objects,
    ...         output_core_dims=[[dim]],
    ...         join="outer",
    ...         dataset_fill_value=np.nan
    ...     )
    ...     result[dim] = new_coord
    ...     return result
    ...

    If your function is not vectorized but can be applied only to core
    dimensions, you can use ``vectorize=True`` to turn into a vectorized
    function. This wraps :py:func:`numpy.vectorize`, so the operation isn't
    terribly fast. Here we'll use it to calculate the distance between
    empirical samples from two probability distributions, using a scipy
    function that needs to be applied to vectors:

    >>> import scipy.stats
    >>> def earth_mover_distance(first_samples, second_samples, dim="ensemble"):
    ...     return apply_ufunc(
    ...         scipy.stats.wasserstein_distance,
    ...         first_samples,
    ...         second_samples,
    ...         input_core_dims=[[dim], [dim]],
    ...         vectorize=True,
    ...     )
    ...

    Most of NumPy's builtin functions already broadcast their inputs
    appropriately for use in ``apply_ufunc``. You may find helper functions such as
    :py:func:`numpy.broadcast_arrays` helpful in writing your function. ``apply_ufunc`` also
    works well with :py:func:`numba.vectorize` and :py:func:`numba.guvectorize`.

    See Also
    --------
    numpy.broadcast_arrays
    numba.vectorize
    numba.guvectorize
    dask.array.apply_gufunc
    xarray.map_blocks

    :ref:`dask.automatic-parallelization`
        User guide describing :py:func:`apply_ufunc` and :py:func:`map_blocks`.

    :doc:`xarray-tutorial:advanced/apply_ufunc/apply_ufunc`
        Advanced Tutorial on applying numpy function using :py:func:`apply_ufunc`

    References
    ----------
    .. [1] https://numpy.org/doc/stable/reference/ufuncs.html
    .. [2] https://numpy.org/doc/stable/reference/c-api/generalized-ufuncs.html
    """
    from xarray.core.dataarray import DataArray
    from xarray.core.groupby import GroupBy
    from xarray.core.variable import Variable

    if input_core_dims is None:
        input_core_dims = ((),) * (len(args))
    elif len(input_core_dims) != len(args):
        raise ValueError(
            f"input_core_dims must be None or a tuple with the length same to "
            f"the number of arguments. "
            f"Given {len(input_core_dims)} input_core_dims: {input_core_dims}, "
            f" but number of args is {len(args)}."
        )

    if kwargs is None:
        kwargs = {}

    signature = _UFuncSignature(input_core_dims, output_core_dims)

    if exclude_dims:
        if not isinstance(exclude_dims, set):
            raise TypeError(
                f"Expected exclude_dims to be a 'set'. Received '{type(exclude_dims).__name__}' instead."
            )
        if not exclude_dims <= signature.all_core_dims:
            raise ValueError(
                f"each dimension in `exclude_dims` must also be a "
                f"core dimension in the function signature. "
                f"Please make {(exclude_dims - signature.all_core_dims)} a core dimension"
            )

    # handle dask_gufunc_kwargs
    if dask == "parallelized":
        if dask_gufunc_kwargs is None:
            dask_gufunc_kwargs = {}
        else:
            dask_gufunc_kwargs = dask_gufunc_kwargs.copy()
        # todo: remove warnings after deprecation cycle
        if meta is not None:
            warnings.warn(
                "``meta`` should be given in the ``dask_gufunc_kwargs`` parameter."
                " It will be removed as direct parameter in a future version.",
                FutureWarning,
                stacklevel=2,
            )
            dask_gufunc_kwargs.setdefault("meta", meta)
        if output_sizes is not None:
            warnings.warn(
                "``output_sizes`` should be given in the ``dask_gufunc_kwargs`` "
                "parameter. It will be removed as direct parameter in a future "
                "version.",
                FutureWarning,
                stacklevel=2,
            )
            dask_gufunc_kwargs.setdefault("output_sizes", output_sizes)

    if kwargs:
        func = functools.partial(func, **kwargs)

    if keep_attrs is None:
        keep_attrs = _get_keep_attrs(default=False)

    if isinstance(keep_attrs, bool):
        keep_attrs = "override" if keep_attrs else "drop"

    variables_vfunc = functools.partial(
        apply_variable_ufunc,
        func,
        signature=signature,
        exclude_dims=exclude_dims,
        keep_attrs=keep_attrs,
        dask=dask,
        vectorize=vectorize,
        output_dtypes=output_dtypes,
        dask_gufunc_kwargs=dask_gufunc_kwargs,
    )

    # feed groupby-apply_ufunc through apply_groupby_func
    if any(isinstance(a, GroupBy) for a in args):
        this_apply = functools.partial(
            apply_ufunc,
            func,
            input_core_dims=input_core_dims,
            output_core_dims=output_core_dims,
            exclude_dims=exclude_dims,
            join=join,
            dataset_join=dataset_join,
            dataset_fill_value=dataset_fill_value,
            keep_attrs=keep_attrs,
            dask=dask,
            vectorize=vectorize,
            output_dtypes=output_dtypes,
            dask_gufunc_kwargs=dask_gufunc_kwargs,
        )
        return apply_groupby_func(this_apply, *args)
    # feed datasets apply_variable_ufunc through apply_dataset_vfunc
    elif any(is_dict_like(a) for a in args):
        return apply_dataset_vfunc(
            variables_vfunc,
            *args,
            signature=signature,
            join=join,
            exclude_dims=exclude_dims,
            dataset_join=dataset_join,
            fill_value=dataset_fill_value,
            keep_attrs=keep_attrs,
            on_missing_core_dim=on_missing_core_dim,
        )
    # feed DataArray apply_variable_ufunc through apply_dataarray_vfunc
    elif any(isinstance(a, DataArray) for a in args):
        return apply_dataarray_vfunc(
            variables_vfunc,
            *args,
            signature=signature,
            join=join,
            exclude_dims=exclude_dims,
            keep_attrs=keep_attrs,
        )
    # feed Variables directly through apply_variable_ufunc
    elif any(isinstance(a, Variable) for a in args):
        return variables_vfunc(*args)
    else:
        # feed anything else through apply_array_ufunc
        return apply_array_ufunc(func, *args, dask=dask)


def cov(
    da_a: T_DataArray,
    da_b: T_DataArray,
    dim: Dims = None,
    ddof: int = 1,
    weights: T_DataArray | None = None,
) -> T_DataArray:
    """
    Compute covariance between two DataArray objects along a shared dimension.

    Parameters
    ----------
    da_a : DataArray
        Array to compute.
    da_b : DataArray
        Array to compute.
    dim : str, iterable of hashable, "..." or None, optional
        The dimension along which the covariance will be computed
    ddof : int, default: 1
        If ddof=1, covariance is normalized by N-1, giving an unbiased estimate,
        else normalization is by N.
    weights : DataArray, optional
        Array of weights.

    Returns
    -------
    covariance : DataArray

    See Also
    --------
    pandas.Series.cov : corresponding pandas function
    xarray.corr : respective function to calculate correlation

    Examples
    --------
    >>> from xarray import DataArray
    >>> da_a = DataArray(
    ...     np.array([[1, 2, 3], [0.1, 0.2, 0.3], [3.2, 0.6, 1.8]]),
    ...     dims=("space", "time"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...         ("time", pd.date_range("2000-01-01", freq="1D", periods=3)),
    ...     ],
    ... )
    >>> da_a
    <xarray.DataArray (space: 3, time: 3)> Size: 72B
    array([[1. , 2. , 3. ],
           [0.1, 0.2, 0.3],
           [3.2, 0.6, 1.8]])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    >>> da_b = DataArray(
    ...     np.array([[0.2, 0.4, 0.6], [15, 10, 5], [3.2, 0.6, 1.8]]),
    ...     dims=("space", "time"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...         ("time", pd.date_range("2000-01-01", freq="1D", periods=3)),
    ...     ],
    ... )
    >>> da_b
    <xarray.DataArray (space: 3, time: 3)> Size: 72B
    array([[ 0.2,  0.4,  0.6],
           [15. , 10. ,  5. ],
           [ 3.2,  0.6,  1.8]])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    >>> xr.cov(da_a, da_b)
    <xarray.DataArray ()> Size: 8B
    array(-3.53055556)
    >>> xr.cov(da_a, da_b, dim="time")
    <xarray.DataArray (space: 3)> Size: 24B
    array([ 0.2       , -0.5       ,  1.69333333])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
    >>> weights = DataArray(
    ...     [4, 2, 1],
    ...     dims=("space"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...     ],
    ... )
    >>> weights
    <xarray.DataArray (space: 3)> Size: 24B
    array([4, 2, 1])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
    >>> xr.cov(da_a, da_b, dim="space", weights=weights)
    <xarray.DataArray (time: 3)> Size: 24B
    array([-4.69346939, -4.49632653, -3.37959184])
    Coordinates:
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    """
    from xarray.core.dataarray import DataArray

    if any(not isinstance(arr, DataArray) for arr in [da_a, da_b]):
        raise TypeError(
            "Only xr.DataArray is supported."
            f"Given {[type(arr) for arr in [da_a, da_b]]}."
        )
    if weights is not None:
        if not isinstance(weights, DataArray):
            raise TypeError(f"Only xr.DataArray is supported. Given {type(weights)}.")
    return _cov_corr(da_a, da_b, weights=weights, dim=dim, ddof=ddof, method="cov")


def corr(
    da_a: T_DataArray,
    da_b: T_DataArray,
    dim: Dims = None,
    weights: T_DataArray | None = None,
) -> T_DataArray:
    """
    Compute the Pearson correlation coefficient between
    two DataArray objects along a shared dimension.

    Parameters
    ----------
    da_a : DataArray
        Array to compute.
    da_b : DataArray
        Array to compute.
    dim : str, iterable of hashable, "..." or None, optional
        The dimension along which the correlation will be computed
    weights : DataArray, optional
        Array of weights.

    Returns
    -------
    correlation: DataArray

    See Also
    --------
    pandas.Series.corr : corresponding pandas function
    xarray.cov : underlying covariance function

    Examples
    --------
    >>> from xarray import DataArray
    >>> da_a = DataArray(
    ...     np.array([[1, 2, 3], [0.1, 0.2, 0.3], [3.2, 0.6, 1.8]]),
    ...     dims=("space", "time"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...         ("time", pd.date_range("2000-01-01", freq="1D", periods=3)),
    ...     ],
    ... )
    >>> da_a
    <xarray.DataArray (space: 3, time: 3)> Size: 72B
    array([[1. , 2. , 3. ],
           [0.1, 0.2, 0.3],
           [3.2, 0.6, 1.8]])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    >>> da_b = DataArray(
    ...     np.array([[0.2, 0.4, 0.6], [15, 10, 5], [3.2, 0.6, 1.8]]),
    ...     dims=("space", "time"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...         ("time", pd.date_range("2000-01-01", freq="1D", periods=3)),
    ...     ],
    ... )
    >>> da_b
    <xarray.DataArray (space: 3, time: 3)> Size: 72B
    array([[ 0.2,  0.4,  0.6],
           [15. , 10. ,  5. ],
           [ 3.2,  0.6,  1.8]])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    >>> xr.corr(da_a, da_b)
    <xarray.DataArray ()> Size: 8B
    array(-0.57087777)
    >>> xr.corr(da_a, da_b, dim="time")
    <xarray.DataArray (space: 3)> Size: 24B
    array([ 1., -1.,  1.])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
    >>> weights = DataArray(
    ...     [4, 2, 1],
    ...     dims=("space"),
    ...     coords=[
    ...         ("space", ["IA", "IL", "IN"]),
    ...     ],
    ... )
    >>> weights
    <xarray.DataArray (space: 3)> Size: 24B
    array([4, 2, 1])
    Coordinates:
      * space    (space) <U2 24B 'IA' 'IL' 'IN'
    >>> xr.corr(da_a, da_b, dim="space", weights=weights)
    <xarray.DataArray (time: 3)> Size: 24B
    array([-0.50240504, -0.83215028, -0.99057446])
    Coordinates:
      * time     (time) datetime64[ns] 24B 2000-01-01 2000-01-02 2000-01-03
    """
    from xarray.core.dataarray import DataArray

    if any(not isinstance(arr, DataArray) for arr in [da_a, da_b]):
        raise TypeError(
            "Only xr.DataArray is supported."
            f"Given {[type(arr) for arr in [da_a, da_b]]}."
        )
    if weights is not None:
        if not isinstance(weights, DataArray):
            raise TypeError(f"Only xr.DataArray is supported. Given {type(weights)}.")
    return _cov_corr(da_a, da_b, weights=weights, dim=dim, method="corr")


def _cov_corr(
    da_a: T_DataArray,
    da_b: T_DataArray,
    weights: T_DataArray | None = None,
    dim: Dims = None,
    ddof: int = 0,
    method: Literal["cov", "corr", None] = None,
) -> T_DataArray:
    """
    Internal method for xr.cov() and xr.corr() so only have to
    sanitize the input arrays once and we don't repeat code.
    """
    # 1. Broadcast the two arrays
    da_a, da_b = align(da_a, da_b, join="inner", copy=False)

    # 2. Ignore the nans
    valid_values = da_a.notnull() & da_b.notnull()
    da_a = da_a.where(valid_values)
    da_b = da_b.where(valid_values)

    # 3. Detrend along the given dim
    if weights is not None:
        demeaned_da_a = da_a - da_a.weighted(weights).mean(dim=dim)
        demeaned_da_b = da_b - da_b.weighted(weights).mean(dim=dim)
    else:
        demeaned_da_a = da_a - da_a.mean(dim=dim)
        demeaned_da_b = da_b - da_b.mean(dim=dim)

    # 4. Compute covariance along the given dim
    # N.B. `skipna=True` is required or auto-covariance is computed incorrectly. E.g.
    # Try xr.cov(da,da) for da = xr.DataArray([[1, 2], [1, np.nan]], dims=["x", "time"])
    if weights is not None:
        cov = (
            (demeaned_da_a.conj() * demeaned_da_b)
            .weighted(weights)
            .mean(dim=dim, skipna=True)
        )
    else:
        cov = (demeaned_da_a.conj() * demeaned_da_b).mean(dim=dim, skipna=True)

    if method == "cov":
        # Adjust covariance for degrees of freedom
        valid_count = valid_values.sum(dim)
        adjust = valid_count / (valid_count - ddof)
        # I think the cast is required because of `T_DataArray` + `T_Xarray` (would be
        # the same with `T_DatasetOrArray`)
        # https://github.com/pydata/xarray/pull/8384#issuecomment-1784228026
        return cast(T_DataArray, cov * adjust)

    else:
        # Compute std and corr
        if weights is not None:
            da_a_std = da_a.weighted(weights).std(dim=dim)
            da_b_std = da_b.weighted(weights).std(dim=dim)
        else:
            da_a_std = da_a.std(dim=dim)
            da_b_std = da_b.std(dim=dim)
        corr = cov / (da_a_std * da_b_std)
        return cast(T_DataArray, corr)


def cross(
    a: DataArray | Variable, b: DataArray | Variable, *, dim: Hashable
) -> DataArray | Variable:
    """
    Compute the cross product of two (arrays of) vectors.

    The cross product of `a` and `b` in :math:`R^3` is a vector
    perpendicular to both `a` and `b`. The vectors in `a` and `b` are
    defined by the values along the dimension `dim` and can have sizes
    1, 2 or 3. Where the size of either `a` or `b` is
    1 or 2, the remaining components of the input vector is assumed to
    be zero and the cross product calculated accordingly. In cases where
    both input vectors have dimension 2, the z-component of the cross
    product is returned.

    Parameters
    ----------
    a, b : DataArray or Variable
        Components of the first and second vector(s).
    dim : hashable
        The dimension along which the cross product will be computed.
        Must be available in both vectors.

    Examples
    --------
    Vector cross-product with 3 dimensions:

    >>> a = xr.DataArray([1, 2, 3])
    >>> b = xr.DataArray([4, 5, 6])
    >>> xr.cross(a, b, dim="dim_0")
    <xarray.DataArray (dim_0: 3)> Size: 24B
    array([-3,  6, -3])
    Dimensions without coordinates: dim_0

    Vector cross-product with 2 dimensions, returns in the perpendicular
    direction:

    >>> a = xr.DataArray([1, 2])
    >>> b = xr.DataArray([4, 5])
    >>> xr.cross(a, b, dim="dim_0")
    <xarray.DataArray ()> Size: 8B
    array(-3)

    Vector cross-product with 3 dimensions but zeros at the last axis
    yields the same results as with 2 dimensions:

    >>> a = xr.DataArray([1, 2, 0])
    >>> b = xr.DataArray([4, 5, 0])
    >>> xr.cross(a, b, dim="dim_0")
    <xarray.DataArray (dim_0: 3)> Size: 24B
    array([ 0,  0, -3])
    Dimensions without coordinates: dim_0

    One vector with dimension 2:

    >>> a = xr.DataArray(
    ...     [1, 2],
    ...     dims=["cartesian"],
    ...     coords=dict(cartesian=(["cartesian"], ["x", "y"])),
    ... )
    >>> b = xr.DataArray(
    ...     [4, 5, 6],
    ...     dims=["cartesian"],
    ...     coords=dict(cartesian=(["cartesian"], ["x", "y", "z"])),
    ... )
    >>> xr.cross(a, b, dim="cartesian")
    <xarray.DataArray (cartesian: 3)> Size: 24B
    array([12, -6, -3])
    Coordinates:
      * cartesian  (cartesian) <U1 12B 'x' 'y' 'z'

    One vector with dimension 2 but coords in other positions:

    >>> a = xr.DataArray(
    ...     [1, 2],
    ...     dims=["cartesian"],
    ...     coords=dict(cartesian=(["cartesian"], ["x", "z"])),
    ... )
    >>> b = xr.DataArray(
    ...     [4, 5, 6],
    ...     dims=["cartesian"],
    ...     coords=dict(cartesian=(["cartesian"], ["x", "y", "z"])),
    ... )
    >>> xr.cross(a, b, dim="cartesian")
    <xarray.DataArray (cartesian: 3)> Size: 24B
    array([-10,   2,   5])
    Coordinates:
      * cartesian  (cartesian) <U1 12B 'x' 'y' 'z'

    Multiple vector cross-products. Note that the direction of the
    cross product vector is defined by the right-hand rule:

    >>> a = xr.DataArray(
    ...     [[1, 2, 3], [4, 5, 6]],
    ...     dims=("time", "cartesian"),
    ...     coords=dict(
    ...         time=(["time"], [0, 1]),
    ...         cartesian=(["cartesian"], ["x", "y", "z"]),
    ...     ),
    ... )
    >>> b = xr.DataArray(
    ...     [[4, 5, 6], [1, 2, 3]],
    ...     dims=("time", "cartesian"),
    ...     coords=dict(
    ...         time=(["time"], [0, 1]),
    ...         cartesian=(["cartesian"], ["x", "y", "z"]),
    ...     ),
    ... )
    >>> xr.cross(a, b, dim="cartesian")
    <xarray.DataArray (time: 2, cartesian: 3)> Size: 48B
    array([[-3,  6, -3],
           [ 3, -6,  3]])
    Coordinates:
      * time       (time) int64 16B 0 1
      * cartesian  (cartesian) <U1 12B 'x' 'y' 'z'

    Cross can be called on Datasets by converting to DataArrays and later
    back to a Dataset:

    >>> ds_a = xr.Dataset(dict(x=("dim_0", [1]), y=("dim_0", [2]), z=("dim_0", [3])))
    >>> ds_b = xr.Dataset(dict(x=("dim_0", [4]), y=("dim_0", [5]), z=("dim_0", [6])))
    >>> c = xr.cross(
    ...     ds_a.to_dataarray("cartesian"),
    ...     ds_b.to_dataarray("cartesian"),
    ...     dim="cartesian",
    ... )
    >>> c.to_dataset(dim="cartesian")
    <xarray.Dataset> Size: 24B
    Dimensions:  (dim_0: 1)
    Dimensions without coordinates: dim_0
    Data variables:
        x        (dim_0) int64 8B -3
        y        (dim_0) int64 8B 6
        z        (dim_0) int64 8B -3

    See Also
    --------
    numpy.cross : Corresponding numpy function
    """

    if dim not in a.dims:
        raise ValueError(f"Dimension {dim!r} not on a")
    elif dim not in b.dims:
        raise ValueError(f"Dimension {dim!r} not on b")

    if not 1 <= a.sizes[dim] <= 3:
        raise ValueError(
            f"The size of {dim!r} on a must be 1, 2, or 3 to be "
            f"compatible with a cross product but is {a.sizes[dim]}"
        )
    elif not 1 <= b.sizes[dim] <= 3:
        raise ValueError(
            f"The size of {dim!r} on b must be 1, 2, or 3 to be "
            f"compatible with a cross product but is {b.sizes[dim]}"
        )

    all_dims = list(dict.fromkeys(a.dims + b.dims))

    if a.sizes[dim] != b.sizes[dim]:
        # Arrays have different sizes. Append zeros where the smaller
        # array is missing a value, zeros will not affect np.cross:

        if (
            not isinstance(a, Variable)  # Only used to make mypy happy.
            and dim in getattr(a, "coords", {})
            and not isinstance(b, Variable)  # Only used to make mypy happy.
            and dim in getattr(b, "coords", {})
        ):
            # If the arrays have coords we know which indexes to fill
            # with zeros:
            a, b = align(
                a,
                b,
                fill_value=0,
                join="outer",
                exclude=set(all_dims) - {dim},
            )
        elif min(a.sizes[dim], b.sizes[dim]) == 2:
            # If the array doesn't have coords we can only infer
            # that it has composite values if the size is at least 2.
            # Once padded, rechunk the padded array because apply_ufunc
            # requires core dimensions not to be chunked:
            if a.sizes[dim] < b.sizes[dim]:
                a = a.pad({dim: (0, 1)}, constant_values=0)
                # TODO: Should pad or apply_ufunc handle correct chunking?
                a = a.chunk({dim: -1}) if is_duck_dask_array(a.data) else a
            else:
                b = b.pad({dim: (0, 1)}, constant_values=0)
                # TODO: Should pad or apply_ufunc handle correct chunking?
                b = b.chunk({dim: -1}) if is_duck_dask_array(b.data) else b
        else:
            raise ValueError(
                f"{dim!r} on {'a' if a.sizes[dim] == 1 else 'b'} is incompatible:"
                " dimensions without coordinates must have have a length of 2 or 3"
            )

    c = apply_ufunc(
        np.cross,
        a,
        b,
        input_core_dims=[[dim], [dim]],
        output_core_dims=[[dim] if a.sizes[dim] == 3 else []],
        dask="parallelized",
        output_dtypes=[np.result_type(a, b)],
    )
    c = c.transpose(*all_dims, missing_dims="ignore")

    return c


@deprecate_dims
def dot(
    *arrays,
    dim: Dims = None,
    **kwargs: Any,
):
    """Generalized dot product for xarray objects. Like ``np.einsum``, but
    provides a simpler interface based on array dimension names.

    Parameters
    ----------
    *arrays : DataArray or Variable
        Arrays to compute.
    dim : str, iterable of hashable, "..." or None, optional
        Which dimensions to sum over. Ellipsis ('...') sums over all dimensions.
        If not specified, then all the common dimensions are summed over.
    **kwargs : dict
        Additional keyword arguments passed to ``numpy.einsum`` or
        ``dask.array.einsum``

    Returns
    -------
    DataArray

    See Also
    --------
    numpy.einsum
    dask.array.einsum
    opt_einsum.contract

    Notes
    -----
    We recommend installing the optional ``opt_einsum`` package, or alternatively passing ``optimize=True``,
    which is passed through to ``np.einsum``, and works for most array backends.

    Examples
    --------
    >>> da_a = xr.DataArray(np.arange(3 * 2).reshape(3, 2), dims=["a", "b"])
    >>> da_b = xr.DataArray(np.arange(3 * 2 * 2).reshape(3, 2, 2), dims=["a", "b", "c"])
    >>> da_c = xr.DataArray(np.arange(2 * 3).reshape(2, 3), dims=["c", "d"])

    >>> da_a
    <xarray.DataArray (a: 3, b: 2)> Size: 48B
    array([[0, 1],
           [2, 3],
           [4, 5]])
    Dimensions without coordinates: a, b

    >>> da_b
    <xarray.DataArray (a: 3, b: 2, c: 2)> Size: 96B
    array([[[ 0,  1],
            [ 2,  3]],
    <BLANKLINE>
           [[ 4,  5],
            [ 6,  7]],
    <BLANKLINE>
           [[ 8,  9],
            [10, 11]]])
    Dimensions without coordinates: a, b, c

    >>> da_c
    <xarray.DataArray (c: 2, d: 3)> Size: 48B
    array([[0, 1, 2],
           [3, 4, 5]])
    Dimensions without coordinates: c, d

    >>> xr.dot(da_a, da_b, dim=["a", "b"])
    <xarray.DataArray (c: 2)> Size: 16B
    array([110, 125])
    Dimensions without coordinates: c

    >>> xr.dot(da_a, da_b, dim=["a"])
    <xarray.DataArray (b: 2, c: 2)> Size: 32B
    array([[40, 46],
           [70, 79]])
    Dimensions without coordinates: b, c

    >>> xr.dot(da_a, da_b, da_c, dim=["b", "c"])
    <xarray.DataArray (a: 3, d: 3)> Size: 72B
    array([[  9,  14,  19],
           [ 93, 150, 207],
           [273, 446, 619]])
    Dimensions without coordinates: a, d

    >>> xr.dot(da_a, da_b)
    <xarray.DataArray (c: 2)> Size: 16B
    array([110, 125])
    Dimensions without coordinates: c

    >>> xr.dot(da_a, da_b, dim=...)
    <xarray.DataArray ()> Size: 8B
    array(235)
    """
    from xarray.core.dataarray import DataArray
    from xarray.core.variable import Variable

    if any(not isinstance(arr, (Variable, DataArray)) for arr in arrays):
        raise TypeError(
            "Only xr.DataArray and xr.Variable are supported."
            f"Given {[type(arr) for arr in arrays]}."
        )

    if len(arrays) == 0:
        raise TypeError("At least one array should be given.")

    common_dims: set[Hashable] = set.intersection(*(set(arr.dims) for arr in arrays))
    all_dims = []
    for arr in arrays:
        all_dims += [d for d in arr.dims if d not in all_dims]

    einsum_axes = "abcdefghijklmnopqrstuvwxyz"
    dim_map = {d: einsum_axes[i] for i, d in enumerate(all_dims)}

    if dim is None:
        # find dimensions that occur more than once
        dim_counts: Counter = Counter()
        for arr in arrays:
            dim_counts.update(arr.dims)
        dim = tuple(d for d, c in dim_counts.items() if c > 1)
    else:
        dim = parse_dims(dim, all_dims=tuple(all_dims))

    dot_dims: set[Hashable] = set(dim)

    # dimensions to be parallelized
    broadcast_dims = common_dims - dot_dims
    input_core_dims = [
        [d for d in arr.dims if d not in broadcast_dims] for arr in arrays
    ]
    output_core_dims = [
        [d for d in all_dims if d not in dot_dims and d not in broadcast_dims]
    ]

    # construct einsum subscripts, such as '...abc,...ab->...c'
    # Note: input_core_dims are always moved to the last position
    subscripts_list = [
        "..." + "".join(dim_map[d] for d in ds) for ds in input_core_dims
    ]
    subscripts = ",".join(subscripts_list)
    subscripts += "->..." + "".join(dim_map[d] for d in output_core_dims[0])

    join = OPTIONS["arithmetic_join"]
    # using "inner" emulates `(a * b).sum()` for all joins (except "exact")
    if join != "exact":
        join = "inner"

    # subscripts should be passed to np.einsum as arg, not as kwargs. We need
    # to construct a partial function for apply_ufunc to work.
    func = functools.partial(duck_array_ops.einsum, subscripts, **kwargs)
    result = apply_ufunc(
        func,
        *arrays,
        input_core_dims=input_core_dims,
        output_core_dims=output_core_dims,
        join=join,
        dask="allowed",
    )
    return result.transpose(*all_dims, missing_dims="ignore")


def where(cond, x, y, keep_attrs=None):
    """Return elements from `x` or `y` depending on `cond`.

    Performs xarray-like broadcasting across input arguments.

    All dimension coordinates on `x` and `y`  must be aligned with each
    other and with `cond`.

    Parameters
    ----------
    cond : scalar, array, Variable, DataArray or Dataset
        When True, return values from `x`, otherwise returns values from `y`.
    x : scalar, array, Variable, DataArray or Dataset
        values to choose from where `cond` is True
    y : scalar, array, Variable, DataArray or Dataset
        values to choose from where `cond` is False
    keep_attrs : bool or str or callable, optional
        How to treat attrs. If True, keep the attrs of `x`.

    Returns
    -------
    Dataset, DataArray, Variable or array
        In priority order: Dataset, DataArray, Variable or array, whichever
        type appears as an input argument.

    Examples
    --------
    >>> x = xr.DataArray(
    ...     0.1 * np.arange(10),
    ...     dims=["lat"],
    ...     coords={"lat": np.arange(10)},
    ...     name="sst",
    ... )
    >>> x
    <xarray.DataArray 'sst' (lat: 10)> Size: 80B
    array([0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    Coordinates:
      * lat      (lat) int64 80B 0 1 2 3 4 5 6 7 8 9

    >>> xr.where(x < 0.5, x, x * 100)
    <xarray.DataArray 'sst' (lat: 10)> Size: 80B
    array([ 0. ,  0.1,  0.2,  0.3,  0.4, 50. , 60. , 70. , 80. , 90. ])
    Coordinates:
      * lat      (lat) int64 80B 0 1 2 3 4 5 6 7 8 9

    >>> y = xr.DataArray(
    ...     0.1 * np.arange(9).reshape(3, 3),
    ...     dims=["lat", "lon"],
    ...     coords={"lat": np.arange(3), "lon": 10 + np.arange(3)},
    ...     name="sst",
    ... )
    >>> y
    <xarray.DataArray 'sst' (lat: 3, lon: 3)> Size: 72B
    array([[0. , 0.1, 0.2],
           [0.3, 0.4, 0.5],
           [0.6, 0.7, 0.8]])
    Coordinates:
      * lat      (lat) int64 24B 0 1 2
      * lon      (lon) int64 24B 10 11 12

    >>> xr.where(y.lat < 1, y, -1)
    <xarray.DataArray (lat: 3, lon: 3)> Size: 72B
    array([[ 0. ,  0.1,  0.2],
           [-1. , -1. , -1. ],
           [-1. , -1. , -1. ]])
    Coordinates:
      * lat      (lat) int64 24B 0 1 2
      * lon      (lon) int64 24B 10 11 12

    >>> cond = xr.DataArray([True, False], dims=["x"])
    >>> x = xr.DataArray([1, 2], dims=["y"])
    >>> xr.where(cond, x, 0)
    <xarray.DataArray (x: 2, y: 2)> Size: 32B
    array([[1, 2],
           [0, 0]])
    Dimensions without coordinates: x, y

    See Also
    --------
    numpy.where : corresponding numpy function
    Dataset.where, DataArray.where :
        equivalent methods
    """
    from xarray.core.dataset import Dataset

    if keep_attrs is None:
        keep_attrs = _get_keep_attrs(default=False)

    # alignment for three arguments is complicated, so don't support it yet
    result = apply_ufunc(
        duck_array_ops.where,
        cond,
        x,
        y,
        join="exact",
        dataset_join="exact",
        dask="allowed",
        keep_attrs=keep_attrs,
    )

    # keep the attributes of x, the second parameter, by default to
    # be consistent with the `where` method of `DataArray` and `Dataset`
    # rebuild the attrs from x at each level of the output, which could be
    # Dataset, DataArray, or Variable, and also handle coords
    if keep_attrs is True and hasattr(result, "attrs"):
        if isinstance(y, Dataset) and not isinstance(x, Dataset):
            # handle special case where x gets promoted to Dataset
            result.attrs = {}
            if getattr(x, "name", None) in result.data_vars:
                result[x.name].attrs = getattr(x, "attrs", {})
        else:
            # otherwise, fill in global attrs and variable attrs (if they exist)
            result.attrs = getattr(x, "attrs", {})
            for v in getattr(result, "data_vars", []):
                result[v].attrs = getattr(getattr(x, v, None), "attrs", {})
        for c in getattr(result, "coords", []):
            # always fill coord attrs of x
            result[c].attrs = getattr(getattr(x, c, None), "attrs", {})

    return result


@overload
def polyval(
    coord: DataArray, coeffs: DataArray, degree_dim: Hashable = "degree"
) -> DataArray: ...


@overload
def polyval(
    coord: DataArray, coeffs: Dataset, degree_dim: Hashable = "degree"
) -> Dataset: ...


@overload
def polyval(
    coord: Dataset, coeffs: DataArray, degree_dim: Hashable = "degree"
) -> Dataset: ...


@overload
def polyval(
    coord: Dataset, coeffs: Dataset, degree_dim: Hashable = "degree"
) -> Dataset: ...


@overload
def polyval(
    coord: Dataset | DataArray,
    coeffs: Dataset | DataArray,
    degree_dim: Hashable = "degree",
) -> Dataset | DataArray: ...


def polyval(
    coord: Dataset | DataArray,
    coeffs: Dataset | DataArray,
    degree_dim: Hashable = "degree",
) -> Dataset | DataArray:
    """Evaluate a polynomial at specific values

    Parameters
    ----------
    coord : DataArray or Dataset
        Values at which to evaluate the polynomial.
    coeffs : DataArray or Dataset
        Coefficients of the polynomial.
    degree_dim : Hashable, default: "degree"
        Name of the polynomial degree dimension in `coeffs`.

    Returns
    -------
    DataArray or Dataset
        Evaluated polynomial.

    See Also
    --------
    xarray.DataArray.polyfit
    numpy.polynomial.polynomial.polyval
    """

    if degree_dim not in coeffs._indexes:
        raise ValueError(
            f"Dimension `{degree_dim}` should be a coordinate variable with labels."
        )
    if not np.issubdtype(coeffs[degree_dim].dtype, np.integer):
        raise ValueError(
            f"Dimension `{degree_dim}` should be of integer dtype. Received {coeffs[degree_dim].dtype} instead."
        )
    max_deg = coeffs[degree_dim].max().item()
    coeffs = coeffs.reindex(
        {degree_dim: np.arange(max_deg + 1)}, fill_value=0, copy=False
    )
    coord = _ensure_numeric(coord)

    # using Horner's method
    # https://en.wikipedia.org/wiki/Horner%27s_method
    res = zeros_like(coord) + coeffs.isel({degree_dim: max_deg}, drop=True)
    for deg in range(max_deg - 1, -1, -1):
        res *= coord
        res += coeffs.isel({degree_dim: deg}, drop=True)

    return res


def _ensure_numeric(data: Dataset | DataArray) -> Dataset | DataArray:
    """Converts all datetime64 variables to float64

    Parameters
    ----------
    data : DataArray or Dataset
        Variables with possible datetime dtypes.

    Returns
    -------
    DataArray or Dataset
        Variables with datetime64 dtypes converted to float64.
    """
    from xarray.core.dataset import Dataset

    def _cfoffset(x: DataArray) -> Any:
        scalar = x.compute().data[0]
        if not is_scalar(scalar):
            # we do not get a scalar back on dask == 2021.04.1
            scalar = scalar.item()
        return type(scalar)(1970, 1, 1)

    def to_floatable(x: DataArray) -> DataArray:
        if x.dtype.kind in "MO":
            # datetimes (CFIndexes are object type)
            offset = (
                np.datetime64("1970-01-01") if x.dtype.kind == "M" else _cfoffset(x)
            )
            return x.copy(
                data=datetime_to_numeric(x.data, offset=offset, datetime_unit="ns"),
            )
        elif x.dtype.kind == "m":
            # timedeltas
            return duck_array_ops.astype(x, dtype=float)
        return x

    if isinstance(data, Dataset):
        return data.map(to_floatable)
    else:
        return to_floatable(data)


def _calc_idxminmax(
    *,
    array,
    func: Callable,
    dim: Hashable | None = None,
    skipna: bool | None = None,
    fill_value: Any = dtypes.NA,
    keep_attrs: bool | None = None,
):
    """Apply common operations for idxmin and idxmax."""
    # This function doesn't make sense for scalars so don't try
    if not array.ndim:
        raise ValueError("This function does not apply for scalars")

    if dim is not None:
        pass  # Use the dim if available
    elif array.ndim == 1:
        # it is okay to guess the dim if there is only 1
        dim = array.dims[0]
    else:
        # The dim is not specified and ambiguous.  Don't guess.
        raise ValueError("Must supply 'dim' argument for multidimensional arrays")

    if dim not in array.dims:
        raise KeyError(
            f"Dimension {dim!r} not found in array dimensions {array.dims!r}"
        )
    if dim not in array.coords:
        raise KeyError(
            f"Dimension {dim!r} is not one of the coordinates {tuple(array.coords.keys())}"
        )

    # These are dtypes with NaN values argmin and argmax can handle
    na_dtypes = "cfO"

    if skipna or (skipna is None and array.dtype.kind in na_dtypes):
        # Need to skip NaN values since argmin and argmax can't handle them
        allna = array.isnull().all(dim)
        array = array.where(~allna, 0)

    # This will run argmin or argmax.
    indx = func(array, dim=dim, axis=None, keep_attrs=keep_attrs, skipna=skipna)

    # Handle chunked arrays (e.g. dask).
    if is_chunked_array(array.data):
        chunkmanager = get_chunked_array_type(array.data)
        chunks = dict(zip(array.dims, array.chunks))
        dask_coord = chunkmanager.from_array(array[dim].data, chunks=chunks[dim])
        data = dask_coord[duck_array_ops.ravel(indx.data)]
        res = indx.copy(data=duck_array_ops.reshape(data, indx.shape))
        # we need to attach back the dim name
        res.name = dim
    else:
        res = array[dim][(indx,)]
        # The dim is gone but we need to remove the corresponding coordinate.
        del res.coords[dim]

    if skipna or (skipna is None and array.dtype.kind in na_dtypes):
        # Put the NaN values back in after removing them
        res = res.where(~allna, fill_value)

    # Copy attributes from argmin/argmax, if any
    res.attrs = indx.attrs

    return res


_T = TypeVar("_T", bound=Union["Dataset", "DataArray"])
_U = TypeVar("_U", bound=Union["Dataset", "DataArray"])
_V = TypeVar("_V", bound=Union["Dataset", "DataArray"])


@overload
def unify_chunks(__obj: _T) -> tuple[_T]: ...


@overload
def unify_chunks(__obj1: _T, __obj2: _U) -> tuple[_T, _U]: ...


@overload
def unify_chunks(__obj1: _T, __obj2: _U, __obj3: _V) -> tuple[_T, _U, _V]: ...


@overload
def unify_chunks(*objects: Dataset | DataArray) -> tuple[Dataset | DataArray, ...]: ...


def unify_chunks(*objects: Dataset | DataArray) -> tuple[Dataset | DataArray, ...]:
    """
    Given any number of Dataset and/or DataArray objects, returns
    new objects with unified chunk size along all chunked dimensions.

    Returns
    -------
    unified (DataArray or Dataset) – Tuple of objects with the same type as
    *objects with consistent chunk sizes for all dask-array variables

    See Also
    --------
    dask.array.core.unify_chunks
    """
    from xarray.core.dataarray import DataArray

    # Convert all objects to datasets
    datasets = [
        obj._to_temp_dataset() if isinstance(obj, DataArray) else obj.copy()
        for obj in objects
    ]

    # Get arguments to pass into dask.array.core.unify_chunks
    unify_chunks_args = []
    sizes: dict[Hashable, int] = {}
    for ds in datasets:
        for v in ds._variables.values():
            if v.chunks is not None:
                # Check that sizes match across different datasets
                for dim, size in v.sizes.items():
                    try:
                        if sizes[dim] != size:
                            raise ValueError(
                                f"Dimension {dim!r} size mismatch: {sizes[dim]} != {size}"
                            )
                    except KeyError:
                        sizes[dim] = size
                unify_chunks_args += [v._data, v._dims]

    # No dask arrays: Return inputs
    if not unify_chunks_args:
        return objects

    chunkmanager = get_chunked_array_type(*[arg for arg in unify_chunks_args])
    _, chunked_data = chunkmanager.unify_chunks(*unify_chunks_args)
    chunked_data_iter = iter(chunked_data)
    out: list[Dataset | DataArray] = []
    for obj, ds in zip(objects, datasets):
        for k, v in ds._variables.items():
            if v.chunks is not None:
                ds._variables[k] = v.copy(data=next(chunked_data_iter))
        out.append(obj._from_temp_dataset(ds) if isinstance(obj, DataArray) else ds)

    return tuple(out)
