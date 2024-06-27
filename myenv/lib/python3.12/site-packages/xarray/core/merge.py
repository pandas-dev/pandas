from __future__ import annotations

from collections import defaultdict
from collections.abc import Hashable, Iterable, Mapping, Sequence, Set
from typing import TYPE_CHECKING, Any, NamedTuple, Optional, Union

import pandas as pd

from xarray.core import dtypes
from xarray.core.alignment import deep_align
from xarray.core.duck_array_ops import lazy_array_equiv
from xarray.core.indexes import (
    Index,
    create_default_index_implicit,
    filter_indexes_from_coords,
    indexes_equal,
)
from xarray.core.utils import Frozen, compat_dict_union, dict_equiv, equivalent
from xarray.core.variable import Variable, as_variable, calculate_dimensions

if TYPE_CHECKING:
    from xarray.core.coordinates import Coordinates
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import CombineAttrsOptions, CompatOptions, JoinOptions

    DimsLike = Union[Hashable, Sequence[Hashable]]
    ArrayLike = Any
    VariableLike = Union[
        ArrayLike,
        tuple[DimsLike, ArrayLike],
        tuple[DimsLike, ArrayLike, Mapping],
        tuple[DimsLike, ArrayLike, Mapping, Mapping],
    ]
    XarrayValue = Union[DataArray, Variable, VariableLike]
    DatasetLike = Union[Dataset, Coordinates, Mapping[Any, XarrayValue]]
    CoercibleValue = Union[XarrayValue, pd.Series, pd.DataFrame]
    CoercibleMapping = Union[Dataset, Mapping[Any, CoercibleValue]]


PANDAS_TYPES = (pd.Series, pd.DataFrame)

_VALID_COMPAT = Frozen(
    {
        "identical": 0,
        "equals": 1,
        "broadcast_equals": 2,
        "minimal": 3,
        "no_conflicts": 4,
        "override": 5,
    }
)


class Context:
    """object carrying the information of a call"""

    def __init__(self, func):
        self.func = func


def broadcast_dimension_size(variables: list[Variable]) -> dict[Hashable, int]:
    """Extract dimension sizes from a dictionary of variables.

    Raises ValueError if any dimensions have different sizes.
    """
    dims: dict[Hashable, int] = {}
    for var in variables:
        for dim, size in zip(var.dims, var.shape):
            if dim in dims and size != dims[dim]:
                raise ValueError(f"index {dim!r} not aligned")
            dims[dim] = size
    return dims


class MergeError(ValueError):
    """Error class for merge failures due to incompatible arguments."""

    # inherits from ValueError for backward compatibility
    # TODO: move this to an xarray.exceptions module?


def unique_variable(
    name: Hashable,
    variables: list[Variable],
    compat: CompatOptions = "broadcast_equals",
    equals: bool | None = None,
) -> Variable:
    """Return the unique variable from a list of variables or raise MergeError.

    Parameters
    ----------
    name : hashable
        Name for this variable.
    variables : list of Variable
        List of Variable objects, all of which go by the same name in different
        inputs.
    compat : {"identical", "equals", "broadcast_equals", "no_conflicts", "override"}, optional
        Type of equality check to use.
    equals : None or bool, optional
        corresponding to result of compat test

    Returns
    -------
    Variable to use in the result.

    Raises
    ------
    MergeError: if any of the variables are not equal.
    """
    out = variables[0]

    if len(variables) == 1 or compat == "override":
        return out

    combine_method = None

    if compat == "minimal":
        compat = "broadcast_equals"

    if compat == "broadcast_equals":
        dim_lengths = broadcast_dimension_size(variables)
        out = out.set_dims(dim_lengths)

    if compat == "no_conflicts":
        combine_method = "fillna"

    if equals is None:
        # first check without comparing values i.e. no computes
        for var in variables[1:]:
            equals = getattr(out, compat)(var, equiv=lazy_array_equiv)
            if equals is not True:
                break

        if equals is None:
            # now compare values with minimum number of computes
            out = out.compute()
            for var in variables[1:]:
                equals = getattr(out, compat)(var)
                if not equals:
                    break

    if not equals:
        raise MergeError(
            f"conflicting values for variable {name!r} on objects to be combined. "
            "You can skip this check by specifying compat='override'."
        )

    if combine_method:
        for var in variables[1:]:
            out = getattr(out, combine_method)(var)

    return out


def _assert_compat_valid(compat):
    if compat not in _VALID_COMPAT:
        raise ValueError(f"compat={compat!r} invalid: must be {set(_VALID_COMPAT)}")


MergeElement = tuple[Variable, Optional[Index]]


def _assert_prioritized_valid(
    grouped: dict[Hashable, list[MergeElement]],
    prioritized: Mapping[Any, MergeElement],
) -> None:
    """Make sure that elements given in prioritized will not corrupt any
    index given in grouped.
    """
    prioritized_names = set(prioritized)
    grouped_by_index: dict[int, list[Hashable]] = defaultdict(list)
    indexes: dict[int, Index] = {}

    for name, elements_list in grouped.items():
        for _, index in elements_list:
            if index is not None:
                grouped_by_index[id(index)].append(name)
                indexes[id(index)] = index

    # An index may be corrupted when the set of its corresponding coordinate name(s)
    # partially overlaps the set of names given in prioritized
    for index_id, index_coord_names in grouped_by_index.items():
        index_names = set(index_coord_names)
        common_names = index_names & prioritized_names
        if common_names and len(common_names) != len(index_names):
            common_names_str = ", ".join(f"{k!r}" for k in common_names)
            index_names_str = ", ".join(f"{k!r}" for k in index_coord_names)
            raise ValueError(
                f"cannot set or update variable(s) {common_names_str}, which would corrupt "
                f"the following index built from coordinates {index_names_str}:\n"
                f"{indexes[index_id]!r}"
            )


def merge_collected(
    grouped: dict[Any, list[MergeElement]],
    prioritized: Mapping[Any, MergeElement] | None = None,
    compat: CompatOptions = "minimal",
    combine_attrs: CombineAttrsOptions = "override",
    equals: dict[Any, bool] | None = None,
) -> tuple[dict[Hashable, Variable], dict[Hashable, Index]]:
    """Merge dicts of variables, while resolving conflicts appropriately.

    Parameters
    ----------
    grouped : mapping
    prioritized : mapping
    compat : str
        Type of equality check to use when checking for conflicts.
    combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts", \
                    "override"} or callable, default: "override"
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
    equals : mapping, optional
        corresponding to result of compat test

    Returns
    -------
    Dict with keys taken by the union of keys on list_of_mappings,
    and Variable values corresponding to those that should be found on the
    merged result.
    """
    if prioritized is None:
        prioritized = {}
    if equals is None:
        equals = {}

    _assert_compat_valid(compat)
    _assert_prioritized_valid(grouped, prioritized)

    merged_vars: dict[Hashable, Variable] = {}
    merged_indexes: dict[Hashable, Index] = {}
    index_cmp_cache: dict[tuple[int, int], bool | None] = {}

    for name, elements_list in grouped.items():
        if name in prioritized:
            variable, index = prioritized[name]
            merged_vars[name] = variable
            if index is not None:
                merged_indexes[name] = index
        else:
            indexed_elements = [
                (variable, index)
                for variable, index in elements_list
                if index is not None
            ]
            if indexed_elements:
                # TODO(shoyer): consider adjusting this logic. Are we really
                # OK throwing away variable without an index in favor of
                # indexed variables, without even checking if values match?
                variable, index = indexed_elements[0]
                for other_var, other_index in indexed_elements[1:]:
                    if not indexes_equal(
                        index, other_index, variable, other_var, index_cmp_cache
                    ):
                        raise MergeError(
                            f"conflicting values/indexes on objects to be combined fo coordinate {name!r}\n"
                            f"first index: {index!r}\nsecond index: {other_index!r}\n"
                            f"first variable: {variable!r}\nsecond variable: {other_var!r}\n"
                        )
                if compat == "identical":
                    for other_variable, _ in indexed_elements[1:]:
                        if not dict_equiv(variable.attrs, other_variable.attrs):
                            raise MergeError(
                                "conflicting attribute values on combined "
                                f"variable {name!r}:\nfirst value: {variable.attrs!r}\nsecond value: {other_variable.attrs!r}"
                            )
                merged_vars[name] = variable
                merged_vars[name].attrs = merge_attrs(
                    [var.attrs for var, _ in indexed_elements],
                    combine_attrs=combine_attrs,
                )
                merged_indexes[name] = index
            else:
                variables = [variable for variable, _ in elements_list]
                try:
                    merged_vars[name] = unique_variable(
                        name, variables, compat, equals.get(name, None)
                    )
                except MergeError:
                    if compat != "minimal":
                        # we need more than "minimal" compatibility (for which
                        # we drop conflicting coordinates)
                        raise

                if name in merged_vars:
                    merged_vars[name].attrs = merge_attrs(
                        [var.attrs for var in variables], combine_attrs=combine_attrs
                    )

    return merged_vars, merged_indexes


def collect_variables_and_indexes(
    list_of_mappings: Iterable[DatasetLike],
    indexes: Mapping[Any, Any] | None = None,
) -> dict[Hashable, list[MergeElement]]:
    """Collect variables and indexes from list of mappings of xarray objects.

    Mappings can be Dataset or Coordinates objects, in which case both
    variables and indexes are extracted from it.

    It can also have values of one of the following types:
    - an xarray.Variable
    - a tuple `(dims, data[, attrs[, encoding]])` that can be converted in
      an xarray.Variable
    - or an xarray.DataArray

    If a mapping of indexes is given, those indexes are assigned to all variables
    with a matching key/name. For dimension variables with no matching index, a
    default (pandas) index is assigned. DataArray indexes that don't match mapping
    keys are also extracted.

    """
    from xarray.core.coordinates import Coordinates
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    if indexes is None:
        indexes = {}

    grouped: dict[Hashable, list[MergeElement]] = defaultdict(list)

    def append(name, variable, index):
        grouped[name].append((variable, index))

    def append_all(variables, indexes):
        for name, variable in variables.items():
            append(name, variable, indexes.get(name))

    for mapping in list_of_mappings:
        if isinstance(mapping, (Coordinates, Dataset)):
            append_all(mapping.variables, mapping.xindexes)
            continue

        for name, variable in mapping.items():
            if isinstance(variable, DataArray):
                coords_ = variable._coords.copy()  # use private API for speed
                indexes_ = dict(variable._indexes)
                # explicitly overwritten variables should take precedence
                coords_.pop(name, None)
                indexes_.pop(name, None)
                append_all(coords_, indexes_)

            variable = as_variable(variable, name=name, auto_convert=False)
            if name in indexes:
                append(name, variable, indexes[name])
            elif variable.dims == (name,):
                idx, idx_vars = create_default_index_implicit(variable)
                append_all(idx_vars, {k: idx for k in idx_vars})
            else:
                append(name, variable, None)

    return grouped


def collect_from_coordinates(
    list_of_coords: list[Coordinates],
) -> dict[Hashable, list[MergeElement]]:
    """Collect variables and indexes to be merged from Coordinate objects."""
    grouped: dict[Hashable, list[MergeElement]] = defaultdict(list)

    for coords in list_of_coords:
        variables = coords.variables
        indexes = coords.xindexes
        for name, variable in variables.items():
            grouped[name].append((variable, indexes.get(name)))

    return grouped


def merge_coordinates_without_align(
    objects: list[Coordinates],
    prioritized: Mapping[Any, MergeElement] | None = None,
    exclude_dims: Set = frozenset(),
    combine_attrs: CombineAttrsOptions = "override",
) -> tuple[dict[Hashable, Variable], dict[Hashable, Index]]:
    """Merge variables/indexes from coordinates without automatic alignments.

    This function is used for merging coordinate from pre-existing xarray
    objects.
    """
    collected = collect_from_coordinates(objects)

    if exclude_dims:
        filtered: dict[Hashable, list[MergeElement]] = {}
        for name, elements in collected.items():
            new_elements = [
                (variable, index)
                for variable, index in elements
                if exclude_dims.isdisjoint(variable.dims)
            ]
            if new_elements:
                filtered[name] = new_elements
    else:
        filtered = collected

    # TODO: indexes should probably be filtered in collected elements
    # before merging them
    merged_coords, merged_indexes = merge_collected(
        filtered, prioritized, combine_attrs=combine_attrs
    )
    merged_indexes = filter_indexes_from_coords(merged_indexes, set(merged_coords))

    return merged_coords, merged_indexes


def determine_coords(
    list_of_mappings: Iterable[DatasetLike],
) -> tuple[set[Hashable], set[Hashable]]:
    """Given a list of dicts with xarray object values, identify coordinates.

    Parameters
    ----------
    list_of_mappings : list of dict or list of Dataset
        Of the same form as the arguments to expand_variable_dicts.

    Returns
    -------
    coord_names : set of variable names
    noncoord_names : set of variable names
        All variable found in the input should appear in either the set of
        coordinate or non-coordinate names.
    """
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    coord_names: set[Hashable] = set()
    noncoord_names: set[Hashable] = set()

    for mapping in list_of_mappings:
        if isinstance(mapping, Dataset):
            coord_names.update(mapping.coords)
            noncoord_names.update(mapping.data_vars)
        else:
            for name, var in mapping.items():
                if isinstance(var, DataArray):
                    coords = set(var._coords)  # use private API for speed
                    # explicitly overwritten variables should take precedence
                    coords.discard(name)
                    coord_names.update(coords)

    return coord_names, noncoord_names


def coerce_pandas_values(objects: Iterable[CoercibleMapping]) -> list[DatasetLike]:
    """Convert pandas values found in a list of labeled objects.

    Parameters
    ----------
    objects : list of Dataset or mapping
        The mappings may contain any sort of objects coercible to
        xarray.Variables as keys, including pandas objects.

    Returns
    -------
    List of Dataset or dictionary objects. Any inputs or values in the inputs
    that were pandas objects have been converted into native xarray objects.
    """
    from xarray.core.coordinates import Coordinates
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    out: list[DatasetLike] = []
    for obj in objects:
        variables: DatasetLike
        if isinstance(obj, (Dataset, Coordinates)):
            variables = obj
        else:
            variables = {}
            if isinstance(obj, PANDAS_TYPES):
                obj = dict(obj.items())
            for k, v in obj.items():
                if isinstance(v, PANDAS_TYPES):
                    v = DataArray(v)
                variables[k] = v
        out.append(variables)
    return out


def _get_priority_vars_and_indexes(
    objects: Sequence[DatasetLike],
    priority_arg: int | None,
    compat: CompatOptions = "equals",
) -> dict[Hashable, MergeElement]:
    """Extract the priority variable from a list of mappings.

    We need this method because in some cases the priority argument itself
    might have conflicting values (e.g., if it is a dict with two DataArray
    values with conflicting coordinate values).

    Parameters
    ----------
    objects : sequence of dict-like of Variable
        Dictionaries in which to find the priority variables.
    priority_arg : int or None
        Integer object whose variable should take priority.
    compat : {"identical", "equals", "broadcast_equals", "no_conflicts", "override"}, optional
        String indicating how to compare non-concatenated variables of the same name for
        potential conflicts. This is passed down to merge.

        - "broadcast_equals": all values must be equal when variables are
          broadcast against each other to ensure common dimensions.
        - "equals": all values and dimensions must be the same.
        - "identical": all values, dimensions and attributes must be the
          same.
        - "no_conflicts": only values which are not null in both datasets
          must be equal. The returned dataset then contains the combination
          of all non-null values.
        - "override": skip comparing and pick variable from first dataset

    Returns
    -------
    A dictionary of variables and associated indexes (if any) to prioritize.
    """
    if priority_arg is None:
        return {}

    collected = collect_variables_and_indexes([objects[priority_arg]])
    variables, indexes = merge_collected(collected, compat=compat)
    grouped: dict[Hashable, MergeElement] = {}
    for name, variable in variables.items():
        grouped[name] = (variable, indexes.get(name))
    return grouped


def merge_coords(
    objects: Iterable[CoercibleMapping],
    compat: CompatOptions = "minimal",
    join: JoinOptions = "outer",
    priority_arg: int | None = None,
    indexes: Mapping[Any, Index] | None = None,
    fill_value: object = dtypes.NA,
) -> tuple[dict[Hashable, Variable], dict[Hashable, Index]]:
    """Merge coordinate variables.

    See merge_core below for argument descriptions. This works similarly to
    merge_core, except everything we don't worry about whether variables are
    coordinates or not.
    """
    _assert_compat_valid(compat)
    coerced = coerce_pandas_values(objects)
    aligned = deep_align(
        coerced, join=join, copy=False, indexes=indexes, fill_value=fill_value
    )
    collected = collect_variables_and_indexes(aligned, indexes=indexes)
    prioritized = _get_priority_vars_and_indexes(aligned, priority_arg, compat=compat)
    variables, out_indexes = merge_collected(collected, prioritized, compat=compat)
    return variables, out_indexes


def merge_attrs(variable_attrs, combine_attrs, context=None):
    """Combine attributes from different variables according to combine_attrs"""
    if not variable_attrs:
        # no attributes to merge
        return None

    if callable(combine_attrs):
        return combine_attrs(variable_attrs, context=context)
    elif combine_attrs == "drop":
        return {}
    elif combine_attrs == "override":
        return dict(variable_attrs[0])
    elif combine_attrs == "no_conflicts":
        result = dict(variable_attrs[0])
        for attrs in variable_attrs[1:]:
            try:
                result = compat_dict_union(result, attrs)
            except ValueError as e:
                raise MergeError(
                    "combine_attrs='no_conflicts', but some values are not "
                    f"the same. Merging {str(result)} with {str(attrs)}"
                ) from e
        return result
    elif combine_attrs == "drop_conflicts":
        result = {}
        dropped_keys = set()
        for attrs in variable_attrs:
            result.update(
                {
                    key: value
                    for key, value in attrs.items()
                    if key not in result and key not in dropped_keys
                }
            )
            result = {
                key: value
                for key, value in result.items()
                if key not in attrs or equivalent(attrs[key], value)
            }
            dropped_keys |= {key for key in attrs if key not in result}
        return result
    elif combine_attrs == "identical":
        result = dict(variable_attrs[0])
        for attrs in variable_attrs[1:]:
            if not dict_equiv(result, attrs):
                raise MergeError(
                    f"combine_attrs='identical', but attrs differ. First is {str(result)} "
                    f", other is {str(attrs)}."
                )
        return result
    else:
        raise ValueError(f"Unrecognised value for combine_attrs={combine_attrs}")


class _MergeResult(NamedTuple):
    variables: dict[Hashable, Variable]
    coord_names: set[Hashable]
    dims: dict[Hashable, int]
    indexes: dict[Hashable, Index]
    attrs: dict[Hashable, Any]


def merge_core(
    objects: Iterable[CoercibleMapping],
    compat: CompatOptions = "broadcast_equals",
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "override",
    priority_arg: int | None = None,
    explicit_coords: Iterable[Hashable] | None = None,
    indexes: Mapping[Any, Any] | None = None,
    fill_value: object = dtypes.NA,
    skip_align_args: list[int] | None = None,
) -> _MergeResult:
    """Core logic for merging labeled objects.

    This is not public API.

    Parameters
    ----------
    objects : list of mapping
        All values must be convertible to labeled arrays.
    compat : {"identical", "equals", "broadcast_equals", "no_conflicts", "override"}, optional
        Compatibility checks to use when merging variables.
    join : {"outer", "inner", "left", "right"}, optional
        How to combine objects with different indexes.
    combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts", \
                     "override"} or callable, default: "override"
        How to combine attributes of objects
    priority_arg : int, optional
        Optional argument in `objects` that takes precedence over the others.
    explicit_coords : set, optional
        An explicit list of variables from `objects` that are coordinates.
    indexes : dict, optional
        Dictionary with values given by xarray.Index objects or anything that
        may be cast to pandas.Index objects.
    fill_value : scalar, optional
        Value to use for newly missing values
    skip_align_args : list of int, optional
        Optional arguments in `objects` that are not included in alignment.

    Returns
    -------
    variables : dict
        Dictionary of Variable objects.
    coord_names : set
        Set of coordinate names.
    dims : dict
        Dictionary mapping from dimension names to sizes.
    attrs : dict
        Dictionary of attributes

    Raises
    ------
    MergeError if the merge cannot be done successfully.
    """
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    _assert_compat_valid(compat)

    objects = list(objects)
    if skip_align_args is None:
        skip_align_args = []

    skip_align_objs = [(pos, objects.pop(pos)) for pos in skip_align_args]

    coerced = coerce_pandas_values(objects)
    aligned = deep_align(
        coerced, join=join, copy=False, indexes=indexes, fill_value=fill_value
    )

    for pos, obj in skip_align_objs:
        aligned.insert(pos, obj)

    collected = collect_variables_and_indexes(aligned, indexes=indexes)
    prioritized = _get_priority_vars_and_indexes(aligned, priority_arg, compat=compat)
    variables, out_indexes = merge_collected(
        collected, prioritized, compat=compat, combine_attrs=combine_attrs
    )

    dims = calculate_dimensions(variables)

    coord_names, noncoord_names = determine_coords(coerced)
    if compat == "minimal":
        # coordinates may be dropped in merged results
        coord_names.intersection_update(variables)
    if explicit_coords is not None:
        coord_names.update(explicit_coords)
    for dim, size in dims.items():
        if dim in variables:
            coord_names.add(dim)
    ambiguous_coords = coord_names.intersection(noncoord_names)
    if ambiguous_coords:
        raise MergeError(
            "unable to determine if these variables should be "
            f"coordinates or not in the merged result: {ambiguous_coords}"
        )

    attrs = merge_attrs(
        [var.attrs for var in coerced if isinstance(var, (Dataset, DataArray))],
        combine_attrs,
    )

    return _MergeResult(variables, coord_names, dims, out_indexes, attrs)


def merge(
    objects: Iterable[DataArray | CoercibleMapping],
    compat: CompatOptions = "no_conflicts",
    join: JoinOptions = "outer",
    fill_value: object = dtypes.NA,
    combine_attrs: CombineAttrsOptions = "override",
) -> Dataset:
    """Merge any number of xarray objects into a single Dataset as variables.

    Parameters
    ----------
    objects : iterable of Dataset or iterable of DataArray or iterable of dict-like
        Merge together all variables from these objects. If any of them are
        DataArray objects, they must have a name.
    compat : {"identical", "equals", "broadcast_equals", "no_conflicts", \
              "override", "minimal"}, default: "no_conflicts"
        String indicating how to compare variables of the same name for
        potential conflicts:

        - "identical": all values, dimensions and attributes must be the
          same.
        - "equals": all values and dimensions must be the same.
        - "broadcast_equals": all values must be equal when variables are
          broadcast against each other to ensure common dimensions.
        - "no_conflicts": only values which are not null in both datasets
          must be equal. The returned dataset then contains the combination
          of all non-null values.
        - "override": skip comparing and pick variable from first dataset
        - "minimal": drop conflicting coordinates

    join : {"outer", "inner", "left", "right", "exact", "override"}, default: "outer"
        String indicating how to combine differing indexes in objects.

        - "outer": use the union of object indexes
        - "inner": use the intersection of object indexes
        - "left": use indexes from the first object with each dimension
        - "right": use indexes from the last object with each dimension
        - "exact": instead of aligning, raise `ValueError` when indexes to be
          aligned are not equal
        - "override": if indexes are of same size, rewrite indexes to be
          those of the first object with that dimension. Indexes for the same
          dimension must have the same size in all objects.

    fill_value : scalar or dict-like, optional
        Value to use for newly missing values. If a dict-like, maps
        variable names to fill values. Use a data array's name to
        refer to its values.
    combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts", \
                     "override"} or callable, default: "override"
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
    Dataset
        Dataset with combined variables from each object.

    Examples
    --------
    >>> x = xr.DataArray(
    ...     [[1.0, 2.0], [3.0, 5.0]],
    ...     dims=("lat", "lon"),
    ...     coords={"lat": [35.0, 40.0], "lon": [100.0, 120.0]},
    ...     name="var1",
    ... )
    >>> y = xr.DataArray(
    ...     [[5.0, 6.0], [7.0, 8.0]],
    ...     dims=("lat", "lon"),
    ...     coords={"lat": [35.0, 42.0], "lon": [100.0, 150.0]},
    ...     name="var2",
    ... )
    >>> z = xr.DataArray(
    ...     [[0.0, 3.0], [4.0, 9.0]],
    ...     dims=("time", "lon"),
    ...     coords={"time": [30.0, 60.0], "lon": [100.0, 150.0]},
    ...     name="var3",
    ... )

    >>> x
    <xarray.DataArray 'var1' (lat: 2, lon: 2)> Size: 32B
    array([[1., 2.],
           [3., 5.]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 40.0
      * lon      (lon) float64 16B 100.0 120.0

    >>> y
    <xarray.DataArray 'var2' (lat: 2, lon: 2)> Size: 32B
    array([[5., 6.],
           [7., 8.]])
    Coordinates:
      * lat      (lat) float64 16B 35.0 42.0
      * lon      (lon) float64 16B 100.0 150.0

    >>> z
    <xarray.DataArray 'var3' (time: 2, lon: 2)> Size: 32B
    array([[0., 3.],
           [4., 9.]])
    Coordinates:
      * time     (time) float64 16B 30.0 60.0
      * lon      (lon) float64 16B 100.0 150.0

    >>> xr.merge([x, y, z])
    <xarray.Dataset> Size: 256B
    Dimensions:  (lat: 3, lon: 3, time: 2)
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 24B 100.0 120.0 150.0
      * time     (time) float64 16B 30.0 60.0
    Data variables:
        var1     (lat, lon) float64 72B 1.0 2.0 nan 3.0 5.0 nan nan nan nan
        var2     (lat, lon) float64 72B 5.0 nan 6.0 nan nan nan 7.0 nan 8.0
        var3     (time, lon) float64 48B 0.0 nan 3.0 4.0 nan 9.0

    >>> xr.merge([x, y, z], compat="identical")
    <xarray.Dataset> Size: 256B
    Dimensions:  (lat: 3, lon: 3, time: 2)
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 24B 100.0 120.0 150.0
      * time     (time) float64 16B 30.0 60.0
    Data variables:
        var1     (lat, lon) float64 72B 1.0 2.0 nan 3.0 5.0 nan nan nan nan
        var2     (lat, lon) float64 72B 5.0 nan 6.0 nan nan nan 7.0 nan 8.0
        var3     (time, lon) float64 48B 0.0 nan 3.0 4.0 nan 9.0

    >>> xr.merge([x, y, z], compat="equals")
    <xarray.Dataset> Size: 256B
    Dimensions:  (lat: 3, lon: 3, time: 2)
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 24B 100.0 120.0 150.0
      * time     (time) float64 16B 30.0 60.0
    Data variables:
        var1     (lat, lon) float64 72B 1.0 2.0 nan 3.0 5.0 nan nan nan nan
        var2     (lat, lon) float64 72B 5.0 nan 6.0 nan nan nan 7.0 nan 8.0
        var3     (time, lon) float64 48B 0.0 nan 3.0 4.0 nan 9.0

    >>> xr.merge([x, y, z], compat="equals", fill_value=-999.0)
    <xarray.Dataset> Size: 256B
    Dimensions:  (lat: 3, lon: 3, time: 2)
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 24B 100.0 120.0 150.0
      * time     (time) float64 16B 30.0 60.0
    Data variables:
        var1     (lat, lon) float64 72B 1.0 2.0 -999.0 3.0 ... -999.0 -999.0 -999.0
        var2     (lat, lon) float64 72B 5.0 -999.0 6.0 -999.0 ... 7.0 -999.0 8.0
        var3     (time, lon) float64 48B 0.0 -999.0 3.0 4.0 -999.0 9.0

    >>> xr.merge([x, y, z], join="override")
    <xarray.Dataset> Size: 144B
    Dimensions:  (lat: 2, lon: 2, time: 2)
    Coordinates:
      * lat      (lat) float64 16B 35.0 40.0
      * lon      (lon) float64 16B 100.0 120.0
      * time     (time) float64 16B 30.0 60.0
    Data variables:
        var1     (lat, lon) float64 32B 1.0 2.0 3.0 5.0
        var2     (lat, lon) float64 32B 5.0 6.0 7.0 8.0
        var3     (time, lon) float64 32B 0.0 3.0 4.0 9.0

    >>> xr.merge([x, y, z], join="inner")
    <xarray.Dataset> Size: 64B
    Dimensions:  (lat: 1, lon: 1, time: 2)
    Coordinates:
      * lat      (lat) float64 8B 35.0
      * lon      (lon) float64 8B 100.0
      * time     (time) float64 16B 30.0 60.0
    Data variables:
        var1     (lat, lon) float64 8B 1.0
        var2     (lat, lon) float64 8B 5.0
        var3     (time, lon) float64 16B 0.0 4.0

    >>> xr.merge([x, y, z], compat="identical", join="inner")
    <xarray.Dataset> Size: 64B
    Dimensions:  (lat: 1, lon: 1, time: 2)
    Coordinates:
      * lat      (lat) float64 8B 35.0
      * lon      (lon) float64 8B 100.0
      * time     (time) float64 16B 30.0 60.0
    Data variables:
        var1     (lat, lon) float64 8B 1.0
        var2     (lat, lon) float64 8B 5.0
        var3     (time, lon) float64 16B 0.0 4.0

    >>> xr.merge([x, y, z], compat="broadcast_equals", join="outer")
    <xarray.Dataset> Size: 256B
    Dimensions:  (lat: 3, lon: 3, time: 2)
    Coordinates:
      * lat      (lat) float64 24B 35.0 40.0 42.0
      * lon      (lon) float64 24B 100.0 120.0 150.0
      * time     (time) float64 16B 30.0 60.0
    Data variables:
        var1     (lat, lon) float64 72B 1.0 2.0 nan 3.0 5.0 nan nan nan nan
        var2     (lat, lon) float64 72B 5.0 nan 6.0 nan nan nan 7.0 nan 8.0
        var3     (time, lon) float64 48B 0.0 nan 3.0 4.0 nan 9.0

    >>> xr.merge([x, y, z], join="exact")
    Traceback (most recent call last):
    ...
    ValueError: cannot align objects with join='exact' where ...

    Raises
    ------
    xarray.MergeError
        If any variables with the same name have conflicting values.

    See also
    --------
    concat
    combine_nested
    combine_by_coords
    """

    from xarray.core.coordinates import Coordinates
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    dict_like_objects = []
    for obj in objects:
        if not isinstance(obj, (DataArray, Dataset, Coordinates, dict)):
            raise TypeError(
                "objects must be an iterable containing only "
                "Dataset(s), DataArray(s), and dictionaries."
            )

        if isinstance(obj, DataArray):
            obj = obj.to_dataset(promote_attrs=True)
        elif isinstance(obj, Coordinates):
            obj = obj.to_dataset()
        dict_like_objects.append(obj)

    merge_result = merge_core(
        dict_like_objects,
        compat,
        join,
        combine_attrs=combine_attrs,
        fill_value=fill_value,
    )
    return Dataset._construct_direct(**merge_result._asdict())


def dataset_merge_method(
    dataset: Dataset,
    other: CoercibleMapping,
    overwrite_vars: Hashable | Iterable[Hashable],
    compat: CompatOptions,
    join: JoinOptions,
    fill_value: Any,
    combine_attrs: CombineAttrsOptions,
) -> _MergeResult:
    """Guts of the Dataset.merge method."""
    # we are locked into supporting overwrite_vars for the Dataset.merge
    # method due for backwards compatibility
    # TODO: consider deprecating it?

    if not isinstance(overwrite_vars, str) and isinstance(overwrite_vars, Iterable):
        overwrite_vars = set(overwrite_vars)
    else:
        overwrite_vars = {overwrite_vars}

    if not overwrite_vars:
        objs = [dataset, other]
        priority_arg = None
    elif overwrite_vars == set(other):
        objs = [dataset, other]
        priority_arg = 1
    else:
        other_overwrite: dict[Hashable, CoercibleValue] = {}
        other_no_overwrite: dict[Hashable, CoercibleValue] = {}
        for k, v in other.items():
            if k in overwrite_vars:
                other_overwrite[k] = v
            else:
                other_no_overwrite[k] = v
        objs = [dataset, other_no_overwrite, other_overwrite]
        priority_arg = 2

    return merge_core(
        objs,
        compat,
        join,
        priority_arg=priority_arg,
        fill_value=fill_value,
        combine_attrs=combine_attrs,
    )


def dataset_update_method(dataset: Dataset, other: CoercibleMapping) -> _MergeResult:
    """Guts of the Dataset.update method.

    This drops a duplicated coordinates from `other` if `other` is not an
    `xarray.Dataset`, e.g., if it's a dict with DataArray values (GH2068,
    GH2180).
    """
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    if not isinstance(other, Dataset):
        other = dict(other)
        for key, value in other.items():
            if isinstance(value, DataArray):
                # drop conflicting coordinates
                coord_names = [
                    c
                    for c in value.coords
                    if c not in value.dims and c in dataset.coords
                ]
                if coord_names:
                    other[key] = value.drop_vars(coord_names)

    return merge_core(
        [dataset, other],
        priority_arg=1,
        indexes=dataset.xindexes,
        combine_attrs="override",
    )
