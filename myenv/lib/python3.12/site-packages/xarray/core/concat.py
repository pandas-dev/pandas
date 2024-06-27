from __future__ import annotations

from collections.abc import Hashable, Iterable
from typing import TYPE_CHECKING, Any, Union, overload

import numpy as np
import pandas as pd

from xarray.core import dtypes, utils
from xarray.core.alignment import align, reindex_variables
from xarray.core.coordinates import Coordinates
from xarray.core.duck_array_ops import lazy_array_equiv
from xarray.core.indexes import Index, PandasIndex
from xarray.core.merge import (
    _VALID_COMPAT,
    collect_variables_and_indexes,
    merge_attrs,
    merge_collected,
)
from xarray.core.types import T_DataArray, T_Dataset, T_Variable
from xarray.core.variable import Variable
from xarray.core.variable import concat as concat_vars

if TYPE_CHECKING:
    from xarray.core.types import (
        CombineAttrsOptions,
        CompatOptions,
        ConcatOptions,
        JoinOptions,
    )

    T_DataVars = Union[ConcatOptions, Iterable[Hashable]]


@overload
def concat(
    objs: Iterable[T_Dataset],
    dim: Hashable | T_Variable | T_DataArray | pd.Index,
    data_vars: T_DataVars = "all",
    coords: ConcatOptions | list[Hashable] = "different",
    compat: CompatOptions = "equals",
    positions: Iterable[Iterable[int]] | None = None,
    fill_value: object = dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "override",
    create_index_for_new_dim: bool = True,
) -> T_Dataset: ...


@overload
def concat(
    objs: Iterable[T_DataArray],
    dim: Hashable | T_Variable | T_DataArray | pd.Index,
    data_vars: T_DataVars = "all",
    coords: ConcatOptions | list[Hashable] = "different",
    compat: CompatOptions = "equals",
    positions: Iterable[Iterable[int]] | None = None,
    fill_value: object = dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "override",
    create_index_for_new_dim: bool = True,
) -> T_DataArray: ...


def concat(
    objs,
    dim,
    data_vars: T_DataVars = "all",
    coords="different",
    compat: CompatOptions = "equals",
    positions=None,
    fill_value=dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "override",
    create_index_for_new_dim: bool = True,
):
    """Concatenate xarray objects along a new or existing dimension.

    Parameters
    ----------
    objs : sequence of Dataset and DataArray
        xarray objects to concatenate together. Each object is expected to
        consist of variables and coordinates with matching shapes except for
        along the concatenated dimension.
    dim : Hashable or Variable or DataArray or pandas.Index
        Name of the dimension to concatenate along. This can either be a new
        dimension name, in which case it is added along axis=0, or an existing
        dimension name, in which case the location of the dimension is
        unchanged. If dimension is provided as a Variable, DataArray or Index, its name
        is used as the dimension to concatenate along and the values are added
        as a coordinate.
    data_vars : {"minimal", "different", "all"} or list of Hashable, optional
        These data variables will be concatenated together:
          * "minimal": Only data variables in which the dimension already
            appears are included.
          * "different": Data variables which are not equal (ignoring
            attributes) across all datasets are also concatenated (as well as
            all for which dimension already appears). Beware: this option may
            load the data payload of data variables into memory if they are not
            already loaded.
          * "all": All data variables will be concatenated.
          * list of dims: The listed data variables will be concatenated, in
            addition to the "minimal" data variables.

        If objects are DataArrays, data_vars must be "all".
    coords : {"minimal", "different", "all"} or list of Hashable, optional
        These coordinate variables will be concatenated together:
          * "minimal": Only coordinates in which the dimension already appears
            are included.
          * "different": Coordinates which are not equal (ignoring attributes)
            across all datasets are also concatenated (as well as all for which
            dimension already appears). Beware: this option may load the data
            payload of coordinate variables into memory if they are not already
            loaded.
          * "all": All coordinate variables will be concatenated, except
            those corresponding to other dimensions.
          * list of Hashable: The listed coordinate variables will be concatenated,
            in addition to the "minimal" coordinates.
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
    positions : None or list of integer arrays, optional
        List of integer arrays which specifies the integer positions to which
        to assign each dataset along the concatenated dimension. If not
        supplied, objects are concatenated in the provided order.
    fill_value : scalar or dict-like, optional
        Value to use for newly missing values. If a dict-like, maps
        variable names to fill values. Use a data array's name to
        refer to its values.
    join : {"outer", "inner", "left", "right", "exact"}, optional
        String indicating how to combine differing indexes
        (excluding dim) in objects

        - "outer": use the union of object indexes
        - "inner": use the intersection of object indexes
        - "left": use indexes from the first object with each dimension
        - "right": use indexes from the last object with each dimension
        - "exact": instead of aligning, raise `ValueError` when indexes to be
          aligned are not equal
        - "override": if indexes are of same size, rewrite indexes to be
          those of the first object with that dimension. Indexes for the same
          dimension must have the same size in all objects.
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
    create_index_for_new_dim : bool, default: True
        Whether to create a new ``PandasIndex`` object when the objects being concatenated contain scalar variables named ``dim``.

    Returns
    -------
    concatenated : type of objs

    See also
    --------
    merge

    Examples
    --------
    >>> da = xr.DataArray(
    ...     np.arange(6).reshape(2, 3), [("x", ["a", "b"]), ("y", [10, 20, 30])]
    ... )
    >>> da
    <xarray.DataArray (x: 2, y: 3)> Size: 48B
    array([[0, 1, 2],
           [3, 4, 5]])
    Coordinates:
      * x        (x) <U1 8B 'a' 'b'
      * y        (y) int64 24B 10 20 30

    >>> xr.concat([da.isel(y=slice(0, 1)), da.isel(y=slice(1, None))], dim="y")
    <xarray.DataArray (x: 2, y: 3)> Size: 48B
    array([[0, 1, 2],
           [3, 4, 5]])
    Coordinates:
      * x        (x) <U1 8B 'a' 'b'
      * y        (y) int64 24B 10 20 30

    >>> xr.concat([da.isel(x=0), da.isel(x=1)], "x")
    <xarray.DataArray (x: 2, y: 3)> Size: 48B
    array([[0, 1, 2],
           [3, 4, 5]])
    Coordinates:
      * x        (x) <U1 8B 'a' 'b'
      * y        (y) int64 24B 10 20 30

    >>> xr.concat([da.isel(x=0), da.isel(x=1)], "new_dim")
    <xarray.DataArray (new_dim: 2, y: 3)> Size: 48B
    array([[0, 1, 2],
           [3, 4, 5]])
    Coordinates:
        x        (new_dim) <U1 8B 'a' 'b'
      * y        (y) int64 24B 10 20 30
    Dimensions without coordinates: new_dim

    >>> xr.concat([da.isel(x=0), da.isel(x=1)], pd.Index([-90, -100], name="new_dim"))
    <xarray.DataArray (new_dim: 2, y: 3)> Size: 48B
    array([[0, 1, 2],
           [3, 4, 5]])
    Coordinates:
        x        (new_dim) <U1 8B 'a' 'b'
      * y        (y) int64 24B 10 20 30
      * new_dim  (new_dim) int64 16B -90 -100

    # Concatenate a scalar variable along a new dimension of the same name with and without creating a new index

    >>> ds = xr.Dataset(coords={"x": 0})
    >>> xr.concat([ds, ds], dim="x")
    <xarray.Dataset> Size: 16B
    Dimensions:  (x: 2)
    Coordinates:
      * x        (x) int64 16B 0 0
    Data variables:
        *empty*

    >>> xr.concat([ds, ds], dim="x").indexes
    Indexes:
        x        Index([0, 0], dtype='int64', name='x')

    >>> xr.concat([ds, ds], dim="x", create_index_for_new_dim=False).indexes
    Indexes:
        *empty*
    """
    # TODO: add ignore_index arguments copied from pandas.concat
    # TODO: support concatenating scalar coordinates even if the concatenated
    # dimension already exists
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    try:
        first_obj, objs = utils.peek_at(objs)
    except StopIteration:
        raise ValueError("must supply at least one object to concatenate")

    if compat not in _VALID_COMPAT:
        raise ValueError(
            f"compat={compat!r} invalid: must be 'broadcast_equals', 'equals', 'identical', 'no_conflicts' or 'override'"
        )

    if isinstance(first_obj, DataArray):
        return _dataarray_concat(
            objs,
            dim=dim,
            data_vars=data_vars,
            coords=coords,
            compat=compat,
            positions=positions,
            fill_value=fill_value,
            join=join,
            combine_attrs=combine_attrs,
            create_index_for_new_dim=create_index_for_new_dim,
        )
    elif isinstance(first_obj, Dataset):
        return _dataset_concat(
            objs,
            dim=dim,
            data_vars=data_vars,
            coords=coords,
            compat=compat,
            positions=positions,
            fill_value=fill_value,
            join=join,
            combine_attrs=combine_attrs,
            create_index_for_new_dim=create_index_for_new_dim,
        )
    else:
        raise TypeError(
            "can only concatenate xarray Dataset and DataArray "
            f"objects, got {type(first_obj)}"
        )


def _calc_concat_dim_index(
    dim_or_data: Hashable | Any,
) -> tuple[Hashable, PandasIndex | None]:
    """Infer the dimension name and 1d index / coordinate variable (if appropriate)
    for concatenating along the new dimension.

    """
    from xarray.core.dataarray import DataArray

    dim: Hashable | None

    if isinstance(dim_or_data, str):
        dim = dim_or_data
        index = None
    else:
        if not isinstance(dim_or_data, (DataArray, Variable)):
            dim = getattr(dim_or_data, "name", None)
            if dim is None:
                dim = "concat_dim"
        else:
            (dim,) = dim_or_data.dims
        coord_dtype = getattr(dim_or_data, "dtype", None)
        index = PandasIndex(dim_or_data, dim, coord_dtype=coord_dtype)

    return dim, index


def _calc_concat_over(datasets, dim, dim_names, data_vars: T_DataVars, coords, compat):
    """
    Determine which dataset variables need to be concatenated in the result,
    """
    # Return values
    concat_over = set()
    equals = {}

    if dim in dim_names:
        concat_over_existing_dim = True
        concat_over.add(dim)
    else:
        concat_over_existing_dim = False

    concat_dim_lengths = []
    for ds in datasets:
        if concat_over_existing_dim:
            if dim not in ds.dims:
                if dim in ds:
                    ds = ds.set_coords(dim)
        concat_over.update(k for k, v in ds.variables.items() if dim in v.dims)
        concat_dim_lengths.append(ds.sizes.get(dim, 1))

    def process_subset_opt(opt, subset):
        if isinstance(opt, str):
            if opt == "different":
                if compat == "override":
                    raise ValueError(
                        f"Cannot specify both {subset}='different' and compat='override'."
                    )
                # all nonindexes that are not the same in each dataset
                for k in getattr(datasets[0], subset):
                    if k not in concat_over:
                        equals[k] = None

                        variables = [
                            ds.variables[k] for ds in datasets if k in ds.variables
                        ]

                        if len(variables) == 1:
                            # coords="different" doesn't make sense when only one object
                            # contains a particular variable.
                            break
                        elif len(variables) != len(datasets) and opt == "different":
                            raise ValueError(
                                f"{k!r} not present in all datasets and coords='different'. "
                                f"Either add {k!r} to datasets where it is missing or "
                                "specify coords='minimal'."
                            )

                        # first check without comparing values i.e. no computes
                        for var in variables[1:]:
                            equals[k] = getattr(variables[0], compat)(
                                var, equiv=lazy_array_equiv
                            )
                            if equals[k] is not True:
                                # exit early if we know these are not equal or that
                                # equality cannot be determined i.e. one or all of
                                # the variables wraps a numpy array
                                break

                        if equals[k] is False:
                            concat_over.add(k)

                        elif equals[k] is None:
                            # Compare the variable of all datasets vs. the one
                            # of the first dataset. Perform the minimum amount of
                            # loads in order to avoid multiple loads from disk
                            # while keeping the RAM footprint low.
                            v_lhs = datasets[0].variables[k].load()
                            # We'll need to know later on if variables are equal.
                            computed = []
                            for ds_rhs in datasets[1:]:
                                v_rhs = ds_rhs.variables[k].compute()
                                computed.append(v_rhs)
                                if not getattr(v_lhs, compat)(v_rhs):
                                    concat_over.add(k)
                                    equals[k] = False
                                    # computed variables are not to be re-computed
                                    # again in the future
                                    for ds, v in zip(datasets[1:], computed):
                                        ds.variables[k].data = v.data
                                    break
                            else:
                                equals[k] = True

            elif opt == "all":
                concat_over.update(
                    set().union(
                        *list(set(getattr(d, subset)) - set(d.dims) for d in datasets)
                    )
                )
            elif opt == "minimal":
                pass
            else:
                raise ValueError(f"unexpected value for {subset}: {opt}")
        else:
            valid_vars = tuple(getattr(datasets[0], subset))
            invalid_vars = [k for k in opt if k not in valid_vars]
            if invalid_vars:
                if subset == "coords":
                    raise ValueError(
                        f"the variables {invalid_vars} in coords are not "
                        f"found in the coordinates of the first dataset {valid_vars}"
                    )
                else:
                    # note: data_vars are not listed in the error message here,
                    # because there may be lots of them
                    raise ValueError(
                        f"the variables {invalid_vars} in data_vars are not "
                        f"found in the data variables of the first dataset"
                    )
            concat_over.update(opt)

    process_subset_opt(data_vars, "data_vars")
    process_subset_opt(coords, "coords")
    return concat_over, equals, concat_dim_lengths


# determine dimensional coordinate names and a dict mapping name to DataArray
def _parse_datasets(
    datasets: list[T_Dataset],
) -> tuple[
    dict[Hashable, Variable],
    dict[Hashable, int],
    set[Hashable],
    set[Hashable],
    list[Hashable],
]:
    dims: set[Hashable] = set()
    all_coord_names: set[Hashable] = set()
    data_vars: set[Hashable] = set()  # list of data_vars
    dim_coords: dict[Hashable, Variable] = {}  # maps dim name to variable
    dims_sizes: dict[Hashable, int] = {}  # shared dimension sizes to expand variables
    variables_order: dict[Hashable, Variable] = {}  # variables in order of appearance

    for ds in datasets:
        dims_sizes.update(ds.sizes)
        all_coord_names.update(ds.coords)
        data_vars.update(ds.data_vars)
        variables_order.update(ds.variables)

        # preserves ordering of dimensions
        for dim in ds.dims:
            if dim in dims:
                continue

            if dim in ds.coords and dim not in dim_coords:
                dim_coords[dim] = ds.coords[dim].variable
        dims = dims | set(ds.dims)

    return dim_coords, dims_sizes, all_coord_names, data_vars, list(variables_order)


def _dataset_concat(
    datasets: list[T_Dataset],
    dim: str | T_Variable | T_DataArray | pd.Index,
    data_vars: T_DataVars,
    coords: str | list[str],
    compat: CompatOptions,
    positions: Iterable[Iterable[int]] | None,
    fill_value: Any = dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "override",
    create_index_for_new_dim: bool = True,
) -> T_Dataset:
    """
    Concatenate a sequence of datasets along a new or existing dimension
    """
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    datasets = list(datasets)

    if not all(isinstance(dataset, Dataset) for dataset in datasets):
        raise TypeError(
            "The elements in the input list need to be either all 'Dataset's or all 'DataArray's"
        )

    if isinstance(dim, DataArray):
        dim_var = dim.variable
    elif isinstance(dim, Variable):
        dim_var = dim
    else:
        dim_var = None

    dim, index = _calc_concat_dim_index(dim)

    # Make sure we're working on a copy (we'll be loading variables)
    datasets = [ds.copy() for ds in datasets]
    datasets = list(
        align(*datasets, join=join, copy=False, exclude=[dim], fill_value=fill_value)
    )

    dim_coords, dims_sizes, coord_names, data_names, vars_order = _parse_datasets(
        datasets
    )
    dim_names = set(dim_coords)

    both_data_and_coords = coord_names & data_names
    if both_data_and_coords:
        raise ValueError(
            f"{both_data_and_coords!r} is a coordinate in some datasets but not others."
        )
    # we don't want the concat dimension in the result dataset yet
    dim_coords.pop(dim, None)
    dims_sizes.pop(dim, None)

    # case where concat dimension is a coordinate or data_var but not a dimension
    if (dim in coord_names or dim in data_names) and dim not in dim_names:
        datasets = [
            ds.expand_dims(dim, create_index_for_new_dim=create_index_for_new_dim)
            for ds in datasets
        ]

    # determine which variables to concatenate
    concat_over, equals, concat_dim_lengths = _calc_concat_over(
        datasets, dim, dim_names, data_vars, coords, compat
    )

    # determine which variables to merge, and then merge them according to compat
    variables_to_merge = (coord_names | data_names) - concat_over

    result_vars = {}
    result_indexes = {}

    if variables_to_merge:
        grouped = {
            k: v
            for k, v in collect_variables_and_indexes(datasets).items()
            if k in variables_to_merge
        }
        merged_vars, merged_indexes = merge_collected(
            grouped, compat=compat, equals=equals
        )
        result_vars.update(merged_vars)
        result_indexes.update(merged_indexes)

    result_vars.update(dim_coords)

    # assign attrs and encoding from first dataset
    result_attrs = merge_attrs([ds.attrs for ds in datasets], combine_attrs)
    result_encoding = datasets[0].encoding

    # check that global attributes are fixed across all datasets if necessary
    if compat == "identical":
        for ds in datasets[1:]:
            if not utils.dict_equiv(ds.attrs, result_attrs):
                raise ValueError("Dataset global attributes not equal.")

    # we've already verified everything is consistent; now, calculate
    # shared dimension sizes so we can expand the necessary variables
    def ensure_common_dims(vars, concat_dim_lengths):
        # ensure each variable with the given name shares the same
        # dimensions and the same shape for all of them except along the
        # concat dimension
        common_dims = tuple(utils.OrderedSet(d for v in vars for d in v.dims))
        if dim not in common_dims:
            common_dims = (dim,) + common_dims
        for var, dim_len in zip(vars, concat_dim_lengths):
            if var.dims != common_dims:
                common_shape = tuple(dims_sizes.get(d, dim_len) for d in common_dims)
                var = var.set_dims(common_dims, common_shape)
            yield var

    # get the indexes to concatenate together, create a PandasIndex
    # for any scalar coordinate variable found with ``name`` matching ``dim``.
    # TODO: depreciate concat a mix of scalar and dimensional indexed coordinates?
    # TODO: (benbovy - explicit indexes): check index types and/or coordinates
    # of all datasets?
    def get_indexes(name):
        for ds in datasets:
            if name in ds._indexes:
                yield ds._indexes[name]
            elif name == dim:
                var = ds._variables[name]
                if not var.dims:
                    data = var.set_dims(dim).values
                    if create_index_for_new_dim:
                        yield PandasIndex(data, dim, coord_dtype=var.dtype)

    # create concatenation index, needed for later reindexing
    file_start_indexes = np.append(0, np.cumsum(concat_dim_lengths))
    concat_index = np.arange(file_start_indexes[-1])
    concat_index_size = concat_index.size
    variable_index_mask = np.ones(concat_index_size, dtype=bool)

    # stack up each variable and/or index to fill-out the dataset (in order)
    # n.b. this loop preserves variable order, needed for groupby.
    ndatasets = len(datasets)
    for name in vars_order:
        if name in concat_over and name not in result_indexes:
            variables = []
            # Initialize the mask to all True then set False if any name is missing in
            # the datasets:
            variable_index_mask.fill(True)
            var_concat_dim_length = []
            for i, ds in enumerate(datasets):
                if name in ds.variables:
                    variables.append(ds[name].variable)
                    var_concat_dim_length.append(concat_dim_lengths[i])
                else:
                    # raise if coordinate not in all datasets
                    if name in coord_names:
                        raise ValueError(
                            f"coordinate {name!r} not present in all datasets."
                        )

                    # Mask out the indexes without the name:
                    start = file_start_indexes[i]
                    end = file_start_indexes[i + 1]
                    variable_index_mask[slice(start, end)] = False

            variable_index = concat_index[variable_index_mask]
            vars = ensure_common_dims(variables, var_concat_dim_length)

            # Try to concatenate the indexes, concatenate the variables when no index
            # is found on all datasets.
            indexes: list[Index] = list(get_indexes(name))
            if indexes:
                if len(indexes) < ndatasets:
                    raise ValueError(
                        f"{name!r} must have either an index or no index in all datasets, "
                        f"found {len(indexes)}/{len(datasets)} datasets with an index."
                    )
                combined_idx = indexes[0].concat(indexes, dim, positions)
                if name in datasets[0]._indexes:
                    idx_vars = datasets[0].xindexes.get_all_coords(name)
                else:
                    # index created from a scalar coordinate
                    idx_vars = {name: datasets[0][name].variable}
                result_indexes.update({k: combined_idx for k in idx_vars})
                combined_idx_vars = combined_idx.create_variables(idx_vars)
                for k, v in combined_idx_vars.items():
                    v.attrs = merge_attrs(
                        [ds.variables[k].attrs for ds in datasets],
                        combine_attrs=combine_attrs,
                    )
                    result_vars[k] = v
            else:
                combined_var = concat_vars(
                    vars, dim, positions, combine_attrs=combine_attrs
                )
                # reindex if variable is not present in all datasets
                if len(variable_index) < concat_index_size:
                    combined_var = reindex_variables(
                        variables={name: combined_var},
                        dim_pos_indexers={
                            dim: pd.Index(variable_index).get_indexer(concat_index)
                        },
                        fill_value=fill_value,
                    )[name]
                result_vars[name] = combined_var

        elif name in result_vars:
            # preserves original variable order
            result_vars[name] = result_vars.pop(name)

    absent_coord_names = coord_names - set(result_vars)
    if absent_coord_names:
        raise ValueError(
            f"Variables {absent_coord_names!r} are coordinates in some datasets but not others."
        )

    result_data_vars = {}
    coord_vars = {}
    for name, result_var in result_vars.items():
        if name in coord_names:
            coord_vars[name] = result_var
        else:
            result_data_vars[name] = result_var

    if index is not None:
        if dim_var is not None:
            index_vars = index.create_variables({dim: dim_var})
        else:
            index_vars = index.create_variables()

        coord_vars[dim] = index_vars[dim]
        result_indexes[dim] = index

    coords_obj = Coordinates(coord_vars, indexes=result_indexes)

    result = type(datasets[0])(result_data_vars, coords=coords_obj, attrs=result_attrs)
    result.encoding = result_encoding

    return result


def _dataarray_concat(
    arrays: Iterable[T_DataArray],
    dim: str | T_Variable | T_DataArray | pd.Index,
    data_vars: T_DataVars,
    coords: str | list[str],
    compat: CompatOptions,
    positions: Iterable[Iterable[int]] | None,
    fill_value: object = dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "override",
    create_index_for_new_dim: bool = True,
) -> T_DataArray:
    from xarray.core.dataarray import DataArray

    arrays = list(arrays)

    if not all(isinstance(array, DataArray) for array in arrays):
        raise TypeError(
            "The elements in the input list need to be either all 'Dataset's or all 'DataArray's"
        )

    if data_vars != "all":
        raise ValueError(
            "data_vars is not a valid argument when concatenating DataArray objects"
        )

    datasets = []
    for n, arr in enumerate(arrays):
        if n == 0:
            name = arr.name
        elif name != arr.name:
            if compat == "identical":
                raise ValueError("array names not identical")
            else:
                arr = arr.rename(name)
        datasets.append(arr._to_temp_dataset())

    ds = _dataset_concat(
        datasets,
        dim,
        data_vars,
        coords,
        compat,
        positions,
        fill_value=fill_value,
        join=join,
        combine_attrs=combine_attrs,
        create_index_for_new_dim=create_index_for_new_dim,
    )

    merged_attrs = merge_attrs([da.attrs for da in arrays], combine_attrs)

    result = arrays[0]._from_temp_dataset(ds, name)
    result.attrs = merged_attrs

    return result
