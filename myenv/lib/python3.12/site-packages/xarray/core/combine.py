from __future__ import annotations

import itertools
from collections import Counter
from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING, Literal, Union

import pandas as pd

from xarray.core import dtypes
from xarray.core.concat import concat
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.merge import merge
from xarray.core.utils import iterate_nested

if TYPE_CHECKING:
    from xarray.core.types import CombineAttrsOptions, CompatOptions, JoinOptions


def _infer_concat_order_from_positions(datasets):
    return dict(_infer_tile_ids_from_nested_list(datasets, ()))


def _infer_tile_ids_from_nested_list(entry, current_pos):
    """
    Given a list of lists (of lists...) of objects, returns a iterator
    which returns a tuple containing the index of each object in the nested
    list structure as the key, and the object. This can then be called by the
    dict constructor to create a dictionary of the objects organised by their
    position in the original nested list.

    Recursively traverses the given structure, while keeping track of the
    current position. Should work for any type of object which isn't a list.

    Parameters
    ----------
    entry : list[list[obj, obj, ...], ...]
        List of lists of arbitrary depth, containing objects in the order
        they are to be concatenated.

    Returns
    -------
    combined_tile_ids : dict[tuple(int, ...), obj]
    """

    if isinstance(entry, list):
        for i, item in enumerate(entry):
            yield from _infer_tile_ids_from_nested_list(item, current_pos + (i,))
    else:
        yield current_pos, entry


def _ensure_same_types(series, dim):
    if series.dtype == object:
        types = set(series.map(type))
        if len(types) > 1:
            try:
                import cftime

                cftimes = any(issubclass(t, cftime.datetime) for t in types)
            except ImportError:
                cftimes = False

            types = ", ".join(t.__name__ for t in types)

            error_msg = (
                f"Cannot combine along dimension '{dim}' with mixed types."
                f" Found: {types}."
            )
            if cftimes:
                error_msg = (
                    f"{error_msg} If importing data directly from a file then "
                    f"setting `use_cftime=True` may fix this issue."
                )

            raise TypeError(error_msg)


def _infer_concat_order_from_coords(datasets):
    concat_dims = []
    tile_ids = [() for ds in datasets]

    # All datasets have same variables because they've been grouped as such
    ds0 = datasets[0]
    for dim in ds0.dims:
        # Check if dim is a coordinate dimension
        if dim in ds0:
            # Need to read coordinate values to do ordering
            indexes = [ds._indexes.get(dim) for ds in datasets]
            if any(index is None for index in indexes):
                raise ValueError(
                    "Every dimension needs a coordinate for "
                    "inferring concatenation order"
                )

            # TODO (benbovy, flexible indexes): support flexible indexes?
            indexes = [index.to_pandas_index() for index in indexes]

            # If dimension coordinate values are same on every dataset then
            # should be leaving this dimension alone (it's just a "bystander")
            if not all(index.equals(indexes[0]) for index in indexes[1:]):
                # Infer order datasets should be arranged in along this dim
                concat_dims.append(dim)

                if all(index.is_monotonic_increasing for index in indexes):
                    ascending = True
                elif all(index.is_monotonic_decreasing for index in indexes):
                    ascending = False
                else:
                    raise ValueError(
                        f"Coordinate variable {dim} is neither "
                        "monotonically increasing nor "
                        "monotonically decreasing on all datasets"
                    )

                # Assume that any two datasets whose coord along dim starts
                # with the same value have the same coord values throughout.
                if any(index.size == 0 for index in indexes):
                    raise ValueError("Cannot handle size zero dimensions")
                first_items = pd.Index([index[0] for index in indexes])

                series = first_items.to_series()

                # ensure series does not contain mixed types, e.g. cftime calendars
                _ensure_same_types(series, dim)

                # Sort datasets along dim
                # We want rank but with identical elements given identical
                # position indices - they should be concatenated along another
                # dimension, not along this one
                rank = series.rank(
                    method="dense", ascending=ascending, numeric_only=False
                )
                order = rank.astype(int).values - 1

                # Append positions along extra dimension to structure which
                # encodes the multi-dimensional concatenation order
                tile_ids = [
                    tile_id + (position,) for tile_id, position in zip(tile_ids, order)
                ]

    if len(datasets) > 1 and not concat_dims:
        raise ValueError(
            "Could not find any dimension coordinates to use to "
            "order the datasets for concatenation"
        )

    combined_ids = dict(zip(tile_ids, datasets))

    return combined_ids, concat_dims


def _check_dimension_depth_tile_ids(combined_tile_ids):
    """
    Check all tuples are the same length, i.e. check that all lists are
    nested to the same depth.
    """
    tile_ids = combined_tile_ids.keys()
    nesting_depths = [len(tile_id) for tile_id in tile_ids]
    if not nesting_depths:
        nesting_depths = [0]
    if set(nesting_depths) != {nesting_depths[0]}:
        raise ValueError(
            "The supplied objects do not form a hypercube because"
            " sub-lists do not have consistent depths"
        )
    # return these just to be reused in _check_shape_tile_ids
    return tile_ids, nesting_depths


def _check_shape_tile_ids(combined_tile_ids):
    """Check all lists along one dimension are same length."""
    tile_ids, nesting_depths = _check_dimension_depth_tile_ids(combined_tile_ids)
    for dim in range(nesting_depths[0]):
        indices_along_dim = [tile_id[dim] for tile_id in tile_ids]
        occurrences = Counter(indices_along_dim)
        if len(set(occurrences.values())) != 1:
            raise ValueError(
                "The supplied objects do not form a hypercube "
                "because sub-lists do not have consistent "
                f"lengths along dimension {dim}"
            )


def _combine_nd(
    combined_ids,
    concat_dims,
    data_vars="all",
    coords="different",
    compat: CompatOptions = "no_conflicts",
    fill_value=dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "drop",
):
    """
    Combines an N-dimensional structure of datasets into one by applying a
    series of either concat and merge operations along each dimension.

    No checks are performed on the consistency of the datasets, concat_dims or
    tile_IDs, because it is assumed that this has already been done.

    Parameters
    ----------
    combined_ids : Dict[Tuple[int, ...]], xarray.Dataset]
        Structure containing all datasets to be concatenated with "tile_IDs" as
        keys, which specify position within the desired final combined result.
    concat_dims : sequence of str
        The dimensions along which the datasets should be concatenated. Must be
        in order, and the length must match the length of the tuples used as
        keys in combined_ids. If the string is a dimension name then concat
        along that dimension, if it is None then merge.

    Returns
    -------
    combined_ds : xarray.Dataset
    """

    example_tile_id = next(iter(combined_ids.keys()))

    n_dims = len(example_tile_id)
    if len(concat_dims) != n_dims:
        raise ValueError(
            f"concat_dims has length {len(concat_dims)} but the datasets "
            f"passed are nested in a {n_dims}-dimensional structure"
        )

    # Each iteration of this loop reduces the length of the tile_ids tuples
    # by one. It always combines along the first dimension, removing the first
    # element of the tuple
    for concat_dim in concat_dims:
        combined_ids = _combine_all_along_first_dim(
            combined_ids,
            dim=concat_dim,
            data_vars=data_vars,
            coords=coords,
            compat=compat,
            fill_value=fill_value,
            join=join,
            combine_attrs=combine_attrs,
        )
    (combined_ds,) = combined_ids.values()
    return combined_ds


def _combine_all_along_first_dim(
    combined_ids,
    dim,
    data_vars,
    coords,
    compat: CompatOptions,
    fill_value=dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "drop",
):
    # Group into lines of datasets which must be combined along dim
    # need to sort by _new_tile_id first for groupby to work
    # TODO: is the sorted need?
    combined_ids = dict(sorted(combined_ids.items(), key=_new_tile_id))
    grouped = itertools.groupby(combined_ids.items(), key=_new_tile_id)

    # Combine all of these datasets along dim
    new_combined_ids = {}
    for new_id, group in grouped:
        combined_ids = dict(sorted(group))
        datasets = combined_ids.values()
        new_combined_ids[new_id] = _combine_1d(
            datasets, dim, compat, data_vars, coords, fill_value, join, combine_attrs
        )
    return new_combined_ids


def _combine_1d(
    datasets,
    concat_dim,
    compat: CompatOptions = "no_conflicts",
    data_vars="all",
    coords="different",
    fill_value=dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "drop",
):
    """
    Applies either concat or merge to 1D list of datasets depending on value
    of concat_dim
    """

    if concat_dim is not None:
        try:
            combined = concat(
                datasets,
                dim=concat_dim,
                data_vars=data_vars,
                coords=coords,
                compat=compat,
                fill_value=fill_value,
                join=join,
                combine_attrs=combine_attrs,
            )
        except ValueError as err:
            if "encountered unexpected variable" in str(err):
                raise ValueError(
                    "These objects cannot be combined using only "
                    "xarray.combine_nested, instead either use "
                    "xarray.combine_by_coords, or do it manually "
                    "with xarray.concat, xarray.merge and "
                    "xarray.align"
                )
            else:
                raise
    else:
        combined = merge(
            datasets,
            compat=compat,
            fill_value=fill_value,
            join=join,
            combine_attrs=combine_attrs,
        )

    return combined


def _new_tile_id(single_id_ds_pair):
    tile_id, ds = single_id_ds_pair
    return tile_id[1:]


def _nested_combine(
    datasets,
    concat_dims,
    compat,
    data_vars,
    coords,
    ids,
    fill_value=dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "drop",
):
    if len(datasets) == 0:
        return Dataset()

    # Arrange datasets for concatenation
    # Use information from the shape of the user input
    if not ids:
        # Determine tile_IDs by structure of input in N-D
        # (i.e. ordering in list-of-lists)
        combined_ids = _infer_concat_order_from_positions(datasets)
    else:
        # Already sorted so just use the ids already passed
        combined_ids = dict(zip(ids, datasets))

    # Check that the inferred shape is combinable
    _check_shape_tile_ids(combined_ids)

    # Apply series of concatenate or merge operations along each dimension
    combined = _combine_nd(
        combined_ids,
        concat_dims,
        compat=compat,
        data_vars=data_vars,
        coords=coords,
        fill_value=fill_value,
        join=join,
        combine_attrs=combine_attrs,
    )
    return combined


# Define type for arbitrarily-nested list of lists recursively:
DATASET_HYPERCUBE = Union[Dataset, Iterable["DATASET_HYPERCUBE"]]


def combine_nested(
    datasets: DATASET_HYPERCUBE,
    concat_dim: str | DataArray | None | Sequence[str | DataArray | pd.Index | None],
    compat: str = "no_conflicts",
    data_vars: str = "all",
    coords: str = "different",
    fill_value: object = dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "drop",
) -> Dataset:
    """
    Explicitly combine an N-dimensional grid of datasets into one by using a
    succession of concat and merge operations along each dimension of the grid.

    Does not sort the supplied datasets under any circumstances, so the
    datasets must be passed in the order you wish them to be concatenated. It
    does align coordinates, but different variables on datasets can cause it to
    fail under some scenarios. In complex cases, you may need to clean up your
    data and use concat/merge explicitly.

    To concatenate along multiple dimensions the datasets must be passed as a
    nested list-of-lists, with a depth equal to the length of ``concat_dims``.
    ``combine_nested`` will concatenate along the top-level list first.

    Useful for combining datasets from a set of nested directories, or for
    collecting the output of a simulation parallelized along multiple
    dimensions.

    Parameters
    ----------
    datasets : list or nested list of Dataset
        Dataset objects to combine.
        If concatenation or merging along more than one dimension is desired,
        then datasets must be supplied in a nested list-of-lists.
    concat_dim : str, or list of str, DataArray, Index or None
        Dimensions along which to concatenate variables, as used by
        :py:func:`xarray.concat`.
        Set ``concat_dim=[..., None, ...]`` explicitly to disable concatenation
        and merge instead along a particular dimension.
        The position of ``None`` in the list specifies the dimension of the
        nested-list input along which to merge.
        Must be the same length as the depth of the list passed to
        ``datasets``.
    compat : {"identical", "equals", "broadcast_equals", \
              "no_conflicts", "override"}, optional
        String indicating how to compare variables of the same name for
        potential merge conflicts:

        - "broadcast_equals": all values must be equal when variables are
          broadcast against each other to ensure common dimensions.
        - "equals": all values and dimensions must be the same.
        - "identical": all values, dimensions and attributes must be the
          same.
        - "no_conflicts": only values which are not null in both datasets
          must be equal. The returned dataset then contains the combination
          of all non-null values.
        - "override": skip comparing and pick variable from first dataset
    data_vars : {"minimal", "different", "all" or list of str}, optional
        Details are in the documentation of concat
    coords : {"minimal", "different", "all" or list of str}, optional
        Details are in the documentation of concat
    fill_value : scalar or dict-like, optional
        Value to use for newly missing values. If a dict-like, maps
        variable names to fill values. Use a data array's name to
        refer to its values.
    join : {"outer", "inner", "left", "right", "exact"}, optional
        String indicating how to combine differing indexes
        (excluding concat_dim) in objects

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
    combined : xarray.Dataset

    Examples
    --------

    A common task is collecting data from a parallelized simulation in which
    each process wrote out to a separate file. A domain which was decomposed
    into 4 parts, 2 each along both the x and y axes, requires organising the
    datasets into a doubly-nested list, e.g:

    >>> x1y1 = xr.Dataset(
    ...     {
    ...         "temperature": (("x", "y"), np.random.randn(2, 2)),
    ...         "precipitation": (("x", "y"), np.random.randn(2, 2)),
    ...     }
    ... )
    >>> x1y1
    <xarray.Dataset> Size: 64B
    Dimensions:        (x: 2, y: 2)
    Dimensions without coordinates: x, y
    Data variables:
        temperature    (x, y) float64 32B 1.764 0.4002 0.9787 2.241
        precipitation  (x, y) float64 32B 1.868 -0.9773 0.9501 -0.1514
    >>> x1y2 = xr.Dataset(
    ...     {
    ...         "temperature": (("x", "y"), np.random.randn(2, 2)),
    ...         "precipitation": (("x", "y"), np.random.randn(2, 2)),
    ...     }
    ... )
    >>> x2y1 = xr.Dataset(
    ...     {
    ...         "temperature": (("x", "y"), np.random.randn(2, 2)),
    ...         "precipitation": (("x", "y"), np.random.randn(2, 2)),
    ...     }
    ... )
    >>> x2y2 = xr.Dataset(
    ...     {
    ...         "temperature": (("x", "y"), np.random.randn(2, 2)),
    ...         "precipitation": (("x", "y"), np.random.randn(2, 2)),
    ...     }
    ... )


    >>> ds_grid = [[x1y1, x1y2], [x2y1, x2y2]]
    >>> combined = xr.combine_nested(ds_grid, concat_dim=["x", "y"])
    >>> combined
    <xarray.Dataset> Size: 256B
    Dimensions:        (x: 4, y: 4)
    Dimensions without coordinates: x, y
    Data variables:
        temperature    (x, y) float64 128B 1.764 0.4002 -0.1032 ... 0.04576 -0.1872
        precipitation  (x, y) float64 128B 1.868 -0.9773 0.761 ... 0.1549 0.3782

    ``combine_nested`` can also be used to explicitly merge datasets with
    different variables. For example if we have 4 datasets, which are divided
    along two times, and contain two different variables, we can pass ``None``
    to ``concat_dim`` to specify the dimension of the nested list over which
    we wish to use ``merge`` instead of ``concat``:

    >>> t1temp = xr.Dataset({"temperature": ("t", np.random.randn(5))})
    >>> t1temp
    <xarray.Dataset> Size: 40B
    Dimensions:      (t: 5)
    Dimensions without coordinates: t
    Data variables:
        temperature  (t) float64 40B -0.8878 -1.981 -0.3479 0.1563 1.23

    >>> t1precip = xr.Dataset({"precipitation": ("t", np.random.randn(5))})
    >>> t1precip
    <xarray.Dataset> Size: 40B
    Dimensions:        (t: 5)
    Dimensions without coordinates: t
    Data variables:
        precipitation  (t) float64 40B 1.202 -0.3873 -0.3023 -1.049 -1.42

    >>> t2temp = xr.Dataset({"temperature": ("t", np.random.randn(5))})
    >>> t2precip = xr.Dataset({"precipitation": ("t", np.random.randn(5))})


    >>> ds_grid = [[t1temp, t1precip], [t2temp, t2precip]]
    >>> combined = xr.combine_nested(ds_grid, concat_dim=["t", None])
    >>> combined
    <xarray.Dataset> Size: 160B
    Dimensions:        (t: 10)
    Dimensions without coordinates: t
    Data variables:
        temperature    (t) float64 80B -0.8878 -1.981 -0.3479 ... -0.4381 -1.253
        precipitation  (t) float64 80B 1.202 -0.3873 -0.3023 ... -0.8955 0.3869

    See also
    --------
    concat
    merge
    """
    mixed_datasets_and_arrays = any(
        isinstance(obj, Dataset) for obj in iterate_nested(datasets)
    ) and any(
        isinstance(obj, DataArray) and obj.name is None
        for obj in iterate_nested(datasets)
    )
    if mixed_datasets_and_arrays:
        raise ValueError("Can't combine datasets with unnamed arrays.")

    if isinstance(concat_dim, (str, DataArray)) or concat_dim is None:
        concat_dim = [concat_dim]

    # The IDs argument tells _nested_combine that datasets aren't yet sorted
    return _nested_combine(
        datasets,
        concat_dims=concat_dim,
        compat=compat,
        data_vars=data_vars,
        coords=coords,
        ids=False,
        fill_value=fill_value,
        join=join,
        combine_attrs=combine_attrs,
    )


def vars_as_keys(ds):
    return tuple(sorted(ds))


def _combine_single_variable_hypercube(
    datasets,
    fill_value=dtypes.NA,
    data_vars="all",
    coords="different",
    compat: CompatOptions = "no_conflicts",
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "no_conflicts",
):
    """
    Attempt to combine a list of Datasets into a hypercube using their
    coordinates.

    All provided Datasets must belong to a single variable, ie. must be
    assigned the same variable name. This precondition is not checked by this
    function, so the caller is assumed to know what it's doing.

    This function is NOT part of the public API.
    """
    if len(datasets) == 0:
        raise ValueError(
            "At least one Dataset is required to resolve variable names "
            "for combined hypercube."
        )

    combined_ids, concat_dims = _infer_concat_order_from_coords(list(datasets))

    if fill_value is None:
        # check that datasets form complete hypercube
        _check_shape_tile_ids(combined_ids)
    else:
        # check only that all datasets have same dimension depth for these
        # vars
        _check_dimension_depth_tile_ids(combined_ids)

    # Concatenate along all of concat_dims one by one to create single ds
    concatenated = _combine_nd(
        combined_ids,
        concat_dims=concat_dims,
        data_vars=data_vars,
        coords=coords,
        compat=compat,
        fill_value=fill_value,
        join=join,
        combine_attrs=combine_attrs,
    )

    # Check the overall coordinates are monotonically increasing
    for dim in concat_dims:
        indexes = concatenated.indexes.get(dim)
        if not (indexes.is_monotonic_increasing or indexes.is_monotonic_decreasing):
            raise ValueError(
                "Resulting object does not have monotonic"
                f" global indexes along dimension {dim}"
            )

    return concatenated


def combine_by_coords(
    data_objects: Iterable[Dataset | DataArray] = [],
    compat: CompatOptions = "no_conflicts",
    data_vars: Literal["all", "minimal", "different"] | list[str] = "all",
    coords: str = "different",
    fill_value: object = dtypes.NA,
    join: JoinOptions = "outer",
    combine_attrs: CombineAttrsOptions = "no_conflicts",
) -> Dataset | DataArray:
    """

    Attempt to auto-magically combine the given datasets (or data arrays)
    into one by using dimension coordinates.

    This function attempts to combine a group of datasets along any number of
    dimensions into a single entity by inspecting coords and metadata and using
    a combination of concat and merge.

    Will attempt to order the datasets such that the values in their dimension
    coordinates are monotonic along all dimensions. If it cannot determine the
    order in which to concatenate the datasets, it will raise a ValueError.
    Non-coordinate dimensions will be ignored, as will any coordinate
    dimensions which do not vary between each dataset.

    Aligns coordinates, but different variables on datasets can cause it
    to fail under some scenarios. In complex cases, you may need to clean up
    your data and use concat/merge explicitly (also see `combine_nested`).

    Works well if, for example, you have N years of data and M data variables,
    and each combination of a distinct time period and set of data variables is
    saved as its own dataset. Also useful for if you have a simulation which is
    parallelized in multiple dimensions, but has global coordinates saved in
    each file specifying the positions of points within the global domain.

    Parameters
    ----------
    data_objects : Iterable of Datasets or DataArrays
        Data objects to combine.

    compat : {"identical", "equals", "broadcast_equals", "no_conflicts", "override"}, optional
        String indicating how to compare variables of the same name for
        potential conflicts:

        - "broadcast_equals": all values must be equal when variables are
          broadcast against each other to ensure common dimensions.
        - "equals": all values and dimensions must be the same.
        - "identical": all values, dimensions and attributes must be the
          same.
        - "no_conflicts": only values which are not null in both datasets
          must be equal. The returned dataset then contains the combination
          of all non-null values.
        - "override": skip comparing and pick variable from first dataset

    data_vars : {"minimal", "different", "all" or list of str}, optional
        These data variables will be concatenated together:

        - "minimal": Only data variables in which the dimension already
          appears are included.
        - "different": Data variables which are not equal (ignoring
          attributes) across all datasets are also concatenated (as well as
          all for which dimension already appears). Beware: this option may
          load the data payload of data variables into memory if they are not
          already loaded.
        - "all": All data variables will be concatenated.
        - list of str: The listed data variables will be concatenated, in
          addition to the "minimal" data variables.

        If objects are DataArrays, `data_vars` must be "all".
    coords : {"minimal", "different", "all"} or list of str, optional
        As per the "data_vars" kwarg, but for coordinate variables.
    fill_value : scalar or dict-like, optional
        Value to use for newly missing values. If a dict-like, maps
        variable names to fill values. Use a data array's name to
        refer to its values. If None, raises a ValueError if
        the passed Datasets do not create a complete hypercube.
    join : {"outer", "inner", "left", "right", "exact"}, optional
        String indicating how to combine differing indexes in objects

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
                     "override"} or callable, default: "no_conflicts"
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
    combined : xarray.Dataset or xarray.DataArray
        Will return a Dataset unless all the inputs are unnamed DataArrays, in which case a
        DataArray will be returned.

    See also
    --------
    concat
    merge
    combine_nested

    Examples
    --------

    Combining two datasets using their common dimension coordinates. Notice
    they are concatenated based on the values in their dimension coordinates,
    not on their position in the list passed to `combine_by_coords`.

    >>> x1 = xr.Dataset(
    ...     {
    ...         "temperature": (("y", "x"), 20 * np.random.rand(6).reshape(2, 3)),
    ...         "precipitation": (("y", "x"), np.random.rand(6).reshape(2, 3)),
    ...     },
    ...     coords={"y": [0, 1], "x": [10, 20, 30]},
    ... )
    >>> x2 = xr.Dataset(
    ...     {
    ...         "temperature": (("y", "x"), 20 * np.random.rand(6).reshape(2, 3)),
    ...         "precipitation": (("y", "x"), np.random.rand(6).reshape(2, 3)),
    ...     },
    ...     coords={"y": [2, 3], "x": [10, 20, 30]},
    ... )
    >>> x3 = xr.Dataset(
    ...     {
    ...         "temperature": (("y", "x"), 20 * np.random.rand(6).reshape(2, 3)),
    ...         "precipitation": (("y", "x"), np.random.rand(6).reshape(2, 3)),
    ...     },
    ...     coords={"y": [2, 3], "x": [40, 50, 60]},
    ... )

    >>> x1
    <xarray.Dataset> Size: 136B
    Dimensions:        (y: 2, x: 3)
    Coordinates:
      * y              (y) int64 16B 0 1
      * x              (x) int64 24B 10 20 30
    Data variables:
        temperature    (y, x) float64 48B 10.98 14.3 12.06 10.9 8.473 12.92
        precipitation  (y, x) float64 48B 0.4376 0.8918 0.9637 0.3834 0.7917 0.5289

    >>> x2
    <xarray.Dataset> Size: 136B
    Dimensions:        (y: 2, x: 3)
    Coordinates:
      * y              (y) int64 16B 2 3
      * x              (x) int64 24B 10 20 30
    Data variables:
        temperature    (y, x) float64 48B 11.36 18.51 1.421 1.743 0.4044 16.65
        precipitation  (y, x) float64 48B 0.7782 0.87 0.9786 0.7992 0.4615 0.7805

    >>> x3
    <xarray.Dataset> Size: 136B
    Dimensions:        (y: 2, x: 3)
    Coordinates:
      * y              (y) int64 16B 2 3
      * x              (x) int64 24B 40 50 60
    Data variables:
        temperature    (y, x) float64 48B 2.365 12.8 2.867 18.89 10.44 8.293
        precipitation  (y, x) float64 48B 0.2646 0.7742 0.4562 0.5684 0.01879 0.6176

    >>> xr.combine_by_coords([x2, x1])
    <xarray.Dataset> Size: 248B
    Dimensions:        (y: 4, x: 3)
    Coordinates:
      * y              (y) int64 32B 0 1 2 3
      * x              (x) int64 24B 10 20 30
    Data variables:
        temperature    (y, x) float64 96B 10.98 14.3 12.06 ... 1.743 0.4044 16.65
        precipitation  (y, x) float64 96B 0.4376 0.8918 0.9637 ... 0.4615 0.7805

    >>> xr.combine_by_coords([x3, x1])
    <xarray.Dataset> Size: 464B
    Dimensions:        (y: 4, x: 6)
    Coordinates:
      * y              (y) int64 32B 0 1 2 3
      * x              (x) int64 48B 10 20 30 40 50 60
    Data variables:
        temperature    (y, x) float64 192B 10.98 14.3 12.06 ... 18.89 10.44 8.293
        precipitation  (y, x) float64 192B 0.4376 0.8918 0.9637 ... 0.01879 0.6176

    >>> xr.combine_by_coords([x3, x1], join="override")
    <xarray.Dataset> Size: 256B
    Dimensions:        (y: 2, x: 6)
    Coordinates:
      * y              (y) int64 16B 0 1
      * x              (x) int64 48B 10 20 30 40 50 60
    Data variables:
        temperature    (y, x) float64 96B 10.98 14.3 12.06 ... 18.89 10.44 8.293
        precipitation  (y, x) float64 96B 0.4376 0.8918 0.9637 ... 0.01879 0.6176

    >>> xr.combine_by_coords([x1, x2, x3])
    <xarray.Dataset> Size: 464B
    Dimensions:        (y: 4, x: 6)
    Coordinates:
      * y              (y) int64 32B 0 1 2 3
      * x              (x) int64 48B 10 20 30 40 50 60
    Data variables:
        temperature    (y, x) float64 192B 10.98 14.3 12.06 ... 18.89 10.44 8.293
        precipitation  (y, x) float64 192B 0.4376 0.8918 0.9637 ... 0.01879 0.6176

    You can also combine DataArray objects, but the behaviour will differ depending on
    whether or not the DataArrays are named. If all DataArrays are named then they will
    be promoted to Datasets before combining, and then the resultant Dataset will be
    returned, e.g.

    >>> named_da1 = xr.DataArray(
    ...     name="a", data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x"
    ... )
    >>> named_da1
    <xarray.DataArray 'a' (x: 2)> Size: 16B
    array([1., 2.])
    Coordinates:
      * x        (x) int64 16B 0 1

    >>> named_da2 = xr.DataArray(
    ...     name="a", data=[3.0, 4.0], coords={"x": [2, 3]}, dims="x"
    ... )
    >>> named_da2
    <xarray.DataArray 'a' (x: 2)> Size: 16B
    array([3., 4.])
    Coordinates:
      * x        (x) int64 16B 2 3

    >>> xr.combine_by_coords([named_da1, named_da2])
    <xarray.Dataset> Size: 64B
    Dimensions:  (x: 4)
    Coordinates:
      * x        (x) int64 32B 0 1 2 3
    Data variables:
        a        (x) float64 32B 1.0 2.0 3.0 4.0

    If all the DataArrays are unnamed, a single DataArray will be returned, e.g.

    >>> unnamed_da1 = xr.DataArray(data=[1.0, 2.0], coords={"x": [0, 1]}, dims="x")
    >>> unnamed_da2 = xr.DataArray(data=[3.0, 4.0], coords={"x": [2, 3]}, dims="x")
    >>> xr.combine_by_coords([unnamed_da1, unnamed_da2])
    <xarray.DataArray (x: 4)> Size: 32B
    array([1., 2., 3., 4.])
    Coordinates:
      * x        (x) int64 32B 0 1 2 3

    Finally, if you attempt to combine a mix of unnamed DataArrays with either named
    DataArrays or Datasets, a ValueError will be raised (as this is an ambiguous operation).
    """

    if not data_objects:
        return Dataset()

    objs_are_unnamed_dataarrays = [
        isinstance(data_object, DataArray) and data_object.name is None
        for data_object in data_objects
    ]
    if any(objs_are_unnamed_dataarrays):
        if all(objs_are_unnamed_dataarrays):
            # Combine into a single larger DataArray
            temp_datasets = [
                unnamed_dataarray._to_temp_dataset()
                for unnamed_dataarray in data_objects
            ]

            combined_temp_dataset = _combine_single_variable_hypercube(
                temp_datasets,
                fill_value=fill_value,
                data_vars=data_vars,
                coords=coords,
                compat=compat,
                join=join,
                combine_attrs=combine_attrs,
            )
            return DataArray()._from_temp_dataset(combined_temp_dataset)
        else:
            # Must be a mix of unnamed dataarrays with either named dataarrays or with datasets
            # Can't combine these as we wouldn't know whether to merge or concatenate the arrays
            raise ValueError(
                "Can't automatically combine unnamed DataArrays with either named DataArrays or Datasets."
            )
    else:
        # Promote any named DataArrays to single-variable Datasets to simplify combining
        data_objects = [
            obj.to_dataset() if isinstance(obj, DataArray) else obj
            for obj in data_objects
        ]

        # Group by data vars
        sorted_datasets = sorted(data_objects, key=vars_as_keys)
        grouped_by_vars = itertools.groupby(sorted_datasets, key=vars_as_keys)

        # Perform the multidimensional combine on each group of data variables
        # before merging back together
        concatenated_grouped_by_data_vars = tuple(
            _combine_single_variable_hypercube(
                tuple(datasets_with_same_vars),
                fill_value=fill_value,
                data_vars=data_vars,
                coords=coords,
                compat=compat,
                join=join,
                combine_attrs=combine_attrs,
            )
            for vars, datasets_with_same_vars in grouped_by_vars
        )

    return merge(
        concatenated_grouped_by_data_vars,
        compat=compat,
        fill_value=fill_value,
        join=join,
        combine_attrs=combine_attrs,
    )
