from __future__ import annotations

import os
from collections.abc import (
    Callable,
    Iterable,
    Mapping,
    Sequence,
)
from functools import partial
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    TypeVar,
    Union,
    cast,
)

from xarray.backends import plugins
from xarray.backends.common import (
    T_PathFileOrDataStore,
    _find_absolute_paths,
    _normalize_path,
)
from xarray.coders import CFDatetimeCoder, CFTimedeltaCoder
from xarray.core import dtypes, indexing
from xarray.core.coordinates import Coordinates
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.datatree import DataTree
from xarray.core.indexes import Index
from xarray.core.treenode import group_subtrees
from xarray.core.types import ReadBuffer
from xarray.core.utils import emit_user_level_warning, is_remote_uri
from xarray.namedarray.daskmanager import DaskManager
from xarray.namedarray.parallelcompat import guess_chunkmanager
from xarray.namedarray.utils import _get_chunk
from xarray.structure.chunks import _maybe_chunk
from xarray.structure.combine import (
    _infer_concat_order_from_positions,
    _nested_combine,
    combine_by_coords,
)
from xarray.util.deprecation_helpers import (
    _COMPAT_DEFAULT,
    _COORDS_DEFAULT,
    _DATA_VARS_DEFAULT,
    _JOIN_DEFAULT,
    CombineKwargDefault,
)

if TYPE_CHECKING:
    try:
        from dask.delayed import Delayed
    except ImportError:
        Delayed = None  # type: ignore[assignment, misc]

    from xarray.backends.common import BackendEntrypoint
    from xarray.core.types import (
        CombineAttrsOptions,
        CompatOptions,
        ErrorOptionsWithWarn,
        JoinOptions,
        NestedSequence,
        T_Chunks,
    )

    T_NetcdfEngine = Literal["netcdf4", "scipy", "h5netcdf"]
    T_Engine = Union[
        T_NetcdfEngine,
        Literal["pydap", "zarr"],  # noqa: PYI051
        type[BackendEntrypoint],
        str,  # no nice typing support for custom backends
        None,
    ]
    T_NetcdfTypes = Literal[
        "NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_64BIT", "NETCDF3_CLASSIC"
    ]

DATAARRAY_NAME = "__xarray_dataarray_name__"
DATAARRAY_VARIABLE = "__xarray_dataarray_variable__"


def _resolve_decoders_kwargs(decode_cf, open_backend_dataset_parameters, **decoders):
    for d in list(decoders):
        if decode_cf is False and d in open_backend_dataset_parameters:
            decoders[d] = False
        if decoders[d] is None:
            decoders.pop(d)
    return decoders


def _get_mtime(filename_or_obj):
    # if passed an actual file path, augment the token with
    # the file modification time
    mtime = None

    try:
        path = os.fspath(filename_or_obj)
    except TypeError:
        path = None

    if path and not is_remote_uri(path):
        mtime = os.path.getmtime(os.path.expanduser(filename_or_obj))

    return mtime


def _protect_dataset_variables_inplace(dataset: Dataset, cache: bool) -> None:
    for name, variable in dataset.variables.items():
        if name not in dataset._indexes:
            # no need to protect IndexVariable objects
            data: indexing.ExplicitlyIndexedNDArrayMixin
            data = indexing.CopyOnWriteArray(variable._data)
            if cache:
                data = indexing.MemoryCachedArray(data)
            variable.data = data


def _protect_datatree_variables_inplace(tree: DataTree, cache: bool) -> None:
    for node in tree.subtree:
        _protect_dataset_variables_inplace(node.dataset, cache)


def _finalize_store(writes, store):
    """Finalize this store by explicitly syncing and closing"""
    del writes  # ensure writing is done first
    store.close()


def delayed_close_after_writes(writes, store):
    import dask

    return dask.delayed(_finalize_store)(writes, store)


def _multi_file_closer(closers):
    for closer in closers:
        closer()


def load_dataset(filename_or_obj: T_PathFileOrDataStore, **kwargs) -> Dataset:
    """Open, load into memory, and close a Dataset from a file or file-like
    object.

    This is a thin wrapper around :py:meth:`~xarray.open_dataset`. It differs
    from `open_dataset` in that it loads the Dataset into memory, closes the
    file, and returns the Dataset. In contrast, `open_dataset` keeps the file
    handle open and lazy loads its contents. All parameters are passed directly
    to `open_dataset`. See that documentation for further details.

    Returns
    -------
    dataset : Dataset
        The newly created Dataset.

    See Also
    --------
    open_dataset
    """
    if "cache" in kwargs:
        raise TypeError("cache has no effect in this context")

    with open_dataset(filename_or_obj, **kwargs) as ds:
        return ds.load()


def load_dataarray(filename_or_obj: T_PathFileOrDataStore, **kwargs) -> DataArray:
    """Open, load into memory, and close a DataArray from a file or file-like
    object containing a single data variable.

    This is a thin wrapper around :py:meth:`~xarray.open_dataarray`. It differs
    from `open_dataarray` in that it loads the Dataset into memory, closes the
    file, and returns the Dataset. In contrast, `open_dataarray` keeps the file
    handle open and lazy loads its contents. All parameters are passed directly
    to `open_dataarray`. See that documentation for further details.

    Returns
    -------
    datarray : DataArray
        The newly created DataArray.

    See Also
    --------
    open_dataarray
    """
    if "cache" in kwargs:
        raise TypeError("cache has no effect in this context")

    with open_dataarray(filename_or_obj, **kwargs) as da:
        return da.load()


def load_datatree(filename_or_obj: T_PathFileOrDataStore, **kwargs) -> DataTree:
    """Open, load into memory, and close a DataTree from a file or file-like
    object.

    This is a thin wrapper around :py:meth:`~xarray.open_datatree`. It differs
    from `open_datatree` in that it loads the DataTree into memory, closes the
    file, and returns the DataTree. In contrast, `open_datatree` keeps the file
    handle open and lazy loads its contents. All parameters are passed directly
    to `open_datatree`. See that documentation for further details.

    Returns
    -------
    datatree : DataTree
        The newly created DataTree.

    See Also
    --------
    open_datatree
    """
    if "cache" in kwargs:
        raise TypeError("cache has no effect in this context")

    with open_datatree(filename_or_obj, **kwargs) as dt:
        return dt.load()


def _chunk_ds(
    backend_ds,
    filename_or_obj,
    engine,
    chunks,
    overwrite_encoded_chunks,
    inline_array,
    chunked_array_type,
    from_array_kwargs,
    **extra_tokens,
):
    chunkmanager = guess_chunkmanager(chunked_array_type)

    # TODO refactor to move this dask-specific logic inside the DaskManager class
    if isinstance(chunkmanager, DaskManager):
        from dask.base import tokenize

        mtime = _get_mtime(filename_or_obj)
        token = tokenize(filename_or_obj, mtime, engine, chunks, **extra_tokens)
        name_prefix = "open_dataset-"
    else:
        # not used
        token = (None,)
        name_prefix = None

    variables = {}
    for name, var in backend_ds.variables.items():
        if var._in_memory:
            variables[name] = var
            continue
        var_chunks = _get_chunk(
            var._data,
            chunks,
            chunkmanager,
            preferred_chunks=var.encoding.get("preferred_chunks", {}),
            dims=var.dims,
        )
        variables[name] = _maybe_chunk(
            name,
            var,
            var_chunks,
            overwrite_encoded_chunks=overwrite_encoded_chunks,
            name_prefix=name_prefix,
            token=token,
            inline_array=inline_array,
            chunked_array_type=chunkmanager,
            from_array_kwargs=from_array_kwargs.copy(),
            just_use_token=True,
        )
    return backend_ds._replace(variables)


def _maybe_create_default_indexes(ds):
    to_index = {
        name: coord.variable
        for name, coord in ds.coords.items()
        if coord.dims == (name,) and name not in ds.xindexes
    }
    return ds.assign_coords(Coordinates(to_index))


def _dataset_from_backend_dataset(
    backend_ds,
    filename_or_obj,
    engine,
    chunks,
    cache,
    overwrite_encoded_chunks,
    inline_array,
    chunked_array_type,
    from_array_kwargs,
    create_default_indexes,
    **extra_tokens,
):
    if not isinstance(chunks, int | dict) and chunks not in {None, "auto"}:
        raise ValueError(
            f"chunks must be an int, dict, 'auto', or None. Instead found {chunks}."
        )

    _protect_dataset_variables_inplace(backend_ds, cache)

    if create_default_indexes:
        ds = _maybe_create_default_indexes(backend_ds)
    else:
        ds = backend_ds

    if chunks is not None:
        ds = _chunk_ds(
            ds,
            filename_or_obj,
            engine,
            chunks,
            overwrite_encoded_chunks,
            inline_array,
            chunked_array_type,
            from_array_kwargs,
            **extra_tokens,
        )

    ds.set_close(backend_ds._close)

    # Ensure source filename always stored in dataset object
    if "source" not in ds.encoding:
        path = getattr(filename_or_obj, "path", filename_or_obj)

        if isinstance(path, str | os.PathLike):
            ds.encoding["source"] = _normalize_path(path)

    return ds


def _datatree_from_backend_datatree(
    backend_tree,
    filename_or_obj,
    engine,
    chunks,
    cache,
    overwrite_encoded_chunks,
    inline_array,
    chunked_array_type,
    from_array_kwargs,
    create_default_indexes,
    **extra_tokens,
):
    if not isinstance(chunks, int | dict) and chunks not in {None, "auto"}:
        raise ValueError(
            f"chunks must be an int, dict, 'auto', or None. Instead found {chunks}."
        )

    _protect_datatree_variables_inplace(backend_tree, cache)
    if create_default_indexes:
        tree = backend_tree.map_over_datasets(_maybe_create_default_indexes)
    else:
        tree = backend_tree
    if chunks is not None:
        tree = DataTree.from_dict(
            {
                path: _chunk_ds(
                    node.dataset,
                    filename_or_obj,
                    engine,
                    chunks,
                    overwrite_encoded_chunks,
                    inline_array,
                    chunked_array_type,
                    from_array_kwargs,
                    node=path,
                    **extra_tokens,
                )
                for path, [node] in group_subtrees(tree)
            },
            name=tree.name,
        )

    if create_default_indexes or chunks is not None:
        for path, [node] in group_subtrees(backend_tree):
            tree[path].set_close(node._close)

    # Ensure source filename always stored in dataset object
    if "source" not in tree.encoding:
        path = getattr(filename_or_obj, "path", filename_or_obj)

        if isinstance(path, str | os.PathLike):
            tree.encoding["source"] = _normalize_path(path)

    return tree


def open_dataset(
    filename_or_obj: T_PathFileOrDataStore,
    *,
    engine: T_Engine = None,
    chunks: T_Chunks = None,
    cache: bool | None = None,
    decode_cf: bool | None = None,
    mask_and_scale: bool | Mapping[str, bool] | None = None,
    decode_times: (
        bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder] | None
    ) = None,
    decode_timedelta: (
        bool | CFTimedeltaCoder | Mapping[str, bool | CFTimedeltaCoder] | None
    ) = None,
    use_cftime: bool | Mapping[str, bool] | None = None,
    concat_characters: bool | Mapping[str, bool] | None = None,
    decode_coords: Literal["coordinates", "all"] | bool | None = None,
    drop_variables: str | Iterable[str] | None = None,
    create_default_indexes: bool = True,
    inline_array: bool = False,
    chunked_array_type: str | None = None,
    from_array_kwargs: dict[str, Any] | None = None,
    backend_kwargs: dict[str, Any] | None = None,
    **kwargs,
) -> Dataset:
    """Open and decode a dataset from a file or file-like object.

    Parameters
    ----------
    filename_or_obj : str, Path, file-like, bytes, memoryview or DataStore
        Strings and Path objects are interpreted as a path to a netCDF file
        or an OpenDAP URL and opened with python-netCDF4, unless the filename
        ends with .gz, in which case the file is gunzipped and opened with
        scipy.io.netcdf (only netCDF3 supported). Bytes, memoryview and
        file-like objects are opened by scipy.io.netcdf (netCDF3) or h5netcdf
        (netCDF4).
    engine : {"netcdf4", "scipy", "pydap", "h5netcdf", "zarr", None}\
        , installed backend \
        or subclass of xarray.backends.BackendEntrypoint, optional
        Engine to use when reading files. If not provided, the default engine
        is chosen based on available dependencies, by default preferring
        "netcdf4" over "h5netcdf" over "scipy" (customizable via
        ``netcdf_engine_order`` in ``xarray.set_options()``). A custom backend
        class (a subclass of ``BackendEntrypoint``) can also be used.
    chunks : int, dict, 'auto' or None, default: None
        If provided, used to load the data into dask arrays.

        - ``chunks="auto"`` will use dask ``auto`` chunking taking into account the
          engine preferred chunks.
        - ``chunks=None`` skips using dask. This uses xarray's internally private
          :ref:`lazy indexing classes <internal design.lazy indexing>`,
          but data is eagerly loaded into memory as numpy arrays when accessed.
          This can be more efficient for smaller arrays or when large arrays are sliced before computation.
        - ``chunks=-1`` loads the data with dask using a single chunk for all arrays.
        - ``chunks={}`` loads the data with dask using the engine's preferred chunk
          size, generally identical to the format's chunk size. If not available, a
          single chunk for all arrays.

        See dask chunking for more details.
    cache : bool, optional
        If True, cache data loaded from the underlying datastore in memory as
        NumPy arrays when accessed to avoid reading from the underlying data-
        store multiple times. Defaults to True unless you specify the `chunks`
        argument to use dask, in which case it defaults to False. Does not
        change the behavior of coordinates corresponding to dimensions, which
        always load their data from disk into a ``pandas.Index``.
    decode_cf : bool, optional
        Whether to decode these variables, assuming they were saved according
        to CF conventions.
    mask_and_scale : bool or dict-like, optional
        If True, replace array values equal to `_FillValue` with NA and scale
        values according to the formula `original_values * scale_factor +
        add_offset`, where `_FillValue`, `scale_factor` and `add_offset` are
        taken from variable attributes (if they exist).  If the `_FillValue` or
        `missing_value` attribute contains multiple values a warning will be
        issued and all array values matching one of the multiple values will
        be replaced by NA. Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_times : bool, CFDatetimeCoder or dict-like, optional
        If True, decode times encoded in the standard NetCDF datetime format
        into datetime objects. Otherwise, use :py:class:`coders.CFDatetimeCoder` or leave them
        encoded as numbers.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_timedelta : bool, CFTimedeltaCoder, or dict-like, optional
        If True, decode variables and coordinates with time units in
        {"days", "hours", "minutes", "seconds", "milliseconds", "microseconds"}
        into timedelta objects. If False, leave them encoded as numbers.
        If None (default), assume the same value of ``decode_times``; if
        ``decode_times`` is a :py:class:`coders.CFDatetimeCoder` instance, this
        takes the form of a :py:class:`coders.CFTimedeltaCoder` instance with a
        matching ``time_unit``.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    use_cftime: bool or dict-like, optional
        Only relevant if encoded dates come from a standard calendar
        (e.g. "gregorian", "proleptic_gregorian", "standard", or not
        specified).  If None (default), attempt to decode times to
        ``np.datetime64[ns]`` objects; if this is not possible, decode times to
        ``cftime.datetime`` objects. If True, always decode times to
        ``cftime.datetime`` objects, regardless of whether or not they can be
        represented using ``np.datetime64[ns]`` objects.  If False, always
        decode times to ``np.datetime64[ns]`` objects; if this is not possible
        raise an error. Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.

        .. deprecated:: 2025.01.1
           Please pass a :py:class:`coders.CFDatetimeCoder` instance initialized with ``use_cftime`` to the ``decode_times`` kwarg instead.

    concat_characters : bool or dict-like, optional
        If True, concatenate along the last dimension of character arrays to
        form string arrays. Dimensions will only be concatenated over (and
        removed) if they have no corresponding variable and if they are only
        used as the last dimension of character arrays.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_coords : bool or {"coordinates", "all"}, optional
        Controls which variables are set as coordinate variables:

        - "coordinates" or True: Set variables referred to in the
          ``'coordinates'`` attribute of the datasets or individual variables
          as coordinate variables.
        - "all": Set variables referred to in  ``'grid_mapping'``, ``'bounds'`` and
          other attributes as coordinate variables.

        Only existing variables can be set as coordinates. Missing variables
        will be silently ignored.
    drop_variables: str or iterable of str, optional
        A variable or list of variables to exclude from being parsed from the
        dataset. This may be useful to drop variables with problems or
        inconsistent values.
    create_default_indexes : bool, default: True
        If True, create pandas indexes for :term:`dimension coordinates <dimension coordinate>`,
        which loads the coordinate data into memory. Set it to False if you want to avoid loading
        data into memory.

        Note that backends can still choose to create other indexes. If you want to control that,
        please refer to the backend's documentation.
    inline_array: bool, default: False
        How to include the array in the dask task graph.
        By default(``inline_array=False``) the array is included in a task by
        itself, and each chunk refers to that task by its key. With
        ``inline_array=True``, Dask will instead inline the array directly
        in the values of the task graph. See :py:func:`dask.array.from_array`.
    chunked_array_type: str, optional
        Which chunked array type to coerce this datasets' arrays to.
        Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
        Experimental API that should not be relied upon.
    from_array_kwargs: dict
        Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
        chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
        For example if :py:func:`dask.array.Array` objects are used for chunking, additional kwargs will be passed
        to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
    backend_kwargs: dict
        Additional keyword arguments passed on to the engine open function,
        equivalent to `**kwargs`.
    **kwargs: dict
        Additional keyword arguments passed on to the engine open function.
        For example:

        - 'group': path to the netCDF4 group in the given file to open given as
          a str,supported by "netcdf4", "h5netcdf", "zarr".
        - 'lock': resource lock to use when reading data from disk. Only
          relevant when using dask or another form of parallelism. By default,
          appropriate locks are chosen to safely read and write files with the
          currently active dask scheduler. Supported by "netcdf4", "h5netcdf",
          "scipy".

        See engine open function for kwargs accepted by each specific engine.

    Returns
    -------
    dataset : Dataset
        The newly created dataset.

    Notes
    -----
    ``open_dataset`` opens the file with read-only access. When you modify
    values of a Dataset, even one linked to files on disk, only the in-memory
    copy you are manipulating in xarray is modified: the original file on disk
    is never touched.

    See Also
    --------
    open_mfdataset
    """

    if cache is None:
        cache = chunks is None

    if backend_kwargs is not None:
        kwargs.update(backend_kwargs)

    if engine is None:
        engine = plugins.guess_engine(filename_or_obj)

    if from_array_kwargs is None:
        from_array_kwargs = {}

    backend = plugins.get_backend(engine)

    decoders = _resolve_decoders_kwargs(
        decode_cf,
        open_backend_dataset_parameters=backend.open_dataset_parameters,
        mask_and_scale=mask_and_scale,
        decode_times=decode_times,
        decode_timedelta=decode_timedelta,
        concat_characters=concat_characters,
        use_cftime=use_cftime,
        decode_coords=decode_coords,
    )

    overwrite_encoded_chunks = kwargs.pop("overwrite_encoded_chunks", None)
    backend_ds = backend.open_dataset(
        filename_or_obj,
        drop_variables=drop_variables,
        **decoders,
        **kwargs,
    )
    ds = _dataset_from_backend_dataset(
        backend_ds,
        filename_or_obj,
        engine,
        chunks,
        cache,
        overwrite_encoded_chunks,
        inline_array,
        chunked_array_type,
        from_array_kwargs,
        drop_variables=drop_variables,
        create_default_indexes=create_default_indexes,
        **decoders,
        **kwargs,
    )
    return ds


def open_dataarray(
    filename_or_obj: T_PathFileOrDataStore,
    *,
    engine: T_Engine = None,
    chunks: T_Chunks = None,
    cache: bool | None = None,
    decode_cf: bool | None = None,
    mask_and_scale: bool | None = None,
    decode_times: (
        bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder] | None
    ) = None,
    decode_timedelta: bool | CFTimedeltaCoder | None = None,
    use_cftime: bool | None = None,
    concat_characters: bool | None = None,
    decode_coords: Literal["coordinates", "all"] | bool | None = None,
    drop_variables: str | Iterable[str] | None = None,
    create_default_indexes: bool = True,
    inline_array: bool = False,
    chunked_array_type: str | None = None,
    from_array_kwargs: dict[str, Any] | None = None,
    backend_kwargs: dict[str, Any] | None = None,
    **kwargs,
) -> DataArray:
    """Open a DataArray from a file or file-like object containing a single
    data variable.

    This is designed to read netCDF files with only one data variable. If
    multiple variables are present then a ValueError is raised.

    Parameters
    ----------
    filename_or_obj : str, Path, file-like, bytes, memoryview or DataStore
        Strings and Path objects are interpreted as a path to a netCDF file
        or an OpenDAP URL and opened with python-netCDF4, unless the filename
        ends with .gz, in which case the file is gunzipped and opened with
        scipy.io.netcdf (only netCDF3 supported). Bytes, memoryview and
        file-like objects are opened by scipy.io.netcdf (netCDF3) or h5netcdf
        (netCDF4).
    engine : {"netcdf4", "scipy", "pydap", "h5netcdf", "zarr", None}\
        , installed backend \
        or subclass of xarray.backends.BackendEntrypoint, optional
        Engine to use when reading files. If not provided, the default engine
        is chosen based on available dependencies, by default preferring
        "netcdf4" over "h5netcdf" over "scipy" (customizable via
        ``netcdf_engine_order`` in ``xarray.set_options()``). A custom backend
        class (a subclass of ``BackendEntrypoint``) can also be used.
    chunks : int, dict, 'auto' or None, default: None
        If provided, used to load the data into dask arrays.

        - ``chunks='auto'`` will use dask ``auto`` chunking taking into account the
          engine preferred chunks.
        - ``chunks=None`` skips using dask. This uses xarray's internally private
          :ref:`lazy indexing classes <internal design.lazy indexing>`,
          but data is eagerly loaded into memory as numpy arrays when accessed.
          This can be more efficient for smaller arrays, though results may vary.
        - ``chunks=-1`` loads the data with dask using a single chunk for all arrays.
        - ``chunks={}`` loads the data with dask using engine preferred chunks if
          exposed by the backend, otherwise with a single chunk for all arrays.

        See dask chunking for more details.

    cache : bool, optional
        If True, cache data loaded from the underlying datastore in memory as
        NumPy arrays when accessed to avoid reading from the underlying data-
        store multiple times. Defaults to True unless you specify the `chunks`
        argument to use dask, in which case it defaults to False. Does not
        change the behavior of coordinates corresponding to dimensions, which
        always load their data from disk into a ``pandas.Index``.
    decode_cf : bool, optional
        Whether to decode these variables, assuming they were saved according
        to CF conventions.
    mask_and_scale : bool, optional
        If True, replace array values equal to `_FillValue` with NA and scale
        values according to the formula `original_values * scale_factor +
        add_offset`, where `_FillValue`, `scale_factor` and `add_offset` are
        taken from variable attributes (if they exist).  If the `_FillValue` or
        `missing_value` attribute contains multiple values a warning will be
        issued and all array values matching one of the multiple values will
        be replaced by NA. This keyword may not be supported by all the backends.
    decode_times : bool, CFDatetimeCoder or dict-like, optional
        If True, decode times encoded in the standard NetCDF datetime format
        into datetime objects. Otherwise, use :py:class:`coders.CFDatetimeCoder` or
        leave them encoded as numbers.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_timedelta : bool, optional
        If True, decode variables and coordinates with time units in
        {"days", "hours", "minutes", "seconds", "milliseconds", "microseconds"}
        into timedelta objects. If False, leave them encoded as numbers.
        If None (default), assume the same value of ``decode_times``; if
        ``decode_times`` is a :py:class:`coders.CFDatetimeCoder` instance, this
        takes the form of a :py:class:`coders.CFTimedeltaCoder` instance with a
        matching ``time_unit``.
        This keyword may not be supported by all the backends.
    use_cftime: bool, optional
        Only relevant if encoded dates come from a standard calendar
        (e.g. "gregorian", "proleptic_gregorian", "standard", or not
        specified).  If None (default), attempt to decode times to
        ``np.datetime64[ns]`` objects; if this is not possible, decode times to
        ``cftime.datetime`` objects. If True, always decode times to
        ``cftime.datetime`` objects, regardless of whether or not they can be
        represented using ``np.datetime64[ns]`` objects.  If False, always
        decode times to ``np.datetime64[ns]`` objects; if this is not possible
        raise an error. This keyword may not be supported by all the backends.

        .. deprecated:: 2025.01.1
           Please pass a :py:class:`coders.CFDatetimeCoder` instance initialized with ``use_cftime`` to the ``decode_times`` kwarg instead.

    concat_characters : bool, optional
        If True, concatenate along the last dimension of character arrays to
        form string arrays. Dimensions will only be concatenated over (and
        removed) if they have no corresponding variable and if they are only
        used as the last dimension of character arrays.
        This keyword may not be supported by all the backends.
    decode_coords : bool or {"coordinates", "all"}, optional
        Controls which variables are set as coordinate variables:

        - "coordinates" or True: Set variables referred to in the
          ``'coordinates'`` attribute of the datasets or individual variables
          as coordinate variables.
        - "all": Set variables referred to in  ``'grid_mapping'``, ``'bounds'`` and
          other attributes as coordinate variables.

        Only existing variables can be set as coordinates. Missing variables
        will be silently ignored.
    drop_variables: str or iterable of str, optional
        A variable or list of variables to exclude from being parsed from the
        dataset. This may be useful to drop variables with problems or
        inconsistent values.
    create_default_indexes : bool, default: True
        If True, create pandas indexes for :term:`dimension coordinates <dimension coordinate>`,
        which loads the coordinate data into memory. Set it to False if you want to avoid loading
        data into memory.

        Note that backends can still choose to create other indexes. If you want to control that,
        please refer to the backend's documentation.
    inline_array: bool, default: False
        How to include the array in the dask task graph.
        By default(``inline_array=False``) the array is included in a task by
        itself, and each chunk refers to that task by its key. With
        ``inline_array=True``, Dask will instead inline the array directly
        in the values of the task graph. See :py:func:`dask.array.from_array`.
    chunked_array_type: str, optional
        Which chunked array type to coerce the underlying data array to.
        Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
        Experimental API that should not be relied upon.
    from_array_kwargs: dict
        Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
        chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
        For example if :py:func:`dask.array.Array` objects are used for chunking, additional kwargs will be passed
        to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
    backend_kwargs: dict
        Additional keyword arguments passed on to the engine open function,
        equivalent to `**kwargs`.
    **kwargs: dict
        Additional keyword arguments passed on to the engine open function.
        For example:

        - 'group': path to the netCDF4 group in the given file to open given as
          a str,supported by "netcdf4", "h5netcdf", "zarr".
        - 'lock': resource lock to use when reading data from disk. Only
          relevant when using dask or another form of parallelism. By default,
          appropriate locks are chosen to safely read and write files with the
          currently active dask scheduler. Supported by "netcdf4", "h5netcdf",
          "scipy".

        See engine open function for kwargs accepted by each specific engine.

    Notes
    -----
    This is designed to be fully compatible with `DataArray.to_netcdf`. Saving
    using `DataArray.to_netcdf` and then loading with this function will
    produce an identical result.

    All parameters are passed directly to `xarray.open_dataset`. See that
    documentation for further details.

    See Also
    --------
    open_dataset
    """

    dataset = open_dataset(
        filename_or_obj,
        decode_cf=decode_cf,
        mask_and_scale=mask_and_scale,
        decode_times=decode_times,
        concat_characters=concat_characters,
        decode_coords=decode_coords,
        engine=engine,
        chunks=chunks,
        cache=cache,
        drop_variables=drop_variables,
        create_default_indexes=create_default_indexes,
        inline_array=inline_array,
        chunked_array_type=chunked_array_type,
        from_array_kwargs=from_array_kwargs,
        backend_kwargs=backend_kwargs,
        use_cftime=use_cftime,
        decode_timedelta=decode_timedelta,
        **kwargs,
    )

    if len(dataset.data_vars) != 1:
        if len(dataset.data_vars) == 0:
            msg = "Given file dataset contains no data variables."
        else:
            msg = (
                "Given file dataset contains more than one data "
                "variable. Please read with xarray.open_dataset and "
                "then select the variable you want."
            )
        raise ValueError(msg)
    else:
        (data_array,) = dataset.data_vars.values()

    data_array.set_close(dataset._close)

    # Reset names if they were changed during saving
    # to ensure that we can 'roundtrip' perfectly
    if DATAARRAY_NAME in dataset.attrs:
        data_array.name = dataset.attrs[DATAARRAY_NAME]
        del dataset.attrs[DATAARRAY_NAME]

    if data_array.name == DATAARRAY_VARIABLE:
        data_array.name = None

    return data_array


def open_datatree(
    filename_or_obj: T_PathFileOrDataStore,
    *,
    engine: T_Engine = None,
    chunks: T_Chunks = None,
    cache: bool | None = None,
    decode_cf: bool | None = None,
    mask_and_scale: bool | Mapping[str, bool] | None = None,
    decode_times: (
        bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder] | None
    ) = None,
    decode_timedelta: (
        bool | CFTimedeltaCoder | Mapping[str, bool | CFTimedeltaCoder] | None
    ) = None,
    use_cftime: bool | Mapping[str, bool] | None = None,
    concat_characters: bool | Mapping[str, bool] | None = None,
    decode_coords: Literal["coordinates", "all"] | bool | None = None,
    drop_variables: str | Iterable[str] | None = None,
    create_default_indexes: bool = True,
    inline_array: bool = False,
    chunked_array_type: str | None = None,
    from_array_kwargs: dict[str, Any] | None = None,
    backend_kwargs: dict[str, Any] | None = None,
    **kwargs,
) -> DataTree:
    """
    Open and decode a DataTree from a file or file-like object, creating one tree node for each group in the file.

    Parameters
    ----------
    filename_or_obj : str, Path, file-like, bytes or DataStore
        Strings and Path objects are interpreted as a path to a netCDF file or
        Zarr store. Bytes and memoryview objects are interpreted as file
        contents.
    engine : {"netcdf4", "h5netcdf", "zarr", None}, \
             installed backend or xarray.backends.BackendEntrypoint, optional
        Engine to use when reading files. If not provided, the default engine
        is chosen based on available dependencies, by default preferring
        "h5netcdf" over "netcdf4" (customizable via ``netcdf_engine_order`` in
        ``xarray.set_options()``). A custom backend class (a subclass of
        ``BackendEntrypoint``) can also be used.
    chunks : int, dict, 'auto' or None, default: None
        If provided, used to load the data into dask arrays.

        - ``chunks="auto"`` will use dask ``auto`` chunking taking into account the
          engine preferred chunks.
        - ``chunks=None`` skips using dask. This uses xarray's internally private
          :ref:`lazy indexing classes <internal design.lazy indexing>`,
          but data is eagerly loaded into memory as numpy arrays when accessed.
          This can be more efficient for smaller arrays, though results may vary.
        - ``chunks=-1`` loads the data with dask using a single chunk for all arrays.
        - ``chunks={}`` loads the data with dask using the engine's preferred chunk
          size, generally identical to the format's chunk size. If not available, a
          single chunk for all arrays.

        See dask chunking for more details.
    cache : bool, optional
        If True, cache data loaded from the underlying datastore in memory as
        NumPy arrays when accessed to avoid reading from the underlying data-
        store multiple times. Defaults to True unless you specify the `chunks`
        argument to use dask, in which case it defaults to False. Does not
        change the behavior of coordinates corresponding to dimensions, which
        always load their data from disk into a ``pandas.Index``.
    decode_cf : bool, optional
        Whether to decode these variables, assuming they were saved according
        to CF conventions.
    mask_and_scale : bool or dict-like, optional
        If True, replace array values equal to `_FillValue` with NA and scale
        values according to the formula `original_values * scale_factor +
        add_offset`, where `_FillValue`, `scale_factor` and `add_offset` are
        taken from variable attributes (if they exist).  If the `_FillValue` or
        `missing_value` attribute contains multiple values a warning will be
        issued and all array values matching one of the multiple values will
        be replaced by NA. Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_times : bool, CFDatetimeCoder or dict-like, optional
        If True, decode times encoded in the standard NetCDF datetime format
        into datetime objects. Otherwise, use :py:class:`coders.CFDatetimeCoder` or
        leave them encoded as numbers.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_timedelta : bool or dict-like, optional
        If True, decode variables and coordinates with time units in
        {"days", "hours", "minutes", "seconds", "milliseconds", "microseconds"}
        into timedelta objects. If False, leave them encoded as numbers.
        If None (default), assume the same value of ``decode_times``; if
        ``decode_times`` is a :py:class:`coders.CFDatetimeCoder` instance, this
        takes the form of a :py:class:`coders.CFTimedeltaCoder` instance with a
        matching ``time_unit``.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    use_cftime: bool or dict-like, optional
        Only relevant if encoded dates come from a standard calendar
        (e.g. "gregorian", "proleptic_gregorian", "standard", or not
        specified).  If None (default), attempt to decode times to
        ``np.datetime64[ns]`` objects; if this is not possible, decode times to
        ``cftime.datetime`` objects. If True, always decode times to
        ``cftime.datetime`` objects, regardless of whether or not they can be
        represented using ``np.datetime64[ns]`` objects.  If False, always
        decode times to ``np.datetime64[ns]`` objects; if this is not possible
        raise an error. Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.

        .. deprecated:: 2025.01.1
           Please pass a :py:class:`coders.CFDatetimeCoder` instance initialized with ``use_cftime`` to the ``decode_times`` kwarg instead.

    concat_characters : bool or dict-like, optional
        If True, concatenate along the last dimension of character arrays to
        form string arrays. Dimensions will only be concatenated over (and
        removed) if they have no corresponding variable and if they are only
        used as the last dimension of character arrays.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_coords : bool or {"coordinates", "all"}, optional
        Controls which variables are set as coordinate variables:

        - "coordinates" or True: Set variables referred to in the
          ``'coordinates'`` attribute of the datasets or individual variables
          as coordinate variables.
        - "all": Set variables referred to in  ``'grid_mapping'``, ``'bounds'`` and
          other attributes as coordinate variables.

        Only existing variables can be set as coordinates. Missing variables
        will be silently ignored.
    drop_variables: str or iterable of str, optional
        A variable or list of variables to exclude from being parsed from the
        dataset. This may be useful to drop variables with problems or
        inconsistent values.
    create_default_indexes : bool, default: True
        If True, create pandas indexes for :term:`dimension coordinates <dimension coordinate>`,
        which loads the coordinate data into memory. Set it to False if you want to avoid loading
        data into memory.

        Note that backends can still choose to create other indexes. If you want to control that,
        please refer to the backend's documentation.
    inline_array: bool, default: False
        How to include the array in the dask task graph.
        By default(``inline_array=False``) the array is included in a task by
        itself, and each chunk refers to that task by its key. With
        ``inline_array=True``, Dask will instead inline the array directly
        in the values of the task graph. See :py:func:`dask.array.from_array`.
    chunked_array_type: str, optional
        Which chunked array type to coerce this datasets' arrays to.
        Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
        Experimental API that should not be relied upon.
    from_array_kwargs: dict
        Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
        chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
        For example if :py:func:`dask.array.Array` objects are used for chunking, additional kwargs will be passed
        to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
    backend_kwargs: dict
        Additional keyword arguments passed on to the engine open function,
        equivalent to `**kwargs`.
    **kwargs: dict
        Additional keyword arguments passed on to the engine open function.
        For example:

        - 'group': path to the group in the given file to open as the root group as
          a str.
        - 'lock': resource lock to use when reading data from disk. Only
          relevant when using dask or another form of parallelism. By default,
          appropriate locks are chosen to safely read and write files with the
          currently active dask scheduler. Supported by "netcdf4", "h5netcdf",
          "scipy".

        See engine open function for kwargs accepted by each specific engine.

    Returns
    -------
    tree : DataTree
        The newly created datatree.

    Notes
    -----
    ``open_datatree`` opens the file with read-only access. When you modify
    values of a DataTree, even one linked to files on disk, only the in-memory
    copy you are manipulating in xarray is modified: the original file on disk
    is never touched.

    See Also
    --------
    xarray.open_groups
    xarray.open_dataset
    """
    if cache is None:
        cache = chunks is None

    if backend_kwargs is not None:
        kwargs.update(backend_kwargs)

    if engine is None:
        engine = plugins.guess_engine(filename_or_obj, must_support_groups=True)

    if from_array_kwargs is None:
        from_array_kwargs = {}

    backend = plugins.get_backend(engine)

    decoders = _resolve_decoders_kwargs(
        decode_cf,
        open_backend_dataset_parameters=backend.open_dataset_parameters,
        mask_and_scale=mask_and_scale,
        decode_times=decode_times,
        decode_timedelta=decode_timedelta,
        concat_characters=concat_characters,
        use_cftime=use_cftime,
        decode_coords=decode_coords,
    )
    overwrite_encoded_chunks = kwargs.pop("overwrite_encoded_chunks", None)

    backend_tree = backend.open_datatree(
        filename_or_obj,
        drop_variables=drop_variables,
        **decoders,
        **kwargs,
    )

    tree = _datatree_from_backend_datatree(
        backend_tree,
        filename_or_obj,
        engine,
        chunks,
        cache,
        overwrite_encoded_chunks,
        inline_array,
        chunked_array_type,
        from_array_kwargs,
        drop_variables=drop_variables,
        create_default_indexes=create_default_indexes,
        **decoders,
        **kwargs,
    )

    return tree


def open_groups(
    filename_or_obj: T_PathFileOrDataStore,
    *,
    engine: T_Engine = None,
    chunks: T_Chunks = None,
    cache: bool | None = None,
    decode_cf: bool | None = None,
    mask_and_scale: bool | Mapping[str, bool] | None = None,
    decode_times: (
        bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder] | None
    ) = None,
    decode_timedelta: (
        bool | CFTimedeltaCoder | Mapping[str, bool | CFTimedeltaCoder] | None
    ) = None,
    use_cftime: bool | Mapping[str, bool] | None = None,
    concat_characters: bool | Mapping[str, bool] | None = None,
    decode_coords: Literal["coordinates", "all"] | bool | None = None,
    drop_variables: str | Iterable[str] | None = None,
    create_default_indexes: bool = True,
    inline_array: bool = False,
    chunked_array_type: str | None = None,
    from_array_kwargs: dict[str, Any] | None = None,
    backend_kwargs: dict[str, Any] | None = None,
    **kwargs,
) -> dict[str, Dataset]:
    """
    Open and decode a file or file-like object, creating a dictionary containing one xarray Dataset for each group in the file.

    Useful for an HDF file ("netcdf4" or "h5netcdf") containing many groups that are not alignable with their parents
    and cannot be opened directly with ``open_datatree``. It is encouraged to use this function to inspect your data,
    then make the necessary changes to make the structure coercible to a `DataTree` object before calling `DataTree.from_dict()` and proceeding with your analysis.

    Parameters
    ----------
    filename_or_obj : str, Path, file-like, bytes, memoryview or DataStore
        Strings and Path objects are interpreted as a path to a netCDF file or
        Zarr store. Bytes and memoryview objects are interpreted as file
        contents.
    engine : {"netcdf4", "h5netcdf", "zarr", None}, \
             installed backend or xarray.backends.BackendEntrypoint, optional
        Engine to use when reading files. If not provided, the default engine
        is chosen based on available dependencies, by default preferring
        "h5netcdf" over "netcdf4" (customizable via ``netcdf_engine_order`` in
        ``xarray.set_options()``). A custom backend class (a subclass of
        ``BackendEntrypoint``) can also be used.
        can also be used.
    chunks : int, dict, 'auto' or None, default: None
        If provided, used to load the data into dask arrays.

        - ``chunks="auto"`` will use dask ``auto`` chunking taking into account the
          engine preferred chunks.
        - ``chunks=None`` skips using dask. This uses xarray's internally private
          :ref:`lazy indexing classes <internal design.lazy indexing>`,
          but data is eagerly loaded into memory as numpy arrays when accessed.
          This can be more efficient for smaller arrays, though results may vary.
        - ``chunks=-1`` loads the data with dask using a single chunk for all arrays.
        - ``chunks={}`` loads the data with dask using the engine's preferred chunk
          size, generally identical to the format's chunk size. If not available, a
          single chunk for all arrays.

        See dask chunking for more details.
    cache : bool, optional
        If True, cache data loaded from the underlying datastore in memory as
        NumPy arrays when accessed to avoid reading from the underlying data-
        store multiple times. Defaults to True unless you specify the `chunks`
        argument to use dask, in which case it defaults to False. Does not
        change the behavior of coordinates corresponding to dimensions, which
        always load their data from disk into a ``pandas.Index``.
    decode_cf : bool, optional
        Whether to decode these variables, assuming they were saved according
        to CF conventions.
    mask_and_scale : bool or dict-like, optional
        If True, replace array values equal to `_FillValue` with NA and scale
        values according to the formula `original_values * scale_factor +
        add_offset`, where `_FillValue`, `scale_factor` and `add_offset` are
        taken from variable attributes (if they exist).  If the `_FillValue` or
        `missing_value` attribute contains multiple values a warning will be
        issued and all array values matching one of the multiple values will
        be replaced by NA. Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_times : bool, CFDatetimeCoder or dict-like, optional
        If True, decode times encoded in the standard NetCDF datetime format
        into datetime objects. Otherwise, use :py:class:`coders.CFDatetimeCoder` or
        leave them encoded as numbers.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_timedelta : bool or dict-like, optional
        If True, decode variables and coordinates with time units in
        {"days", "hours", "minutes", "seconds", "milliseconds", "microseconds"}
        into timedelta objects. If False, leave them encoded as numbers.
        If None (default), assume the same value of ``decode_times``; if
        ``decode_times`` is a :py:class:`coders.CFDatetimeCoder` instance, this
        takes the form of a :py:class:`coders.CFTimedeltaCoder` instance with a
        matching ``time_unit``.
        This keyword may not be supported by all the backends.
    use_cftime: bool or dict-like, optional
        Only relevant if encoded dates come from a standard calendar
        (e.g. "gregorian", "proleptic_gregorian", "standard", or not
        specified).  If None (default), attempt to decode times to
        ``np.datetime64[ns]`` objects; if this is not possible, decode times to
        ``cftime.datetime`` objects. If True, always decode times to
        ``cftime.datetime`` objects, regardless of whether or not they can be
        represented using ``np.datetime64[ns]`` objects.  If False, always
        decode times to ``np.datetime64[ns]`` objects; if this is not possible
        raise an error. Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.

        .. deprecated:: 2025.01.1
           Please pass a :py:class:`coders.CFDatetimeCoder` instance initialized with ``use_cftime`` to the ``decode_times`` kwarg instead.

    concat_characters : bool or dict-like, optional
        If True, concatenate along the last dimension of character arrays to
        form string arrays. Dimensions will only be concatenated over (and
        removed) if they have no corresponding variable and if they are only
        used as the last dimension of character arrays.
        Pass a mapping, e.g. ``{"my_variable": False}``,
        to toggle this feature per-variable individually.
        This keyword may not be supported by all the backends.
    decode_coords : bool or {"coordinates", "all"}, optional
        Controls which variables are set as coordinate variables:

        - "coordinates" or True: Set variables referred to in the
          ``'coordinates'`` attribute of the datasets or individual variables
          as coordinate variables.
        - "all": Set variables referred to in  ``'grid_mapping'``, ``'bounds'`` and
          other attributes as coordinate variables.

        Only existing variables can be set as coordinates. Missing variables
        will be silently ignored.
    drop_variables: str or iterable of str, optional
        A variable or list of variables to exclude from being parsed from the
        dataset. This may be useful to drop variables with problems or
        inconsistent values.
    create_default_indexes : bool, default: True
        If True, create pandas indexes for :term:`dimension coordinates <dimension coordinate>`,
        which loads the coordinate data into memory. Set it to False if you want to avoid loading
        data into memory.

        Note that backends can still choose to create other indexes. If you want to control that,
        please refer to the backend's documentation.
    inline_array: bool, default: False
        How to include the array in the dask task graph.
        By default(``inline_array=False``) the array is included in a task by
        itself, and each chunk refers to that task by its key. With
        ``inline_array=True``, Dask will instead inline the array directly
        in the values of the task graph. See :py:func:`dask.array.from_array`.
    chunked_array_type: str, optional
        Which chunked array type to coerce this datasets' arrays to.
        Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEnetryPoint` system.
        Experimental API that should not be relied upon.
    from_array_kwargs: dict
        Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
        chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
        For example if :py:func:`dask.array.Array` objects are used for chunking, additional kwargs will be passed
        to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
    backend_kwargs: dict
        Additional keyword arguments passed on to the engine open function,
        equivalent to `**kwargs`.
    **kwargs: dict
        Additional keyword arguments passed on to the engine open function.
        For example:

        - 'group': path to the group in the given file to open as the root group as
          a str.
        - 'lock': resource lock to use when reading data from disk. Only
          relevant when using dask or another form of parallelism. By default,
          appropriate locks are chosen to safely read and write files with the
          currently active dask scheduler. Supported by "netcdf4", "h5netcdf",
          "scipy".

        See engine open function for kwargs accepted by each specific engine.

    Returns
    -------
    groups : dict of str to xarray.Dataset
        The groups as Dataset objects

    Notes
    -----
    ``open_groups`` opens the file with read-only access. When you modify
    values of a Dataset, even one linked to files on disk, only the in-memory
    copy you are manipulating in xarray is modified: the original file on disk
    is never touched.

    See Also
    --------
    xarray.open_datatree
    xarray.open_dataset
    xarray.DataTree.from_dict
    """
    if cache is None:
        cache = chunks is None

    if backend_kwargs is not None:
        kwargs.update(backend_kwargs)

    if engine is None:
        engine = plugins.guess_engine(filename_or_obj, must_support_groups=True)

    if from_array_kwargs is None:
        from_array_kwargs = {}

    backend = plugins.get_backend(engine)

    decoders = _resolve_decoders_kwargs(
        decode_cf,
        open_backend_dataset_parameters=(),
        mask_and_scale=mask_and_scale,
        decode_times=decode_times,
        decode_timedelta=decode_timedelta,
        concat_characters=concat_characters,
        use_cftime=use_cftime,
        decode_coords=decode_coords,
    )
    overwrite_encoded_chunks = kwargs.pop("overwrite_encoded_chunks", None)

    backend_groups = backend.open_groups_as_dict(
        filename_or_obj,
        drop_variables=drop_variables,
        **decoders,
        **kwargs,
    )

    groups = {
        name: _dataset_from_backend_dataset(
            backend_ds,
            filename_or_obj,
            engine,
            chunks,
            cache,
            overwrite_encoded_chunks,
            inline_array,
            chunked_array_type,
            from_array_kwargs,
            drop_variables=drop_variables,
            create_default_indexes=create_default_indexes,
            **decoders,
            **kwargs,
        )
        for name, backend_ds in backend_groups.items()
    }

    return groups


_FLike = TypeVar("_FLike", bound=Union[str, ReadBuffer])


def _remove_path(
    paths: NestedSequence[_FLike], paths_to_remove: set[_FLike]
) -> NestedSequence[_FLike]:
    # Initialize an empty list to store the result
    result: list[Union[_FLike, NestedSequence[_FLike]]] = []

    for item in paths:
        if isinstance(item, list):
            # If the current item is a list, recursively call remove_elements on it
            nested_result = _remove_path(item, paths_to_remove)
            if nested_result:  # Only add non-empty lists to avoid adding empty lists
                result.append(nested_result)
        elif item not in paths_to_remove:
            # Add the item to the result if it is not in the set of elements to remove
            result.append(item)

    return result


def open_mfdataset(
    paths: (
        str | os.PathLike | ReadBuffer | NestedSequence[str | os.PathLike | ReadBuffer]
    ),
    chunks: T_Chunks = None,
    concat_dim: (
        str
        | DataArray
        | Index
        | Sequence[str]
        | Sequence[DataArray]
        | Sequence[Index]
        | None
    ) = None,
    compat: CompatOptions | CombineKwargDefault = _COMPAT_DEFAULT,
    preprocess: Callable[[Dataset], Dataset] | None = None,
    engine: T_Engine = None,
    data_vars: (
        Literal["all", "minimal", "different"] | None | list[str] | CombineKwargDefault
    ) = _DATA_VARS_DEFAULT,
    coords=_COORDS_DEFAULT,
    combine: Literal["by_coords", "nested"] = "by_coords",
    parallel: bool = False,
    join: JoinOptions | CombineKwargDefault = _JOIN_DEFAULT,
    attrs_file: str | os.PathLike | None = None,
    combine_attrs: CombineAttrsOptions = "override",
    errors: ErrorOptionsWithWarn = "raise",
    **kwargs,
) -> Dataset:
    """Open multiple files as a single dataset.

    If combine='by_coords' then the function ``combine_by_coords`` is used to combine
    the datasets into one before returning the result, and if combine='nested' then
    ``combine_nested`` is used. The filepaths must be structured according to which
    combining function is used, the details of which are given in the documentation for
    ``combine_by_coords`` and ``combine_nested``. By default ``combine='by_coords'``
    will be used. Requires dask to be installed. See documentation for
    details on dask [1]_. Global attributes from the ``attrs_file`` are used
    for the combined dataset.

    Parameters
    ----------
    paths : str or nested sequence of paths
        Either a string glob in the form ``"path/to/my/files/*.nc"`` or an explicit list of
        files to open. Paths can be given as strings or as pathlib Paths. If
        concatenation along more than one dimension is desired, then ``paths`` must be a
        nested list-of-lists (see ``combine_nested`` for details). (A string glob will
        be expanded to a 1-dimensional list.)
    chunks : int, dict, 'auto' or None, optional
        Dictionary with keys given by dimension names and values given by chunk sizes.
        In general, these should divide the dimensions of each dataset. If int, chunk
        each dimension by ``chunks``. By default, chunks will be chosen to match the
        chunks on disk. This may impact performance: please see the full documentation
        for more details [2]_. This argument is evaluated on a per-file basis, so chunk
        sizes that span multiple files will be ignored.
    concat_dim : str, DataArray, Index or a Sequence of these or None, optional
        Dimensions to concatenate files along.  You only need to provide this argument
        if ``combine='nested'``, and if any of the dimensions along which you want to
        concatenate is not a dimension in the original datasets, e.g., if you want to
        stack a collection of 2D arrays along a third dimension. Set
        ``concat_dim=[..., None, ...]`` explicitly to disable concatenation along a
        particular dimension. Default is None, which for a 1D list of filepaths is
        equivalent to opening the files separately and then merging them with
        ``xarray.merge``.
    combine : {"by_coords", "nested"}, optional
        Whether ``xarray.combine_by_coords`` or ``xarray.combine_nested`` is used to
        combine all the data. Default is to use ``xarray.combine_by_coords``.
    compat : {"identical", "equals", "broadcast_equals", \
              "no_conflicts", "override"}, default: "no_conflicts"
        String indicating how to compare variables of the same name for
        potential conflicts when merging:

         * "broadcast_equals": all values must be equal when variables are
           broadcast against each other to ensure common dimensions.
         * "equals": all values and dimensions must be the same.
         * "identical": all values, dimensions and attributes must be the
           same.
         * "no_conflicts": only values which are not null in both datasets
           must be equal. The returned dataset then contains the combination
           of all non-null values.
         * "override": skip comparing and pick variable from first dataset

    preprocess : callable, optional
        If provided, call this function on each dataset prior to concatenation.
        You can find the file-name from which each dataset was loaded in
        ``ds.encoding["source"]``.
    engine : {"netcdf4", "scipy", "pydap", "h5netcdf", "zarr", None}\
        , installed backend \
        or subclass of xarray.backends.BackendEntrypoint, optional
        Engine to use when reading files. If not provided, the default engine
        is chosen based on available dependencies, by default preferring
        "netcdf4" over "h5netcdf" over "scipy" (customizable via
        ``netcdf_engine_order`` in ``xarray.set_options()``). A custom backend
        class (a subclass of ``BackendEntrypoint``) can also be used.
    data_vars : {"minimal", "different", "all", None} or list of str, default: "all"
        These data variables will be concatenated together:
          * "minimal": Only data variables in which the dimension already
            appears are included.
          * "different": Data variables which are not equal (ignoring
            attributes) across all datasets are also concatenated (as well as
            all for which dimension already appears). Beware: this option may
            load the data payload of data variables into memory if they are not
            already loaded.
          * "all": All data variables will be concatenated.
          * None: Means ``"all"`` if ``concat_dim`` is not present in any of
            the ``objs``, and ``"minimal"`` if ``concat_dim`` is present
            in any of ``objs``.
          * list of str: The listed data variables will be concatenated, in
            addition to the "minimal" data variables.
    coords : {"minimal", "different", "all"} or list of str, default: "different"
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
         * list of str: The listed coordinate variables will be concatenated,
           in addition the "minimal" coordinates.
    parallel : bool, default: False
        If True, the open and preprocess steps of this function will be
        performed in parallel using ``dask.delayed``. Default is False.
    join : {"outer", "inner", "left", "right", "exact", "override"}, default: "outer"
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
    attrs_file : str or path-like, optional
        Path of the file used to read global attributes from.
        By default global attributes are read from the first file provided,
        with wildcard matches sorted by filename.
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
    errors : {"raise", "warn", "ignore"}, default: "raise"
        String indicating how to handle errors in opening dataset.

        - "raise": invalid dataset will raise an exception.
        - "warn": a warning will be issued for each invalid dataset.
        - "ignore": invalid dataset will be ignored.
    **kwargs : optional
        Additional arguments passed on to :py:func:`xarray.open_dataset`. For an
        overview of some of the possible options, see the documentation of
        :py:func:`xarray.open_dataset`

    Returns
    -------
    xarray.Dataset

    Notes
    -----
    ``open_mfdataset`` opens files with read-only access. When you modify values
    of a Dataset, even one linked to files on disk, only the in-memory copy you
    are manipulating in xarray is modified: the original file on disk is never
    touched.

    See Also
    --------
    combine_by_coords
    combine_nested
    open_dataset

    Examples
    --------
    A user might want to pass additional arguments into ``preprocess`` when
    applying some operation to many individual files that are being opened. One route
    to do this is through the use of ``functools.partial``.

    >>> from functools import partial
    >>> def _preprocess(x, lon_bnds, lat_bnds):
    ...     return x.sel(lon=slice(*lon_bnds), lat=slice(*lat_bnds))
    ...
    >>> lon_bnds, lat_bnds = (-110, -105), (40, 45)
    >>> partial_func = partial(_preprocess, lon_bnds=lon_bnds, lat_bnds=lat_bnds)
    >>> ds = xr.open_mfdataset(
    ...     "file_*.nc", concat_dim="time", preprocess=partial_func
    ... )  # doctest: +SKIP

    It is also possible to use any argument to ``open_dataset`` together
    with ``open_mfdataset``, such as for example ``drop_variables``:

    >>> ds = xr.open_mfdataset(
    ...     "file.nc", drop_variables=["varname_1", "varname_2"]  # any list of vars
    ... )  # doctest: +SKIP

    References
    ----------

    .. [1] https://docs.xarray.dev/en/stable/dask.html
    .. [2] https://docs.xarray.dev/en/stable/dask.html#chunking-and-performance
    """
    paths = _find_absolute_paths(paths, engine=engine, **kwargs)

    if not paths:
        raise OSError("no files to open")

    paths1d: list[str | ReadBuffer]
    if combine == "nested":
        if isinstance(concat_dim, str | DataArray) or concat_dim is None:
            concat_dim = [concat_dim]  # type: ignore[assignment]

        # This creates a flat list which is easier to iterate over, whilst
        # encoding the originally-supplied structure as "ids".
        # The "ids" are not used at all if combine='by_coords`.
        combined_ids_paths = _infer_concat_order_from_positions(paths)
        ids, paths1d = (
            list(combined_ids_paths.keys()),
            list(combined_ids_paths.values()),
        )
    elif concat_dim is not None:
        raise ValueError(
            "When combine='by_coords', passing a value for `concat_dim` has no "
            "effect. To manually combine along a specific dimension you should "
            "instead specify combine='nested' along with a value for `concat_dim`.",
        )
    else:
        paths1d = paths  # type: ignore[assignment]

    open_kwargs = dict(engine=engine, chunks=chunks or {}, **kwargs)

    if parallel:
        import dask

        # wrap the open_dataset, getattr, and preprocess with delayed
        open_ = dask.delayed(open_dataset)
        getattr_ = dask.delayed(getattr)
        if preprocess is not None:
            preprocess = dask.delayed(preprocess)
    else:
        open_ = open_dataset
        getattr_ = getattr

    if errors not in ("raise", "warn", "ignore"):
        raise ValueError(
            f"'errors' must be 'raise', 'warn' or 'ignore', got '{errors}'"
        )

    datasets = []
    invalid_paths = set()
    for p in paths1d:
        try:
            ds = open_(p, **open_kwargs)
            datasets.append(ds)
        except Exception as e:
            if errors == "raise":
                raise
            elif errors == "warn":
                emit_user_level_warning(f"Could not open {p} due to {e}. Ignoring.")
            # remove invalid paths
            invalid_paths.add(p)

    if invalid_paths:
        paths = _remove_path(paths, invalid_paths)
        if combine == "nested":
            # Create new ids and paths based on removed items
            combined_ids_paths = _infer_concat_order_from_positions(paths)
            ids = list(combined_ids_paths.keys())

    closers = [getattr_(ds, "_close") for ds in datasets]
    if preprocess is not None:
        datasets = [preprocess(ds) for ds in datasets]

    if parallel:
        # calling compute here will return the datasets/file_objs lists,
        # the underlying datasets will still be stored as dask arrays
        datasets, closers = dask.compute(datasets, closers)

    # Combine all datasets, closing them in case of a ValueError
    try:
        if combine == "nested":
            # Combined nested list by successive concat and merge operations
            # along each dimension, using structure given by "ids"
            combined = _nested_combine(
                datasets,
                concat_dims=concat_dim,
                compat=compat,
                data_vars=data_vars,
                coords=coords,
                ids=ids,
                join=join,
                combine_attrs=combine_attrs,
                fill_value=dtypes.NA,
            )
        elif combine == "by_coords":
            # Redo ordering from coordinates, ignoring how they were ordered
            # previously
            combined = combine_by_coords(
                datasets,
                compat=compat,
                data_vars=data_vars,
                coords=coords,
                join=join,
                combine_attrs=combine_attrs,
            )
        else:
            raise ValueError(
                f"{combine} is an invalid option for the keyword argument ``combine``"
            )
    except ValueError:
        for ds in datasets:
            ds.close()
        raise

    combined.set_close(partial(_multi_file_closer, closers))

    # read global attributes from the attrs_file or from the first dataset
    if attrs_file is not None:
        if isinstance(attrs_file, os.PathLike):
            attrs_file = cast(str, os.fspath(attrs_file))
        combined.attrs = datasets[paths1d.index(attrs_file)].attrs

    return combined
