from __future__ import annotations

import importlib.util
import os
from collections.abc import (
    Callable,
    Hashable,
    Iterable,
    Mapping,
    MutableMapping,
    Sequence,
)
from functools import partial
from io import IOBase
from itertools import starmap
from numbers import Number
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    TypeVar,
    Union,
    cast,
    overload,
)

import numpy as np

from xarray import backends, conventions
from xarray.backends import plugins
from xarray.backends.common import (
    AbstractDataStore,
    AbstractWritableDataStore,
    ArrayWriter,
    BytesIOProxy,
    T_PathFileOrDataStore,
    _find_absolute_paths,
    _normalize_path,
)
from xarray.backends.locks import get_dask_scheduler
from xarray.coders import CFDatetimeCoder, CFTimedeltaCoder
from xarray.core import dtypes, indexing
from xarray.core.coordinates import Coordinates
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.datatree import DataTree
from xarray.core.indexes import Index
from xarray.core.treenode import group_subtrees
from xarray.core.types import NetcdfWriteModes, ReadBuffer, ZarrWriteModes
from xarray.core.utils import emit_user_level_warning, is_remote_uri
from xarray.namedarray.daskmanager import DaskManager
from xarray.namedarray.parallelcompat import guess_chunkmanager
from xarray.structure.chunks import _get_chunk, _maybe_chunk
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
        ReadBuffer,
        T_Chunks,
        ZarrStoreLike,
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


def get_default_netcdf_write_engine(
    format: T_NetcdfTypes | None,
    to_fileobject_or_memoryview: bool,
) -> Literal["netcdf4", "h5netcdf", "scipy"]:
    """Return the default netCDF library to use for writing a netCDF file."""
    module_names = {
        "netcdf4": "netCDF4",
        "scipy": "scipy",
        "h5netcdf": "h5netcdf",
    }

    candidates = list(plugins.NETCDF_BACKENDS_ORDER)

    if format is not None:
        if format.upper().startswith("NETCDF3"):
            candidates.remove("h5netcdf")
        elif format.upper().startswith("NETCDF4"):
            candidates.remove("scipy")
        else:
            raise ValueError(f"unexpected {format=}")

    if to_fileobject_or_memoryview:
        candidates.remove("netcdf4")

    for engine in candidates:
        module_name = module_names[engine]
        if importlib.util.find_spec(module_name) is not None:
            return cast(Literal["netcdf4", "h5netcdf", "scipy"], engine)

    format_str = f" with {format=}" if format is not None else ""
    libraries = ", ".join(module_names[c] for c in candidates)
    raise ValueError(
        f"cannot write NetCDF files{format_str} because none of the suitable "
        f"backend libraries ({libraries}) are installed"
    )


def _validate_dataset_names(dataset: Dataset) -> None:
    """DataArray.name and Dataset keys must be a string or None"""

    def check_name(name: Hashable):
        if isinstance(name, str):
            if not name:
                raise ValueError(
                    f"Invalid name {name!r} for DataArray or Dataset key: "
                    "string must be length 1 or greater for "
                    "serialization to netCDF or zarr files"
                )
        elif name is not None:
            raise TypeError(
                f"Invalid name {name!r} for DataArray or Dataset key: "
                "must be either a string or None for serialization to netCDF "
                "or zarr files"
            )

    for k in dataset.variables:
        check_name(k)


def _validate_attrs(dataset, engine, invalid_netcdf=False):
    """`attrs` must have a string key and a value which is either: a number,
    a string, an ndarray, a list/tuple of numbers/strings, or a numpy.bool_.

    Notes
    -----
    A numpy.bool_ is only allowed when using the h5netcdf engine with
    `invalid_netcdf=True`.
    """

    valid_types = (str, Number, np.ndarray, np.number, list, tuple, bytes)
    if invalid_netcdf and engine == "h5netcdf":
        valid_types += (np.bool_,)

    def check_attr(name, value, valid_types):
        if isinstance(name, str):
            if not name:
                raise ValueError(
                    f"Invalid name for attr {name!r}: string must be "
                    "length 1 or greater for serialization to "
                    "netCDF files"
                )
        else:
            raise TypeError(
                f"Invalid name for attr: {name!r} must be a string for "
                "serialization to netCDF files"
            )

        if not isinstance(value, valid_types):
            raise TypeError(
                f"Invalid value for attr {name!r}: {value!r}. For serialization to "
                "netCDF files, its value must be of one of the following types: "
                f"{', '.join([vtype.__name__ for vtype in valid_types])}"
            )

        if isinstance(value, bytes) and engine == "h5netcdf":
            try:
                value.decode("utf-8")
            except UnicodeDecodeError as e:
                raise ValueError(
                    f"Invalid value provided for attribute '{name!r}': {value!r}. "
                    "Only binary data derived from UTF-8 encoded strings is allowed "
                    f"for the '{engine}' engine. Consider using the 'netcdf4' engine."
                ) from e

            if b"\x00" in value:
                raise ValueError(
                    f"Invalid value provided for attribute '{name!r}': {value!r}. "
                    f"Null characters are not permitted for the '{engine}' engine. "
                    "Consider using the 'netcdf4' engine."
                )

    # Check attrs on the dataset itself
    for k, v in dataset.attrs.items():
        check_attr(k, v, valid_types)

    # Check attrs on each variable within the dataset
    for variable in dataset.variables.values():
        for k, v in variable.attrs.items():
            check_attr(k, v, valid_types)


def _sanitize_unlimited_dims(dataset, unlimited_dims):
    msg_origin = "unlimited_dims-kwarg"
    if unlimited_dims is None:
        unlimited_dims = dataset.encoding.get("unlimited_dims", None)
        msg_origin = "dataset.encoding"
    if unlimited_dims is not None:
        if isinstance(unlimited_dims, str) or not isinstance(unlimited_dims, Iterable):
            unlimited_dims = [unlimited_dims]
        else:
            unlimited_dims = list(unlimited_dims)
        dataset_dims = set(dataset.dims)
        unlimited_dims = set(unlimited_dims)
        if undeclared_dims := (unlimited_dims - dataset_dims):
            msg = (
                f"Unlimited dimension(s) {undeclared_dims!r} declared in {msg_origin!r}, "
                f"but not part of current dataset dimensions. "
                f"Consider removing {undeclared_dims!r} from {msg_origin!r}."
            )
            if msg_origin == "unlimited_dims-kwarg":
                raise ValueError(msg)
            else:
                emit_user_level_warning(msg)
        return unlimited_dims


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
        var_chunks = _get_chunk(var, chunks, chunkmanager)
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
    decode_times: bool
    | CFDatetimeCoder
    | Mapping[str, bool | CFDatetimeCoder]
    | None = None,
    decode_timedelta: bool
    | CFTimedeltaCoder
    | Mapping[str, bool | CFTimedeltaCoder]
    | None = None,
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
        is chosen based on available dependencies, with a preference for
        "netcdf4". A custom backend class (a subclass of ``BackendEntrypoint``)
        can also be used.
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
    decode_times: bool
    | CFDatetimeCoder
    | Mapping[str, bool | CFDatetimeCoder]
    | None = None,
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
    """Open an DataArray from a file or file-like object containing a single
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
        is chosen based on available dependencies, with a preference for
        "netcdf4".
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

    See also
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
    decode_times: bool
    | CFDatetimeCoder
    | Mapping[str, bool | CFDatetimeCoder]
    | None = None,
    decode_timedelta: bool
    | CFTimedeltaCoder
    | Mapping[str, bool | CFTimedeltaCoder]
    | None = None,
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
        is chosen based on available dependencies, with a preference for
        "netcdf4". A custom backend class (a subclass of ``BackendEntrypoint``)
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
    decode_times: bool
    | CFDatetimeCoder
    | Mapping[str, bool | CFDatetimeCoder]
    | None = None,
    decode_timedelta: bool
    | CFTimedeltaCoder
    | Mapping[str, bool | CFTimedeltaCoder]
    | None = None,
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
        is chosen based on available dependencies, with a preference for
        "netcdf4". A custom backend class (a subclass of ``BackendEntrypoint``)
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
        engine = plugins.guess_engine(filename_or_obj)

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
    paths: str
    | os.PathLike
    | ReadBuffer
    | NestedSequence[str | os.PathLike | ReadBuffer],
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
    data_vars: Literal["all", "minimal", "different"]
    | None
    | list[str]
    | CombineKwargDefault = _DATA_VARS_DEFAULT,
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
        each dimension by ``chunks``. By default, chunks will be chosen to load entire
        input files into memory at once. This has a major impact on performance: please
        see the full documentation for more details [2]_. This argument is evaluated
        on a per-file basis, so chunk sizes that span multiple files will be ignored.
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
        is chosen based on available dependencies, with a preference for
        "netcdf4".
    data_vars : {"minimal", "different", "all"} or list of str, default: "all"
        These data variables will be concatenated together:
          * "minimal": Only data variables in which the dimension already
            appears are included.
          * "different": Data variables which are not equal (ignoring
            attributes) across all datasets are also concatenated (as well as
            all for which dimension already appears). Beware: this option may
            load the data payload of data variables into memory if they are not
            already loaded.
          * "all": All data variables will be concatenated.
          * list of str: The listed data variables will be concatenated, in
            addition to the "minimal" data variables.
    coords : {"minimal", "different", "all"} or list of str, optional
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


WRITEABLE_STORES: dict[T_NetcdfEngine, Callable] = {
    "netcdf4": backends.NetCDF4DataStore.open,
    "scipy": backends.ScipyDataStore,
    "h5netcdf": backends.H5NetCDFStore.open,
}


def get_writable_netcdf_store(
    target,
    engine: T_NetcdfEngine,
    *,
    format: T_NetcdfTypes | None,
    mode: NetcdfWriteModes,
    autoclose: bool,
    invalid_netcdf: bool,
    auto_complex: bool | None,
) -> AbstractWritableDataStore:
    """Create a store for writing to a netCDF file."""
    try:
        store_open = WRITEABLE_STORES[engine]
    except KeyError as err:
        raise ValueError(f"unrecognized engine for to_netcdf: {engine!r}") from err

    if format is not None:
        format = format.upper()  # type: ignore[assignment]

    kwargs = dict(autoclose=True) if autoclose else {}
    if invalid_netcdf:
        if engine == "h5netcdf":
            kwargs["invalid_netcdf"] = invalid_netcdf
        else:
            raise ValueError(
                f"unrecognized option 'invalid_netcdf' for engine {engine}"
            )
    if auto_complex is not None:
        kwargs["auto_complex"] = auto_complex

    return store_open(target, mode=mode, format=format, **kwargs)


# multifile=True returns writer and datastore
@overload
def to_netcdf(
    dataset: Dataset,
    path_or_file: str | os.PathLike | None = None,
    mode: NetcdfWriteModes = "w",
    format: T_NetcdfTypes | None = None,
    group: str | None = None,
    engine: T_NetcdfEngine | None = None,
    encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
    unlimited_dims: Iterable[Hashable] | None = None,
    compute: bool = True,
    *,
    multifile: Literal[True],
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> tuple[ArrayWriter, AbstractDataStore]: ...


# path=None writes to bytes or memoryview, depending on store
@overload
def to_netcdf(
    dataset: Dataset,
    path_or_file: None = None,
    mode: NetcdfWriteModes = "w",
    format: T_NetcdfTypes | None = None,
    group: str | None = None,
    engine: T_NetcdfEngine | None = None,
    encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
    unlimited_dims: Iterable[Hashable] | None = None,
    compute: bool = True,
    multifile: Literal[False] = False,
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> memoryview: ...


# compute=False returns dask.Delayed
@overload
def to_netcdf(
    dataset: Dataset,
    path_or_file: str | os.PathLike,
    mode: NetcdfWriteModes = "w",
    format: T_NetcdfTypes | None = None,
    group: str | None = None,
    engine: T_NetcdfEngine | None = None,
    encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
    unlimited_dims: Iterable[Hashable] | None = None,
    *,
    compute: Literal[False],
    multifile: Literal[False] = False,
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> Delayed: ...


# default return None
@overload
def to_netcdf(
    dataset: Dataset,
    path_or_file: str | os.PathLike | IOBase,
    mode: NetcdfWriteModes = "w",
    format: T_NetcdfTypes | None = None,
    group: str | None = None,
    engine: T_NetcdfEngine | None = None,
    encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
    unlimited_dims: Iterable[Hashable] | None = None,
    compute: Literal[True] = True,
    multifile: Literal[False] = False,
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> None: ...


# if compute cannot be evaluated at type check time
# we may get back either Delayed or None
@overload
def to_netcdf(
    dataset: Dataset,
    path_or_file: str | os.PathLike,
    mode: NetcdfWriteModes = "w",
    format: T_NetcdfTypes | None = None,
    group: str | None = None,
    engine: T_NetcdfEngine | None = None,
    encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
    unlimited_dims: Iterable[Hashable] | None = None,
    compute: bool = False,
    multifile: Literal[False] = False,
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> Delayed | None: ...


# if multifile cannot be evaluated at type check time
# we may get back either writer and datastore or Delayed or None
@overload
def to_netcdf(
    dataset: Dataset,
    path_or_file: str | os.PathLike,
    mode: NetcdfWriteModes = "w",
    format: T_NetcdfTypes | None = None,
    group: str | None = None,
    engine: T_NetcdfEngine | None = None,
    encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
    unlimited_dims: Iterable[Hashable] | None = None,
    compute: bool = False,
    multifile: bool = False,
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> tuple[ArrayWriter, AbstractDataStore] | Delayed | None: ...


# Any
@overload
def to_netcdf(
    dataset: Dataset,
    path_or_file: str | os.PathLike | IOBase | None,
    mode: NetcdfWriteModes = "w",
    format: T_NetcdfTypes | None = None,
    group: str | None = None,
    engine: T_NetcdfEngine | None = None,
    encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
    unlimited_dims: Iterable[Hashable] | None = None,
    compute: bool = False,
    multifile: bool = False,
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> tuple[ArrayWriter, AbstractDataStore] | memoryview | Delayed | None: ...


def to_netcdf(
    dataset: Dataset,
    path_or_file: str | os.PathLike | IOBase | None = None,
    mode: NetcdfWriteModes = "w",
    format: T_NetcdfTypes | None = None,
    group: str | None = None,
    engine: T_NetcdfEngine | None = None,
    encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
    unlimited_dims: Iterable[Hashable] | None = None,
    compute: bool = True,
    multifile: bool = False,
    invalid_netcdf: bool = False,
    auto_complex: bool | None = None,
) -> tuple[ArrayWriter, AbstractDataStore] | memoryview | Delayed | None:
    """This function creates an appropriate datastore for writing a dataset to
    disk as a netCDF file

    See `Dataset.to_netcdf` for full API docs.

    The ``multifile`` argument is only for the private use of save_mfdataset.
    """
    if encoding is None:
        encoding = {}

    path_or_file = _normalize_path(path_or_file)

    if engine is None:
        to_fileobject_or_memoryview = not isinstance(path_or_file, str)
        engine = get_default_netcdf_write_engine(format, to_fileobject_or_memoryview)

    # validate Dataset keys, DataArray names, and attr keys/values
    _validate_dataset_names(dataset)
    _validate_attrs(dataset, engine, invalid_netcdf)
    # sanitize unlimited_dims
    unlimited_dims = _sanitize_unlimited_dims(dataset, unlimited_dims)

    # handle scheduler specific logic
    scheduler = get_dask_scheduler()
    have_chunks = any(v.chunks is not None for v in dataset.variables.values())

    autoclose = have_chunks and scheduler in ["distributed", "multiprocessing"]
    if autoclose and engine == "scipy":
        raise NotImplementedError(
            f"Writing netCDF files with the {engine} backend "
            f"is not currently supported with dask's {scheduler} scheduler"
        )

    if path_or_file is None:
        if not compute:
            raise NotImplementedError(
                "to_netcdf() with compute=False is not yet implemented when "
                "returning a memoryview"
            )
        target = BytesIOProxy()
    else:
        target = path_or_file  # type: ignore[assignment]

    store = get_writable_netcdf_store(
        target,
        engine,
        mode=mode,
        format=format,
        autoclose=autoclose,
        invalid_netcdf=invalid_netcdf,
        auto_complex=auto_complex,
    )
    if group is not None:
        store = store.get_child_store(group)

    writer = ArrayWriter()

    # TODO: figure out how to refactor this logic (here and in save_mfdataset)
    # to avoid this mess of conditionals
    try:
        # TODO: allow this work (setting up the file for writing array data)
        # to be parallelized with dask
        dump_to_store(
            dataset, store, writer, encoding=encoding, unlimited_dims=unlimited_dims
        )
        if autoclose:
            store.close()

        if multifile:
            return writer, store

        writes = writer.sync(compute=compute)

    finally:
        if not multifile:
            if compute:
                store.close()
            else:
                store.sync()

    if path_or_file is None:
        assert isinstance(target, BytesIOProxy)  # created in this function
        return target.getbuffer()

    if not compute:
        return delayed_close_after_writes(writes, store)

    return None


def dump_to_store(
    dataset, store, writer=None, encoder=None, encoding=None, unlimited_dims=None
):
    """Store dataset contents to a backends.*DataStore object."""
    if writer is None:
        writer = ArrayWriter()

    if encoding is None:
        encoding = {}

    variables, attrs = conventions.encode_dataset_coordinates(dataset)

    check_encoding = set()
    for k, enc in encoding.items():
        # no need to shallow copy the variable again; that already happened
        # in encode_dataset_coordinates
        variables[k].encoding = enc
        check_encoding.add(k)

    if encoder:
        variables, attrs = encoder(variables, attrs)

    store.store(variables, attrs, check_encoding, writer, unlimited_dims=unlimited_dims)


def save_mfdataset(
    datasets,
    paths,
    mode="w",
    format=None,
    groups=None,
    engine=None,
    compute=True,
    **kwargs,
):
    """Write multiple datasets to disk as netCDF files simultaneously.

    This function is intended for use with datasets consisting of dask.array
    objects, in which case it can write the multiple datasets to disk
    simultaneously using a shared thread pool.

    When not using dask, it is no different than calling ``to_netcdf``
    repeatedly.

    Parameters
    ----------
    datasets : list of Dataset
        List of datasets to save.
    paths : list of str or list of path-like objects
        List of paths to which to save each corresponding dataset.
    mode : {"w", "a"}, optional
        Write ("w") or append ("a") mode. If mode="w", any existing file at
        these locations will be overwritten.
    format : {"NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_64BIT", \
              "NETCDF3_CLASSIC"}, optional
        File format for the resulting netCDF file:

        * NETCDF4: Data is stored in an HDF5 file, using netCDF4 API
          features.
        * NETCDF4_CLASSIC: Data is stored in an HDF5 file, using only
          netCDF 3 compatible API features.
        * NETCDF3_64BIT: 64-bit offset version of the netCDF 3 file format,
          which fully supports 2+ GB files, but is only compatible with
          clients linked against netCDF version 3.6.0 or later.
        * NETCDF3_CLASSIC: The classic netCDF 3 file format. It does not
          handle 2+ GB files very well.

        All formats are supported by the netCDF4-python library.
        scipy.io.netcdf only supports the last two formats.

        The default format is NETCDF4 if you are saving a file to disk and
        have the netCDF4-python library available. Otherwise, xarray falls
        back to using scipy to write netCDF files and defaults to the
        NETCDF3_64BIT format (scipy does not support netCDF4).
    groups : list of str, optional
        Paths to the netCDF4 group in each corresponding file to which to save
        datasets (only works for format="NETCDF4"). The groups will be created
        if necessary.
    engine : {"netcdf4", "scipy", "h5netcdf"}, optional
        Engine to use when writing netCDF files. If not provided, the
        default engine is chosen based on available dependencies, with a
        preference for "netcdf4" if writing to a file on disk.
        See `Dataset.to_netcdf` for additional information.
    compute : bool
        If true compute immediately, otherwise return a
        ``dask.delayed.Delayed`` object that can be computed later.
    **kwargs : dict, optional
        Additional arguments are passed along to ``to_netcdf``.

    Examples
    --------
    Save a dataset into one netCDF per year of data:

    >>> ds = xr.Dataset(
    ...     {"a": ("time", np.linspace(0, 1, 48))},
    ...     coords={"time": pd.date_range("2010-01-01", freq="ME", periods=48)},
    ... )
    >>> ds
    <xarray.Dataset> Size: 768B
    Dimensions:  (time: 48)
    Coordinates:
      * time     (time) datetime64[ns] 384B 2010-01-31 2010-02-28 ... 2013-12-31
    Data variables:
        a        (time) float64 384B 0.0 0.02128 0.04255 ... 0.9574 0.9787 1.0
    >>> years, datasets = zip(*ds.groupby("time.year"))
    >>> paths = [f"{y}.nc" for y in years]
    >>> xr.save_mfdataset(datasets, paths)
    """
    if mode == "w" and len(set(paths)) < len(paths):
        raise ValueError(
            "cannot use mode='w' when writing multiple datasets to the same path"
        )

    for obj in datasets:
        if not isinstance(obj, Dataset):
            raise TypeError(
                "save_mfdataset only supports writing Dataset "
                f"objects, received type {type(obj)}"
            )

    if groups is None:
        groups = [None] * len(datasets)

    if len({len(datasets), len(paths), len(groups)}) > 1:
        raise ValueError(
            "must supply lists of the same length for the "
            "datasets, paths and groups arguments to "
            "save_mfdataset"
        )

    writers, stores = zip(
        *[
            to_netcdf(
                ds,
                path,
                mode,
                format,
                group,
                engine,
                compute=compute,
                multifile=True,
                **kwargs,
            )
            for ds, path, group in zip(datasets, paths, groups, strict=True)
        ],
        strict=True,
    )

    try:
        writes = [w.sync(compute=compute) for w in writers]
    finally:
        for store in stores:
            if compute:
                store.close()
            else:
                store.sync()

    if not compute:
        import dask

        return dask.delayed(
            list(starmap(delayed_close_after_writes, zip(writes, stores, strict=True)))
        )


def get_writable_zarr_store(
    store: ZarrStoreLike | None = None,
    *,
    chunk_store: MutableMapping | str | os.PathLike | None = None,
    mode: ZarrWriteModes | None = None,
    synchronizer=None,
    group: str | None = None,
    consolidated: bool | None = None,
    append_dim: Hashable | None = None,
    region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None = None,
    safe_chunks: bool = True,
    align_chunks: bool = False,
    storage_options: dict[str, str] | None = None,
    zarr_version: int | None = None,
    zarr_format: int | None = None,
    write_empty_chunks: bool | None = None,
) -> backends.ZarrStore:
    """Create a store for writing to Zarr."""
    from xarray.backends.zarr import _choose_default_mode, _get_mappers

    kwargs, mapper, chunk_mapper = _get_mappers(
        storage_options=storage_options, store=store, chunk_store=chunk_store
    )
    mode = _choose_default_mode(mode=mode, append_dim=append_dim, region=region)

    if mode == "r+":
        already_consolidated = consolidated
        consolidate_on_close = False
    else:
        already_consolidated = False
        consolidate_on_close = consolidated or consolidated is None

    return backends.ZarrStore.open_group(
        store=mapper,
        mode=mode,
        synchronizer=synchronizer,
        group=group,
        consolidated=already_consolidated,
        consolidate_on_close=consolidate_on_close,
        chunk_store=chunk_mapper,
        append_dim=append_dim,
        write_region=region,
        safe_chunks=safe_chunks,
        align_chunks=align_chunks,
        zarr_version=zarr_version,
        zarr_format=zarr_format,
        write_empty=write_empty_chunks,
        **kwargs,
    )


# compute=True returns ZarrStore
@overload
def to_zarr(
    dataset: Dataset,
    store: ZarrStoreLike | None = None,
    chunk_store: MutableMapping | str | os.PathLike | None = None,
    mode: ZarrWriteModes | None = None,
    synchronizer=None,
    group: str | None = None,
    encoding: Mapping | None = None,
    *,
    compute: Literal[True] = True,
    consolidated: bool | None = None,
    append_dim: Hashable | None = None,
    region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None = None,
    safe_chunks: bool = True,
    align_chunks: bool = False,
    storage_options: dict[str, str] | None = None,
    zarr_version: int | None = None,
    write_empty_chunks: bool | None = None,
    chunkmanager_store_kwargs: dict[str, Any] | None = None,
) -> backends.ZarrStore: ...


# compute=False returns dask.Delayed
@overload
def to_zarr(
    dataset: Dataset,
    store: ZarrStoreLike | None = None,
    chunk_store: MutableMapping | str | os.PathLike | None = None,
    mode: ZarrWriteModes | None = None,
    synchronizer=None,
    group: str | None = None,
    encoding: Mapping | None = None,
    *,
    compute: Literal[False],
    consolidated: bool | None = None,
    append_dim: Hashable | None = None,
    region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None = None,
    safe_chunks: bool = True,
    align_chunks: bool = False,
    storage_options: dict[str, str] | None = None,
    zarr_version: int | None = None,
    write_empty_chunks: bool | None = None,
    chunkmanager_store_kwargs: dict[str, Any] | None = None,
) -> Delayed: ...


def to_zarr(
    dataset: Dataset,
    store: ZarrStoreLike | None = None,
    chunk_store: MutableMapping | str | os.PathLike | None = None,
    mode: ZarrWriteModes | None = None,
    synchronizer=None,
    group: str | None = None,
    encoding: Mapping | None = None,
    *,
    compute: bool = True,
    consolidated: bool | None = None,
    append_dim: Hashable | None = None,
    region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None = None,
    safe_chunks: bool = True,
    align_chunks: bool = False,
    storage_options: dict[str, str] | None = None,
    zarr_version: int | None = None,
    zarr_format: int | None = None,
    write_empty_chunks: bool | None = None,
    chunkmanager_store_kwargs: dict[str, Any] | None = None,
) -> backends.ZarrStore | Delayed:
    """This function creates an appropriate datastore for writing a dataset to
    a zarr ztore

    See `Dataset.to_zarr` for full API docs.
    """

    # validate Dataset keys, DataArray names
    _validate_dataset_names(dataset)

    # Load empty arrays to avoid bug saving zero length dimensions (Issue #5741)
    # TODO: delete when min dask>=2023.12.1
    # https://github.com/dask/dask/pull/10506
    for v in dataset.variables.values():
        if v.size == 0:
            v.load()

    if encoding is None:
        encoding = {}

    zstore = get_writable_zarr_store(
        store,
        chunk_store=chunk_store,
        mode=mode,
        synchronizer=synchronizer,
        group=group,
        consolidated=consolidated,
        append_dim=append_dim,
        region=region,
        safe_chunks=safe_chunks,
        align_chunks=align_chunks,
        storage_options=storage_options,
        zarr_version=zarr_version,
        zarr_format=zarr_format,
        write_empty_chunks=write_empty_chunks,
    )

    dataset = zstore._validate_and_autodetect_region(dataset)
    zstore._validate_encoding(encoding)

    writer = ArrayWriter()

    # TODO: figure out how to properly handle unlimited_dims
    try:
        dump_to_store(dataset, zstore, writer, encoding=encoding)
        writes = writer.sync(
            compute=compute, chunkmanager_store_kwargs=chunkmanager_store_kwargs
        )
    finally:
        if compute:
            zstore.close()

    if not compute:
        return delayed_close_after_writes(writes, zstore)

    return zstore
