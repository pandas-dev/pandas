from __future__ import annotations

import base64
import json
import os
import struct
from collections.abc import Hashable, Iterable, Mapping
from typing import TYPE_CHECKING, Any, Literal, Self, cast

import numpy as np
import pandas as pd

from xarray import coding, conventions
from xarray.backends.chunks import grid_rechunk, validate_grid_chunks_alignment
from xarray.backends.common import (
    BACKEND_ENTRYPOINTS,
    AbstractWritableDataStore,
    BackendArray,
    BackendEntrypoint,
    T_PathFileOrDataStore,
    _encode_variable_name,
    _normalize_path,
    datatree_from_dict_with_io_cleanup,
    ensure_dtype_not_object,
)
from xarray.backends.store import StoreBackendEntrypoint
from xarray.core import indexing
from xarray.core.treenode import NodePath
from xarray.core.types import ZarrWriteModes
from xarray.core.utils import (
    FrozenDict,
    HiddenKeyDict,
    attempt_import,
    close_on_error,
    emit_user_level_warning,
)
from xarray.core.variable import Variable
from xarray.namedarray.parallelcompat import guess_chunkmanager
from xarray.namedarray.pycompat import integer_types
from xarray.namedarray.utils import module_available

if TYPE_CHECKING:
    from xarray.core.dataset import Dataset
    from xarray.core.datatree import DataTree
    from xarray.core.types import ZarrArray, ZarrGroup


def _get_mappers(*, storage_options, store, chunk_store):
    # expand str and path-like arguments
    store = _normalize_path(store)
    chunk_store = _normalize_path(chunk_store)

    kwargs = {}
    if storage_options is None:
        mapper = store
        chunk_mapper = chunk_store
    else:
        if not isinstance(store, str):
            raise ValueError(
                f"store must be a string to use storage_options. Got {type(store)}"
            )

        if _zarr_v3():
            kwargs["storage_options"] = storage_options
            mapper = store
            chunk_mapper = chunk_store
        else:
            from fsspec import get_mapper

            mapper = get_mapper(store, **storage_options)
            if chunk_store is not None:
                chunk_mapper = get_mapper(chunk_store, **storage_options)
            else:
                chunk_mapper = chunk_store
    return kwargs, mapper, chunk_mapper


def _choose_default_mode(
    *,
    mode: ZarrWriteModes | None,
    append_dim: Hashable | None,
    region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None,
) -> ZarrWriteModes:
    if mode is None:
        if append_dim is not None:
            mode = "a"
        elif region is not None:
            mode = "r+"
        else:
            mode = "w-"

    if mode not in ["a", "a-"] and append_dim is not None:
        raise ValueError("cannot set append_dim unless mode='a' or mode=None")

    if mode not in ["a", "a-", "r+"] and region is not None:
        raise ValueError(
            "cannot set region unless mode='a', mode='a-', mode='r+' or mode=None"
        )

    if mode not in ["w", "w-", "a", "a-", "r+"]:
        raise ValueError(
            "The only supported options for mode are 'w', "
            f"'w-', 'a', 'a-', and 'r+', but mode={mode!r}"
        )
    return mode


def _zarr_v3() -> bool:
    return module_available("zarr", minversion="3")


# need some special secret attributes to tell us the dimensions
DIMENSION_KEY = "_ARRAY_DIMENSIONS"
ZarrFormat = Literal[2, 3]


class FillValueCoder:
    """Handle custom logic to safely encode and decode fill values in Zarr.
    Possibly redundant with logic in xarray/coding/variables.py but needs to be
    isolated from NetCDF-specific logic.
    """

    @classmethod
    def encode(cls, value: int | float | str | bytes, dtype: np.dtype[Any]) -> Any:
        if dtype.kind in "S":
            # byte string, this implies that 'value' must also be `bytes` dtype.
            assert isinstance(value, bytes)
            return base64.standard_b64encode(value).decode()
        elif dtype.kind in "b":
            # boolean
            return bool(value)
        elif dtype.kind in "iu":
            # todo: do we want to check for decimals?
            return int(value)
        elif dtype.kind in "f":
            return base64.standard_b64encode(struct.pack("<d", float(value))).decode()
        elif dtype.kind in "U":
            return str(value)
        else:
            raise ValueError(f"Failed to encode fill_value. Unsupported dtype {dtype}")

    @classmethod
    def decode(cls, value: int | float | str | bytes, dtype: str | np.dtype[Any]):
        if dtype == "string":
            # zarr V3 string type
            return str(value)
        elif dtype == "bytes":
            # zarr V3 bytes type
            assert isinstance(value, str | bytes)
            return base64.standard_b64decode(value)
        np_dtype = np.dtype(dtype)
        if np_dtype.kind in "f":
            assert isinstance(value, str | bytes)
            return struct.unpack("<d", base64.standard_b64decode(value))[0]
        elif np_dtype.kind in "b":
            return bool(value)
        elif np_dtype.kind in "iu":
            return int(value)
        else:
            raise ValueError(f"Failed to decode fill_value. Unsupported dtype {dtype}")


def encode_zarr_attr_value(value):
    """
    Encode a attribute value as something that can be serialized as json

    Many xarray datasets / variables have numpy arrays and values. This
    function handles encoding / decoding of such items.

    ndarray -> list
    scalar array -> scalar
    other -> other (no change)
    """
    if isinstance(value, np.ndarray):
        encoded = value.tolist()
    elif isinstance(value, np.generic):
        encoded = value.item()
    else:
        encoded = value
    return encoded


def has_zarr_async_index() -> bool:
    try:
        import zarr

        return hasattr(zarr.AsyncArray, "oindex")
    except (ImportError, AttributeError):
        return False


class ZarrArrayWrapper(BackendArray):
    __slots__ = ("_array", "dtype", "shape")

    def __init__(self, zarr_array):
        # some callers attempt to evaluate an array if an `array` property exists on the object.
        # we prefix with _ to avoid this inference.

        # TODO type hint this?
        self._array = zarr_array
        self.shape = self._array.shape

        # preserve vlen string object dtype (GH 7328)
        if (
            not _zarr_v3()
            and self._array.filters is not None
            and any(filt.codec_id == "vlen-utf8" for filt in self._array.filters)
        ):
            dtype = coding.strings.create_vlen_dtype(str)
        else:
            dtype = self._array.dtype

        self.dtype = dtype

    def get_array(self):
        return self._array

    def _oindex(self, key):
        return self._array.oindex[key]

    def _vindex(self, key):
        return self._array.vindex[key]

    def _getitem(self, key):
        return self._array[key]

    async def _async_getitem(self, key):
        if not _zarr_v3():
            raise NotImplementedError(
                "For lazy basic async indexing with zarr, zarr-python=>v3.0.0 is required"
            )

        async_array = self._array._async_array
        return await async_array.getitem(key)

    async def _async_oindex(self, key):
        if not has_zarr_async_index():
            raise NotImplementedError(
                "For lazy orthogonal async indexing with zarr, zarr-python=>v3.1.2 is required"
            )

        async_array = self._array._async_array
        return await async_array.oindex.getitem(key)

    async def _async_vindex(self, key):
        if not has_zarr_async_index():
            raise NotImplementedError(
                "For lazy vectorized async indexing with zarr, zarr-python=>v3.1.2 is required"
            )

        async_array = self._array._async_array
        return await async_array.vindex.getitem(key)

    def __getitem__(self, key):
        array = self._array
        if isinstance(key, indexing.BasicIndexer):
            method = self._getitem
        elif isinstance(key, indexing.VectorizedIndexer):
            method = self._vindex
        elif isinstance(key, indexing.OuterIndexer):
            method = self._oindex
        return indexing.explicit_indexing_adapter(
            key, array.shape, indexing.IndexingSupport.VECTORIZED, method
        )

        # if self.ndim == 0:
        # could possibly have a work-around for 0d data here

    async def async_getitem(self, key):
        array = self._array
        if isinstance(key, indexing.BasicIndexer):
            method = self._async_getitem
        elif isinstance(key, indexing.VectorizedIndexer):
            method = self._async_vindex
        elif isinstance(key, indexing.OuterIndexer):
            method = self._async_oindex
        return await indexing.async_explicit_indexing_adapter(
            key, array.shape, indexing.IndexingSupport.VECTORIZED, method
        )


def _determine_zarr_chunks(enc_chunks, var_chunks, ndim, name):
    """
    Given encoding chunks (possibly None or []) and variable chunks
    (possibly None or []).
    """

    # zarr chunk spec:
    # chunks : int or tuple of ints, optional
    #   Chunk shape. If not provided, will be guessed from shape and dtype.

    # if there are no chunks in encoding and the variable data is a numpy
    # array, then we let zarr use its own heuristics to pick the chunks
    if not var_chunks and not enc_chunks:
        return None

    # if there are no chunks in encoding but there are dask chunks, we try to
    # use the same chunks in zarr
    # However, zarr chunks needs to be uniform for each array
    # https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html#chunks
    # while dask chunks can be variable sized
    # https://dask.pydata.org/en/latest/array-design.html#chunks
    if var_chunks and not enc_chunks:
        if any(len(set(chunks[:-1])) > 1 for chunks in var_chunks):
            raise ValueError(
                "Zarr requires uniform chunk sizes except for final chunk. "
                f"Variable named {name!r} has incompatible dask chunks: {var_chunks!r}. "
                "Consider rechunking using `chunk()`."
            )
        if any((chunks[0] < chunks[-1]) for chunks in var_chunks):
            raise ValueError(
                "Final chunk of Zarr array must be the same size or smaller "
                f"than the first. Variable named {name!r} has incompatible Dask chunks {var_chunks!r}."
                "Consider either rechunking using `chunk()` or instead deleting "
                "or modifying `encoding['chunks']`."
            )
        # return the first chunk for each dimension
        return tuple(chunk[0] for chunk in var_chunks)

    # From here on, we are dealing with user-specified chunks in encoding
    # zarr allows chunks to be an integer, in which case it uses the same chunk
    # size on each dimension.
    # Here we re-implement this expansion ourselves. That makes the logic of
    # checking chunk compatibility easier

    if isinstance(enc_chunks, integer_types):
        enc_chunks_tuple = ndim * (enc_chunks,)
    else:
        enc_chunks_tuple = tuple(enc_chunks)

    if len(enc_chunks_tuple) != ndim:
        # throw away encoding chunks, start over
        return _determine_zarr_chunks(
            None,
            var_chunks,
            ndim,
            name,
        )

    for x in enc_chunks_tuple:
        if not isinstance(x, int):
            raise TypeError(
                "zarr chunk sizes specified in `encoding['chunks']` "
                "must be an int or a tuple of ints. "
                f"Instead found encoding['chunks']={enc_chunks_tuple!r} "
                f"for variable named {name!r}."
            )

    # if there are chunks in encoding and the variable data is a numpy array,
    # we use the specified chunks
    if not var_chunks:
        return enc_chunks_tuple

    return enc_chunks_tuple


def _get_zarr_dims_and_attrs(zarr_obj, dimension_key, try_nczarr):
    # Zarr V3 explicitly stores the dimension names in the metadata
    try:
        # if this exists, we are looking at a Zarr V3 array
        # convert None to empty tuple
        dimensions = zarr_obj.metadata.dimension_names or ()
    except AttributeError:
        # continue to old code path
        pass
    else:
        attributes = dict(zarr_obj.attrs)
        if len(zarr_obj.shape) != len(dimensions):
            raise KeyError(
                "Zarr object is missing the `dimension_names` metadata which is "
                "required for xarray to determine variable dimensions."
            )
        return dimensions, attributes

    # Zarr arrays do not have dimensions. To get around this problem, we add
    # an attribute that specifies the dimension. We have to hide this attribute
    # when we send the attributes to the user.
    # zarr_obj can be either a zarr group or zarr array
    try:
        # Xarray-Zarr
        dimensions = zarr_obj.attrs[dimension_key]
    except KeyError as e:
        if not try_nczarr:
            raise KeyError(
                f"Zarr object is missing the attribute `{dimension_key}`, which is "
                "required for xarray to determine variable dimensions."
            ) from e

        # NCZarr defines dimensions through metadata in .zarray
        zarray_path = os.path.join(zarr_obj.path, ".zarray")
        if _zarr_v3():
            import asyncio

            zarray_str = asyncio.run(zarr_obj.store.get(zarray_path)).to_bytes()
        else:
            zarray_str = zarr_obj.store.get(zarray_path)
        zarray = json.loads(zarray_str)
        try:
            # NCZarr uses Fully Qualified Names
            dimensions = [
                os.path.basename(dim) for dim in zarray["_NCZARR_ARRAY"]["dimrefs"]
            ]
        except KeyError as e:
            raise KeyError(
                f"Zarr object is missing the attribute `{dimension_key}` and the NCZarr metadata, "
                "which are required for xarray to determine variable dimensions."
            ) from e

    nc_attrs = [attr for attr in zarr_obj.attrs if attr.lower().startswith("_nc")]
    attributes = HiddenKeyDict(zarr_obj.attrs, [dimension_key] + nc_attrs)
    return dimensions, attributes


def extract_zarr_variable_encoding(
    variable,
    raise_on_invalid=False,
    name=None,
    *,
    zarr_format: ZarrFormat,
):
    """
    Extract zarr encoding dictionary from xarray Variable

    Parameters
    ----------
    variable : Variable
    raise_on_invalid : bool, optional
    name: str | Hashable, optional
    zarr_format: Literal[2,3]
    Returns
    -------
    encoding : dict
        Zarr encoding for `variable`
    """

    encoding = variable.encoding.copy()

    safe_to_drop = {"source", "original_shape", "preferred_chunks"}
    valid_encodings = {
        "chunks",
        "shards",
        "compressor",  # TODO: delete when min zarr >=3
        "compressors",
        "filters",
        "serializer",
        "cache_metadata",
        "write_empty_chunks",
        "chunk_key_encoding",
    }
    if zarr_format == 3:
        valid_encodings.add("fill_value")

    for k in safe_to_drop:
        if k in encoding:
            del encoding[k]

    if raise_on_invalid:
        invalid = [k for k in encoding if k not in valid_encodings]
        if "fill_value" in invalid and zarr_format == 2:
            msg = " Use `_FillValue` to set the Zarr array `fill_value`"
        else:
            msg = ""

        if invalid:
            raise ValueError(
                f"unexpected encoding parameters for zarr backend:  {invalid!r}." + msg
            )
    else:
        for k in list(encoding):
            if k not in valid_encodings:
                del encoding[k]

    chunks = _determine_zarr_chunks(
        enc_chunks=encoding.get("chunks"),
        var_chunks=variable.chunks,
        ndim=variable.ndim,
        name=name,
    )
    if _zarr_v3() and chunks is None:
        chunks = "auto"
    encoding["chunks"] = chunks
    return encoding


# Function below is copied from conventions.encode_cf_variable.
# The only change is to raise an error for object dtypes.
def encode_zarr_variable(var, needs_copy=True, name=None):
    """
    Converts an Variable into an Variable which follows some
    of the CF conventions:

        - Nans are masked using _FillValue (or the deprecated missing_value)
        - Rescaling via: scale_factor and add_offset
        - datetimes are converted to the CF 'units since time' format
        - dtype encodings are enforced.

    Parameters
    ----------
    var : Variable
        A variable holding un-encoded data.

    Returns
    -------
    out : Variable
        A variable which has been encoded as described above.
    """

    var = conventions.encode_cf_variable(var, name=name)
    var = ensure_dtype_not_object(var, name=name)

    # zarr allows unicode, but not variable-length strings, so it's both
    # simpler and more compact to always encode as UTF-8 explicitly.
    # TODO: allow toggling this explicitly via dtype in encoding.
    # TODO: revisit this now that Zarr _does_ allow variable-length strings
    coder = coding.strings.EncodedStringCoder(allows_unicode=True)
    var = coder.encode(var, name=name)
    var = coding.strings.ensure_fixed_length_bytes(var)

    return var


def _validate_datatypes_for_zarr_append(vname, existing_var, new_var):
    """If variable exists in the store, confirm dtype of the data to append is compatible with
    existing dtype.
    """
    if (
        np.issubdtype(new_var.dtype, np.number)
        or np.issubdtype(new_var.dtype, np.datetime64)
        or np.issubdtype(new_var.dtype, np.bool_)
        or new_var.dtype == object
        or (new_var.dtype.kind in ("S", "U") and existing_var.dtype == object)
    ):
        # We can skip dtype equality checks under two conditions: (1) if the var to append is
        # new to the dataset, because in this case there is no existing var to compare it to;
        # or (2) if var to append's dtype is known to be easy-to-append, because in this case
        # we can be confident appending won't cause problems. Examples of dtypes which are not
        # easy-to-append include length-specified strings of type `|S*` or `<U*` (where * is a
        # positive integer character length). For these dtypes, appending dissimilar lengths
        # can result in truncation of appended data. Therefore, variables which already exist
        # in the dataset, and with dtypes which are not known to be easy-to-append, necessitate
        # exact dtype equality, as checked below.
        pass
    elif new_var.dtype != existing_var.dtype:
        raise ValueError(
            f"Mismatched dtypes for variable {vname} between Zarr store on disk "
            f"and dataset to append. Store has dtype {existing_var.dtype} but "
            f"dataset to append has dtype {new_var.dtype}."
        )


def _validate_and_transpose_existing_dims(
    var_name, new_var, existing_var, region, append_dim
):
    if new_var.dims != existing_var.dims:
        if set(existing_var.dims) == set(new_var.dims):
            new_var = new_var.transpose(*existing_var.dims)
        else:
            raise ValueError(
                f"variable {var_name!r} already exists with different "
                f"dimension names {existing_var.dims} != "
                f"{new_var.dims}, but changing variable "
                f"dimensions is not supported by to_zarr()."
            )

    existing_sizes = {}
    for dim, size in existing_var.sizes.items():
        if region is not None and dim in region:
            start, stop, stride = region[dim].indices(size)
            assert stride == 1  # region was already validated
            size = stop - start
        if dim != append_dim:
            existing_sizes[dim] = size

    new_sizes = {dim: size for dim, size in new_var.sizes.items() if dim != append_dim}
    if existing_sizes != new_sizes:
        raise ValueError(
            f"variable {var_name!r} already exists with different "
            f"dimension sizes: {existing_sizes} != {new_sizes}. "
            f"to_zarr() only supports changing dimension sizes when "
            f"explicitly appending, but append_dim={append_dim!r}. "
            f"If you are attempting to write to a subset of the "
            f"existing store without changing dimension sizes, "
            f"consider using the region argument in to_zarr()."
        )

    return new_var


def _put_attrs(zarr_obj, attrs):
    """Raise a more informative error message for invalid attrs."""
    try:
        zarr_obj.attrs.put(attrs)
    except TypeError as e:
        raise TypeError("Invalid attribute in Dataset.attrs.") from e
    return zarr_obj


class ZarrStore(AbstractWritableDataStore):
    """Store for reading and writing data via zarr"""

    __slots__ = (
        "_align_chunks",
        "_append_dim",
        "_cache_members",
        "_close_store_on_close",
        "_consolidate_on_close",
        "_group",
        "_members",
        "_mode",
        "_read_only",
        "_safe_chunks",
        "_synchronizer",
        "_use_zarr_fill_value_as_mask",
        "_write_empty",
        "_write_region",
        "zarr_group",
    )

    @classmethod
    def open_store(
        cls,
        store,
        mode: ZarrWriteModes = "r",
        synchronizer=None,
        group=None,
        consolidated=False,
        consolidate_on_close=False,
        chunk_store=None,
        storage_options=None,
        append_dim=None,
        write_region=None,
        safe_chunks=True,
        align_chunks=False,
        zarr_version=None,
        zarr_format=None,
        use_zarr_fill_value_as_mask=None,
        write_empty: bool | None = None,
        cache_members: bool = True,
    ):
        (
            zarr_group,
            consolidate_on_close,
            close_store_on_close,
            use_zarr_fill_value_as_mask,
        ) = _get_open_params(
            store=store,
            mode=mode,
            synchronizer=synchronizer,
            group=group,
            consolidated=consolidated,
            consolidate_on_close=consolidate_on_close,
            chunk_store=chunk_store,
            storage_options=storage_options,
            zarr_version=zarr_version,
            use_zarr_fill_value_as_mask=use_zarr_fill_value_as_mask,
            zarr_format=zarr_format,
        )

        from zarr import Group

        group_members: dict[str, Group] = {}
        group_paths = list(_iter_zarr_groups(zarr_group, parent=group))
        for path in group_paths:
            if path == group:
                group_members[path] = zarr_group
            else:
                rel_path = path.removeprefix(f"{group}/")
                group_members[path] = zarr_group[rel_path.removeprefix("/")]

        out = {
            group: cls(
                group_store,
                mode,
                consolidate_on_close,
                append_dim,
                write_region,
                safe_chunks,
                write_empty,
                close_store_on_close,
                use_zarr_fill_value_as_mask,
                align_chunks=align_chunks,
                cache_members=cache_members,
            )
            for group, group_store in group_members.items()
        }
        return out

    @classmethod
    def open_group(
        cls,
        store,
        mode: ZarrWriteModes = "r",
        synchronizer=None,
        group=None,
        consolidated=False,
        consolidate_on_close=False,
        chunk_store=None,
        storage_options=None,
        append_dim=None,
        write_region=None,
        safe_chunks=True,
        align_chunks=False,
        zarr_version=None,
        zarr_format=None,
        use_zarr_fill_value_as_mask=None,
        write_empty: bool | None = None,
        cache_members: bool = True,
    ):
        (
            zarr_group,
            consolidate_on_close,
            close_store_on_close,
            use_zarr_fill_value_as_mask,
        ) = _get_open_params(
            store=store,
            mode=mode,
            synchronizer=synchronizer,
            group=group,
            consolidated=consolidated,
            consolidate_on_close=consolidate_on_close,
            chunk_store=chunk_store,
            storage_options=storage_options,
            zarr_version=zarr_version,
            use_zarr_fill_value_as_mask=use_zarr_fill_value_as_mask,
            zarr_format=zarr_format,
        )

        return cls(
            zarr_group,
            mode,
            consolidate_on_close,
            append_dim,
            write_region,
            safe_chunks,
            write_empty,
            close_store_on_close,
            use_zarr_fill_value_as_mask,
            align_chunks=align_chunks,
            cache_members=cache_members,
        )

    def __init__(
        self,
        zarr_group,
        mode=None,
        consolidate_on_close=False,
        append_dim=None,
        write_region=None,
        safe_chunks=True,
        write_empty: bool | None = None,
        close_store_on_close: bool = False,
        use_zarr_fill_value_as_mask=None,
        align_chunks: bool = False,
        cache_members: bool = True,
    ):
        if align_chunks:
            # Disabled the safe_chunks validations if the alignment is going to be applied
            safe_chunks = False

        self.zarr_group = zarr_group
        self._read_only = self.zarr_group.read_only
        self._synchronizer = self.zarr_group.synchronizer
        self._group = self.zarr_group.path
        self._mode = mode
        self._consolidate_on_close = consolidate_on_close
        self._append_dim = append_dim
        self._write_region = write_region
        self._align_chunks = align_chunks
        self._safe_chunks = safe_chunks
        self._write_empty = write_empty
        self._close_store_on_close = close_store_on_close
        self._use_zarr_fill_value_as_mask = use_zarr_fill_value_as_mask
        self._cache_members: bool = cache_members
        self._members: dict[str, ZarrArray | ZarrGroup] = {}

        if self._cache_members:
            # initialize the cache
            # this cache is created here and never updated.
            # If the `ZarrStore` instance creates a new zarr array, or if an external process
            # removes an existing zarr array, then the cache will be invalid.
            # We use this cache only to record any pre-existing arrays when the group was opened
            # create a new ZarrStore instance if you want to
            # capture the current state of the zarr group, or create a ZarrStore with
            # `cache_members` set to `False` to disable this cache and instead fetch members
            # on demand.
            self._members = self._fetch_members()

    def get_child_store(self, group: str) -> Self:
        zarr_group = self.zarr_group.require_group(group)
        return type(self)(
            zarr_group=zarr_group,
            mode=self._mode,
            consolidate_on_close=self._consolidate_on_close,
            append_dim=self._append_dim,
            write_region=self._write_region,
            safe_chunks=self._safe_chunks,
            write_empty=self._write_empty,
            close_store_on_close=self._close_store_on_close,
            use_zarr_fill_value_as_mask=self._use_zarr_fill_value_as_mask,
            align_chunks=self._align_chunks,
            cache_members=self._cache_members,
        )

    @property
    def members(self) -> dict[str, ZarrArray | ZarrGroup]:
        """
        Model the arrays and groups contained in self.zarr_group as a dict. If `self._cache_members`
        is true, the dict is cached. Otherwise, it is retrieved from storage.
        """
        if not self._cache_members:
            return self._fetch_members()
        else:
            return self._members

    def _fetch_members(self) -> dict[str, ZarrArray | ZarrGroup]:
        """
        Get the arrays and groups defined in the zarr group modelled by this Store
        """
        import zarr

        if zarr.__version__ >= "3":
            return dict(self.zarr_group.members())
        else:
            return dict(self.zarr_group.items())

    def array_keys(self) -> tuple[str, ...]:
        from zarr import Array as ZarrArray

        return tuple(
            key for (key, node) in self.members.items() if isinstance(node, ZarrArray)
        )

    def arrays(self) -> tuple[tuple[str, ZarrArray], ...]:
        from zarr import Array as ZarrArray

        return tuple(
            (key, node)
            for (key, node) in self.members.items()
            if isinstance(node, ZarrArray)
        )

    @property
    def ds(self):
        # TODO: consider deprecating this in favor of zarr_group
        return self.zarr_group

    def open_store_variable(self, name):
        zarr_array = self.members[name]
        data = indexing.LazilyIndexedArray(ZarrArrayWrapper(zarr_array))
        try_nczarr = self._mode == "r"
        dimensions, attributes = _get_zarr_dims_and_attrs(
            zarr_array, DIMENSION_KEY, try_nczarr
        )
        attributes = dict(attributes)

        encoding = {
            "chunks": zarr_array.chunks,
            "preferred_chunks": dict(zip(dimensions, zarr_array.chunks, strict=True)),
        }

        if _zarr_v3():
            encoding.update(
                {
                    "compressors": zarr_array.compressors,
                    "filters": zarr_array.filters,
                    "shards": zarr_array.shards,
                }
            )
            if self.zarr_group.metadata.zarr_format == 3:
                encoding.update({"serializer": zarr_array.serializer})
        else:
            encoding.update(
                {
                    "compressor": zarr_array.compressor,
                    "filters": zarr_array.filters,
                }
            )

        if self._use_zarr_fill_value_as_mask:
            # Setting this attribute triggers CF decoding for missing values
            # by interpreting Zarr's fill_value to mean the same as netCDF's _FillValue
            if zarr_array.fill_value is not None:
                attributes["_FillValue"] = zarr_array.fill_value
        elif "_FillValue" in attributes:
            attributes["_FillValue"] = FillValueCoder.decode(
                attributes["_FillValue"], zarr_array.dtype
            )

        return Variable(dimensions, data, attributes, encoding)

    def get_variables(self):
        return FrozenDict((k, self.open_store_variable(k)) for k in self.array_keys())

    def get_attrs(self):
        return {
            k: v
            for k, v in self.zarr_group.attrs.asdict().items()
            if not k.lower().startswith("_nc")
        }

    def get_dimensions(self):
        try_nczarr = self._mode == "r"
        dimensions = {}
        for _k, v in self.arrays():
            dim_names, _ = _get_zarr_dims_and_attrs(v, DIMENSION_KEY, try_nczarr)
            for d, s in zip(dim_names, v.shape, strict=True):
                if d in dimensions and dimensions[d] != s:
                    raise ValueError(
                        f"found conflicting lengths for dimension {d} "
                        f"({s} != {dimensions[d]})"
                    )
                dimensions[d] = s
        return dimensions

    def set_dimensions(self, variables, unlimited_dims=None):
        if unlimited_dims is not None:
            raise NotImplementedError(
                "Zarr backend doesn't know how to handle unlimited dimensions"
            )

    def set_attributes(self, attributes):
        _put_attrs(self.zarr_group, attributes)

    def encode_variable(self, variable, name=None):
        variable = encode_zarr_variable(variable, name=name)
        return variable

    def encode_attribute(self, a):
        return encode_zarr_attr_value(a)

    def store(
        self,
        variables,
        attributes,
        check_encoding_set=frozenset(),
        writer=None,
        unlimited_dims=None,
    ):
        """
        Top level method for putting data on this store, this method:
          - encodes variables/attributes
          - sets dimensions
          - sets variables

        Parameters
        ----------
        variables : dict-like
            Dictionary of key/value (variable name / xr.Variable) pairs
        attributes : dict-like
            Dictionary of key/value (attribute name / attribute) pairs
        check_encoding_set : list-like
            List of variables that should be checked for invalid encoding
            values
        writer : ArrayWriter
        unlimited_dims : list-like
            List of dimension names that should be treated as unlimited
            dimensions.
            dimension on which the zarray will be appended
            only needed in append mode
        """
        if TYPE_CHECKING:
            import zarr
        else:
            zarr = attempt_import("zarr")

        if self._mode == "w":
            # always overwrite, so we don't care about existing names,
            # and consistency of encoding
            new_variable_names = set(variables)
            existing_keys = {}
            existing_variable_names = {}
        else:
            existing_keys = self.array_keys()
            existing_variable_names = {
                vn for vn in variables if _encode_variable_name(vn) in existing_keys
            }
            new_variable_names = set(variables) - existing_variable_names

        if self._mode == "r+" and (
            new_names := [k for k in variables if k not in existing_keys]
        ):
            raise ValueError(
                f"dataset contains non-pre-existing variables {new_names!r}, "
                "which is not allowed in ``xarray.Dataset.to_zarr()`` with "
                "``mode='r+'``. To allow writing new variables, set ``mode='a'``."
            )

        if self._append_dim is not None and self._append_dim not in existing_keys:
            # For dimensions without coordinate values, we must parse
            # the _ARRAY_DIMENSIONS attribute on *all* arrays to check if it
            # is a valid existing dimension name.
            # TODO: This `get_dimensions` method also does shape checking
            # which isn't strictly necessary for our check.
            existing_dims = self.get_dimensions()
            if self._append_dim not in existing_dims:
                raise ValueError(
                    f"append_dim={self._append_dim!r} does not match any existing "
                    f"dataset dimensions {existing_dims}"
                )

        variables_encoded, attributes = self.encode(
            {vn: variables[vn] for vn in new_variable_names}, attributes
        )

        if existing_variable_names:
            # We make sure that values to be appended are encoded *exactly*
            # as the current values in the store.
            # To do so, we decode variables directly to access the proper encoding,
            # without going via xarray.Dataset to avoid needing to load
            # index variables into memory.
            existing_vars, _, _ = conventions.decode_cf_variables(
                variables={
                    k: self.open_store_variable(name=k) for k in existing_variable_names
                },
                # attributes = {} since we don't care about parsing the global
                # "coordinates" attribute
                attributes={},
            )
            # Modified variables must use the same encoding as the store.
            vars_with_encoding = {}
            for vn in existing_variable_names:
                _validate_datatypes_for_zarr_append(
                    vn, existing_vars[vn], variables[vn]
                )
                vars_with_encoding[vn] = variables[vn].copy(deep=False)
                vars_with_encoding[vn].encoding = existing_vars[vn].encoding
            vars_with_encoding, _ = self.encode(vars_with_encoding, {})
            variables_encoded.update(vars_with_encoding)

            for var_name in existing_variable_names:
                variables_encoded[var_name] = _validate_and_transpose_existing_dims(
                    var_name,
                    variables_encoded[var_name],
                    existing_vars[var_name],
                    self._write_region,
                    self._append_dim,
                )

        if self._mode not in ["r", "r+"]:
            self.set_attributes(attributes)
            self.set_dimensions(variables_encoded, unlimited_dims=unlimited_dims)

        # if we are appending to an append_dim, only write either
        # - new variables not already present, OR
        # - variables with the append_dim in their dimensions
        # We do NOT overwrite other variables.
        if self._mode == "a-" and self._append_dim is not None:
            variables_to_set = {
                k: v
                for k, v in variables_encoded.items()
                if (k not in existing_variable_names) or (self._append_dim in v.dims)
            }
        else:
            variables_to_set = variables_encoded

        self.set_variables(
            variables_to_set, check_encoding_set, writer, unlimited_dims=unlimited_dims
        )
        if self._consolidate_on_close:
            kwargs = {}
            if _zarr_v3():
                kwargs["zarr_format"] = self.zarr_group.metadata.zarr_format
            zarr.consolidate_metadata(self.zarr_group.store, **kwargs)

    def _open_existing_array(self, *, name) -> ZarrArray:
        import zarr
        from zarr import Array as ZarrArray

        # TODO: if mode="a", consider overriding the existing variable
        # metadata. This would need some case work properly with region
        # and append_dim.
        if self._write_empty is not None:
            # Write to zarr_group.chunk_store instead of zarr_group.store
            # See https://github.com/pydata/xarray/pull/8326#discussion_r1365311316 for a longer explanation
            #    The open_consolidated() enforces a mode of r or r+
            #    (and to_zarr with region provided enforces a read mode of r+),
            #    and this function makes sure the resulting Group has a store of type ConsolidatedMetadataStore
            #    and a 'normal Store subtype for chunk_store.
            #    The exact type depends on if a local path was used, or a URL of some sort,
            #    but the point is that it's not a read-only ConsolidatedMetadataStore.
            #    It is safe to write chunk data to the chunk_store because no metadata would be changed by
            #    to_zarr with the region parameter:
            #     - Because the write mode is enforced to be r+, no new variables can be added to the store
            #       (this is also checked and enforced in xarray.backends.api.py::to_zarr()).
            #     - Existing variables already have their attrs included in the consolidated metadata file.
            #     - The size of dimensions can not be expanded, that would require a call using `append_dim`
            #        which is mutually exclusive with `region`
            empty: dict[str, bool] | dict[str, dict[str, bool]]
            if _zarr_v3():
                empty = dict(config={"write_empty_chunks": self._write_empty})
            else:
                empty = dict(write_empty_chunks=self._write_empty)

            zarr_array = zarr.open(
                store=(
                    self.zarr_group.store if _zarr_v3() else self.zarr_group.chunk_store
                ),
                # TODO: see if zarr should normalize these strings.
                path="/".join([self.zarr_group.name.rstrip("/"), name]).lstrip("/"),
                **empty,
            )
        else:
            zarr_array = self.zarr_group[name]

        return cast(ZarrArray, zarr_array)

    def _create_new_array(
        self, *, name, shape, dtype, fill_value, encoding, attrs
    ) -> ZarrArray:
        if coding.strings.check_vlen_dtype(dtype) is str:
            dtype = str

        if self._write_empty is not None:
            if (
                "write_empty_chunks" in encoding
                and encoding["write_empty_chunks"] != self._write_empty
            ):
                raise ValueError(
                    'Differing "write_empty_chunks" values in encoding and parameters'
                    f'Got {encoding["write_empty_chunks"] = } and {self._write_empty = }'
                )
            else:
                encoding["write_empty_chunks"] = self._write_empty

        if _zarr_v3():
            # zarr v3 deprecated origin and write_empty_chunks
            # instead preferring to pass them via the config argument
            encoding["config"] = {}
            for c in ("write_empty_chunks", "order"):
                if c in encoding:
                    encoding["config"][c] = encoding.pop(c)

        zarr_array = self.zarr_group.create(
            name,
            shape=shape,
            dtype=dtype,
            fill_value=fill_value,
            **encoding,
        )
        zarr_array = _put_attrs(zarr_array, attrs)
        return zarr_array

    def set_variables(
        self,
        variables: dict[str, Variable],
        check_encoding_set,
        writer,
        unlimited_dims=None,
    ):
        """
        This provides a centralized method to set the variables on the data
        store.

        Parameters
        ----------
        variables : dict-like
            Dictionary of key/value (variable name / xr.Variable) pairs
        check_encoding_set : list-like
            List of variables that should be checked for invalid encoding
            values
        writer
        unlimited_dims : list-like
            List of dimension names that should be treated as unlimited
            dimensions.
        """

        existing_keys = self.array_keys()
        is_zarr_v3_format = _zarr_v3() and self.zarr_group.metadata.zarr_format == 3

        for vn, v in variables.items():
            name = _encode_variable_name(vn)
            attrs = v.attrs.copy()
            dims = v.dims
            dtype = v.dtype
            shape = v.shape

            if self._use_zarr_fill_value_as_mask:
                fill_value = attrs.pop("_FillValue", None)
            else:
                fill_value = v.encoding.pop("fill_value", None)
                if "_FillValue" in attrs:
                    # replace with encoded fill value
                    fv = attrs.pop("_FillValue")
                    if fv is not None:
                        attrs["_FillValue"] = FillValueCoder.encode(fv, dtype)

            # _FillValue is never a valid encoding for Zarr
            # TODO: refactor this logic so we don't need to check this here
            if "_FillValue" in v.encoding:
                if v.encoding.get("_FillValue") is not None:
                    raise ValueError("Zarr does not support _FillValue in encoding.")
                else:
                    del v.encoding["_FillValue"]

            zarr_shape = None
            write_region = self._write_region if self._write_region is not None else {}
            write_region = {dim: write_region.get(dim, slice(None)) for dim in dims}

            if self._mode != "w" and name in existing_keys:
                # existing variable
                zarr_array = self._open_existing_array(name=name)
                if self._append_dim is not None and self._append_dim in dims:
                    # resize existing variable
                    append_axis = dims.index(self._append_dim)
                    assert write_region[self._append_dim] == slice(None)
                    write_region[self._append_dim] = slice(
                        zarr_array.shape[append_axis], None
                    )

                    new_shape = (
                        zarr_array.shape[:append_axis]
                        + (zarr_array.shape[append_axis] + v.shape[append_axis],)
                        + zarr_array.shape[append_axis + 1 :]
                    )
                    zarr_array.resize(new_shape)

                zarr_shape = zarr_array.shape
            region = tuple(write_region[dim] for dim in dims)

            # We need to do this for both new and existing variables to ensure we're not
            # writing to a partial chunk, even though we don't use the `encoding` value
            # when writing to an existing variable. See
            # https://github.com/pydata/xarray/issues/8371 for details.
            # Note: Ideally there should be two functions, one for validating the chunks and
            # another one for extracting the encoding.
            encoding = extract_zarr_variable_encoding(
                v,
                raise_on_invalid=vn in check_encoding_set,
                name=vn,
                zarr_format=3 if is_zarr_v3_format else 2,
            )

            if self._align_chunks and isinstance(encoding["chunks"], tuple):
                v = grid_rechunk(
                    v=v,
                    enc_chunks=encoding["chunks"],
                    region=region,
                )

            if self._safe_chunks and isinstance(encoding["chunks"], tuple):
                # the hard case
                # DESIGN CHOICE: do not allow multiple dask chunks on a single zarr chunk
                # this avoids the need to get involved in zarr synchronization / locking
                # From zarr docs:
                #  "If each worker in a parallel computation is writing to a
                #   separate region of the array, and if region boundaries are perfectly aligned
                #   with chunk boundaries, then no synchronization is required."
                # TODO: incorporate synchronizer to allow writes from multiple dask
                # threads
                shape = zarr_shape or v.shape
                validate_grid_chunks_alignment(
                    nd_var_chunks=v.chunks,
                    enc_chunks=encoding["chunks"],
                    region=region,
                    allow_partial_chunks=self._mode != "r+",
                    name=name,
                    backend_shape=shape,
                )

            if self._mode == "w" or name not in existing_keys:
                # new variable
                encoded_attrs = {k: self.encode_attribute(v) for k, v in attrs.items()}
                # the magic for storing the hidden dimension data
                if is_zarr_v3_format:
                    encoding["dimension_names"] = dims
                else:
                    encoded_attrs[DIMENSION_KEY] = dims

                encoding["overwrite"] = self._mode == "w"

                zarr_array = self._create_new_array(
                    name=name,
                    dtype=dtype,
                    shape=shape,
                    fill_value=fill_value,
                    encoding=encoding,
                    attrs=encoded_attrs,
                )

            writer.add(v.data, zarr_array, region)

    def sync(self) -> None:
        pass

    def close(self) -> None:
        if self._close_store_on_close:
            self.zarr_group.store.close()

    def _auto_detect_regions(self, ds, region):
        for dim, val in region.items():
            if val != "auto":
                continue

            if dim not in ds._variables:
                # unindexed dimension
                region[dim] = slice(0, ds.sizes[dim])
                continue

            variable = conventions.decode_cf_variable(
                dim, self.open_store_variable(dim).compute()
            )
            assert variable.dims == (dim,)
            index = pd.Index(variable.data)
            idxs = index.get_indexer(ds[dim].data)
            if (idxs == -1).any():
                raise KeyError(
                    f"Not all values of coordinate '{dim}' in the new array were"
                    " found in the original store. Writing to a zarr region slice"
                    " requires that no dimensions or metadata are changed by the write."
                )

            if (np.diff(idxs) != 1).any():
                raise ValueError(
                    f"The auto-detected region of coordinate '{dim}' for writing new data"
                    " to the original store had non-contiguous indices. Writing to a zarr"
                    " region slice requires that the new data constitute a contiguous subset"
                    " of the original store."
                )
            region[dim] = slice(idxs[0], idxs[-1] + 1)
        return region

    def _validate_and_autodetect_region(self, ds: Dataset) -> Dataset:
        if self._write_region is None:
            return ds

        region = self._write_region

        if region == "auto":
            region = dict.fromkeys(ds.dims, "auto")

        if not isinstance(region, dict):
            raise TypeError(f"``region`` must be a dict, got {type(region)}")
        if any(v == "auto" for v in region.values()):
            if self._mode not in ["r+", "a"]:
                raise ValueError(
                    f"``mode`` must be 'r+' or 'a' when using ``region='auto'``, got {self._mode!r}"
                )
            region = self._auto_detect_regions(ds, region)

        # validate before attempting to auto-detect since the auto-detection
        # should always return a valid slice.
        for k, v in region.items():
            if k not in ds.dims:
                raise ValueError(
                    f"all keys in ``region`` are not in Dataset dimensions, got "
                    f"{list(region)} and {list(ds.dims)}"
                )
            if not isinstance(v, slice):
                raise TypeError(
                    "all values in ``region`` must be slice objects, got "
                    f"region={region}"
                )
            if v.step not in {1, None}:
                raise ValueError(
                    "step on all slices in ``region`` must be 1 or None, got "
                    f"region={region}"
                )

        non_matching_vars = [
            k for k, v in ds.variables.items() if not set(region).intersection(v.dims)
        ]
        if non_matching_vars:
            raise ValueError(
                f"when setting `region` explicitly in to_zarr(), all "
                f"variables in the dataset to write must have at least "
                f"one dimension in common with the region's dimensions "
                f"{list(region.keys())}, but that is not "
                f"the case for some variables here. To drop these variables "
                f"from this dataset before exporting to zarr, write: "
                f".drop_vars({non_matching_vars!r})"
            )

        if self._append_dim is not None and self._append_dim in region:
            raise ValueError(
                f"cannot list the same dimension in both ``append_dim`` and "
                f"``region`` with to_zarr(), got {self._append_dim} in both"
            )

        self._write_region = region

        # can't modify indexes with region writes
        return ds.drop_vars(ds.indexes)

    def _validate_encoding(self, encoding) -> None:
        if encoding and self._mode in ["a", "a-", "r+"]:
            existing_var_names = self.array_keys()
            for var_name in existing_var_names:
                if var_name in encoding:
                    raise ValueError(
                        f"variable {var_name!r} already exists, but encoding was provided"
                    )


def open_zarr(
    store,
    group=None,
    synchronizer=None,
    chunks="auto",
    decode_cf=True,
    mask_and_scale=True,
    decode_times=True,
    concat_characters=True,
    decode_coords=True,
    drop_variables=None,
    consolidated=None,
    overwrite_encoded_chunks=False,
    chunk_store=None,
    storage_options=None,
    decode_timedelta=None,
    use_cftime=None,
    zarr_version=None,
    zarr_format=None,
    use_zarr_fill_value_as_mask=None,
    chunked_array_type: str | None = None,
    from_array_kwargs: dict[str, Any] | None = None,
    create_default_indexes=True,
    **kwargs,
):
    """Load and decode a dataset from a Zarr store.

    The `store` object should be a valid store for a Zarr group. `store`
    variables must contain dimension metadata encoded in the
    `_ARRAY_DIMENSIONS` attribute or must have NCZarr format.

    Parameters
    ----------
    store : MutableMapping or str
        A MutableMapping where a Zarr Group has been stored or a path to a
        directory in file system where a Zarr DirectoryStore has been stored.
    synchronizer : object, optional
        Array synchronizer provided to zarr
    group : str, optional
        Group path. (a.k.a. `path` in zarr terminology.)
    chunks : int, dict, 'auto' or None, default: 'auto'
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
    overwrite_encoded_chunks : bool, optional
        Whether to drop the zarr chunks encoded for each variable when a
        dataset is loaded with specified chunk sizes (default: False)
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
        be replaced by NA.
    decode_times : bool, optional
        If True, decode times encoded in the standard NetCDF datetime format
        into datetime objects. Otherwise, leave them encoded as numbers.
    concat_characters : bool, optional
        If True, concatenate along the last dimension of character arrays to
        form string arrays. Dimensions will only be concatenated over (and
        removed) if they have no corresponding variable and if they are only
        used as the last dimension of character arrays.
    decode_coords : bool, optional
        If True, decode the 'coordinates' attribute to identify coordinates in
        the resulting dataset.
    drop_variables : str or iterable, optional
        A variable or list of variables to exclude from being parsed from the
        dataset. This may be useful to drop variables with problems or
        inconsistent values.
    consolidated : bool, optional
        Whether to open the store using zarr's consolidated metadata
        capability. Only works for stores that have already been consolidated.
        By default (`consolidate=None`), attempts to read consolidated metadata,
        falling back to read non-consolidated metadata if that fails.

        When the experimental ``zarr_version=3``, ``consolidated`` must be
        either be ``None`` or ``False``.
    chunk_store : MutableMapping, optional
        A separate Zarr store only for chunk data.
    storage_options : dict, optional
        Any additional parameters for the storage backend (ignored for local
        paths).
    decode_timedelta : bool, optional
        If True, decode variables and coordinates with time units in
        {'days', 'hours', 'minutes', 'seconds', 'milliseconds', 'microseconds'}
        into timedelta objects. If False, leave them encoded as numbers.
        If None (default), assume the same value of decode_time.
    use_cftime : bool, optional
        Only relevant if encoded dates come from a standard calendar
        (e.g. "gregorian", "proleptic_gregorian", "standard", or not
        specified).  If None (default), attempt to decode times to
        ``np.datetime64[ns]`` objects; if this is not possible, decode times to
        ``cftime.datetime`` objects. If True, always decode times to
        ``cftime.datetime`` objects, regardless of whether or not they can be
        represented using ``np.datetime64[ns]`` objects.  If False, always
        decode times to ``np.datetime64[ns]`` objects; if this is not possible
        raise an error.
    zarr_version : int or None, optional

        .. deprecated:: 2024.9.1
           Use ``zarr_format`` instead.

    zarr_format : int or None, optional
        The desired zarr format to target (currently 2 or 3). The default
        of None will attempt to determine the zarr version from ``store`` when
        possible, otherwise defaulting to the default version used by
        the zarr-python library installed.
    use_zarr_fill_value_as_mask : bool, optional
        If True, use the zarr Array ``fill_value`` to mask the data, the same as done
        for NetCDF data with ``_FillValue`` or ``missing_value`` attributes. If False,
        the ``fill_value`` is ignored and the data are not masked. If None, this defaults
        to True for ``zarr_version=2`` and False for ``zarr_version=3``.
    chunked_array_type: str, optional
        Which chunked array type to coerce this datasets' arrays to.
        Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEntryPoint` system.
        Experimental API that should not be relied upon.
    from_array_kwargs: dict, optional
        Additional keyword arguments passed on to the ``ChunkManagerEntrypoint.from_array`` method used to create
        chunked arrays, via whichever chunk manager is specified through the ``chunked_array_type`` kwarg.
        Defaults to ``{'manager': 'dask'}``, meaning additional kwargs will be passed eventually to
        :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
    create_default_indexes : bool, default: True
        If True, create pandas indexes for :term:`dimension coordinates <dimension coordinate>`,
        which loads the coordinate data into memory. Set it to False if you want to avoid loading
        data into memory.

        Note that backends can still choose to create other indexes. If you want to control that,
        please refer to the backend's documentation.

    Returns
    -------
    dataset : Dataset
        The newly created dataset.

    See Also
    --------
    open_dataset
    open_mfdataset

    References
    ----------
    https://zarr.readthedocs.io/
    """
    from xarray.backends.api import open_dataset

    if from_array_kwargs is None:
        from_array_kwargs = {}

    if chunks == "auto":
        try:
            guess_chunkmanager(
                chunked_array_type
            )  # attempt to import that parallel backend

            chunks = {}
        except (ValueError, ImportError):
            chunks = None

    if kwargs:
        raise TypeError(
            "open_zarr() got unexpected keyword arguments " + ",".join(kwargs.keys())
        )

    backend_kwargs = {
        "synchronizer": synchronizer,
        "consolidated": consolidated,
        "overwrite_encoded_chunks": overwrite_encoded_chunks,
        "chunk_store": chunk_store,
        "storage_options": storage_options,
        "zarr_version": zarr_version,
        "zarr_format": zarr_format,
    }

    ds = open_dataset(
        filename_or_obj=store,
        group=group,
        decode_cf=decode_cf,
        mask_and_scale=mask_and_scale,
        decode_times=decode_times,
        concat_characters=concat_characters,
        decode_coords=decode_coords,
        engine="zarr",
        chunks=chunks,
        drop_variables=drop_variables,
        create_default_indexes=create_default_indexes,
        chunked_array_type=chunked_array_type,
        from_array_kwargs=from_array_kwargs,
        backend_kwargs=backend_kwargs,
        decode_timedelta=decode_timedelta,
        use_cftime=use_cftime,
        zarr_version=zarr_version,
        use_zarr_fill_value_as_mask=use_zarr_fill_value_as_mask,
    )
    return ds


class ZarrBackendEntrypoint(BackendEntrypoint):
    """
    Backend for ".zarr" files based on the zarr package.

    For more information about the underlying library, visit:
    https://zarr.readthedocs.io/en/stable

    See Also
    --------
    backends.ZarrStore
    """

    description = "Open zarr files (.zarr) using zarr in Xarray"
    url = "https://docs.xarray.dev/en/stable/generated/xarray.backends.ZarrBackendEntrypoint.html"

    def guess_can_open(self, filename_or_obj: T_PathFileOrDataStore) -> bool:
        if isinstance(filename_or_obj, str | os.PathLike):
            # allow a trailing slash to account for an autocomplete
            # adding it.
            _, ext = os.path.splitext(str(filename_or_obj).rstrip("/"))
            return ext in [".zarr"]

        return False

    def open_dataset(
        self,
        filename_or_obj: T_PathFileOrDataStore,
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime=None,
        decode_timedelta=None,
        group=None,
        mode="r",
        synchronizer=None,
        consolidated=None,
        chunk_store=None,
        storage_options=None,
        zarr_version=None,
        zarr_format=None,
        store=None,
        engine=None,
        use_zarr_fill_value_as_mask=None,
        cache_members: bool = True,
    ) -> Dataset:
        filename_or_obj = _normalize_path(filename_or_obj)
        if not store:
            store = ZarrStore.open_group(
                filename_or_obj,
                group=group,
                mode=mode,
                synchronizer=synchronizer,
                consolidated=consolidated,
                consolidate_on_close=False,
                chunk_store=chunk_store,
                storage_options=storage_options,
                zarr_version=zarr_version,
                use_zarr_fill_value_as_mask=None,
                zarr_format=zarr_format,
                cache_members=cache_members,
            )

        store_entrypoint = StoreBackendEntrypoint()
        with close_on_error(store):
            ds = store_entrypoint.open_dataset(
                store,
                mask_and_scale=mask_and_scale,
                decode_times=decode_times,
                concat_characters=concat_characters,
                decode_coords=decode_coords,
                drop_variables=drop_variables,
                use_cftime=use_cftime,
                decode_timedelta=decode_timedelta,
            )
        return ds

    def open_datatree(
        self,
        filename_or_obj: T_PathFileOrDataStore,
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime=None,
        decode_timedelta=None,
        group: str | None = None,
        mode="r",
        synchronizer=None,
        consolidated=None,
        chunk_store=None,
        storage_options=None,
        zarr_version=None,
        zarr_format=None,
    ) -> DataTree:
        filename_or_obj = _normalize_path(filename_or_obj)
        groups_dict = self.open_groups_as_dict(
            filename_or_obj=filename_or_obj,
            mask_and_scale=mask_and_scale,
            decode_times=decode_times,
            concat_characters=concat_characters,
            decode_coords=decode_coords,
            drop_variables=drop_variables,
            use_cftime=use_cftime,
            decode_timedelta=decode_timedelta,
            group=group,
            mode=mode,
            synchronizer=synchronizer,
            consolidated=consolidated,
            chunk_store=chunk_store,
            storage_options=storage_options,
            zarr_version=zarr_version,
            zarr_format=zarr_format,
        )

        return datatree_from_dict_with_io_cleanup(groups_dict)

    def open_groups_as_dict(
        self,
        filename_or_obj: T_PathFileOrDataStore,
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime=None,
        decode_timedelta=None,
        group: str | None = None,
        mode="r",
        synchronizer=None,
        consolidated=None,
        chunk_store=None,
        storage_options=None,
        zarr_version=None,
        zarr_format=None,
    ) -> dict[str, Dataset]:
        filename_or_obj = _normalize_path(filename_or_obj)

        # Check for a group and make it a parent if it exists
        if group:
            parent = str(NodePath("/") / NodePath(group))
        else:
            parent = str(NodePath("/"))

        stores = ZarrStore.open_store(
            filename_or_obj,
            group=parent,
            mode=mode,
            synchronizer=synchronizer,
            consolidated=consolidated,
            consolidate_on_close=False,
            chunk_store=chunk_store,
            storage_options=storage_options,
            zarr_version=zarr_version,
            zarr_format=zarr_format,
        )

        groups_dict = {}
        for path_group, store in stores.items():
            store_entrypoint = StoreBackendEntrypoint()

            with close_on_error(store):
                group_ds = store_entrypoint.open_dataset(
                    store,
                    mask_and_scale=mask_and_scale,
                    decode_times=decode_times,
                    concat_characters=concat_characters,
                    decode_coords=decode_coords,
                    drop_variables=drop_variables,
                    use_cftime=use_cftime,
                    decode_timedelta=decode_timedelta,
                )
            if group:
                group_name = str(NodePath(path_group).relative_to(parent))
            else:
                group_name = str(NodePath(path_group))
            groups_dict[group_name] = group_ds
        return groups_dict


def _iter_zarr_groups(root: ZarrGroup, parent: str = "/") -> Iterable[str]:
    parent_nodepath = NodePath(parent)
    yield str(parent_nodepath)
    for path, group in root.groups():
        gpath = parent_nodepath / path
        yield from _iter_zarr_groups(group, parent=str(gpath))


def _get_open_params(
    store,
    mode,
    synchronizer,
    group,
    consolidated,
    consolidate_on_close,
    chunk_store,
    storage_options,
    zarr_version,
    use_zarr_fill_value_as_mask,
    zarr_format,
):
    if TYPE_CHECKING:
        import zarr
    else:
        zarr = attempt_import("zarr")

    # zarr doesn't support pathlib.Path objects yet. zarr-python#601
    if isinstance(store, os.PathLike):
        store = os.fspath(store)

    open_kwargs = dict(
        # mode='a-' is a handcrafted xarray specialty
        mode="a" if mode == "a-" else mode,
        synchronizer=synchronizer,
        path=group,
    )
    open_kwargs["storage_options"] = storage_options

    zarr_format = _handle_zarr_version_or_format(
        zarr_version=zarr_version, zarr_format=zarr_format
    )

    if _zarr_v3():
        open_kwargs["zarr_format"] = zarr_format
    else:
        open_kwargs["zarr_version"] = zarr_format

    if chunk_store is not None:
        open_kwargs["chunk_store"] = chunk_store
        if consolidated is None:
            consolidated = False

    if _zarr_v3():
        # TODO: replace AssertionError after https://github.com/zarr-developers/zarr-python/issues/2821 is resolved
        missing_exc = AssertionError
    else:
        missing_exc = zarr.errors.GroupNotFoundError

    if _zarr_v3():
        # zarr 3.0.8 and earlier did not support this property - it was effectively assumed true
        if not getattr(store, "supports_consolidated_metadata", True):
            consolidated = consolidate_on_close = False

    if consolidated in [None, True]:
        # open the root of the store, in case there is metadata consolidated there
        group = open_kwargs.pop("path")

        if consolidated:
            # TODO: an option to pass the metadata_key keyword
            zarr_root_group = zarr.open_consolidated(store, **open_kwargs)
        elif consolidated is None:
            # same but with more error handling in case no consolidated metadata found
            try:
                zarr_root_group = zarr.open_consolidated(store, **open_kwargs)
            except (ValueError, KeyError):
                # ValueError in zarr-python 3.x, KeyError in 2.x.
                try:
                    zarr_root_group = zarr.open_group(store, **open_kwargs)
                    emit_user_level_warning(
                        "Failed to open Zarr store with consolidated metadata, "
                        "but successfully read with non-consolidated metadata. "
                        "This is typically much slower for opening a dataset. "
                        "To silence this warning, consider:\n"
                        "1. Consolidating metadata in this existing store with "
                        "zarr.consolidate_metadata().\n"
                        "2. Explicitly setting consolidated=False, to avoid trying "
                        "to read consolidate metadata, or\n"
                        "3. Explicitly setting consolidated=True, to raise an "
                        "error in this case instead of falling back to try "
                        "reading non-consolidated metadata.",
                        RuntimeWarning,
                    )
                except missing_exc as err:
                    raise FileNotFoundError(
                        f"No such file or directory: '{store}'"
                    ) from err

        # but the user should still receive a DataTree whose root is the group they asked for
        if group and group != "/":
            zarr_group = zarr_root_group[group.removeprefix("/")]
        else:
            zarr_group = zarr_root_group
    else:
        if _zarr_v3():
            # we have determined that we don't want to use consolidated metadata
            # so we set that to False to avoid trying to read it
            open_kwargs["use_consolidated"] = False
        zarr_group = zarr.open_group(store, **open_kwargs)

    close_store_on_close = zarr_group.store is not store

    # we use this to determine how to handle fill_value
    is_zarr_v3_format = _zarr_v3() and zarr_group.metadata.zarr_format == 3
    if use_zarr_fill_value_as_mask is None:
        if is_zarr_v3_format:
            # for new data, we use a better default
            use_zarr_fill_value_as_mask = False
        else:
            # this was the default for v2 and should apply to most existing Zarr data
            use_zarr_fill_value_as_mask = True

    return (
        zarr_group,
        consolidate_on_close,
        close_store_on_close,
        use_zarr_fill_value_as_mask,
    )


def _handle_zarr_version_or_format(
    *, zarr_version: ZarrFormat | None, zarr_format: ZarrFormat | None
) -> ZarrFormat | None:
    """handle the deprecated zarr_version kwarg and return zarr_format"""
    if (
        zarr_format is not None
        and zarr_version is not None
        and zarr_format != zarr_version
    ):
        raise ValueError(
            f"zarr_format {zarr_format} does not match zarr_version {zarr_version}, please only set one"
        )
    if zarr_version is not None:
        emit_user_level_warning(
            "zarr_version is deprecated, use zarr_format", FutureWarning
        )
        return zarr_version
    return zarr_format


BACKEND_ENTRYPOINTS["zarr"] = ("zarr", ZarrBackendEntrypoint)
