from __future__ import annotations

import functools
import operator
import os
from collections.abc import Iterable
from contextlib import suppress
from dataclasses import dataclass
from io import IOBase
from typing import TYPE_CHECKING, Any, Self

import numpy as np

from xarray.backends.common import (
    BACKEND_ENTRYPOINTS,
    BackendArray,
    BackendEntrypoint,
    BytesIOProxy,
    T_PathFileOrDataStore,
    WritableCFDataStore,
    _normalize_path,
    collect_ancestor_dimensions,
    datatree_from_dict_with_io_cleanup,
    find_root_and_group,
    robust_getitem,
)
from xarray.backends.file_manager import (
    CachingFileManager,
    DummyFileManager,
    PickleableFileManager,
)
from xarray.backends.locks import (
    HDF5_LOCK,
    NETCDFC_LOCK,
    combine_locks,
    ensure_lock,
    get_write_lock,
)
from xarray.backends.netcdf3 import encode_nc3_attr_value, encode_nc3_variable
from xarray.backends.store import StoreBackendEntrypoint
from xarray.coding.strings import (
    CharacterArrayCoder,
    EncodedStringCoder,
    create_vlen_dtype,
    is_unicode_dtype,
)
from xarray.coding.variables import pop_to
from xarray.core import indexing
from xarray.core.utils import (
    FrozenDict,
    close_on_error,
    is_remote_uri,
    strip_uri_params,
    try_read_magic_number_from_path,
)
from xarray.core.variable import Variable

if TYPE_CHECKING:
    import netCDF4
    from h5netcdf.core import EnumType as h5EnumType
    from netCDF4 import EnumType as ncEnumType

    from xarray.core.dataset import Dataset
    from xarray.core.datatree import DataTree

# This lookup table maps from dtype.byteorder to a readable endian
# string used by netCDF4.
_endian_lookup = {"=": "native", ">": "big", "<": "little", "|": "native"}

NETCDF4_PYTHON_LOCK = combine_locks([NETCDFC_LOCK, HDF5_LOCK])


class BaseNetCDF4Array(BackendArray):
    __slots__ = ("datastore", "dtype", "shape", "variable_name")

    def __init__(self, variable_name, datastore):
        self.datastore = datastore
        self.variable_name = variable_name

        array = self.get_array()
        self.shape = array.shape

        dtype = array.dtype
        if dtype is str:
            # use object dtype (with additional vlen string metadata) because that's
            # the only way in numpy to represent variable length strings and to
            # check vlen string dtype in further steps
            # it also prevents automatic string concatenation via
            # conventions.decode_cf_variable
            dtype = create_vlen_dtype(str)
        self.dtype = dtype

    def __setitem__(self, key, value):
        with self.datastore.lock:
            data = self.get_array(needs_lock=False)
            data[key] = value
            if self.datastore.autoclose:
                self.datastore.close(needs_lock=False)

    def get_array(self, needs_lock=True):
        raise NotImplementedError("Virtual Method")


class NetCDF4ArrayWrapper(BaseNetCDF4Array):
    __slots__ = ()

    def get_array(self, needs_lock=True):
        ds = self.datastore._acquire(needs_lock)
        variable = ds.variables[self.variable_name]
        variable.set_auto_maskandscale(False)
        # only added in netCDF4-python v1.2.8
        with suppress(AttributeError):
            variable.set_auto_chartostring(False)
        return variable

    def __getitem__(self, key):
        return indexing.explicit_indexing_adapter(
            key, self.shape, indexing.IndexingSupport.OUTER, self._getitem
        )

    def _getitem(self, key):
        if self.datastore.is_remote:  # pragma: no cover
            getitem = functools.partial(robust_getitem, catch=RuntimeError)
        else:
            getitem = operator.getitem

        try:
            with self.datastore.lock:
                original_array = self.get_array(needs_lock=False)
                array = getitem(original_array, key)
        except IndexError as err:
            # Catch IndexError in netCDF4 and return a more informative
            # error message.  This is most often called when an unsorted
            # indexer is used before the data is loaded from disk.
            msg = (
                "The indexing operation you are attempting to perform "
                "is not valid on netCDF4.Variable object. Try loading "
                "your data into memory first by calling .load()."
            )
            raise IndexError(msg) from err
        return array


def _encode_nc4_variable(var, name=None):
    for coder in [
        EncodedStringCoder(allows_unicode=True),
        CharacterArrayCoder(),
    ]:
        var = coder.encode(var, name=name)
    return var


def _check_encoding_dtype_is_vlen_string(dtype):
    if dtype is not str:
        raise AssertionError(  # pragma: no cover
            f"unexpected dtype encoding {dtype!r}. This shouldn't happen: please "
            "file a bug report at github.com/pydata/xarray"
        )


def _get_datatype(
    var, nc_format="NETCDF4", raise_on_invalid_encoding=False
) -> np.dtype:
    if nc_format == "NETCDF4":
        return _nc4_dtype(var)
    if "dtype" in var.encoding:
        encoded_dtype = var.encoding["dtype"]
        _check_encoding_dtype_is_vlen_string(encoded_dtype)
        if raise_on_invalid_encoding:
            raise ValueError(
                "encoding dtype=str for vlen strings is only supported "
                "with format='NETCDF4'."
            )
    return var.dtype


def _nc4_dtype(var):
    if "dtype" in var.encoding:
        dtype = var.encoding.pop("dtype")
        _check_encoding_dtype_is_vlen_string(dtype)
    elif is_unicode_dtype(var.dtype):
        dtype = str
    elif var.dtype.kind in ["i", "u", "f", "c", "S"]:
        dtype = var.dtype
    else:
        raise ValueError(f"unsupported dtype for netCDF4 variable: {var.dtype}")
    return dtype


def _netcdf4_create_group(dataset, name):
    return dataset.createGroup(name)


def _nc4_require_group(ds, group, mode, create_group=_netcdf4_create_group):
    if group in {None, "", "/"}:
        # use the root group
        return ds
    else:
        # make sure it's a string
        if not isinstance(group, str):
            raise ValueError("group must be a string or None")
        # support path-like syntax
        path = group.strip("/").split("/")
        for key in path:
            try:
                ds = ds.groups[key]
            except KeyError as e:
                if mode != "r":
                    ds = create_group(ds, key)
                else:
                    # wrap error to provide slightly more helpful message
                    raise OSError(f"group not found: {key}", e) from e
        return ds


def _ensure_no_forward_slash_in_name(name):
    if "/" in name:
        raise ValueError(
            f"Forward slashes '/' are not allowed in variable and dimension names (got {name!r}). "
            "Forward slashes are used as hierarchy-separators for "
            "HDF5-based files ('netcdf4'/'h5netcdf')."
        )


def _ensure_fill_value_valid(data, attributes):
    # work around for netCDF4/scipy issue where _FillValue has the wrong type:
    # https://github.com/Unidata/netcdf4-python/issues/271
    if data.dtype.kind == "S" and "_FillValue" in attributes:
        attributes["_FillValue"] = np.bytes_(attributes["_FillValue"])


def _force_native_endianness(var):
    # possible values for byteorder are:
    #     =    native
    #     <    little-endian
    #     >    big-endian
    #     |    not applicable
    # Below we check if the data type is not native or NA
    if var.dtype.byteorder not in ["=", "|"]:
        # if endianness is specified explicitly, convert to the native type
        data = var.data.astype(var.dtype.newbyteorder("="))
        var = Variable(var.dims, data, var.attrs, var.encoding)
        # if endian exists, remove it from the encoding.
        var.encoding.pop("endian", None)
    # check to see if encoding has a value for endian its 'native'
    if var.encoding.get("endian", "native") != "native":
        raise NotImplementedError(
            "Attempt to write non-native endian type, "
            "this is not supported by the netCDF4 "
            "python library."
        )
    return var


def _extract_nc4_variable_encoding(
    variable: Variable,
    raise_on_invalid=False,
    lsd_okay=True,
    h5py_okay=False,
    backend="netCDF4",
    unlimited_dims=None,
) -> dict[str, Any]:
    if unlimited_dims is None:
        unlimited_dims = ()

    encoding = variable.encoding.copy()

    safe_to_drop = {"source", "original_shape"}
    valid_encodings = {
        "zlib",
        "complevel",
        "fletcher32",
        "contiguous",
        "chunksizes",
        "shuffle",
        "_FillValue",
        "dtype",
        "compression",
        "significant_digits",
        "quantize_mode",
        "blosc_shuffle",
        "szip_coding",
        "szip_pixels_per_block",
        "endian",
    }
    if lsd_okay:
        valid_encodings.add("least_significant_digit")
    if h5py_okay:
        valid_encodings.add("compression_opts")

    if not raise_on_invalid and encoding.get("chunksizes") is not None:
        # It's possible to get encoded chunksizes larger than a dimension size
        # if the original file had an unlimited dimension. This is problematic
        # if the new file no longer has an unlimited dimension.
        chunksizes = encoding["chunksizes"]
        chunks_too_big = any(
            c > d and dim not in unlimited_dims
            for c, d, dim in zip(
                chunksizes, variable.shape, variable.dims, strict=False
            )
        )
        has_original_shape = "original_shape" in encoding
        changed_shape = (
            has_original_shape and encoding.get("original_shape") != variable.shape
        )
        if chunks_too_big or changed_shape:
            del encoding["chunksizes"]

    var_has_unlim_dim = any(dim in unlimited_dims for dim in variable.dims)
    if not raise_on_invalid and var_has_unlim_dim and "contiguous" in encoding.keys():
        del encoding["contiguous"]

    for k in safe_to_drop:
        if k in encoding:
            del encoding[k]

    if raise_on_invalid:
        invalid = [k for k in encoding if k not in valid_encodings]
        if invalid:
            raise ValueError(
                f"unexpected encoding parameters for {backend!r} backend: {invalid!r}. Valid "
                f"encodings are: {valid_encodings!r}"
            )
    else:
        for k in list(encoding):
            if k not in valid_encodings:
                del encoding[k]

    return encoding


def _is_list_of_strings(value) -> bool:
    arr = np.asarray(value)
    return arr.dtype.kind in ["U", "S"] and arr.size > 1


def _build_and_get_enum(
    store, var_name: str, dtype: np.dtype, enum_name: str, enum_dict: dict[str, int]
) -> ncEnumType | h5EnumType:
    """
    Add or get the netCDF4 Enum based on the dtype in encoding.
    The return type should be ``netCDF4.EnumType``,
    but we avoid importing netCDF4 globally for performances.
    """
    if enum_name not in store.ds.enumtypes:
        create_func = (
            store.ds.createEnumType
            if isinstance(store, NetCDF4DataStore)
            else store.ds.create_enumtype
        )
        return create_func(
            dtype,
            enum_name,
            enum_dict,
        )
    datatype = store.ds.enumtypes[enum_name]
    if datatype.enum_dict != enum_dict:
        error_msg = (
            f"Cannot save variable `{var_name}` because an enum"
            f" `{enum_name}` already exists in the Dataset but has"
            " a different definition. To fix this error, make sure"
            " all variables have a uniquely named enum in their"
            " `encoding['dtype'].metadata` or, if they should share"
            " the same enum type, make sure the enums are identical."
        )
        raise ValueError(error_msg)
    return datatype


@dataclass
class _Thunk:
    """Pickleable equivalent of `lambda: value`."""

    value: Any

    def __call__(self):
        return self.value


@dataclass
class _CloseWithCopy:
    """Wrapper around netCDF4's esoteric interface for writing in-memory data."""

    proxy: BytesIOProxy
    nc4_dataset: netCDF4.Dataset

    def __call__(self):
        value = self.nc4_dataset.close()
        self.proxy.getvalue = _Thunk(value)


class NetCDF4DataStore(WritableCFDataStore):
    """Store for reading and writing data via the Python-NetCDF4 library.

    This store supports NetCDF3, NetCDF4 and OpenDAP datasets.
    """

    __slots__ = (
        "_filename",
        "_group",
        "_manager",
        "_mode",
        "autoclose",
        "format",
        "is_remote",
        "lock",
    )

    def __init__(
        self, manager, group=None, mode=None, lock=NETCDF4_PYTHON_LOCK, autoclose=False
    ):
        import netCDF4

        if isinstance(manager, netCDF4.Dataset):
            if group is None:
                root, group = find_root_and_group(manager)
            else:
                if type(manager) is not netCDF4.Dataset:
                    raise ValueError(
                        "must supply a root netCDF4.Dataset if the group "
                        "argument is provided"
                    )
                root = manager
            manager = DummyFileManager(root, lock=NETCDF4_PYTHON_LOCK)

        self._manager = manager
        self._group = group
        self._mode = mode
        self.format = self.ds.data_model
        self._filename = self.ds.filepath()
        self.is_remote = is_remote_uri(self._filename)
        self.lock = ensure_lock(lock)
        self.autoclose = autoclose

    def get_child_store(self, group: str) -> Self:
        if self._group is not None:
            group = os.path.join(self._group, group)
        return type(self)(
            self._manager,
            group=group,
            mode=self._mode,
            lock=self.lock,
            autoclose=self.autoclose,
        )

    @classmethod
    def open(
        cls,
        filename,
        mode="r",
        format="NETCDF4",
        group=None,
        clobber=True,
        diskless=False,
        persist=False,
        auto_complex=None,
        lock=None,
        lock_maker=None,
        autoclose=False,
    ):
        import netCDF4

        if isinstance(filename, os.PathLike):
            filename = os.fspath(filename)

        if isinstance(filename, IOBase):
            raise TypeError(
                f"file objects are not supported by the netCDF4 backend: {filename}"
            )

        if not isinstance(filename, str | bytes | memoryview | BytesIOProxy):
            raise TypeError(f"invalid filename for netCDF4 backend: {filename}")

        if format is None:
            format = "NETCDF4"

        if lock is None:
            if mode == "r":
                if isinstance(filename, str) and is_remote_uri(filename):
                    lock = NETCDFC_LOCK
                else:
                    lock = NETCDF4_PYTHON_LOCK
            else:
                if format is None or format.startswith("NETCDF4"):
                    lock = NETCDF4_PYTHON_LOCK
                else:
                    lock = NETCDFC_LOCK

                if isinstance(filename, str):
                    lock = combine_locks([lock, get_write_lock(filename)])

        kwargs = dict(
            clobber=clobber,
            diskless=diskless,
            persist=persist,
            format=format,
        )
        if auto_complex is not None:
            kwargs["auto_complex"] = auto_complex

        if isinstance(filename, BytesIOProxy):
            assert mode == "w"
            # Size hint used for creating netCDF3 files. Per the documentation
            # for nc__create(), the special value NC_SIZEHINT_DEFAULT (which is
            # the value 0), lets the netcdf library choose a suitable initial
            # size.
            memory = 0
            kwargs["diskless"] = False
            nc4_dataset = netCDF4.Dataset(
                "<xarray-in-memory-write>", mode=mode, memory=memory, **kwargs
            )
            close = _CloseWithCopy(filename, nc4_dataset)
            manager = DummyFileManager(nc4_dataset, close=close, lock=lock)

        elif isinstance(filename, bytes | memoryview):
            assert mode == "r"
            kwargs["memory"] = filename
            manager = PickleableFileManager(
                netCDF4.Dataset,
                "<xarray-in-memory-read>",
                mode=mode,
                kwargs=kwargs,
                lock=lock,
            )
        else:
            manager = CachingFileManager(
                netCDF4.Dataset, filename, mode=mode, kwargs=kwargs, lock=lock
            )
        return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)

    def _acquire(self, needs_lock=True):
        with self._manager.acquire_context(needs_lock) as root:
            ds = _nc4_require_group(root, self._group, self._mode)
        return ds

    @property
    def ds(self):
        return self._acquire()

    def open_store_variable(self, name: str, var):
        import netCDF4

        dimensions = var.dimensions
        attributes = {k: var.getncattr(k) for k in var.ncattrs()}
        data = indexing.LazilyIndexedArray(NetCDF4ArrayWrapper(name, self))
        encoding: dict[str, Any] = {}
        if isinstance(var.datatype, netCDF4.EnumType):
            encoding["dtype"] = np.dtype(
                data.dtype,
                metadata={
                    "enum": var.datatype.enum_dict,
                    "enum_name": var.datatype.name,
                },
            )
        else:
            encoding["dtype"] = var.dtype
        _ensure_fill_value_valid(data, attributes)
        # netCDF4 specific encoding; save _FillValue for later
        filters = var.filters()
        if filters is not None:
            encoding.update(filters)
        chunking = var.chunking()
        if chunking is not None:
            if chunking == "contiguous":
                encoding["contiguous"] = True
                encoding["chunksizes"] = None
            else:
                encoding["contiguous"] = False
                encoding["chunksizes"] = tuple(chunking)
                encoding["preferred_chunks"] = dict(
                    zip(var.dimensions, chunking, strict=True)
                )
        # TODO: figure out how to round-trip "endian-ness" without raising
        # warnings from netCDF4
        # encoding['endian'] = var.endian()
        pop_to(attributes, encoding, "least_significant_digit")
        # save source so __repr__ can detect if it's local or not
        encoding["source"] = self._filename
        encoding["original_shape"] = data.shape

        return Variable(dimensions, data, attributes, encoding)

    def get_variables(self):
        return FrozenDict(
            (k, self.open_store_variable(k, v)) for k, v in self.ds.variables.items()
        )

    def get_attrs(self):
        return FrozenDict((k, self.ds.getncattr(k)) for k in self.ds.ncattrs())

    def get_dimensions(self):
        return FrozenDict((k, len(v)) for k, v in self.ds.dimensions.items())

    def get_parent_dimensions(self):
        return FrozenDict(collect_ancestor_dimensions(self.ds))

    def get_encoding(self):
        return {
            "unlimited_dims": {
                k for k, v in self.ds.dimensions.items() if v.isunlimited()
            }
        }

    def set_dimension(self, name, length, is_unlimited=False):
        _ensure_no_forward_slash_in_name(name)
        dim_length = length if not is_unlimited else None
        self.ds.createDimension(name, size=dim_length)

    def set_attribute(self, key, value):
        if self.format != "NETCDF4":
            value = encode_nc3_attr_value(value)
        if _is_list_of_strings(value):
            # encode as NC_STRING if attr is list of strings
            self.ds.setncattr_string(key, value)
        else:
            self.ds.setncattr(key, value)

    def encode_variable(self, variable, name=None):
        variable = _force_native_endianness(variable)
        if self.format == "NETCDF4":
            variable = _encode_nc4_variable(variable, name=name)
        else:
            variable = encode_nc3_variable(variable, name=name)
        return variable

    def prepare_variable(
        self, name, variable: Variable, check_encoding=False, unlimited_dims=None
    ):
        _ensure_no_forward_slash_in_name(name)
        attrs = variable.attrs.copy()
        fill_value = attrs.pop("_FillValue", None)
        datatype: np.dtype | ncEnumType | h5EnumType
        datatype = _get_datatype(
            variable, self.format, raise_on_invalid_encoding=check_encoding
        )
        # check enum metadata and use netCDF4.EnumType
        if (
            (meta := np.dtype(datatype).metadata)
            and (e_name := meta.get("enum_name"))
            and (e_dict := meta.get("enum"))
        ):
            datatype = _build_and_get_enum(self, name, datatype, e_name, e_dict)
        encoding = _extract_nc4_variable_encoding(
            variable, raise_on_invalid=check_encoding, unlimited_dims=unlimited_dims
        )
        if name in self.ds.variables:
            nc4_var = self.ds.variables[name]
        else:
            default_args = dict(
                varname=name,
                datatype=datatype,
                dimensions=variable.dims,
                zlib=False,
                complevel=4,
                shuffle=True,
                fletcher32=False,
                contiguous=False,
                chunksizes=None,
                endian="native",
                least_significant_digit=None,
                fill_value=fill_value,
            )
            default_args.update(encoding)
            default_args.pop("_FillValue", None)
            nc4_var = self.ds.createVariable(**default_args)

        nc4_var.setncatts(attrs)

        target = NetCDF4ArrayWrapper(name, self)

        return target, variable.data

    def sync(self):
        self.ds.sync()

    def close(self, **kwargs):
        self._manager.close(**kwargs)


class NetCDF4BackendEntrypoint(BackendEntrypoint):
    """
    Backend for netCDF files based on the netCDF4 package.

    It can open ".nc", ".nc4", ".cdf" files and will be chosen
    as default for these files.

    Additionally it can open valid HDF5 files, see
    https://h5netcdf.org/#invalid-netcdf-files for more info.
    It will not be detected as valid backend for such files, so make
    sure to specify ``engine="netcdf4"`` in ``open_dataset``.

    For more information about the underlying library, visit:
    https://unidata.github.io/netcdf4-python

    See Also
    --------
    backends.NetCDF4DataStore
    backends.H5netcdfBackendEntrypoint
    backends.ScipyBackendEntrypoint
    """

    description = (
        "Open netCDF (.nc, .nc4 and .cdf) and most HDF5 files using netCDF4 in Xarray"
    )
    url = "https://docs.xarray.dev/en/stable/generated/xarray.backends.NetCDF4BackendEntrypoint.html"
    supports_groups = True

    def guess_can_open(self, filename_or_obj: T_PathFileOrDataStore) -> bool:
        # Helper to check if magic number is netCDF or HDF5
        def _is_netcdf_magic(magic: bytes) -> bool:
            return magic.startswith((b"CDF", b"\211HDF\r\n\032\n"))

        # Helper to check if extension is netCDF
        def _has_netcdf_ext(path: str | os.PathLike, is_remote: bool = False) -> bool:
            path = str(path).rstrip("/")
            # For remote URIs, strip query parameters and fragments
            if is_remote:
                path = strip_uri_params(path)
            _, ext = os.path.splitext(path)
            return ext in {".nc", ".nc4", ".cdf"}

        if isinstance(filename_or_obj, str):
            if is_remote_uri(filename_or_obj):
                # For remote URIs, check extension (accounting for query params/fragments)
                # Remote netcdf-c can handle both regular URLs and DAP URLs
                if _has_netcdf_ext(filename_or_obj, is_remote=True):
                    return True
                elif "zarr" in filename_or_obj.lower():
                    return False
                # return true for non-zarr URLs so we don't have a breaking change for people relying on this
                # netcdf backend guessing true for all remote sources.
                # TODO: emit a warning here about deprecation of this behavior
                # https://github.com/pydata/xarray/pull/10931
                return True

        if isinstance(filename_or_obj, str | os.PathLike):
            # For local paths, check magic number first, then extension
            magic_number = try_read_magic_number_from_path(filename_or_obj)
            if magic_number is not None:
                return _is_netcdf_magic(magic_number)
            # No magic number available, fallback to extension
            return _has_netcdf_ext(filename_or_obj)

        if isinstance(filename_or_obj, bytes | memoryview):
            return _is_netcdf_magic(bytes(filename_or_obj[:8]))

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
        format="NETCDF4",
        clobber=True,
        diskless=False,
        persist=False,
        auto_complex=None,
        lock=None,
        autoclose=False,
    ) -> Dataset:
        filename_or_obj = _normalize_path(filename_or_obj)
        store = NetCDF4DataStore.open(
            filename_or_obj,
            mode=mode,
            format=format,
            group=group,
            clobber=clobber,
            diskless=diskless,
            persist=persist,
            auto_complex=auto_complex,
            lock=lock,
            autoclose=autoclose,
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
        format="NETCDF4",
        clobber=True,
        diskless=False,
        persist=False,
        auto_complex=None,
        lock=None,
        autoclose=False,
        **kwargs,
    ) -> DataTree:
        groups_dict = self.open_groups_as_dict(
            filename_or_obj,
            mask_and_scale=mask_and_scale,
            decode_times=decode_times,
            concat_characters=concat_characters,
            decode_coords=decode_coords,
            drop_variables=drop_variables,
            use_cftime=use_cftime,
            decode_timedelta=decode_timedelta,
            group=group,
            format=format,
            clobber=clobber,
            diskless=diskless,
            persist=persist,
            auto_complex=auto_complex,
            lock=lock,
            autoclose=autoclose,
            **kwargs,
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
        format="NETCDF4",
        clobber=True,
        diskless=False,
        persist=False,
        auto_complex=None,
        lock=None,
        autoclose=False,
        **kwargs,
    ) -> dict[str, Dataset]:
        from xarray.backends.common import _iter_nc_groups
        from xarray.core.treenode import NodePath

        filename_or_obj = _normalize_path(filename_or_obj)
        store = NetCDF4DataStore.open(
            filename_or_obj,
            group=group,
            format=format,
            clobber=clobber,
            diskless=diskless,
            persist=persist,
            auto_complex=auto_complex,
            lock=lock,
            autoclose=autoclose,
        )

        # Check for a group and make it a parent if it exists
        if group:
            parent = NodePath("/") / NodePath(group)
        else:
            parent = NodePath("/")

        manager = store._manager
        groups_dict = {}
        for path_group in _iter_nc_groups(store.ds, parent=parent):
            group_store = NetCDF4DataStore(manager, group=path_group, **kwargs)
            store_entrypoint = StoreBackendEntrypoint()
            with close_on_error(group_store):
                group_ds = store_entrypoint.open_dataset(
                    group_store,
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


BACKEND_ENTRYPOINTS["netcdf4"] = ("netCDF4", NetCDF4BackendEntrypoint)
