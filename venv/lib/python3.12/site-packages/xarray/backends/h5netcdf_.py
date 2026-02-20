from __future__ import annotations

import functools
import io
import os
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any, Self

import numpy as np
from packaging.version import Version

from xarray.backends.common import (
    BACKEND_ENTRYPOINTS,
    BackendEntrypoint,
    BytesIOProxy,
    T_PathFileOrDataStore,
    WritableCFDataStore,
    _normalize_path,
    _open_remote_file,
    collect_ancestor_dimensions,
    datatree_from_dict_with_io_cleanup,
    find_root_and_group,
)
from xarray.backends.file_manager import (
    CachingFileManager,
    DummyFileManager,
    FileManager,
    PickleableFileManager,
)
from xarray.backends.locks import HDF5_LOCK, combine_locks, ensure_lock, get_write_lock
from xarray.backends.netcdf3 import encode_nc3_attr_value, encode_nc3_variable
from xarray.backends.netCDF4_ import (
    BaseNetCDF4Array,
    _build_and_get_enum,
    _encode_nc4_variable,
    _ensure_no_forward_slash_in_name,
    _extract_nc4_variable_encoding,
    _get_datatype,
    _nc4_require_group,
)
from xarray.backends.store import StoreBackendEntrypoint
from xarray.core import indexing
from xarray.core.utils import (
    FrozenDict,
    emit_user_level_warning,
    is_remote_uri,
    read_magic_number_from_file,
    try_read_magic_number_from_file_or_path,
)
from xarray.core.variable import Variable

if TYPE_CHECKING:
    import h5netcdf

    from xarray.backends.common import AbstractDataStore
    from xarray.core.dataset import Dataset
    from xarray.core.datatree import DataTree
    from xarray.core.types import ReadBuffer


class H5NetCDFArrayWrapper(BaseNetCDF4Array):
    def get_array(self, needs_lock=True):
        ds = self.datastore._acquire(needs_lock)
        return ds.variables[self.variable_name]

    def __getitem__(self, key):
        return indexing.explicit_indexing_adapter(
            key, self.shape, indexing.IndexingSupport.OUTER_1VECTOR, self._getitem
        )

    def _getitem(self, key):
        with self.datastore.lock:
            array = self.get_array(needs_lock=False)
            return array[key]


def _read_attributes(h5netcdf_var):
    # GH451
    # to ensure conventions decoding works properly on Python 3, decode all
    # bytes attributes to strings
    attrs = {}
    for k, v in h5netcdf_var.attrs.items():
        if k not in ["_FillValue", "missing_value"] and isinstance(v, bytes):
            try:
                v = v.decode("utf-8")
            except UnicodeDecodeError:
                emit_user_level_warning(
                    f"'utf-8' codec can't decode bytes for attribute "
                    f"{k!r} of h5netcdf object {h5netcdf_var.name!r}, "
                    f"returning bytes undecoded.",
                    UnicodeWarning,
                )
        attrs[k] = v
    return attrs


_extract_h5nc_encoding = functools.partial(
    _extract_nc4_variable_encoding,
    lsd_okay=False,
    h5py_okay=True,
    backend="h5netcdf",
    unlimited_dims=None,
)


def _h5netcdf_create_group(dataset, name):
    return dataset.create_group(name)


class H5NetCDFStore(WritableCFDataStore):
    """Store for reading and writing data via h5netcdf"""

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
        self,
        manager: FileManager | h5netcdf.File | h5netcdf.Group,
        group=None,
        mode=None,
        format="NETCDF4",
        lock=HDF5_LOCK,
        autoclose=False,
    ):
        import h5netcdf

        if isinstance(manager, h5netcdf.File | h5netcdf.Group):
            if group is None:
                root, group = find_root_and_group(manager)
            else:
                if type(manager) is not h5netcdf.File:
                    raise ValueError(
                        "must supply a h5netcdf.File if the group argument is provided"
                    )
                root = manager
            manager = DummyFileManager(root)

        self._manager = manager
        self._group = group
        self._mode = mode
        self.format = format or "NETCDF4"
        # todo: utilizing find_root_and_group seems a bit clunky
        #  making filename available on h5netcdf.Group seems better
        self._filename = find_root_and_group(self.ds)[0].filename
        self.is_remote = is_remote_uri(self._filename)
        self.lock = ensure_lock(lock)
        self.autoclose = autoclose

    def get_child_store(self, group: str) -> Self:
        if self.format == "NETCDF4_CLASSIC":
            raise ValueError("Cannot create sub-groups in `NETCDF4_CLASSIC` format.")

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
        lock=None,
        autoclose=False,
        invalid_netcdf=None,
        phony_dims=None,
        decode_vlen_strings=True,
        driver=None,
        driver_kwds=None,
        storage_options: dict[str, Any] | None = None,
    ):
        import h5netcdf

        if isinstance(filename, str) and is_remote_uri(filename) and driver is None:
            mode_ = "rb" if mode == "r" else mode
            filename = _open_remote_file(
                filename, mode=mode_, storage_options=storage_options
            )

        if isinstance(filename, BytesIOProxy):
            source = filename
            filename = io.BytesIO()
            source.getvalue = filename.getbuffer

        if isinstance(filename, io.IOBase) and mode == "r":
            magic_number = read_magic_number_from_file(filename)
            if not magic_number.startswith(b"\211HDF\r\n\032\n"):
                raise ValueError(
                    f"{magic_number!r} is not the signature of a valid netCDF4 file"
                )

        if format is None:
            format = "NETCDF4"

        if format not in ["NETCDF4", "NETCDF4_CLASSIC"]:
            raise ValueError(f"invalid format for h5netcdf backend: {format}")

        kwargs = {
            "invalid_netcdf": invalid_netcdf,
            "decode_vlen_strings": decode_vlen_strings,
            "driver": driver,
        }
        if driver_kwds is not None:
            kwargs.update(driver_kwds)
        if phony_dims is not None:
            kwargs["phony_dims"] = phony_dims
        if Version(h5netcdf.__version__) > Version("1.6.4"):
            kwargs["format"] = format
        elif format == "NETCDF4_CLASSIC":
            raise ValueError(
                "h5netcdf >= 1.7.0 is required to save output in NETCDF4_CLASSIC format."
            )

        if lock is None:
            if mode == "r":
                lock = HDF5_LOCK
            else:
                lock = combine_locks([HDF5_LOCK, get_write_lock(filename)])

        manager_cls = (
            CachingFileManager
            if isinstance(filename, str) and not is_remote_uri(filename)
            else PickleableFileManager
        )
        manager = manager_cls(h5netcdf.File, filename, mode=mode, kwargs=kwargs)

        return cls(
            manager,
            group=group,
            format=format,
            mode=mode,
            lock=lock,
            autoclose=autoclose,
        )

    def _acquire(self, needs_lock=True):
        with self._manager.acquire_context(needs_lock) as root:
            ds = _nc4_require_group(
                root, self._group, self._mode, create_group=_h5netcdf_create_group
            )
        return ds

    @property
    def ds(self):
        return self._acquire()

    def open_store_variable(self, name, var):
        import h5netcdf.core
        import h5py

        dimensions = var.dimensions
        data = indexing.LazilyIndexedArray(H5NetCDFArrayWrapper(name, self))
        attrs = _read_attributes(var)

        # netCDF4 specific encoding
        encoding = {
            "chunksizes": var.chunks,
            "fletcher32": var.fletcher32,
            "shuffle": var.shuffle,
        }
        if var.chunks:
            encoding["preferred_chunks"] = dict(
                zip(var.dimensions, var.chunks, strict=True)
            )
        # Convert h5py-style compression options to NetCDF4-Python
        # style, if possible
        if var.compression == "gzip":
            encoding["zlib"] = True
            encoding["complevel"] = var.compression_opts
        elif var.compression is not None:
            encoding["compression"] = var.compression
            encoding["compression_opts"] = var.compression_opts

        # save source so __repr__ can detect if it's local or not
        encoding["source"] = self._filename
        encoding["original_shape"] = data.shape

        vlen_dtype = h5py.check_dtype(vlen=var.dtype)
        if vlen_dtype is str:
            encoding["dtype"] = str
        elif vlen_dtype is not None:  # pragma: no cover
            # xarray doesn't support writing arbitrary vlen dtypes yet.
            pass
        # just check if datatype is available and create dtype
        # this check can be removed if h5netcdf >= 1.4.0 for any environment
        elif (datatype := getattr(var, "datatype", None)) and isinstance(
            datatype, h5netcdf.core.EnumType
        ):
            encoding["dtype"] = np.dtype(
                data.dtype,
                metadata={
                    "enum": datatype.enum_dict,
                    "enum_name": datatype.name,
                },
            )
        else:
            encoding["dtype"] = var.dtype

        return Variable(dimensions, data, attrs, encoding)

    def get_variables(self):
        return FrozenDict(
            (k, self.open_store_variable(k, v)) for k, v in self.ds.variables.items()
        )

    def get_attrs(self):
        return FrozenDict(_read_attributes(self.ds))

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
        if is_unlimited:
            self.ds.dimensions[name] = None
            self.ds.resize_dimension(name, length)
        else:
            self.ds.dimensions[name] = length

    def set_attribute(self, key, value):
        if self.format == "NETCDF4_CLASSIC":
            value = encode_nc3_attr_value(value)
        self.ds.attrs[key] = value

    def encode_variable(self, variable, name=None):
        if self.format == "NETCDF4_CLASSIC":
            return encode_nc3_variable(variable, name=name)
        else:
            return _encode_nc4_variable(variable, name=name)

    def prepare_variable(
        self, name, variable, check_encoding=False, unlimited_dims=None
    ):
        import h5py

        _ensure_no_forward_slash_in_name(name)
        attrs = variable.attrs.copy()
        dtype = _get_datatype(
            variable, nc_format=self.format, raise_on_invalid_encoding=check_encoding
        )

        fillvalue = attrs.pop("_FillValue", None)

        if dtype is str:
            dtype = h5py.special_dtype(vlen=str)

        # check enum metadata and use h5netcdf.core.EnumType
        if (
            hasattr(self.ds, "enumtypes")
            and (meta := np.dtype(dtype).metadata)
            and (e_name := meta.get("enum_name"))
            and (e_dict := meta.get("enum"))
        ):
            dtype = _build_and_get_enum(self, name, dtype, e_name, e_dict)
        encoding = _extract_h5nc_encoding(variable, raise_on_invalid=check_encoding)
        kwargs = {}

        # Convert from NetCDF4-Python style compression settings to h5py style
        # If both styles are used together, h5py takes precedence
        # If set_encoding=True, raise ValueError in case of mismatch
        if encoding.pop("zlib", False):
            if check_encoding and encoding.get("compression") not in (None, "gzip"):
                raise ValueError("'zlib' and 'compression' encodings mismatch")
            encoding.setdefault("compression", "gzip")

        if (
            check_encoding
            and "complevel" in encoding
            and "compression_opts" in encoding
            and encoding["complevel"] != encoding["compression_opts"]
        ):
            raise ValueError("'complevel' and 'compression_opts' encodings mismatch")
        complevel = encoding.pop("complevel", 0)
        if complevel != 0:
            encoding.setdefault("compression_opts", complevel)

        encoding["chunks"] = encoding.pop("chunksizes", None)

        # Do not apply compression, filters or chunking to scalars.
        if variable.shape:
            for key in [
                "compression",
                "compression_opts",
                "shuffle",
                "chunks",
                "fletcher32",
            ]:
                if key in encoding:
                    kwargs[key] = encoding[key]
        if name not in self.ds:
            nc4_var = self.ds.create_variable(
                name,
                dtype=dtype,
                dimensions=variable.dims,
                fillvalue=fillvalue,
                **kwargs,
            )
        else:
            nc4_var = self.ds[name]

        for k, v in attrs.items():
            if self.format == "NETCDF4_CLASSIC":
                v = encode_nc3_attr_value(v)
            nc4_var.attrs[k] = v

        target = H5NetCDFArrayWrapper(name, self)

        return target, variable.data

    def sync(self):
        self.ds.sync()

    def close(self, **kwargs):
        self._manager.close(**kwargs)


def _check_phony_dims(phony_dims):
    emit_phony_dims_warning = False
    if phony_dims is None:
        emit_phony_dims_warning = True
        phony_dims = "access"
    return emit_phony_dims_warning, phony_dims


def _emit_phony_dims_warning():
    emit_user_level_warning(
        "The 'phony_dims' kwarg now defaults to 'access'. "
        "Previously 'phony_dims=None' would raise an error. "
        "For full netcdf equivalence please use phony_dims='sort'.",
        UserWarning,
    )


def _normalize_filename_or_obj(
    filename_or_obj: T_PathFileOrDataStore,
) -> str | ReadBuffer | AbstractDataStore:
    if isinstance(filename_or_obj, bytes | memoryview):
        return io.BytesIO(filename_or_obj)
    else:
        return _normalize_path(filename_or_obj)


class H5netcdfBackendEntrypoint(BackendEntrypoint):
    """
    Backend for netCDF files based on the h5netcdf package.

    It can open ".nc", ".nc4", ".cdf" files but will only be
    selected as the default if the "netcdf4" engine is not available.

    Additionally it can open valid HDF5 files, see
    https://h5netcdf.org/#invalid-netcdf-files for more info.
    It will not be detected as valid backend for such files, so make
    sure to specify ``engine="h5netcdf"`` in ``open_dataset``.

    For more information about the underlying library, visit:
    https://h5netcdf.org

    See Also
    --------
    backends.H5NetCDFStore
    backends.NetCDF4BackendEntrypoint
    backends.ScipyBackendEntrypoint
    """

    description = (
        "Open netCDF (.nc, .nc4 and .cdf) and most HDF5 files using h5netcdf in Xarray"
    )
    url = "https://docs.xarray.dev/en/stable/generated/xarray.backends.H5netcdfBackendEntrypoint.html"
    supports_groups = True

    def guess_can_open(self, filename_or_obj: T_PathFileOrDataStore) -> bool:
        from xarray.core.utils import is_remote_uri

        filename_or_obj = _normalize_filename_or_obj(filename_or_obj)

        # Try to read magic number for local files only
        is_remote = isinstance(filename_or_obj, str) and is_remote_uri(filename_or_obj)
        if not is_remote:
            magic_number = try_read_magic_number_from_file_or_path(filename_or_obj)
            if magic_number is not None:
                return magic_number.startswith(b"\211HDF\r\n\032\n")

        if isinstance(filename_or_obj, str | os.PathLike):
            _, ext = os.path.splitext(filename_or_obj)
            return ext in {".nc", ".nc4", ".cdf"}

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
        format="NETCDF4",
        group=None,
        lock=None,
        invalid_netcdf=None,
        phony_dims=None,
        decode_vlen_strings=True,
        driver=None,
        driver_kwds=None,
        storage_options: dict[str, Any] | None = None,
    ) -> Dataset:
        # Keep this message for some versions
        # remove and set phony_dims="access" above
        emit_phony_dims_warning, phony_dims = _check_phony_dims(phony_dims)

        filename_or_obj = _normalize_filename_or_obj(filename_or_obj)
        store = H5NetCDFStore.open(
            filename_or_obj,
            format=format,
            group=group,
            lock=lock,
            invalid_netcdf=invalid_netcdf,
            phony_dims=phony_dims,
            decode_vlen_strings=decode_vlen_strings,
            driver=driver,
            driver_kwds=driver_kwds,
            storage_options=storage_options,
        )

        store_entrypoint = StoreBackendEntrypoint()

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

        # only warn if phony_dims exist in file
        # remove together with the above check
        # after some versions
        if store.ds._root._phony_dim_count > 0 and emit_phony_dims_warning:
            _emit_phony_dims_warning()

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
        format="NETCDF4",
        group: str | None = None,
        lock=None,
        invalid_netcdf=None,
        phony_dims=None,
        decode_vlen_strings=True,
        driver=None,
        driver_kwds=None,
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
            format=format,
            group=group,
            lock=lock,
            invalid_netcdf=invalid_netcdf,
            phony_dims=phony_dims,
            decode_vlen_strings=decode_vlen_strings,
            driver=driver,
            driver_kwds=driver_kwds,
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
        format="NETCDF4",
        group: str | None = None,
        lock=None,
        invalid_netcdf=None,
        phony_dims=None,
        decode_vlen_strings=True,
        driver=None,
        driver_kwds=None,
        **kwargs,
    ) -> dict[str, Dataset]:
        from xarray.backends.common import _iter_nc_groups
        from xarray.core.treenode import NodePath
        from xarray.core.utils import close_on_error

        # Keep this message for some versions
        # remove and set phony_dims="access" above
        emit_phony_dims_warning, phony_dims = _check_phony_dims(phony_dims)

        filename_or_obj = _normalize_filename_or_obj(filename_or_obj)
        store = H5NetCDFStore.open(
            filename_or_obj,
            format=format,
            group=group,
            lock=lock,
            invalid_netcdf=invalid_netcdf,
            phony_dims=phony_dims,
            decode_vlen_strings=decode_vlen_strings,
            driver=driver,
            driver_kwds=driver_kwds,
        )

        # Check for a group and make it a parent if it exists
        if group:
            parent = NodePath("/") / NodePath(group)
        else:
            parent = NodePath("/")

        manager = store._manager
        groups_dict = {}
        for path_group in _iter_nc_groups(store.ds, parent=parent):
            group_store = H5NetCDFStore(manager, group=path_group, **kwargs)
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

        # only warn if phony_dims exist in file
        # remove together with the above check
        # after some versions
        if store.ds._root._phony_dim_count > 0 and emit_phony_dims_warning:
            _emit_phony_dims_warning()

        return groups_dict


BACKEND_ENTRYPOINTS["h5netcdf"] = ("h5netcdf", H5netcdfBackendEntrypoint)
