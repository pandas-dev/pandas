from __future__ import annotations

import gzip
import io
import os
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

import numpy as np

from xarray.backends.common import (
    BACKEND_ENTRYPOINTS,
    BackendArray,
    BackendEntrypoint,
    WritableCFDataStore,
    _normalize_path,
)
from xarray.backends.file_manager import CachingFileManager, DummyFileManager
from xarray.backends.locks import ensure_lock, get_write_lock
from xarray.backends.netcdf3 import (
    encode_nc3_attr_value,
    encode_nc3_variable,
    is_valid_nc3_name,
)
from xarray.backends.store import StoreBackendEntrypoint
from xarray.core import indexing
from xarray.core.utils import (
    Frozen,
    FrozenDict,
    close_on_error,
    module_available,
    try_read_magic_number_from_file_or_path,
)
from xarray.core.variable import Variable

if TYPE_CHECKING:
    from io import BufferedIOBase

    from xarray.backends.common import AbstractDataStore
    from xarray.core.dataset import Dataset


HAS_NUMPY_2_0 = module_available("numpy", minversion="2.0.0.dev0")


def _decode_string(s):
    if isinstance(s, bytes):
        return s.decode("utf-8", "replace")
    return s


def _decode_attrs(d):
    # don't decode _FillValue from bytes -> unicode, because we want to ensure
    # that its type matches the data exactly
    return {k: v if k == "_FillValue" else _decode_string(v) for (k, v) in d.items()}


class ScipyArrayWrapper(BackendArray):
    def __init__(self, variable_name, datastore):
        self.datastore = datastore
        self.variable_name = variable_name
        array = self.get_variable().data
        self.shape = array.shape
        self.dtype = np.dtype(array.dtype.kind + str(array.dtype.itemsize))

    def get_variable(self, needs_lock=True):
        ds = self.datastore._manager.acquire(needs_lock)
        return ds.variables[self.variable_name]

    def _getitem(self, key):
        with self.datastore.lock:
            data = self.get_variable(needs_lock=False).data
            return data[key]

    def __getitem__(self, key):
        data = indexing.explicit_indexing_adapter(
            key, self.shape, indexing.IndexingSupport.OUTER_1VECTOR, self._getitem
        )
        # Copy data if the source file is mmapped. This makes things consistent
        # with the netCDF4 library by ensuring we can safely read arrays even
        # after closing associated files.
        copy = self.datastore.ds.use_mmap

        # adapt handling of copy-kwarg to numpy 2.0
        # see https://github.com/numpy/numpy/issues/25916
        # and https://github.com/numpy/numpy/pull/25922
        copy = None if HAS_NUMPY_2_0 and copy is False else copy

        return np.array(data, dtype=self.dtype, copy=copy)

    def __setitem__(self, key, value):
        with self.datastore.lock:
            data = self.get_variable(needs_lock=False)
            try:
                data[key] = value
            except TypeError:
                if key is Ellipsis:
                    # workaround for GH: scipy/scipy#6880
                    data[:] = value
                else:
                    raise


def _open_scipy_netcdf(filename, mode, mmap, version):
    import scipy.io

    # if the string ends with .gz, then gunzip and open as netcdf file
    if isinstance(filename, str) and filename.endswith(".gz"):
        try:
            return scipy.io.netcdf_file(
                gzip.open(filename), mode=mode, mmap=mmap, version=version
            )
        except TypeError as e:
            # TODO: gzipped loading only works with NetCDF3 files.
            errmsg = e.args[0]
            if "is not a valid NetCDF 3 file" in errmsg:
                raise ValueError("gzipped file loading only supports NetCDF 3 files.")
            else:
                raise

    if isinstance(filename, bytes) and filename.startswith(b"CDF"):
        # it's a NetCDF3 bytestring
        filename = io.BytesIO(filename)

    try:
        return scipy.io.netcdf_file(filename, mode=mode, mmap=mmap, version=version)
    except TypeError as e:  # netcdf3 message is obscure in this case
        errmsg = e.args[0]
        if "is not a valid NetCDF 3 file" in errmsg:
            msg = """
            If this is a NetCDF4 file, you may need to install the
            netcdf4 library, e.g.,

            $ pip install netcdf4
            """
            errmsg += msg
            raise TypeError(errmsg)
        else:
            raise


class ScipyDataStore(WritableCFDataStore):
    """Store for reading and writing data via scipy.io.netcdf.

    This store has the advantage of being able to be initialized with a
    StringIO object, allow for serialization without writing to disk.

    It only supports the NetCDF3 file-format.
    """

    def __init__(
        self, filename_or_obj, mode="r", format=None, group=None, mmap=None, lock=None
    ):
        if group is not None:
            raise ValueError("cannot save to a group with the scipy.io.netcdf backend")

        if format is None or format == "NETCDF3_64BIT":
            version = 2
        elif format == "NETCDF3_CLASSIC":
            version = 1
        else:
            raise ValueError(f"invalid format for scipy.io.netcdf backend: {format!r}")

        if lock is None and mode != "r" and isinstance(filename_or_obj, str):
            lock = get_write_lock(filename_or_obj)

        self.lock = ensure_lock(lock)

        if isinstance(filename_or_obj, str):
            manager = CachingFileManager(
                _open_scipy_netcdf,
                filename_or_obj,
                mode=mode,
                lock=lock,
                kwargs=dict(mmap=mmap, version=version),
            )
        else:
            scipy_dataset = _open_scipy_netcdf(
                filename_or_obj, mode=mode, mmap=mmap, version=version
            )
            manager = DummyFileManager(scipy_dataset)

        self._manager = manager

    @property
    def ds(self):
        return self._manager.acquire()

    def open_store_variable(self, name, var):
        return Variable(
            var.dimensions,
            ScipyArrayWrapper(name, self),
            _decode_attrs(var._attributes),
        )

    def get_variables(self):
        return FrozenDict(
            (k, self.open_store_variable(k, v)) for k, v in self.ds.variables.items()
        )

    def get_attrs(self):
        return Frozen(_decode_attrs(self.ds._attributes))

    def get_dimensions(self):
        return Frozen(self.ds.dimensions)

    def get_encoding(self):
        return {
            "unlimited_dims": {k for k, v in self.ds.dimensions.items() if v is None}
        }

    def set_dimension(self, name, length, is_unlimited=False):
        if name in self.ds.dimensions:
            raise ValueError(
                f"{type(self).__name__} does not support modifying dimensions"
            )
        dim_length = length if not is_unlimited else None
        self.ds.createDimension(name, dim_length)

    def _validate_attr_key(self, key):
        if not is_valid_nc3_name(key):
            raise ValueError("Not a valid attribute name")

    def set_attribute(self, key, value):
        self._validate_attr_key(key)
        value = encode_nc3_attr_value(value)
        setattr(self.ds, key, value)

    def encode_variable(self, variable):
        variable = encode_nc3_variable(variable)
        return variable

    def prepare_variable(
        self, name, variable, check_encoding=False, unlimited_dims=None
    ):
        if (
            check_encoding
            and variable.encoding
            and variable.encoding != {"_FillValue": None}
        ):
            raise ValueError(
                f"unexpected encoding for scipy backend: {list(variable.encoding)}"
            )

        data = variable.data
        # nb. this still creates a numpy array in all memory, even though we
        # don't write the data yet; scipy.io.netcdf does not not support
        # incremental writes.
        if name not in self.ds.variables:
            self.ds.createVariable(name, data.dtype, variable.dims)
        scipy_var = self.ds.variables[name]
        for k, v in variable.attrs.items():
            self._validate_attr_key(k)
            setattr(scipy_var, k, v)

        target = ScipyArrayWrapper(name, self)

        return target, data

    def sync(self):
        self.ds.sync()

    def close(self):
        self._manager.close()


class ScipyBackendEntrypoint(BackendEntrypoint):
    """
    Backend for netCDF files based on the scipy package.

    It can open ".nc", ".nc4", ".cdf" and ".gz" files but will only be
    selected as the default if the "netcdf4" and "h5netcdf" engines are
    not available. It has the advantage that is is a lightweight engine
    that has no system requirements (unlike netcdf4 and h5netcdf).

    Additionally it can open gizp compressed (".gz") files.

    For more information about the underlying library, visit:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.netcdf_file.html

    See Also
    --------
    backends.ScipyDataStore
    backends.NetCDF4BackendEntrypoint
    backends.H5netcdfBackendEntrypoint
    """

    description = "Open netCDF files (.nc, .nc4, .cdf and .gz) using scipy in Xarray"
    url = "https://docs.xarray.dev/en/stable/generated/xarray.backends.ScipyBackendEntrypoint.html"

    def guess_can_open(
        self,
        filename_or_obj: str | os.PathLike[Any] | BufferedIOBase | AbstractDataStore,
    ) -> bool:
        magic_number = try_read_magic_number_from_file_or_path(filename_or_obj)
        if magic_number is not None and magic_number.startswith(b"\x1f\x8b"):
            with gzip.open(filename_or_obj) as f:  # type: ignore[arg-type]
                magic_number = try_read_magic_number_from_file_or_path(f)
        if magic_number is not None:
            return magic_number.startswith(b"CDF")

        if isinstance(filename_or_obj, (str, os.PathLike)):
            _, ext = os.path.splitext(filename_or_obj)
            return ext in {".nc", ".nc4", ".cdf", ".gz"}

        return False

    def open_dataset(  # type: ignore[override]  # allow LSP violation, not supporting **kwargs
        self,
        filename_or_obj: str | os.PathLike[Any] | BufferedIOBase | AbstractDataStore,
        *,
        mask_and_scale=True,
        decode_times=True,
        concat_characters=True,
        decode_coords=True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime=None,
        decode_timedelta=None,
        mode="r",
        format=None,
        group=None,
        mmap=None,
        lock=None,
    ) -> Dataset:
        filename_or_obj = _normalize_path(filename_or_obj)
        store = ScipyDataStore(
            filename_or_obj, mode=mode, format=format, group=group, mmap=mmap, lock=lock
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


BACKEND_ENTRYPOINTS["scipy"] = ("scipy", ScipyBackendEntrypoint)
