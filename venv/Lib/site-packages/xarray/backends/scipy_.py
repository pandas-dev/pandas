from __future__ import annotations

import gzip
import io
import os
from collections.abc import Iterable, Mapping
from typing import IO, TYPE_CHECKING, Any, Literal, TypeVar, overload

import numpy as np

from xarray.backends.common import (
    BACKEND_ENTRYPOINTS,
    BackendArray,
    BackendEntrypoint,
    BytesIOProxy,
    T_PathFileOrDataStore,
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
from xarray.core.types import Lock
from xarray.core.utils import (
    Frozen,
    FrozenDict,
    close_on_error,
    module_available,
    try_read_magic_number_from_file_or_path,
)
from xarray.core.variable import Variable

if TYPE_CHECKING:
    import scipy.io

    from xarray.backends.common import AbstractDataStore
    from xarray.backends.file_manager import FileManager
    from xarray.core.dataset import Dataset
    from xarray.core.types import ReadBuffer

T = TypeVar("T")
K = TypeVar("K")
V = TypeVar("V")

HAS_NUMPY_2_0 = module_available("numpy", minversion="2.0.0.dev0")


@overload
def _decode_string(s: bytes) -> str: ...


@overload
def _decode_string(s: T) -> T: ...


def _decode_string(s: bytes | T) -> str | T:
    if isinstance(s, bytes):
        return s.decode("utf-8", "replace")
    return s


def _decode_attrs(d: Mapping[K, V]) -> dict[K, V]:
    # don't decode _FillValue from bytes -> unicode, because we want to ensure
    # that its type matches the data exactly
    return {k: v if k == "_FillValue" else _decode_string(v) for (k, v) in d.items()}


class ScipyArrayWrapper(BackendArray):
    datastore: ScipyDataStore
    variable_name: str
    shape: tuple[int, ...]
    dtype: np.dtype

    def __init__(self, variable_name: str, datastore: ScipyDataStore) -> None:
        self.datastore = datastore
        self.variable_name = variable_name
        array = self.get_variable().data
        self.shape = array.shape
        self.dtype = np.dtype(array.dtype.kind + str(array.dtype.itemsize))

    def get_variable(self, needs_lock: bool = True) -> scipy.io.netcdf_variable:
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
        copy: bool | None = self.datastore.ds.use_mmap

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


# This is a dirty workaround to allow pickling of the flush_only_netcdf_file class.
# https://stackoverflow.com/questions/72766345/attributeerror-cant-pickle-local-object-in-multiprocessing
# TODO: Remove this after upstreaming the fixes to scipy.
class _PickleWorkaround:
    flush_only_netcdf_file: type[scipy.io.netcdf_file]

    @classmethod
    def add_cls(cls, new_class: type[Any]) -> None:
        setattr(cls, new_class.__name__, new_class)
        new_class.__qualname__ = cls.__qualname__ + "." + new_class.__name__


def _open_scipy_netcdf(
    filename: str | os.PathLike[Any] | IO[bytes],
    mode: Literal["r", "w", "a"],
    mmap: bool | None,
    version: Literal[1, 2],
    flush_only: bool = False,
) -> scipy.io.netcdf_file:
    import scipy.io

    # TODO: Remove this after upstreaming these fixes.
    class flush_only_netcdf_file(scipy.io.netcdf_file):
        # scipy.io.netcdf_file.close() incorrectly closes file objects that
        # were passed in as constructor arguments:
        # https://github.com/scipy/scipy/issues/13905

        # Instead of closing such files, only call flush(), which is
        # equivalent as long as the netcdf_file object is not mmapped.
        # This suffices to keep BytesIO objects open long enough to read
        # their contents from to_netcdf(), but underlying files still get
        # closed when the netcdf_file is garbage collected (via __del__),
        # and will need to be fixed upstream in scipy.
        def close(self):
            if hasattr(self, "fp") and not self.fp.closed:
                self.flush()
                self.fp.seek(0)  # allow file to be read again

        def __del__(self):
            # Remove the __del__ method, which in scipy is aliased to close().
            # These files need to be closed explicitly by xarray.
            pass

    _PickleWorkaround.add_cls(flush_only_netcdf_file)

    netcdf_file = (
        _PickleWorkaround.flush_only_netcdf_file if flush_only else scipy.io.netcdf_file
    )

    # if the string ends with .gz, then gunzip and open as netcdf file
    if isinstance(filename, str) and filename.endswith(".gz"):
        try:
            return netcdf_file(
                gzip.open(filename),  # type: ignore[arg-type]  # not sure if gzip issue or scipy-stubs issue
                mode=mode,
                mmap=mmap,
                version=version,
            )
        except TypeError as e:
            # TODO: gzipped loading only works with NetCDF3 files.
            errmsg = e.args[0]
            if "is not a valid NetCDF 3 file" in errmsg:
                raise ValueError(
                    "gzipped file loading only supports NetCDF 3 files."
                ) from e
            else:
                raise

    try:
        return netcdf_file(filename, mode=mode, mmap=mmap, version=version)
    except TypeError as e:  # netcdf3 message is obscure in this case
        errmsg = e.args[0]
        if "is not a valid NetCDF 3 file" in errmsg:
            msg = """
            If this is a NetCDF4 file, you may need to install the
            netcdf4 library, e.g.,

            $ pip install netcdf4
            """
            errmsg += msg
            raise TypeError(errmsg) from e
        else:
            raise


class ScipyDataStore(WritableCFDataStore):
    """Store for reading and writing data via scipy.io.netcdf_file.

    This store has the advantage of being able to be initialized with a
    StringIO object, allow for serialization without writing to disk.

    It only supports the NetCDF3 file-format.
    """

    lock: Lock
    _manager: FileManager[scipy.io.netcdf_file]

    def __init__(
        self,
        filename_or_obj: T_PathFileOrDataStore,
        mode: Literal["r", "w", "a"] = "r",
        format: str | None = None,
        group: str | None = None,
        mmap: bool | None = None,
        lock: Lock | Literal[False] | None = None,
    ) -> None:
        if group is not None:
            raise ValueError("cannot save to a group with the scipy.io.netcdf backend")

        version: Literal[1, 2]
        if format is None or format == "NETCDF3_64BIT":
            version = 2
        elif format == "NETCDF3_CLASSIC":
            version = 1
        else:
            raise ValueError(f"invalid format for scipy.io.netcdf backend: {format!r}")

        if lock is None and mode != "r" and isinstance(filename_or_obj, str):
            lock = get_write_lock(filename_or_obj)

        self.lock = ensure_lock(lock)

        if isinstance(filename_or_obj, BytesIOProxy):
            source = filename_or_obj
            filename_or_obj = io.BytesIO()
            source.getvalue = filename_or_obj.getbuffer

        manager: FileManager
        if isinstance(filename_or_obj, str):  # path
            manager = CachingFileManager(
                _open_scipy_netcdf,
                filename_or_obj,
                mode=mode,
                lock=lock,
                kwargs=dict(mmap=mmap, version=version),
            )
        elif hasattr(filename_or_obj, "seek"):  # file object
            # Note: checking for .seek matches the check for file objects
            # in scipy.io.netcdf_file
            scipy_dataset = _open_scipy_netcdf(
                filename_or_obj,  # type: ignore[arg-type]  # unsupported cases are caught above
                mode=mode,
                mmap=mmap,
                version=version,
                flush_only=True,
            )
            assert not scipy_dataset.use_mmap  # no mmap for file objects
            manager = DummyFileManager(scipy_dataset)
        else:
            raise ValueError(
                f"cannot open {filename_or_obj=} with scipy.io.netcdf_file"
            )

        self._manager = manager

    @property
    def ds(self) -> scipy.io.netcdf_file:
        return self._manager.acquire()

    def open_store_variable(self, name: str, var: scipy.io.netcdf_variable) -> Variable:
        return Variable(
            var.dimensions,
            indexing.LazilyIndexedArray(ScipyArrayWrapper(name, self)),
            _decode_attrs(var._attributes),  # type: ignore[attr-defined]  # using private attribute
        )

    def get_variables(self) -> Frozen[str, Variable]:
        return FrozenDict(
            (k, self.open_store_variable(k, v)) for k, v in self.ds.variables.items()
        )

    def get_attrs(self) -> Frozen[str, Any]:
        return Frozen(_decode_attrs(self.ds._attributes))  # type: ignore[attr-defined]  # using private attribute

    def get_dimensions(self) -> Frozen[str, int | None]:
        return Frozen(self.ds.dimensions)

    def get_encoding(self) -> dict[Literal["unlimited_dims"], set[str]]:
        return {
            "unlimited_dims": {k for k, v in self.ds.dimensions.items() if v is None}
        }

    def set_dimension(self, name: str, length: int, is_unlimited: bool = False) -> None:
        if name in self.ds.dimensions:
            raise ValueError(
                f"{type(self).__name__} does not support modifying dimensions"
            )
        dim_length = length if not is_unlimited else None
        self.ds.createDimension(name, dim_length)

    def _validate_attr_key(self, key: Any) -> None:
        if not is_valid_nc3_name(key):
            raise ValueError("Not a valid attribute name")

    def set_attribute(self, key: str, value: Any) -> None:
        self._validate_attr_key(key)
        value = encode_nc3_attr_value(value)
        setattr(self.ds, key, value)

    def encode_variable(self, variable: Variable, name: str | None = None) -> Variable:
        variable = encode_nc3_variable(variable, name=name)
        return variable

    def prepare_variable(
        self,
        name: str,
        variable: Variable,
        check_encoding: bool = False,
        unlimited_dims: set[str] | None = None,
    ) -> tuple[ScipyArrayWrapper, Any]:
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
        # don't write the data yet; scipy.io.netcdf does not support incremental
        # writes.
        if name not in self.ds.variables:
            self.ds.createVariable(name, data.dtype, [str(v) for v in variable.dims])
        scipy_var = self.ds.variables[name]
        for k, v in variable.attrs.items():
            self._validate_attr_key(k)
            setattr(scipy_var, k, v)

        target = ScipyArrayWrapper(name, self)

        return target, data

    def sync(self) -> None:
        self.ds.sync()

    def close(self) -> None:
        self._manager.close()


def _normalize_filename_or_obj(
    filename_or_obj: T_PathFileOrDataStore,
) -> str | ReadBuffer | AbstractDataStore:
    if isinstance(filename_or_obj, bytes | memoryview):
        return io.BytesIO(filename_or_obj)
    else:
        return _normalize_path(filename_or_obj)


class ScipyBackendEntrypoint(BackendEntrypoint):
    """
    Backend for netCDF files based on the scipy package.

    It can open ".nc", ".cdf", and "nc..gz" files but will only be
    selected as the default if the "netcdf4" and "h5netcdf" engines are
    not available. It has the advantage that is is a lightweight engine
    that has no system requirements (unlike netcdf4 and h5netcdf).

    Additionally it can open gzip compressed (".gz") files.

    For more information about the underlying library, visit:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.netcdf_file.html

    See Also
    --------
    backends.ScipyDataStore
    backends.NetCDF4BackendEntrypoint
    backends.H5netcdfBackendEntrypoint
    """

    description = "Open netCDF files (.nc, .cdf and .nc.gz) using scipy in Xarray"
    url = "https://docs.xarray.dev/en/stable/generated/xarray.backends.ScipyBackendEntrypoint.html"

    def guess_can_open(
        self,
        filename_or_obj: T_PathFileOrDataStore,
    ) -> bool:
        from xarray.core.utils import is_remote_uri

        filename_or_obj = _normalize_filename_or_obj(filename_or_obj)

        # scipy can only handle local files - check this before trying to read magic number
        if isinstance(filename_or_obj, str) and is_remote_uri(filename_or_obj):
            return False

        magic_number = try_read_magic_number_from_file_or_path(filename_or_obj)
        if magic_number is not None and magic_number.startswith(b"\x1f\x8b"):
            with gzip.open(filename_or_obj) as f:  # type: ignore[arg-type]
                magic_number = try_read_magic_number_from_file_or_path(f)
        if magic_number is not None:
            return magic_number.startswith(b"CDF")

        if isinstance(filename_or_obj, str | os.PathLike):
            from pathlib import Path

            suffix = "".join(Path(filename_or_obj).suffixes)
            return suffix in {".nc", ".cdf", ".nc.gz"}

        return False

    def open_dataset(
        self,
        filename_or_obj: T_PathFileOrDataStore,
        *,
        mask_and_scale: bool = True,
        decode_times: bool = True,
        concat_characters: bool = True,
        decode_coords: bool = True,
        drop_variables: str | Iterable[str] | None = None,
        use_cftime: bool | None = None,
        decode_timedelta: bool | None = None,
        mode: Literal["r", "w", "a"] = "r",
        format: str | None = None,
        group: str | None = None,
        mmap: bool | None = None,
        lock: Lock | Literal[False] | None = None,
    ) -> Dataset:
        filename_or_obj = _normalize_filename_or_obj(filename_or_obj)
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
