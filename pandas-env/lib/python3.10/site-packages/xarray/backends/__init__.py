"""Backend objects for saving and loading data

DataStores provide a uniform interface for saving and loading data in different
formats. They should not be used directly, but rather through Dataset objects.
"""

from xarray.backends.common import AbstractDataStore, BackendArray, BackendEntrypoint
from xarray.backends.file_manager import (
    CachingFileManager,
    DummyFileManager,
    FileManager,
)
from xarray.backends.h5netcdf_ import H5netcdfBackendEntrypoint, H5NetCDFStore
from xarray.backends.memory import InMemoryDataStore
from xarray.backends.netCDF4_ import NetCDF4BackendEntrypoint, NetCDF4DataStore
from xarray.backends.plugins import list_engines, refresh_engines
from xarray.backends.pydap_ import PydapBackendEntrypoint, PydapDataStore
from xarray.backends.scipy_ import ScipyBackendEntrypoint, ScipyDataStore
from xarray.backends.store import StoreBackendEntrypoint
from xarray.backends.zarr import ZarrBackendEntrypoint, ZarrStore

__all__ = [
    "AbstractDataStore",
    "BackendArray",
    "BackendEntrypoint",
    "FileManager",
    "CachingFileManager",
    "DummyFileManager",
    "InMemoryDataStore",
    "NetCDF4DataStore",
    "PydapDataStore",
    "ScipyDataStore",
    "H5NetCDFStore",
    "ZarrStore",
    "H5netcdfBackendEntrypoint",
    "NetCDF4BackendEntrypoint",
    "PydapBackendEntrypoint",
    "ScipyBackendEntrypoint",
    "StoreBackendEntrypoint",
    "ZarrBackendEntrypoint",
    "list_engines",
    "refresh_engines",
]
