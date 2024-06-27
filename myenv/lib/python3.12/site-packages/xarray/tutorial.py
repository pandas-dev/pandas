"""
Useful for:

* users learning xarray
* building tutorials in the documentation.

"""

from __future__ import annotations

import os
import pathlib
from typing import TYPE_CHECKING

import numpy as np

from xarray.backends.api import open_dataset as _open_dataset
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset

if TYPE_CHECKING:
    from xarray.backends.api import T_Engine


_default_cache_dir_name = "xarray_tutorial_data"
base_url = "https://github.com/pydata/xarray-data"
version = "master"


def _construct_cache_dir(path):
    import pooch

    if isinstance(path, os.PathLike):
        path = os.fspath(path)
    elif path is None:
        path = pooch.os_cache(_default_cache_dir_name)

    return path


external_urls = {}  # type: dict
file_formats = {
    "air_temperature": 3,
    "air_temperature_gradient": 4,
    "ASE_ice_velocity": 4,
    "basin_mask": 4,
    "ersstv5": 4,
    "rasm": 3,
    "ROMS_example": 4,
    "tiny": 3,
    "eraint_uvz": 3,
}


def _check_netcdf_engine_installed(name):
    version = file_formats.get(name)
    if version == 3:
        try:
            import scipy  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either scipy or "
                    "netCDF4 to be installed."
                )
    if version == 4:
        try:
            import h5netcdf  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either h5netcdf "
                    "or netCDF4 to be installed."
                )


# idea borrowed from Seaborn
def open_dataset(
    name: str,
    cache: bool = True,
    cache_dir: None | str | os.PathLike = None,
    *,
    engine: T_Engine = None,
    **kws,
) -> Dataset:
    """
    Open a dataset from the online repository (requires internet).

    If a local copy is found then always use that to avoid network traffic.

    Available datasets:

    * ``"air_temperature"``: NCEP reanalysis subset
    * ``"air_temperature_gradient"``: NCEP reanalysis subset with approximate x,y gradients
    * ``"basin_mask"``: Dataset with ocean basins marked using integers
    * ``"ASE_ice_velocity"``: MEaSUREs InSAR-Based Ice Velocity of the Amundsen Sea Embayment, Antarctica, Version 1
    * ``"rasm"``: Output of the Regional Arctic System Model (RASM)
    * ``"ROMS_example"``: Regional Ocean Model System (ROMS) output
    * ``"tiny"``: small synthetic dataset with a 1D data variable
    * ``"era5-2mt-2019-03-uk.grib"``: ERA5 temperature data over the UK
    * ``"eraint_uvz"``: data from ERA-Interim reanalysis, monthly averages of upper level data
    * ``"ersstv5"``: NOAA's Extended Reconstructed Sea Surface Temperature monthly averages

    Parameters
    ----------
    name : str
        Name of the file containing the dataset.
        e.g. 'air_temperature'
    cache_dir : path-like, optional
        The directory in which to search for and write cached data.
    cache : bool, optional
        If True, then cache data locally for use on subsequent calls
    **kws : dict, optional
        Passed to xarray.open_dataset

    See Also
    --------
    tutorial.load_dataset
    open_dataset
    load_dataset
    """
    try:
        import pooch
    except ImportError as e:
        raise ImportError(
            "tutorial.open_dataset depends on pooch to download and manage datasets."
            " To proceed please install pooch."
        ) from e

    logger = pooch.get_logger()
    logger.setLevel("WARNING")

    cache_dir = _construct_cache_dir(cache_dir)
    if name in external_urls:
        url = external_urls[name]
    else:
        path = pathlib.Path(name)
        if not path.suffix:
            # process the name
            default_extension = ".nc"
            if engine is None:
                _check_netcdf_engine_installed(name)
            path = path.with_suffix(default_extension)
        elif path.suffix == ".grib":
            if engine is None:
                engine = "cfgrib"
                try:
                    import cfgrib  # noqa
                except ImportError as e:
                    raise ImportError(
                        "Reading this tutorial dataset requires the cfgrib package."
                    ) from e

        url = f"{base_url}/raw/{version}/{path.name}"

    # retrieve the file
    filepath = pooch.retrieve(url=url, known_hash=None, path=cache_dir)
    ds = _open_dataset(filepath, engine=engine, **kws)
    if not cache:
        ds = ds.load()
        pathlib.Path(filepath).unlink()

    return ds


def load_dataset(*args, **kwargs) -> Dataset:
    """
    Open, load into memory, and close a dataset from the online repository
    (requires internet).

    If a local copy is found then always use that to avoid network traffic.

    Available datasets:

    * ``"air_temperature"``: NCEP reanalysis subset
    * ``"air_temperature_gradient"``: NCEP reanalysis subset with approximate x,y gradients
    * ``"basin_mask"``: Dataset with ocean basins marked using integers
    * ``"rasm"``: Output of the Regional Arctic System Model (RASM)
    * ``"ROMS_example"``: Regional Ocean Model System (ROMS) output
    * ``"tiny"``: small synthetic dataset with a 1D data variable
    * ``"era5-2mt-2019-03-uk.grib"``: ERA5 temperature data over the UK
    * ``"eraint_uvz"``: data from ERA-Interim reanalysis, monthly averages of upper level data
    * ``"ersstv5"``: NOAA's Extended Reconstructed Sea Surface Temperature monthly averages

    Parameters
    ----------
    name : str
        Name of the file containing the dataset.
        e.g. 'air_temperature'
    cache_dir : path-like, optional
        The directory in which to search for and write cached data.
    cache : bool, optional
        If True, then cache data locally for use on subsequent calls
    **kws : dict, optional
        Passed to xarray.open_dataset

    See Also
    --------
    tutorial.open_dataset
    open_dataset
    load_dataset
    """
    with open_dataset(*args, **kwargs) as ds:
        return ds.load()


def scatter_example_dataset(*, seed: None | int = None) -> Dataset:
    """
    Create an example dataset.

    Parameters
    ----------
    seed : int, optional
        Seed for the random number generation.
    """
    rng = np.random.default_rng(seed)
    A = DataArray(
        np.zeros([3, 11, 4, 4]),
        dims=["x", "y", "z", "w"],
        coords={
            "x": np.arange(3),
            "y": np.linspace(0, 1, 11),
            "z": np.arange(4),
            "w": 0.1 * rng.standard_normal(4),
        },
    )
    B = 0.1 * A.x**2 + A.y**2.5 + 0.1 * A.z * A.w
    A = -0.1 * A.x + A.y / (5 + A.z) + A.w
    ds = Dataset({"A": A, "B": B})
    ds["w"] = ["one", "two", "three", "five"]

    ds.x.attrs["units"] = "xunits"
    ds.y.attrs["units"] = "yunits"
    ds.z.attrs["units"] = "zunits"
    ds.w.attrs["units"] = "wunits"

    ds.A.attrs["units"] = "Aunits"
    ds.B.attrs["units"] = "Bunits"

    return ds
