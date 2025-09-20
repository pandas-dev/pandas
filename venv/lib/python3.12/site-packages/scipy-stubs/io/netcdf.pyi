# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from ._netcdf import netcdf_file as _netcdf_file, netcdf_variable as _netcdf_variable

__all__ = ["netcdf_file", "netcdf_variable"]

@deprecated("will be removed in SciPy v2.0.0")
class netcdf_file(_netcdf_file): ...

@deprecated("will be removed in SciPy v2.0.0")
class netcdf_variable(_netcdf_variable): ...
