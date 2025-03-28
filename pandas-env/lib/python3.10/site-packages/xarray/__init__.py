from importlib.metadata import version as _version

from xarray import groupers, testing, tutorial
from xarray.backends.api import (
    load_dataarray,
    load_dataset,
    open_dataarray,
    open_dataset,
    open_mfdataset,
    save_mfdataset,
)
from xarray.backends.zarr import open_zarr
from xarray.coding.cftime_offsets import cftime_range, date_range, date_range_like
from xarray.coding.cftimeindex import CFTimeIndex
from xarray.coding.frequencies import infer_freq
from xarray.conventions import SerializationWarning, decode_cf
from xarray.core.alignment import align, broadcast
from xarray.core.combine import combine_by_coords, combine_nested
from xarray.core.common import ALL_DIMS, full_like, ones_like, zeros_like
from xarray.core.computation import (
    apply_ufunc,
    corr,
    cov,
    cross,
    dot,
    polyval,
    unify_chunks,
    where,
)
from xarray.core.concat import concat
from xarray.core.coordinates import Coordinates
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.extensions import (
    register_dataarray_accessor,
    register_dataset_accessor,
)
from xarray.core.indexes import Index
from xarray.core.indexing import IndexSelResult
from xarray.core.merge import Context, MergeError, merge
from xarray.core.options import get_options, set_options
from xarray.core.parallel import map_blocks
from xarray.core.variable import IndexVariable, Variable, as_variable
from xarray.namedarray.core import NamedArray
from xarray.util.print_versions import show_versions

try:
    __version__ = _version("xarray")
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "9999"

# A hardcoded __all__ variable is necessary to appease
# `mypy --strict` running in projects that import xarray.
__all__ = (
    # Sub-packages
    "groupers",
    "testing",
    "tutorial",
    # Top-level functions
    "align",
    "apply_ufunc",
    "as_variable",
    "broadcast",
    "cftime_range",
    "combine_by_coords",
    "combine_nested",
    "concat",
    "date_range",
    "date_range_like",
    "decode_cf",
    "dot",
    "cov",
    "corr",
    "cross",
    "full_like",
    "get_options",
    "infer_freq",
    "load_dataarray",
    "load_dataset",
    "map_blocks",
    "merge",
    "ones_like",
    "open_dataarray",
    "open_dataset",
    "open_mfdataset",
    "open_zarr",
    "polyval",
    "register_dataarray_accessor",
    "register_dataset_accessor",
    "save_mfdataset",
    "set_options",
    "show_versions",
    "unify_chunks",
    "where",
    "zeros_like",
    # Classes
    "CFTimeIndex",
    "Context",
    "Coordinates",
    "DataArray",
    "Dataset",
    "Index",
    "IndexSelResult",
    "IndexVariable",
    "Variable",
    "NamedArray",
    # Exceptions
    "MergeError",
    "SerializationWarning",
    # Constants
    "__version__",
    "ALL_DIMS",
)
