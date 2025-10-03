"""
Support for pandas ExtensionArray in dask.dataframe.

See :ref:`extensionarrays` for more.
"""

from __future__ import annotations

from dask.dataframe.accessor import (
    register_dataframe_accessor,
    register_index_accessor,
    register_series_accessor,
)
from dask.utils import Dispatch

make_array_nonempty = Dispatch("make_array_nonempty")
make_scalar = Dispatch("make_scalar")


__all__ = [
    "make_array_nonempty",
    "make_scalar",
    "register_dataframe_accessor",
    "register_index_accessor",
    "register_series_accessor",
]
