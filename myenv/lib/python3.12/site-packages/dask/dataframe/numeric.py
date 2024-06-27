from __future__ import annotations

import pandas as pd
from pandas.api.types import is_scalar as pd_is_scalar

from dask.array import Array
from dask.dataframe.core import Series
from dask.delayed import delayed
from dask.utils import derived_from

__all__ = ("to_numeric",)


@derived_from(pd, ua_args=["downcast"])
def to_numeric(arg, errors="raise", meta=None):
    """
    Return type depends on input. Delayed if scalar, otherwise same as input.
    For errors, only "raise" and "coerce" are allowed.
    """
    if errors not in ("raise", "coerce"):
        raise ValueError("invalid error value specified")

    is_series = isinstance(arg, Series)
    is_array = isinstance(arg, Array)
    is_scalar = pd_is_scalar(arg)

    if not any([is_series, is_array, is_scalar]):
        raise TypeError(
            "arg must be a list, tuple, dask.array.Array, or dask.dataframe.Series"
        )

    if meta is not None:
        if is_scalar:
            raise KeyError("``meta`` is not allowed when input is a scalar.")
    else:
        if is_series or is_array:
            meta = pd.to_numeric(arg._meta)

    if is_series:
        return arg.map_partitions(
            pd.to_numeric,
            token=arg._name + "-to_numeric",
            meta=meta,
            enforce_metadata=False,
            errors=errors,
        )
    if is_array:
        return arg.map_blocks(
            pd.to_numeric,
            name=arg._name + "-to_numeric",
            meta=meta,
            errors=errors,
        )
    if is_scalar:
        return delayed(pd.to_numeric, pure=True)(arg, errors=errors)
