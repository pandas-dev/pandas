"""support pyarrow compatibility across versions"""

from __future__ import annotations

import sys
from typing import Any

from pandas.util.version import Version

PYARROW_MIN_VERSION = "13.0.0"
try:
    import pyarrow as pa

    _palv = Version(Version(pa.__version__).base_version)
    pa_version_under14p0 = _palv < Version("14.0.0")
    pa_version_under14p1 = _palv < Version("14.0.1")
    pa_version_under15p0 = _palv < Version("15.0.0")
    pa_version_under16p0 = _palv < Version("16.0.0")
    pa_version_under17p0 = _palv < Version("17.0.0")
    pa_version_under18p0 = _palv < Version("18.0.0")
    pa_version_under19p0 = _palv < Version("19.0.0")
    pa_version_under20p0 = _palv < Version("20.0.0")
    pa_version_under21p0 = _palv < Version("21.0.0")
    pa_version_under22p0 = _palv < Version("22.0.0")
    HAS_PYARROW = _palv >= Version(PYARROW_MIN_VERSION)
except ImportError:
    pa_version_under14p0 = True
    pa_version_under14p1 = True
    pa_version_under15p0 = True
    pa_version_under16p0 = True
    pa_version_under17p0 = True
    pa_version_under18p0 = True
    pa_version_under19p0 = True
    pa_version_under20p0 = True
    pa_version_under21p0 = True
    pa_version_under22p0 = True
    HAS_PYARROW = False


def _safe_fill_null(
    arr: pa.Array | pa.ChunkedArray, fill_value: Any
) -> pa.Array | pa.ChunkedArray:
    """
    Safe wrapper for pyarrow.compute.fill_null with fallback for Windows + pyarrow 21.

    pyarrow 21.0.0 on Windows has a bug in fill_null that incorrectly fills null values.
    This function uses a fallback implementation for that specific case, otherwise uses
    the standard pyarrow.compute.fill_null.

    Parameters
    ----------
    arr : pyarrow.Array | pyarrow.ChunkedArray
        Input array with potential null values.
    fill_value : Any
        Value to fill nulls with.

    Returns
    -------
    pyarrow.Array | pyarrow.ChunkedArray
        Array with nulls filled with fill_value.
    """
    import pyarrow.compute as pc

    is_windows = sys.platform in ["win32", "cygwin"]
    use_fallback = (
        HAS_PYARROW and is_windows and not pa_version_under21p0 and pa_version_under22p0
    )
    if not use_fallback or isinstance(fill_value, (pa.Array, pa.ChunkedArray)):
        return pc.fill_null(arr, fill_value)

    fill_scalar = pa.scalar(fill_value, type=arr.type)

    if pa.types.is_duration(arr.type):

        def fill_null_duration(arr: pa.Array, fill_scalar: pa.Scalar) -> pa.Array:
            mask = pc.is_null(arr)
            zero_duration = pa.scalar(0, type=arr.type)
            arr_zeroed = pc.if_else(mask, zero_duration, arr)
            return pc.if_else(mask, fill_scalar, arr_zeroed)

        if isinstance(arr, pa.ChunkedArray):
            return pa.chunked_array(
                [fill_null_duration(chunk, fill_scalar) for chunk in arr.chunks]
            )
        return fill_null_duration(arr, fill_scalar)

    if isinstance(arr, pa.ChunkedArray):
        return pa.chunked_array(
            [pc.if_else(pc.is_null(chunk), fill_scalar, chunk) for chunk in arr.chunks]
        )
    return pc.if_else(pc.is_null(arr), fill_scalar, arr)
