"""support pyarrow compatibility across versions"""

from __future__ import annotations

import sys

import numpy as np

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


def _safe_fill_null(arr: pa.Array, fill_value: object) -> pa.Array:
    """
    Safe wrapper for pyarrow.compute.fill_null with fallback for Windows + pyarrow 21.

    pyarrow 21.0.0 on Windows has a bug in fill_null that incorrectly fills null values.
    This function uses a fallback implementation for that specific case, otherwise uses
    the standard pyarrow.compute.fill_null.

    Parameters
    ----------
    arr : pyarrow.Array
        Input array with potential null values.
    fill_value : object
        Value to fill nulls with.

    Returns
    -------
    pyarrow.Array
        Array with nulls filled with fill_value.
    """
    import pyarrow as pa
    import pyarrow.compute as pc

    # Convert ChunkedArray to Array for consistent handling
    if isinstance(arr, pa.ChunkedArray):
        arr = arr.combine_chunks()

    is_windows = sys.platform in ["win32", "cygwin"]
    use_fallback = (
        HAS_PYARROW and is_windows and not pa_version_under21p0 and pa_version_under22p0
    )

    if use_fallback:
        null_mask = pc.is_null(arr).to_numpy(zero_copy_only=False)
        np_arr = arr.to_numpy(zero_copy_only=False, writable=True)
        np.putmask(np_arr, null_mask, fill_value)
        return pa.array(np_arr, type=arr.type)
    else:
        return pc.fill_null(arr, fill_value)
