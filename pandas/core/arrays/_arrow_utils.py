from distutils.version import LooseVersion

import numpy as np

try:
    import pyarrow

    _PYARROW_INSTALLED = True
except ImportError:
    _PYARROW_INSTALLED = False
    pyarrow = None


if _PYARROW_INSTALLED:
    _pyarrow_version_ge_015 = LooseVersion(pyarrow.__version__) >= LooseVersion("0.15")
else:
    _pyarrow_version_ge_015 = False


def pyarrow_array_to_numpy_and_mask(arr, dtype):
    """
    Convert a primitive pyarrow.Array to a numpy array and boolean mask based
    on the buffers of the Array.

    Parameters
    ----------
    arr : pyarrow.Array
    dtype : numpy.dtype

    Returns
    -------
    (data, mask)
        Tuple of two numpy arrays with the raw data (with specified dtype) and
        a boolean mask (validity mask, so False means missing)
    """
    buflist = arr.buffers()
    data = np.frombuffer(buflist[-1], dtype=dtype)[arr.offset : arr.offset + len(arr)]
    bitmask = buflist[0]
    if bitmask is not None:
        mask = pyarrow.BooleanArray.from_buffers(
            pyarrow.bool_(), len(arr), [None, bitmask]
        )
        mask = np.asarray(mask)
    else:
        mask = np.ones(len(arr), dtype=bool)
    return data, mask
