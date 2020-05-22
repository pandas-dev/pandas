from distutils.version import LooseVersion
import json

import numpy as np
import pyarrow

from pandas.core.arrays.interval import _VALID_CLOSED

_pyarrow_version_ge_015 = LooseVersion(pyarrow.__version__) >= LooseVersion("0.15")


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
    data = np.frombuffer(buflist[1], dtype=dtype)[arr.offset : arr.offset + len(arr)]
    bitmask = buflist[0]
    if bitmask is not None:
        mask = pyarrow.BooleanArray.from_buffers(
            pyarrow.bool_(), len(arr), [None, bitmask]
        )
        mask = np.asarray(mask)
    else:
        mask = np.ones(len(arr), dtype=bool)
    return data, mask


if _pyarrow_version_ge_015:
    # the pyarrow extension types are only available for pyarrow 0.15+

    class ArrowPeriodType(pyarrow.ExtensionType):
        def __init__(self, freq):
            # attributes need to be set first before calling
            # super init (as that calls serialize)
            self._freq = freq
            pyarrow.ExtensionType.__init__(self, pyarrow.int64(), "pandas.period")

        @property
        def freq(self):
            return self._freq

        def __arrow_ext_serialize__(self):
            metadata = {"freq": self.freq}
            return json.dumps(metadata).encode()

        @classmethod
        def __arrow_ext_deserialize__(cls, storage_type, serialized):
            metadata = json.loads(serialized.decode())
            return ArrowPeriodType(metadata["freq"])

        def __eq__(self, other):
            if isinstance(other, pyarrow.BaseExtensionType):
                return type(self) == type(other) and self.freq == other.freq
            else:
                return NotImplemented

        def __hash__(self):
            return hash((str(self), self.freq))

        def to_pandas_dtype(self):
            import pandas as pd

            return pd.PeriodDtype(freq=self.freq)

    # register the type with a dummy instance
    _period_type = ArrowPeriodType("D")
    pyarrow.register_extension_type(_period_type)

    class ArrowIntervalType(pyarrow.ExtensionType):
        def __init__(self, subtype, closed):
            # attributes need to be set first before calling
            # super init (as that calls serialize)
            assert closed in _VALID_CLOSED
            self._closed = closed
            if not isinstance(subtype, pyarrow.DataType):
                subtype = pyarrow.type_for_alias(str(subtype))
            self._subtype = subtype

            storage_type = pyarrow.struct([("left", subtype), ("right", subtype)])
            pyarrow.ExtensionType.__init__(self, storage_type, "pandas.interval")

        @property
        def subtype(self):
            return self._subtype

        @property
        def closed(self):
            return self._closed

        def __arrow_ext_serialize__(self):
            metadata = {"subtype": str(self.subtype), "closed": self.closed}
            return json.dumps(metadata).encode()

        @classmethod
        def __arrow_ext_deserialize__(cls, storage_type, serialized):
            metadata = json.loads(serialized.decode())
            subtype = pyarrow.type_for_alias(metadata["subtype"])
            closed = metadata["closed"]
            return ArrowIntervalType(subtype, closed)

        def __eq__(self, other):
            if isinstance(other, pyarrow.BaseExtensionType):
                return (
                    type(self) == type(other)
                    and self.subtype == other.subtype
                    and self.closed == other.closed
                )
            else:
                return NotImplemented

        def __hash__(self):
            return hash((str(self), str(self.subtype), self.closed))

        def to_pandas_dtype(self):
            import pandas as pd

            return pd.IntervalDtype(self.subtype.to_pandas_dtype())

    # register the type with a dummy instance
    _interval_type = ArrowIntervalType(pyarrow.int64(), "left")
    pyarrow.register_extension_type(_interval_type)
