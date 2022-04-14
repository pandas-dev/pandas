from typing import (
    Optional,
    Tuple,
)

import numpy as np
from packaging import version

from pandas.core.exchange.dataframe_protocol import (
    Buffer,
    DlpackDeviceType,
)

_NUMPY_DLPACK = version.parse("1.22.0")


class PandasBuffer(Buffer):
    """
    Data in the buffer is guaranteed to be contiguous in memory.
    """

    def __init__(self, x: np.ndarray, allow_copy: bool = True) -> None:
        """
        Handle only regular columns (= numpy arrays) for now.
        """
        if not x.strides == (x.dtype.itemsize,):
            # The protocol does not support strided buffers, so a copy is
            # necessary. If that's not allowed, we need to raise an exception.
            if allow_copy:
                x = x.copy()
            else:
                raise RuntimeError(
                    "Exports cannot be zero-copy in the case "
                    "of a non-contiguous buffer"
                )

        # Store the numpy array in which the data resides as a private
        # attribute, so we can use it to retrieve the public attributes
        self._x = x

    @property
    def bufsize(self) -> int:
        """
        Buffer size in bytes.
        """
        return self._x.size * self._x.dtype.itemsize

    @property
    def ptr(self) -> int:
        """
        Pointer to start of the buffer as an integer.
        """
        return self._x.__array_interface__["data"][0]

    if version.parse(np.__version__) >= _NUMPY_DLPACK:

        def __dlpack__(self):
            """
            Represent this structure as DLPack interface.
            """
            return self._x.__dlpack__()

    else:

        def __dlpack__(self):
            """
            Represent this structure as DLPack interface.
            """
            raise NotImplementedError("__dlpack__")

    def __dlpack_device__(self) -> Tuple[DlpackDeviceType, Optional[int]]:
        """
        Device type and device ID for where the data in the buffer resides.
        """
        return (DlpackDeviceType.CPU, None)

    def __repr__(self) -> str:
        return (
            "PandasBuffer("
            + str(
                {
                    "bufsize": self.bufsize,
                    "ptr": self.ptr,
                    "device": self.__dlpack_device__()[0].name,
                }
            )
            + ")"
        )
