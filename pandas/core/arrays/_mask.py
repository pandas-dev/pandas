import numpy as np


class NAMask():
    """Generic class which can be used to represent missing data.

    Will use bitarray if available; otherwise will use numpy."""

    def __init__(self, mask):
        """
        Parameters
        ----------
        mask : numpy array
            Mask of missing values.
        """

        self._has_bitarray = False
        try:
            import bitarray
            globals()['bitarray'] = bitarray
            self._has_bitarray = True
            self._data = self._numpy_to_bitarray(mask)
        except (ImportError, ModuleNotFoundError):
            self._data = mask.astype(bool, copy=False)

    def _numpy_to_bitarray(self, arr):
        bit_arr = bitarray()
        bit_arr.pack(arr.astype(bool, copy=False))

    def _bitarray_to_numpy(self, arr):
        return np.fromstring(arr.unpack(), dtype=bool)

    def __getitem__(self, item):
        if self._has_bitarray:
            raise NotImplementedError

        return self._data[item]

    def __setitem__(self, key, value):
        if self._has_bitarray:
            raise NotImplementedError

        self._data[key] = value

    def __array__(self):
        if self._has_bitarray:
            raise NotImplementedError

        return self._data

    def __iter__(self):
        for i in range(len(self._data)):
            yield self._data[i]

    def __invert__(self):
        if self._has_bitarray:
            raise NotImplementedError

        return type(self)(~self._data)

    def __or__(self, other):
        if self._has_bitarray:
            raise NotImplementedError

        return self._data.__or__(other)

    def __ior__(self, other):
        if self._has_bitarray:
            raise NotImplementedError

        return self._data | other

    @property
    def nbytes(self):
        if self._has_bitarray:
            return self._data.buffer_info()[1]

        return self._data.nbytes

    @property
    def size(self):
        if self._has_bitarray:
            raise NotImplementedError

        return self._data.size        

    def astype(self, dtype, copy=False):
        if self._has_bitarray:
            raise NotImplementedError

        return self._data.astype(dtype, copy=copy)

    def any(self):
        return self._data.any()

    def copy(self):
        return type(self)(self._data.copy())

    def sum(self):
        if self._has_bitarray:
            raise NotImplementedError

        return self._data.sum()
