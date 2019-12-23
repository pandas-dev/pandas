from typing import TYPE_CHECKING

import numpy as np

from pandas.core.dtypes.common import is_integer
from pandas.core.dtypes.missing import isna, notna

from pandas.core.algorithms import take
from pandas.core.arrays import ExtensionArray, ExtensionOpsMixin
import pandas.core.common as com
from pandas.core.indexers import check_bool_array_indexer

if TYPE_CHECKING:
    from pandas._typing import Scalar


class BaseMaskedArray(ExtensionArray, ExtensionOpsMixin):
    """
    Base class for masked arrays (which use _data and _mask to store the data).

    numpy based
    """

    _data: np.ndarray
    _mask: np.ndarray

    # The value used to fill '_data' to avoid upcasting
    _internal_fill_value: "Scalar"

    def __getitem__(self, item):
        if is_integer(item):
            if self._mask[item]:
                return self.dtype.na_value
            return self._data[item]

        elif com.is_bool_indexer(item):
            item = check_bool_array_indexer(self, item)

        return type(self)(self._data[item], self._mask[item])

    def __iter__(self):
        for i in range(len(self)):
            if self._mask[i]:
                yield self.dtype.na_value
            else:
                yield self._data[i]

    def __len__(self) -> int:
        return len(self._data)

    __array_priority__ = 1000  # higher than ndarray so ops dispatch to us

    def __array__(self, dtype=None):
        """
        the array interface, return my values
        We return an object array here to preserve our scalar values
        """
        return self.to_numpy(dtype=dtype)

    def __arrow_array__(self, type=None):
        """
        Convert myself into a pyarrow Array.
        """
        import pyarrow as pa

        return pa.array(self._data, mask=self._mask, type=type)

    @property
    def _hasna(self) -> bool:
        # Note: this is expensive right now! The hope is that we can
        # make this faster by having an optional mask, but not have to change
        # source code using it..
        return self._mask.any()

    def isna(self):
        return self._mask

    @property
    def _na_value(self):
        return self._dtype.na_value

    @property
    def nbytes(self):
        return self._data.nbytes + self._mask.nbytes

    @classmethod
    def _concat_same_type(cls, to_concat):
        data = np.concatenate([x._data for x in to_concat])
        mask = np.concatenate([x._mask for x in to_concat])
        return cls(data, mask)

    def take(self, indexer, allow_fill=False, fill_value=None):
        # we always fill with 1 internally
        # to avoid upcasting
        data_fill_value = self._internal_fill_value if isna(fill_value) else fill_value
        result = take(
            self._data, indexer, fill_value=data_fill_value, allow_fill=allow_fill
        )

        mask = take(self._mask, indexer, fill_value=True, allow_fill=allow_fill)

        # if we are filling
        # we only fill where the indexer is null
        # not existing missing values
        # TODO(jreback) what if we have a non-na float as a fill value?
        if allow_fill and notna(fill_value):
            fill_mask = np.asarray(indexer) == -1
            result[fill_mask] = fill_value
            mask = mask ^ fill_mask

        return type(self)(result, mask, copy=False)

    def copy(self):
        data, mask = self._data, self._mask
        data = data.copy()
        mask = mask.copy()
        return type(self)(data, mask, copy=False)
