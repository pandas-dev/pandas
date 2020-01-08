from typing import TYPE_CHECKING

import numpy as np

from pandas._libs import lib, missing as libmissing

from pandas.core.dtypes.common import is_integer, is_object_dtype, is_string_dtype
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

    def to_numpy(
        self, dtype=None, copy=False, na_value: "Scalar" = lib.no_default,
    ):
        """
        Convert to a NumPy Array.

        By default converts to an object-dtype NumPy array. Specify the `dtype` and
        `na_value` keywords to customize the conversion.

        Parameters
        ----------
        dtype : dtype, default object
            The numpy dtype to convert to.
        copy : bool, default False
            Whether to ensure that the returned value is a not a view on
            the array. Note that ``copy=False`` does not *ensure* that
            ``to_numpy()`` is no-copy. Rather, ``copy=True`` ensure that
            a copy is made, even if not strictly necessary. This is typically
            only possible when no missing values are present and `dtype`
            is the equivalent numpy dtype.
        na_value : scalar, optional
             Scalar missing value indicator to use in numpy array. Defaults
             to the native missing value indicator of this array (pd.NA).

        Returns
        -------
        numpy.ndarray

        Examples
        --------
        An object-dtype is the default result

        >>> a = pd.array([True, False, pd.NA], dtype="boolean")
        >>> a.to_numpy()
        array([True, False, NA], dtype=object)

        When no missing values are present, an equivalent dtype can be used.

        >>> pd.array([True, False], dtype="boolean").to_numpy(dtype="bool")
        array([ True, False])
        >>> pd.array([1, 2], dtype="Int64").to_numpy("int64")
        array([1, 2])

        However, requesting such dtype will raise a ValueError if
        missing values are present and the default missing value :attr:`NA`
        is used.

        >>> a = pd.array([True, False, pd.NA], dtype="boolean")
        >>> a
        <BooleanArray>
        [True, False, NA]
        Length: 3, dtype: boolean

        >>> a.to_numpy(dtype="bool")
        Traceback (most recent call last):
        ...
        ValueError: cannot convert to bool numpy array in presence of missing values

        Specify a valid `na_value` instead

        >>> a.to_numpy(dtype="bool", na_value=False)
        array([ True, False, False])
        """
        if na_value is lib.no_default:
            na_value = libmissing.NA
        if dtype is None:
            dtype = object
        if self._hasna:
            if (
                not (is_object_dtype(dtype) or is_string_dtype(dtype))
                and na_value is libmissing.NA
            ):
                raise ValueError(
                    f"cannot convert to '{dtype}'-dtype NumPy array "
                    "with missing values. Specify an appropriate 'na_value' "
                    "for this dtype."
                )
            # don't pass copy to astype -> always need a copy since we are mutating
            data = self._data.astype(dtype)
            data[self._mask] = na_value
        else:
            data = self._data.astype(dtype, copy=copy)
        return data

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
        return self.dtype.na_value

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
