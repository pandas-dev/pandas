from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.arrays.datetimelike import DatelikeOps, DatetimeLikeArrayMixin
from pandas.core.arrays.datetimes import sequence_to_dt64ns
from pandas.core.dtypes.common import (
    is_integer_dtype,
    is_datetime64_dtype,
    is_object_dtype,
    is_string_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.generic import ABCSeries, ABCIndexClass
from pandas.core.dtypes.dtypes import DateDtype
from pandas.core.construction import array
from pandas._libs.tslibs import Timestamp
from pandas._libs.tslibs.conversion import DT64NS_DTYPE
from pandas._libs import tslib, lib
from pandas.core.arrays._mixins import _T

import numpy as np

D_DATETIME_DTYPE = "datetime64[D]"
INTEGER_BACKEND = "i8"
VALID_TYPES = {INTEGER_BACKEND, "datetime64[ns]", D_DATETIME_DTYPE, "object"}


def _to_date_values(values, copy=False):
    data, _, _ = sequence_to_dt64ns(values, copy=copy)
    return data.astype(D_DATETIME_DTYPE)


class DateArray(DatetimeLikeArrayMixin, DatelikeOps):
    """
    Pandas ExtensionArray for date (year, month, day only) data.

    Parameters
    ----------
    values : Series, Index, DateArray, ndarray
        The date data.
    copy : bool, default False
        Whether to copy the underlying array of values.

    Attributes
    ----------
    None

    Methods
    -------
    None
    """

    freq = "D"

    def __init__(self, values, copy=False):
        if isinstance(values, (ABCSeries, ABCIndexClass)):
            values = values._values

        if isinstance(values, type(self)):
            values = values._data

        if not isinstance(values, np.ndarray):
            msg = (
                f"Unexpected type '{type(values).__name__}'. 'values' must be "
                "a DateArray ndarray, or Series or Index containing one of"
                " those."
            )
            raise ValueError(msg)

        if not self._is_compatible_dtype(values.dtype):
            msg = (
                f"The dtype of 'values' is incorrect. Must be one of {VALID_TYPES}."
                f" Got {values.dtype} instead."
            )
            raise ValueError(msg)

        if values.dtype == INTEGER_BACKEND:
            values = values.view(D_DATETIME_DTYPE)
        elif values.dtype != "datetime64[D]":
            values = _to_date_values(values, copy)

        if copy:
            values = values.copy()

        self._data = values

    @staticmethod
    def _is_compatible_dtype(dtype):
        return (
            is_integer_dtype(dtype)
            or is_object_dtype(dtype)
            or is_datetime64_dtype(dtype)
            or dtype == "datetime64[D]"
        )

    @classmethod
    def _simple_new(cls, values, **kwargs):
        assert isinstance(values, np.ndarray)
        if values.dtype == INTEGER_BACKEND:
            values = values.view(D_DATETIME_DTYPE)

        result = object.__new__(cls)
        result._data = values
        return result

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        """
        Construct a new ExtensionArray from a sequence of scalars.

        Parameters
        ----------
        scalars : Sequence
            Each element will be an instance of the scalar type for this
            array, ``cls.dtype.type``.
        dtype : dtype, optional
            Construct for this particular dtype. This should be a Dtype
            compatible with the ExtensionArray.
        copy : bool, default False
            If True, copy the underlying data.

        Returns
        -------
        DateArray
        """
        if (
            isinstance(scalars, np.ndarray)
            and lib.infer_dtype(scalars, skipna=True) == "integer"
        ):
            values = scalars.astype(INTEGER_BACKEND)
        elif is_integer_dtype(scalars):
            values = scalars._data
        else:
            values = _to_date_values(scalars, copy)
        return cls._simple_new(values)

    def _from_backing_data(self: _T, arr: np.ndarray) -> _T:
        return type(self)(arr)

    @property
    def dtype(self) -> ExtensionDtype:
        return DateDtype()

    def __iter__(self):
        for date_data in self._data:
            yield date_data

    @property
    def _box_func(self):
        # TODO Implement Datestamp of a similar form in cython
        return lambda x: Timestamp(x, freq="D", tz="utc")

    @property
    def asi8(self) -> np.ndarray:
        return self._data.view(INTEGER_BACKEND)

    @property
    def as_datetime_i8(self) -> np.ndarray:
        return self._data.astype(DT64NS_DTYPE).view(INTEGER_BACKEND)

    @property
    def date(self):
        timestamps = self.as_datetime_i8
        return tslib.ints_to_pydatetime(timestamps, box="date")

    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)
        if isinstance(dtype, type(self.dtype)):
            if copy:
                return self.copy()
            return self
        if is_datetime64_dtype(dtype):
            return array(self._data, dtype=DT64NS_DTYPE)
        if is_object_dtype(dtype):
            return self._box_values(self.as_datetime_i8)
        if is_string_dtype(dtype):
            return array(self._format_native_types())
        return super().astype(dtype, copy)

    def _format_native_types(self, na_rep="NaT", date_format=None):
        return tslib.format_array_from_datetime(
            self.as_datetime_i8, tz="utc", format="%Y-%m-%d", na_rep=na_rep
        )

    def __len__(self):
        return len(self._data)
