from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.arrays.datetimelike import DatelikeOps, DatetimeLikeArrayMixin
from pandas.core.arrays.datetimes import sequence_to_dt64ns
from pandas.core.dtypes.common import is_integer_dtype, is_datetime64_dtype, is_object_dtype
from pandas.core.dtypes.generic import ABCSeries, ABCIndexClass
from pandas.core.dtypes.dtypes import DateDtype
from pandas.core.construction import array
from pandas._libs.tslibs import Timestamp, NaT
from pandas._libs.tslibs.conversion import NS_DTYPE
from pandas._libs import tslib

import numpy as np

D_DATETIME_DTYPE = "datetime64[D]"
INTEGER_BACKEND = "i8"

def _to_date_values(values, copy=False):
    data, _, _ = sequence_to_dt64ns(values, copy=copy)
    return data.astype(D_DATETIME_DTYPE)

class DateArray(DatetimeLikeArrayMixin, DatelikeOps):
    """
    Pandas ExtensionArray for date (year, month, day only) data.

     .. warning::

       DateArray is currently experimental, and its API may change
       without warning. In particular, :attr:`DateArray.dtype` is
       expected to change to always be an instance of an ``ExtensionDtype``
       subclass.

    Parameters
    ----------
    values : Series, Index, DateArray, ndarray
        The date data.
    freq : str or Offset, optional
    dtype : pd.DateDtype

    copy : bool, default False
        Whether to copy the underlying array of values.

    Attributes
    ----------
    None

    Methods
    -------
    None
    """

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


        if values.dtype == INTEGER_BACKEND:
            values = values.view(D_DATETIME_DTYPE)
        else:
            values = _to_date_values(values, copy)

        if copy:
            values = values.copy()

        self._data = values

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
        if is_integer_dtype(scalars):
            values = scalars._data
        else:
            values = _to_date_values(scalars, copy)
        return cls._simple_new(values)

    @property
    def dtype(self) -> ExtensionDtype:
        return DateDtype()

    @property
    def freq(self):
        return "D"

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
        return self._data.astype("datetime64[ns]").view(INTEGER_BACKEND)

    @property
    def date(self):
        timestamps = self.as_datetime_i8
        return tslib.ints_to_pydatetime(timestamps, box="date")

    def astype(self, dtype, copy=True):
        if is_datetime64_dtype(dtype):
            return array(self._data, dtype="datetime64[ns]")
        if is_object_dtype(dtype):
            return self._box_values(self.as_datetime_i8)
        return super().astype(dtype, copy)

    def _format_native_types(self, na_rep="NaT", date_format=None):
        from pandas.io.formats.format import _get_format_datetime64_from_values

        fmt = _get_format_datetime64_from_values(self, date_format)

        return tslib.format_array_from_datetime(
            self.as_datetime_i8, tz="utc", format=fmt, na_rep=na_rep
        )

    def __len__(self):
        return len(self._data)