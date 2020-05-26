from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.arrays.datetimelike import DatelikeOps, DatetimeLikeArrayMixin
from pandas.core.arrays.datetimes import sequence_to_dt64ns
from pandas.core.dtypes.generic import ABCSeries, ABCIndexClass
from pandas.core.dtypes.dtypes import DateDtype
from pandas._libs.tslibs import Timestamp, NaT
from pandas._libs import tslib

import numpy as np

D_DATETIME_DTYPE = "datetime64[D]"


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
        print("initing", values)
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

        if copy:
            values = values.copy()

        self._data = values
        print("backend numpy values:", self._data)
        print("self representation", self)
        print("done initing")

    @classmethod
    def _simple_new(cls, values, **kwargs):
        assert isinstance(values, np.ndarray)
        if values.dtype == "i8":
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
        data, _, _ = sequence_to_dt64ns(scalars, copy=copy)
        return cls._simple_new(data.astype(D_DATETIME_DTYPE))

    # def __repr__(self):

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
        def test(x: np.int64):
            return x

        return test

    @property
    def date(self):
        timestamps = self.asi8
        return tslib.ints_to_pydatetime(timestamps, box="date")


    # def astype(self, dtype, copy=True):
    #     print(dtype)
    #     print("test")
    #     exit()

    def __len__(self):
        return len(self._data)

    # TODO Add month name

    # TODO Add day name