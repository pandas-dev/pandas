from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.common import is_object_dtype, is_string_dtype, \
    is_datetime64_dtype, is_integer_dtype, is_float_dtype, _NS_DTYPE
from pandas.core.arrays.base import ExtensionArray
from pandas.core.arrays.datetimelike import DatelikeOps, DatetimeLikeArrayMixin
from pandas.core.dtypes.generic import ABCSeries, ABCIndexClass
from pandas.core.dtypes.dtypes import register_extension_dtype, PandasExtensionDtype
from datetime import datetime
from typing import Type
from pandas._libs import lib
from pandas._libs.tslibs import Timestamp, NaT

import numpy as np

D_DATETIME_DTYPE = "datetime64[D]"

@register_extension_dtype
class Date64Dtype(PandasExtensionDtype):
    """
    An ExtensionDtype to hold a single date.

    The attributes name & type are set when subclasses are created.
    """

    _date_aliases = {"date", "date64"}
    _unit = "D"
    _numpy_dtype = np.datetime64


    @property
    def name(self) -> str:
        """
        The alias for DateDtype is ``'string'``.
        """
        return "date"

    @property
    def type(self):
        return Timestamp

    @property
    def na_value(self):
        return NaT

    def __repr__(self):
        return type(self)

    @property
    def kind(self):
        return self.type.kind

    @property
    def itemsize(self):
        """ Return the number of bytes in this dtype """
        return self.numpy_dtype.itemsize

    @classmethod
    def construct_from_string(cls, string: str):
        if string in cls._date_aliases:
            return cls()
        return super().construct_from_string(string)


    @classmethod
    def construct_array_type(cls):
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return Date64Array

    # TODO make from arrow


class Date64Array(DatetimeLikeArrayMixin, DatelikeOps):
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

        self._data = values.astype("datetime64[D]")
        print("backend numpy values:", self._data)
        print("self representation", self)
        print("done initing")

    @classmethod
    def _simple_new(cls, values, freq="D", dtype=D_DATETIME_DTYPE):
        if values.dtype == "i8":
            values = values.view(D_DATETIME_DTYPE)

        freq = "D"
        dtype = D_DATETIME_DTYPE

        result = object.__new__(cls)
        result._data = values
        result._freq = freq
        result._dtype = dtype
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
        Date64Array
        """
        return Date64Array(scalars, copy=copy)

    @property
    def dtype(self) -> ExtensionDtype:
        return Date64Dtype()

    @property
    def freq(self):
        return "D"

    @property
    def _box_func(self):
        def test(x):
            return Timestamp(x, freq="D", tz="utc")

        return test


    # def astype(self, dtype, copy=True):
    #     print(dtype)
    #     print("test")
    #     exit()

    def __len__(self):
        return len(self._data)

    # TODO Add month name

    # TODO Add day name