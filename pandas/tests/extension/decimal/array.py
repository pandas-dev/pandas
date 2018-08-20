import decimal
import numbers
import sys

import numpy as np

import pandas as pd
from pandas.core.arrays import (ExtensionArray,
                                ExtensionScalarOpsMixin)
from pandas.core.dtypes.base import ExtensionDtype


class DecimalDtype(ExtensionDtype):
    type = decimal.Decimal
    name = 'decimal'
    na_value = decimal.Decimal('NaN')

    def __init__(self, context=None):
        self.context = context or decimal.getcontext()

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.context == other.context
        return super(DecimalDtype, self).__eq__(other)

    def __repr__(self):
        return 'DecimalDtype(context={})'.format(self.context)

    @classmethod
    def construct_array_type(cls):
        """Return the array type associated with this dtype

        Returns
        -------
        type
        """
        return DecimalArray

    @classmethod
    def construct_from_string(cls, string):
        if string == cls.name:
            return cls()
        else:
            raise TypeError("Cannot construct a '{}' from "
                            "'{}'".format(cls, string))


class DecimalArray(ExtensionArray, ExtensionScalarOpsMixin):

    def __init__(self, values, dtype=None, copy=False, context=None):
        for val in values:
            if not isinstance(val, decimal.Decimal):
                raise TypeError("All values must be of type " +
                                str(decimal.Decimal))
        values = np.asarray(values, dtype=object)

        self._data = values
        # Some aliases for common attribute names to ensure pandas supports
        # these
        self._items = self.data = self._data
        # those aliases are currently not working due to assumptions
        # in internal code (GH-20735)
        # self._values = self.values = self.data
        self._dtype = DecimalDtype(context)

    @property
    def dtype(self):
        return self._dtype

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        return cls(scalars)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls(values)

    def __getitem__(self, item):
        if isinstance(item, numbers.Integral):
            return self._data[item]
        else:
            return type(self)(self._data[item])

    def take(self, indexer, allow_fill=False, fill_value=None):
        from pandas.api.extensions import take

        data = self._data
        if allow_fill and fill_value is None:
            fill_value = self.dtype.na_value

        result = take(data, indexer, fill_value=fill_value,
                      allow_fill=allow_fill)
        return self._from_sequence(result)

    def copy(self, deep=False):
        if deep:
            return type(self)(self._data.copy())
        return type(self)(self)

    def astype(self, dtype, copy=True):
        if isinstance(dtype, type(self.dtype)):
            return type(self)(self._data, context=dtype.context)
        return super(DecimalArray, self).astype(dtype, copy)

    def __setitem__(self, key, value):
        if pd.api.types.is_list_like(value):
            value = [decimal.Decimal(v) for v in value]
        else:
            value = decimal.Decimal(value)
        self._data[key] = value

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        return 'DecimalArray({!r})'.format(self._data)

    @property
    def nbytes(self):
        n = len(self)
        if n:
            return n * sys.getsizeof(self[0])
        return 0

    def isna(self):
        return np.array([x.is_nan() for x in self._data], dtype=bool)

    @property
    def _na_value(self):
        return decimal.Decimal('NaN')

    @classmethod
    def _concat_same_type(cls, to_concat):
        return cls(np.concatenate([x._data for x in to_concat]))


DecimalArray._add_arithmetic_ops()
DecimalArray._add_comparison_ops()
