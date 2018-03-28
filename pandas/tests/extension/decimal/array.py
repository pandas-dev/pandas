import decimal
import numbers
import random
import sys

import numpy as np

import pandas as pd
from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.common import _ensure_platform_int


class DecimalDtype(ExtensionDtype):
    type = decimal.Decimal
    name = 'decimal'

    @classmethod
    def construct_from_string(cls, string):
        if string == cls.name:
            return cls()
        else:
            raise TypeError("Cannot construct a '{}' from "
                            "'{}'".format(cls, string))


class DecimalArray(ExtensionArray):
    dtype = DecimalDtype()

    def __init__(self, values):
        values = np.asarray(values, dtype=object)

        self.values = values

    @classmethod
    def _constructor_from_sequence(cls, scalars):
        return cls(scalars)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls(values)

    def __getitem__(self, item):
        if isinstance(item, numbers.Integral):
            return self.values[item]
        else:
            return type(self)(self.values[item])

    def copy(self, deep=False):
        if deep:
            return type(self)(self.values.copy())
        return type(self)(self)

    def __setitem__(self, key, value):
        if pd.api.types.is_list_like(value):
            value = [decimal.Decimal(v) for v in value]
        else:
            value = decimal.Decimal(value)
        self.values[key] = value

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return repr(self.values)

    @property
    def nbytes(self):
        n = len(self)
        if n:
            return n * sys.getsizeof(self[0])
        return 0

    def isna(self):
        return np.array([x.is_nan() for x in self.values])

    def take(self, indexer, allow_fill=True, fill_value=None):
        indexer = np.asarray(indexer)
        mask = indexer == -1

        indexer = _ensure_platform_int(indexer)
        out = self.values.take(indexer)
        out[mask] = self._na_value

        return type(self)(out)

    @property
    def _na_value(self):
        return decimal.Decimal('NaN')

    @classmethod
    def _concat_same_type(cls, to_concat):
        return cls(np.concatenate([x.values for x in to_concat]))


def make_data():
    return [decimal.Decimal(random.random()) for _ in range(100)]
