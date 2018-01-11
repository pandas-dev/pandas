import sys

import itertools
import numpy as np
import pandas as pd
from pandas.core.extensions import ExtensionArray, ExtensionDtype


class MyDtypeType(type):
    pass


class MyDtype(ExtensionDtype):
    name = 'str'
    type = MyDtypeType
    base = None


class MyArray(ExtensionArray):
    can_hold_na = True

    def __init__(self, values):
        if not all(isinstance(v, str) for v in values):
            raise TypeError("Strings Only")
        self.data = list(values)

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.data[item]

    @property
    def dtype(self):
        return MyDtype()

    @property
    def nbytes(self):
        return sum(map(sys.getsizeof, self.data))

    @property
    def ndim(self):
        return 1

    @property
    def shape(self):
        return len(self.data)

    def take(self, indexer, allow_fill=True, fill_value=None):
        return list(np.array(self.data).take(indexer))

    take_nd = take

    def copy(self, deep=False):
        return type(self)(self.data.copy(deep=deep))

    def isna(self):
        return [x is None or x == '' for x in self.data]

    @classmethod
    def concat_same_type(cls, to_concat):
        return cls(list(itertools.chain(to_concat)))

    def formatting_values(self):
        return np.array(self.data)

    def get_values(self):
        return np.array(self.data)

    def to_dense(self):
        return np.array(self.data)



def test_series_constructor():
    result = pd.Series(MyArray(['a', 'b']))
    assert isinstance(result.dtype, MyDtype)
