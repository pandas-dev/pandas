import itertools
import json
import random
import string
import sys

import numpy as np
import pytest

from pandas.core.extensions import ExtensionArray, ExtensionDtype

from .base import BaseArrayTests


class JSONType(ExtensionDtype):
    name = 'json'
    base = None
    kind = 'O'


class JSONArray(ExtensionArray):
    dtype = JSONType()
    fill_value = []
    can_hold_na = True

    def __init__(self, data):
        if isinstance(data, str):
            data = json.loads(data)
        elif isinstance(data, type(self)):
            data = data.data
        assert isinstance(data, list), "'data' must be a list of records."
        self.data = data

    def __getitem__(self, item):
        if isinstance(item, slice):
            result = self.data[item]
        else:
            result = [self.data[i] for i in item]
        return type(self)(result)

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    @property
    def nbytes(self):
        return sum(sys.getsizeof(x) for x in self)

    def isna(self):
        return np.array(x == [] for x in self)

    def take(self, indexer, allow_fill=True, fill_value=None):
        return type(self)(self[indexer])

    take_nd = take

    def formatting_values(self):
        return np.array(self.data)

    def get_values(self):
        return np.array(self.data)

    def slice(self, slicer):
        return self[slicer]

    @classmethod
    def concat_same_type(cls, to_concat):
        return cls(list(itertools.chain(to_concat)))

    def copy(self, deep=False):
        data = self.data
        if deep:
            data = self.data.copy()
        return type(self)(data)


@pytest.fixture
def test_data():
    choices = list(string.ascii_letters) + list(range(100))
    data = [dict([random.choices(choices, k=2)])
            for _ in range(100)]
    return JSONArray(data)


class TestJSONArray(BaseArrayTests):
    pass
