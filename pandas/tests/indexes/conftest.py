from operator import attrgetter

import numpy as np
import pytest

import pandas as pd
from pandas.core.indexes.api import Index, MultiIndex, PeriodIndex, RangeIndex
import pandas.util.testing as tm

indices_list = [tm.makeUnicodeIndex(100),
                tm.makeStringIndex(100),
                tm.makeDateIndex(100),
                tm.makePeriodIndex(100),
                tm.makeTimedeltaIndex(100),
                tm.makeIntIndex(100),
                tm.makeUIntIndex(100),
                tm.makeRangeIndex(100),
                tm.makeFloatIndex(100),
                Index([True, False]),
                tm.makeCategoricalIndex(100),
                tm.makeIntervalIndex(100),
                Index([]),
                MultiIndex.from_tuples(zip(
                    ['foo', 'bar', 'baz'], [1, 2, 3])),
                Index([0, 0, 1, 1, 2, 2])]


@pytest.fixture(params=indices_list, ids=lambda x: type(x).__name__)
def indices(request):
    return request.param


@pytest.fixture(params=[1, np.array(1, dtype=np.int64)])
def one(request):
    # zero-dim integer array behaves like an integer
    return request.param


zeros = [box([0] * 5, dtype=dtype)
         for box in [pd.Index, np.array]
         for dtype in [np.int64, np.uint64, np.float64]]
zeros.extend([np.array(0, dtype=dtype)
              for dtype in [np.int64, np.uint64, np.float64]])
zeros.extend([0, 0.0])


@pytest.fixture(params=zeros)
def zero(request):
    # For testing division by (or of) zero for Index with length 5, this
    # gives several scalar-zeros and length-5 vector-zeros
    return request.param


def _get_subclasses(cls):
    for subclass in cls.__subclasses__():
        yield from _get_subclasses(subclass)
        yield subclass


all_indexes_inc_abc = [Index] + list(set(_get_subclasses(Index)))
all_indexes_inc_abc_sorted = sorted(all_indexes_inc_abc,
                                    key=attrgetter('__name__'))
all_indexes = [index for index in all_indexes_inc_abc_sorted
               if getattr(pd, index.__name__, None) is not None]


@pytest.fixture(params=all_indexes)
def all_index_types(request):
    """
    A Fixture for all indexes types. Index and subclasses in pandas namespace.
    """
    return request.param


@pytest.fixture
def all_index_empty(all_index_types):
    """
    A Fixture for empty instances of all indexes types in pandas namespace.
    """
    cls = all_index_types
    if issubclass(cls, RangeIndex):
        return cls(0, name='foo')
    elif issubclass(cls, MultiIndex):
        return cls.from_arrays([[], []], names=['foo', 'bar'])
    elif issubclass(cls, PeriodIndex):
        return cls([], freq='M', name='foo')
    else:
        return cls([], name='foo')
