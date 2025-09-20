import inspect
import warnings
import copy
import pickle

import numpy as np

from hypothesis import given, example, settings
from hypothesis.strategies import sampled_from

from pytest import raises

from ..array import ArrayIndex
from ..ndindex import ndindex
from ..booleanarray import BooleanArray
from ..integer import Integer
from ..ellipsis import ellipsis
from ..integerarray import IntegerArray
from ..slice import Slice
from ..tuple import Tuple
from .helpers import ndindices, check_same, assert_equal

@example(None)
@example([1, 2])
@given(ndindices)
def test_eq(idx):
    index = ndindex(idx)
    new = type(index)(*index.args)

    assert_equal(new.raw, index.raw)

    def check_eq(a, b):
        assert a == b
        assert b == a
        assert not (a != b)
        assert not (b != a)

    def check_neq(a, b):
        assert a != b
        assert b != a
        assert not (a == b)
        assert not (b == a)

    check_eq(new, index)
    check_eq(new.raw, index)
    check_eq(new, index.raw)

    check_eq(index.raw, index)
    assert hash(new) == hash(index)
    check_neq(index, 'a')

    try:
        h = hash(idx)
    except TypeError:
        pass
    else:
        assert hash(index) == h

    try:
        h = hash(index.raw)
    except TypeError:
        pass
    else:
        assert hash(index) == h

def test_eq_array_raises():
    index = ndindex([1, 2, 3])
    with raises(TypeError):
        np.equal(index.raw, index)
    with raises(TypeError):
        np.array_equal(index.raw, index)

def test_eq_explicit():
    assert Integer(0) != False
    assert Integer(1) != True
    assert Integer(0) != IntegerArray(0)
    assert IntegerArray([0, 1]) != [False, True]
    assert IntegerArray([0, 1]) == [0, 1]
    assert BooleanArray([False, True]) != [0, 1]
    assert BooleanArray([False, True]) == [False, True]

@example((np.array([1, 2]), 0))
@example([1, 2, 3])
@given(ndindices)
def test_ndindex(idx):
    index = ndindex(idx)
    assert index == idx
    assert ndindex[idx] == index

    def assert_raw_eq(idx, index):
        if isinstance(idx, (list, bool, np.bool_)):
            assert isinstance(index, ArrayIndex)
            assert index.dtype in [np.intp, np.bool_]
            assert_equal(index.raw, np.asarray(idx, dtype=index.dtype))
        elif isinstance(idx, np.integer):
            assert type(index) is Integer
            assert_equal(index.raw, int(idx))
        elif isinstance(idx, tuple):
            assert type(index.raw) is tuple
            assert len(idx) == len(index.raw)
            assert index.args == index.raw
            for i, j in zip(index.args, idx):
                assert_raw_eq(j, i)
        else:
            assert_equal(index.raw, idx)

    assert_raw_eq(idx, index)
    assert ndindex(index.raw) == index

def test_ndindex_invalid():
    a = np.arange(10)
    for idx in [1.0, [1.0], np.array([1.0]), np.array([1], dtype=object),
                np.array([])]:
        check_same(a, idx)

    # Older versions of NumPy gives a deprecation warning for this index. We
    # are not going to allow indices that give deprecation warnings in
    # ndindex.
    with warnings.catch_warnings(record=True) as r:
        # Make sure no warnings are emitted from ndindex()
        warnings.simplefilter("error")
        # Newer numpy versions raise ValueError with this index (although
        # perhaps they shouldn't)
        raises((IndexError, ValueError), lambda: ndindex([1, []]))
    assert not r

def test_ndindex_ellipsis():
    raises(IndexError, lambda: ndindex(ellipsis))

def test_signature():
    sig = inspect.signature(Integer)
    assert sig.parameters.keys() == {'idx'}


@example(([0, 1],))
@example((IntegerArray([], (0, 1)),))
@example(IntegerArray([], (0, 1)))
@example((1, ..., slice(1, 2)))
# eval can sometimes be slower than the default deadline of 200ms for large
# array indices
@settings(deadline=None)
@given(ndindices)
def test_repr_str(idx):
    # The repr form should be re-creatable
    index = ndindex(idx)
    d = {}
    exec("from ndindex import *", d)
    assert eval(repr(index), d) == idx

    # Str may not be re-creatable. Just test that it doesn't give an exception.
    str(index)

# _Tuple does not serialize properly with protocols 0 and 1. Support could
# probably be added if this is necessary.
LOWEST_SUPPORTED_PROTOCOL = 2
protocols = ["copy", "deepcopy"] + list(range(LOWEST_SUPPORTED_PROTOCOL, pickle.HIGHEST_PROTOCOL + 1))

@given(ndindices, sampled_from(protocols))
def test_serialization(idx, protocol):
    index = ndindex(idx)

    def serialize(index):
        if protocol == "copy":
            return copy.copy(index)
        elif protocol == "deepcopy":
            return copy.deepcopy(index)
        else:
            return pickle.loads(pickle.dumps(index, protocol=protocol))

    roundtripped = serialize(index)
    assert type(roundtripped) is type(index)
    assert roundtripped == index
    assert_equal(roundtripped.raw, index.raw)
    assert_equal(roundtripped.args, index.args)

    if isinstance(index, Slice):
        assert index._reduced == roundtripped._reduced == False
        s = index.reduce()
        assert s._reduced == True
        roundtripped_s = serialize(s)
        assert roundtripped_s._reduced == True

    if isinstance(index, Tuple):
        assert all([i._reduced == False for i in index.args if isinstance(i, Slice)])
        t = index.reduce()
        assert all([i._reduced == True for i in t.args if isinstance(i, Slice)])
        roundtripped_t = serialize(t)
        assert all([i._reduced == True for i in roundtripped_t.args if isinstance(i, Slice)])
