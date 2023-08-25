import itertools
import pickle

import numpy as np
import pytest

from pandas._libs.arrays import BitMaskArray

import pandas._testing as tm


@pytest.mark.parametrize(
    "array,expected",
    [
        (np.array([False, False]), bytes([0x0])),
        (np.array([True, False]), bytes([0x1])),
        (np.array([False, True]), bytes([0x2])),
        (np.array([True, True]), bytes([0x3])),
        (np.array([True, False] * 8), bytes([0x55, 0x55])),
    ],
)
def test_constructor_ndarray(array, expected):
    bma = BitMaskArray(array)
    assert bma.bytes == expected
    assert not bma.parent
    assert bma.array_shape == array.shape


@pytest.mark.parametrize(
    "parent,expected",
    [
        (BitMaskArray(np.array([False, False])), bytes([0x0])),
        (BitMaskArray(np.array([True, False])), bytes([0x1])),
        (BitMaskArray(np.array([False, True])), bytes([0x2])),
        (BitMaskArray(np.array([True, True])), bytes([0x3])),
        (BitMaskArray(np.array([True, False] * 8)), bytes([0x55, 0x55])),
    ],
)
def test_constructor_bitmap(parent, expected):
    bma = BitMaskArray(parent)
    assert bma.bytes == expected
    assert bma.parent is parent
    assert bma.array_shape == parent.shape


def test_len():
    bma = BitMaskArray(np.array([True, False, False]))
    assert len(bma) == 3


def test_repr_no_parent():
    bma = BitMaskArray(np.array([True, False, False]))
    result = repr(bma)
    assert "parent: None" in result
    assert "shape: (3,)" in result
    assert "data: b'\\x01'" in result


def test_repr_parent():
    parent = BitMaskArray(np.array([False, False, True]))
    bma = BitMaskArray(parent)
    result = repr(bma)
    parent_id = hex(id(parent))
    assert f"parent: <pandas._libs.arrays.BitMaskArray object at {parent_id}" in result
    assert "shape: (3,)" in result
    assert "data: b'\\x04'" in result


@pytest.mark.parametrize(
    "input_data",
    [
        pytest.param([[True]], id="identity_case"),
        pytest.param([[True], [True]], id="base_case"),
        pytest.param([[True], [False]], id="base_case_2"),
        pytest.param([[True], [False] * 7], id="single_byte_boundary_end"),
        pytest.param([[True], [False] * 8], id="multi_byte_non_boundary"),
        pytest.param(
            [[True] * 4, [False] * 4, [True] * 6, [False] * 2],
            id="multi_byte_boundary_end",
        ),
    ],
)
def test_concatenate(input_data):
    masks = [BitMaskArray(np.array(x)) for x in input_data]

    result = BitMaskArray.concatenate(masks, axis=0)
    expected = BitMaskArray(np.array(list(itertools.chain.from_iterable(input_data))))

    assert result.bytes == expected.bytes


def test_concatenate_raises_not_axis0():
    with pytest.raises(NotImplementedError, match="only implemented for axis=0"):
        BitMaskArray.concatenate([], axis=1)


@pytest.mark.parametrize(
    "indexer,expected",
    [
        (0, True),
        (1, False),
    ],
)
def test_getitem_scalar(indexer, expected):
    bma = BitMaskArray(np.array([True, False, True]))
    result = bma[indexer]

    assert result == expected


def test_getitem_null_slice():
    bma = BitMaskArray(np.array([True, False, True]))
    result = bma[:]

    assert result.array_shape == bma.array_shape
    assert not result.parent

    assert result.bytes[0] & 1 == 1
    assert (result.bytes[0] >> 1) & 1 == 0
    assert (result.bytes[0] >> 2) & 1 == 1


@pytest.mark.parametrize(
    "indexer,expected",
    [
        ([0, 1], np.array([True, False])),
        (np.array([2, 1]), np.array([True, False])),
        (slice(1, 2), np.array([False])),
    ],
)
def test_getitem_numpy_fallback(indexer, expected):
    bma = BitMaskArray(np.array([True, False, True]))
    result = bma[indexer]

    tm.assert_numpy_array_equal(result, expected)


def test_setitem_scalar():
    bma = BitMaskArray(np.array([True, False, True]))

    bma[0] = False
    assert not bma[0]

    bma[:] = True
    assert bma[0] and bma[1] and bma[2]

    bma[np.array([False, False, True])] = False
    assert bma[0] and bma[1] and not bma[2]

    bma[[False, True, False]] = False
    assert bma[0] and not bma[1] and not bma[2]


def test_setitem_array():
    bma = BitMaskArray(np.array([True, False, True]))

    bma[:] = [False, True, False]
    assert not bma[0] and bma[1] and not bma[2]

    bma[:] = np.array([True, False, True])
    assert bma[0] and not bma[1] and bma[2]


def test_invert():
    result1 = ~BitMaskArray(np.array([True, False]))
    assert (result1.bytes[0] & 0x1) == 0
    assert ((result1.bytes[0] >> 1) & 0x1) == 1

    result2 = ~BitMaskArray(np.array([False, True]))
    assert (result2.bytes[0] & 0x1) == 1
    assert ((result2.bytes[0] >> 1) & 0x1) == 0


@pytest.mark.parametrize("rhs_as_bitmask", [True, False])
@pytest.mark.parametrize(
    "lhs,rhs,expected",
    [
        ([True], [True], [True]),
        ([True], [False], [False]),
        ([False], [False], [False]),
        ([True] * 10, [True] * 10, [True] * 10),
        ([False] * 10, [True] * 10, [False] * 10),
    ],
)
def test_and(rhs_as_bitmask, lhs, rhs, expected):
    bma1 = BitMaskArray(np.array(lhs))

    if rhs_as_bitmask:
        bma2 = BitMaskArray(np.array(rhs))
    else:
        bma2 = np.array(rhs)

    expected = np.array(expected)
    result = bma1 & bma2
    assert (result == expected).all()


@pytest.mark.parametrize("rhs_as_bitmask", [True, False])
@pytest.mark.parametrize(
    "lhs,rhs,expected",
    [
        ([True], [True], [True]),
        ([True], [False], [True]),
        ([False], [False], [False]),
        ([True] * 10, [True] * 10, [True] * 10),
        ([False] * 10, [True] * 10, [True] * 10),
    ],
)
def test_or(rhs_as_bitmask, lhs, rhs, expected):
    bma1 = BitMaskArray(np.array(lhs))

    if rhs_as_bitmask:
        bma2 = BitMaskArray(np.array(rhs))
    else:
        bma2 = np.array(rhs)

    expected = np.array(expected)
    result = bma1 | bma2
    assert (result == expected).all()


@pytest.mark.parametrize("rhs_as_bitmask", [True, False])
@pytest.mark.parametrize(
    "lhs,rhs,expected",
    [
        ([True], [True], [False]),
        ([True], [False], [True]),
        ([False], [False], [False]),
        ([True] * 10, [True] * 10, [False] * 10),
        ([False] * 10, [True] * 10, [True] * 10),
    ],
)
def test_xor(rhs_as_bitmask, lhs, rhs, expected):
    bma1 = BitMaskArray(np.array(lhs))

    if rhs_as_bitmask:
        bma2 = BitMaskArray(np.array(rhs))
    else:
        bma2 = np.array(rhs)

    expected = np.array(expected)
    result = bma1 ^ bma2
    assert (result == expected).all()


def test_pickle():
    parent = BitMaskArray(np.array([True, False, True]))
    child = BitMaskArray(parent)

    result_child = pickle.loads(pickle.dumps(child))

    assert result_child.shape == child.shape
    assert result_child.bytes == child.bytes

    assert result_child.parent.shape == parent.shape
    assert result_child.parent.bytes == parent.bytes
    assert not result_child.parent.parent


def test_iter():
    bma = BitMaskArray(np.array([True, False, True]))
    itr = iter(bma)

    assert next(itr) is True
    assert next(itr) is False
    assert next(itr) is True

    with pytest.raises(StopIteration, match=""):
        next(itr)


@pytest.mark.parametrize(
    "data,expected",
    [
        (np.array([], dtype=bool), 0),
        (np.array([True, False, True]), 3),
        (np.array([True] * 8), 8),
        (np.array([True] * 8 + [False]), 9),
    ],
)
def test_size(data, expected):
    bma = BitMaskArray(data)
    result = bma.size
    assert result == expected


@pytest.mark.parametrize(
    "data,expected",
    [
        (np.array([], dtype=bool), 0),
        (np.array([True, False, True]), 1),
        (np.array([True] * 8), 1),
        (np.array([True] * 8 + [False]), 2),
    ],
)
def test_nbytes(data, expected):
    bma = BitMaskArray(data)
    result = bma.nbytes
    assert result == expected


@pytest.mark.parametrize(
    "data",
    [
        np.array([True, False]),
        np.array([True, False]).reshape(2, -1),
        np.array([True, False]).reshape(-1, 2),
    ],
)
def test_shape(data):
    bma = BitMaskArray(data)
    assert bma.array_shape == data.shape


@pytest.mark.parametrize(
    "data,expected",
    [
        (np.array([], dtype=bool), False),
        (np.array([True]), True),
        (np.array([False]), False),
        (np.array([True] * 8 + [False]), True),
    ],
)
def test_any(data, expected):
    bma = BitMaskArray(data)
    assert bma.any() == expected


@pytest.mark.parametrize(
    "data,expected",
    [
        (np.array([], dtype=bool), True),
        (np.array([True]), True),
        (np.array([False]), False),
        (np.array([True] * 8 + [False]), False),
    ],
)
def test_all(data, expected):
    bma = BitMaskArray(data)
    assert bma.all() == expected


@pytest.mark.parametrize(
    "data,expected",
    [
        (np.array([], dtype=bool), 0),
        (np.array([True]), 1),
        (np.array([False]), 0),
        (np.array([True] * 8 + [False]), 8),
    ],
)
def test_sum(data, expected):
    bma = BitMaskArray(data)
    assert bma.sum() == expected


def test_take1d():
    bma = BitMaskArray(np.array([True, False, True, False]))

    result1 = bma.take_1d(np.array([0]), axis=0)
    assert (result1.bytes[0] & 0x1) == 1

    result2 = bma.take_1d(np.array([1]), axis=0)
    assert (result2.bytes[0] & 0x1) == 0

    result3 = bma.take_1d(np.array([0, 1]), axis=0)
    assert (result3.bytes[0] & 0x1) == 1
    assert ((result3.bytes[0] >> 1) & 0x1) == 0

    result4 = bma.take_1d(np.array([0, 0]), axis=0)
    assert (result4.bytes[0] & 0x1) == 1
    assert ((result4.bytes[0] >> 1) & 0x1) == 1

    result5 = bma.take_1d(np.array([3, 2, 1, 0]), axis=0)
    assert (result5.bytes[0] & 0x1) == 0
    assert ((result5.bytes[0] >> 1) & 0x1) == 1
    assert ((result5.bytes[0] >> 2) & 0x1) == 0
    assert ((result5.bytes[0] >> 3) & 0x1) == 1


def test_take1d_raises_not_axis0():
    bma = BitMaskArray(np.array([True, False, True]))
    with pytest.raises(NotImplementedError, match="only implemented for axis=0"):
        bma.take_1d(np.array([1]), axis=1)


def test_take_1d_raises_empty_indices():
    bma = BitMaskArray(np.array([True, False, True]))
    with pytest.raises(NotImplementedError, match="does not support empty takes"):
        bma.take_1d(np.array([], dtype="int64"), axis=0)


def test_take_1d_raises_negative_indices():
    bma = BitMaskArray(np.array([True, False, True]))
    with pytest.raises(NotImplementedError, match="does not support negative indexing"):
        bma.take_1d(np.array([-1], dtype="int64"), axis=0)


def test_copy():
    old_bma = BitMaskArray(np.array([True, False, True, False]))
    bma = old_bma.copy()

    assert bma.bytes == old_bma.bytes
    assert bma.shape == old_bma.shape
    assert not bma.parent


@pytest.mark.parametrize(
    "data",
    [
        np.array([], dtype=bool),
        np.array([True] * 100, dtype=bool),
        np.array([[True, False], [True, False], [True, True], [False, False]]),
        np.array([[True, False, True, False], [True, True, False, False]]),
    ],
)
def test_to_numpy(data):
    bma = BitMaskArray(data)

    result = bma.to_numpy()
    tm.assert_numpy_array_equal(result, data)
