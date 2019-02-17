# -*- coding: utf-8 -*-
import numpy as np
import pytest

from pandas.util import testing as tm


@pytest.fixture(params=['numpy', 'arrow', 'mask'])
def mask_dtype(request):
    """ dtype type """
    if request.param == 'numpy':
        from pandas.core.arrays.mask._numpy import NumpyMaskDtype
        return NumpyMaskDtype
    elif request.param == 'arrow':
        pytest.importorskip('pyarrow', minversion="0.10.0")
        from pandas.core.arrays.mask._pyarrow import ArrowMaskDtype
        return ArrowMaskDtype
    elif request.param == 'mask':
        from pandas.core.arrays.mask import get_mask_array_type
        return type(get_mask_array_type().dtype)


@pytest.fixture
def mask_type(mask_dtype):
    """ array type """
    return mask_dtype.construct_array_type()


@pytest.fixture
def mask(mask_type):
    """ array object """
    return mask_type._from_sequence([1, 0, 1])


def test_construction(mask_type):
    expected = np.array([1, 0, 1], dtype=bool)

    # list
    result = np.array(mask_type._from_sequence([1, 0, 1]))
    tm.assert_numpy_array_equal(result, expected)

    # array
    result = np.array(mask_type._from_sequence(np.array([1, 0, 1])))
    tm.assert_numpy_array_equal(result, expected)

    result = np.array(mask_type._from_sequence(
        np.array([1, 0, 1], dtype=bool)))
    tm.assert_numpy_array_equal(result, expected)


def test_str(mask):

    result = repr(mask)
    expected = '<{}>\n[True, False, True]\nLength: 3, dtype: {}'.format(
        mask.__class__.__name__, mask.dtype)
    assert result == expected


def test_indexing(mask):

    # getitem
    assert mask[0]
    assert not mask[1]
    assert mask[2]

    # slice
    assert (mask[:] == mask).all()
    assert (mask[[0, 1]] == mask._from_sequence([1, 0])).all()

    # setitem
    mask[0] = False
    assert not mask[0]
    mask[[0, 1]] = [1, 1]
    assert mask.all()


def test_ops(mask):

    mask2 = mask._from_sequence([0, 0, 0])
    assert not mask.all()
    assert mask.any()
    assert (mask2 | mask == mask).all()
    assert (mask2 & mask == mask2).any()

    assert (~mask2).all()

    # inplace
    mask2 |= mask
    assert (mask2 == mask._from_sequence([1, 0, 1])).all()

    mask2 &= np.array([0, 0, 0], dtype=bool)
    assert (mask2 == mask._from_sequence([0, 0, 0])).all()


def test_functions(mask):

    assert mask.sum() == 2

    mask2 = mask.copy()
    assert mask2 is not mask
    assert (mask2 == mask).all()

    assert mask.size == len(mask)


def test_dtype(mask_dtype):
    m = mask_dtype()
    assert m == m
    assert m == mask_dtype()
    assert hash(m) is not None
