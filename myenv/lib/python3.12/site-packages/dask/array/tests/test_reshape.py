from __future__ import annotations

import numpy as np
import pytest

import dask.array as da
from dask.array.reshape import contract_tuple, expand_tuple, reshape_rechunk
from dask.array.utils import assert_eq


@pytest.mark.parametrize(
    "inshape,outshape,prechunks,inchunks,outchunks",
    [
        ((4,), (4,), ((2, 2),), ((2, 2),), ((2, 2),)),
        ((4,), (2, 2), ((2, 2),), ((2, 2),), ((1, 1), (2,))),
        ((4,), (4, 1), ((2, 2),), ((2, 2),), ((2, 2), (1,))),
        ((4,), (1, 4), ((2, 2),), ((2, 2),), ((1,), (2, 2))),
        ((1, 4), (4,), ((1,), (2, 2)), ((1,), (2, 2)), ((2, 2),)),
        ((4, 1), (4,), ((2, 2), (1,)), ((2, 2), (1,)), ((2, 2),)),
        (
            (4, 1, 4),
            (4, 4),
            ((2, 2), (1,), (2, 2)),
            ((2, 2), (1,), (2, 2)),
            ((2, 2), (2, 2)),
        ),
        ((4, 4), (4, 1, 4), ((2, 2), (2, 2)), ((2, 2), (2, 2)), ((2, 2), (1,), (2, 2))),
        ((2, 2), (4,), ((2,), (2,)), ((2,), (2,)), ((4,),)),
        ((2, 2), (4,), ((1, 1), (2,)), ((1, 1), (2,)), ((2, 2),)),
        ((2, 2), (4,), ((2,), (1, 1)), ((1, 1), (2,)), ((2, 2),)),
        (
            (64,),
            (4, 4, 4),
            ((8, 8, 8, 8, 8, 8, 8, 8),),
            ((16, 16, 16, 16),),
            ((1, 1, 1, 1), (4,), (4,)),
        ),
        ((64,), (4, 4, 4), ((32, 32),), ((32, 32),), ((2, 2), (4,), (4,))),
        ((64,), (4, 4, 4), ((16, 48),), ((16, 48),), ((1, 3), (4,), (4,))),
        ((64,), (4, 4, 4), ((20, 44),), ((16, 48),), ((1, 3), (4,), (4,))),
        (
            (64, 4),
            (8, 8, 4),
            ((16, 16, 16, 16), (2, 2)),
            ((16, 16, 16, 16), (2, 2)),
            ((2, 2, 2, 2), (8,), (2, 2)),
        ),
    ],
)
def test_reshape_rechunk(inshape, outshape, prechunks, inchunks, outchunks):
    result_in, result_out = reshape_rechunk(inshape, outshape, prechunks)
    assert result_in == inchunks
    assert result_out == outchunks
    assert np.prod(list(map(len, result_in))) == np.prod(list(map(len, result_out)))


def test_expand_tuple():
    assert expand_tuple((2, 4), 2) == (1, 1, 2, 2)
    assert expand_tuple((2, 4), 3) == (1, 1, 1, 1, 2)
    assert expand_tuple((3, 4), 2) == (1, 2, 2, 2)
    assert expand_tuple((7, 4), 3) == (2, 2, 3, 1, 1, 2)


def test_contract_tuple():
    assert contract_tuple((1, 1, 2, 3, 1), 2) == (2, 2, 2, 2)
    assert contract_tuple((1, 1, 2, 5, 1), 2) == (2, 2, 4, 2)
    assert contract_tuple((2, 4), 2) == (2, 4)
    assert contract_tuple((2, 4), 3) == (6,)


def test_reshape_unknown_sizes():
    a = np.random.random((10, 6, 6))
    A = da.from_array(a, chunks=(5, 2, 3))

    a2 = a.reshape((60, -1))
    A2 = A.reshape((60, -1))

    assert A2.shape == (60, 6)
    assert_eq(A2, a2)

    with pytest.raises(ValueError):
        a.reshape((60, -1, -1))
    with pytest.raises(ValueError):
        A.reshape((60, -1, -1))


@pytest.mark.parametrize(
    "inshape, inchunks, outshape, outchunks",
    [
        # (2, 3, 4) -> (6, 4)
        ((2, 3, 4), ((1, 1), (1, 2), (2, 2)), (6, 4), ((1, 2, 1, 2), (2, 2))),
        # (1, 2, 3, 4) -> (12, 4)
        ((1, 2, 3, 4), ((1,), (1, 1), (1, 2), (2, 2)), (6, 4), ((1, 2, 1, 2), (2, 2))),
        # (2, 2, 3, 4) -> (12, 4)
        (
            (2, 2, 3, 4),
            ((1, 1), (1, 1), (1, 2), (2, 2)),
            (12, 4),
            ((1, 2, 1, 2, 1, 2, 1, 2), (2, 2)),
        ),
        # (2, 2, 3, 4) -> (4, 3, 4)
        (
            (2, 2, 3, 4),
            ((1, 1), (1, 1), (1, 2), (2, 2)),
            (4, 3, 4),
            ((1, 1, 1, 1), (1, 2), (2, 2)),
        ),
        # (2, 2, 3, 4) -> (4, 3, 4)
        ((2, 2, 3, 4), ((1, 1), (2,), (1, 2), (4,)), (4, 3, 4), ((2, 2), (1, 2), (4,))),
        # (2, 3, 4) -> (24,).
        ((2, 3, 4), ((1, 1), (1, 1, 1), (2, 2)), (24,), ((2,) * 12,)),
        # (2, 3, 4) -> (2, 12)
        ((2, 3, 4), ((1, 1), (1, 1, 1), (4,)), (2, 12), ((1, 1), (4,) * 3)),
    ],
)
def test_reshape_all_chunked_no_merge(inshape, inchunks, outshape, outchunks):
    # https://github.com/dask/dask/issues/5544#issuecomment-712280433
    # When the early axes are completely chunked then we are just moving blocks
    # and can avoid any rechunking. The result inchunks are the same as the
    # input chunks.
    base = np.arange(np.prod(inshape)).reshape(inshape)
    a = da.from_array(base, chunks=inchunks)

    # test directly
    inchunks2, outchunks2 = reshape_rechunk(a.shape, outshape, inchunks)
    assert inchunks2 == inchunks
    assert outchunks2 == outchunks

    # and via reshape
    result = a.reshape(outshape)
    assert result.chunks == outchunks
    assert_eq(result, base.reshape(outshape))


@pytest.mark.parametrize(
    "inshape, inchunks, expected_inchunks, outshape, outchunks",
    [
        # (2, 3, 4) -> (24,). This does merge, since the second dim isn't fully chunked!
        ((2, 3, 4), ((1, 1), (1, 2), (2, 2)), ((1, 1), (3,), (4,)), (24,), ((12, 12),)),
    ],
)
def test_reshape_all_not_chunked_merge(
    inshape, inchunks, expected_inchunks, outshape, outchunks
):
    base = np.arange(np.prod(inshape)).reshape(inshape)
    a = da.from_array(base, chunks=inchunks)

    # test directly
    inchunks2, outchunks2 = reshape_rechunk(a.shape, outshape, inchunks)
    assert inchunks2 == expected_inchunks
    assert outchunks2 == outchunks

    # and via reshape
    result = a.reshape(outshape)
    assert result.chunks == outchunks
    assert_eq(result, base.reshape(outshape))


@pytest.mark.parametrize(
    "inshape, inchunks, outshape, outchunks",
    [
        # (2, 3, 4) -> (6, 4)
        ((2, 3, 4), ((2,), (1, 2), (2, 2)), (6, 4), ((1, 2, 1, 2), (2, 2))),
        # (1, 2, 3, 4) -> (12, 4)
        ((1, 2, 3, 4), ((1,), (2,), (1, 2), (2, 2)), (6, 4), ((1, 2, 1, 2), (2, 2))),
        # (2, 2, 3, 4) -> (12, 4) (3 cases)
        (
            (2, 2, 3, 4),
            ((1, 1), (2,), (1, 2), (2, 2)),
            (12, 4),
            ((1, 2, 1, 2, 1, 2, 1, 2), (2, 2)),
        ),
        (
            (2, 2, 3, 4),
            ((2,), (1, 1), (1, 2), (2, 2)),
            (12, 4),
            ((1, 2, 1, 2, 1, 2, 1, 2), (2, 2)),
        ),
        (
            (2, 2, 3, 4),
            ((2,), (2,), (1, 2), (2, 2)),
            (12, 4),
            ((1, 2, 1, 2, 1, 2, 1, 2), (2, 2)),
        ),
        # (2, 2, 3, 4) -> (4, 3, 4)
        # TODO: I'm confused about the behavior in this case.
        # (
        #     (2, 2, 3, 4),
        #     ((2,), (2,), (1, 2), (2, 2)),
        #     (4, 3, 4),
        #     ((1, 1, 1, 1), (1, 2), (2, 2)),
        # ),
        # (2, 2, 3, 4) -> (4, 3, 4)
        ((2, 2, 3, 4), ((2,), (2,), (1, 2), (4,)), (4, 3, 4), ((2, 2), (1, 2), (4,))),
    ],
)
def test_reshape_merge_chunks(inshape, inchunks, outshape, outchunks):
    # https://github.com/dask/dask/issues/5544#issuecomment-712280433
    # When the early axes are completely chunked then we are just moving blocks
    # and can avoid any rechunking. The outchunks will always be ...
    base = np.arange(np.prod(inshape)).reshape(inshape)
    a = da.from_array(base, chunks=inchunks)

    # and via reshape
    result = a.reshape(outshape, merge_chunks=False)
    assert result.chunks == outchunks
    assert_eq(result, base.reshape(outshape))

    assert result.chunks != a.reshape(outshape).chunks
