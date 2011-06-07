import itertools

from numpy import nan
import numpy as np

from pandas.core.index import Index
from pandas.core.common import _ensure_index
import pandas.lib.tseries as tseries

class _Block(object):

    def __init__(self, values, columns):
        if values.ndim == 1:
            values = np.atleast_2d(values).T
        self.values = values
        self.columns = _ensure_index(columns)

    def __len__(self):
        return len(self.values)

    @property
    def dtype(self):
        return self.values.dtype

    def reindex(self, new_columns):
        pass

def _interleave(blocks, columns):
    dtype = _interleaved_dtype(blocks)
    columns = _ensure_index(columns)

    result = np.empty((len(blocks[0]), len(columns)), dtype=dtype)
    result.fill(nan)

    for block in blocks:
        indexer, mask = columns.get_indexer(block.columns)

        if mask.all():
            result[:, indexer] = block.values
        else:
            result[:, indexer[mask]] = block.values[mask]

    return result

def _consolidate(blocks):
    """
    Merge blocks having same dtype
    """
    get_dtype = lambda x: x.dtype
    grouper = itertools.groupby(sorted(blocks, key=get_dtype),
                                lambda x: x.dtype)

    new_blocks = []
    for dtype, group_blocks in grouper:
        group_blocks = list(group_blocks)
        new_values = np.hstack([b.values for b in group_blocks])
        new_columns = np.concatenate([b.columns for b in group_blocks])
        new_blocks.append(_Block(new_values, new_columns))

    return new_blocks

# def _merge_blocks(blocks):
#     return np.hstack([b.values for b in blocks])

def _xs(blocks, i, copy=True):
    if copy:
        return np.concatenate([b[i] for b in blocks])
    else:
        if len(blocks) == 1:
            return blocks[0].values[i]
        else:
            raise Exception('cannot get view with mixed-type data')

def _interleaved_dtype(blocks):
    for block in blocks:
        if not issubclass(block.dtype.type, np.floating):
            return object
    return np.float64

class _MixedTypeData(object):

    def __init__(self, floats, objects, float_cols, object_cols):
        pass

    def as_matrix(self, columns=None):
        if columns is None:
            if self.nblocks == 0:
                return self.blocks[0]

            return self.hstack(self.blocks)
        else:
            pass

class _BlockedData(object):

    def __init__(self, columns, blocks, block_cols):

        self.columns = columns
        self.blocks = blocks
        self.block_cols = block_cols
        self.block_widths = np.asarray([len(cols) for cols in block_cols])

    @property
    def nblocks(self):
        return len(self.blocks)

    def as_matrix(self, columns=None):
        if columns is None:
            if self.nblocks == 0:
                return self.blocks.values[0]
            return _interleave(self.blocks, self.columns)
        else:
            pass

    def xs(self, i, copy=True):
        return np.concatenate([b[i] for b in blocks])

    def get(self, col):
        for block in self.blocks:
            pass

    def set(self, key, value):
        pass

    def _chunk_index(self, col):
        pass

import unittest
class TestBlockOperations(unittest.TestCase):

    def test_interleave(self):
        pass

    def test_consolidate(self):
        pass

    def test_xs(self):
        pass

if __name__ == '__main__':
    floats = np.repeat(np.atleast_2d(np.arange(3.)), 10, axis=0)
    objects = np.empty((10, 2), dtype=object)
    objects[:, 0] = 'foo'
    objects[:, 1] = 'bar'

    float_cols = Index(['a', 'c', 'e'])
    object_cols = Index(['b', 'd'])
    columns = Index(sorted(float_cols + object_cols))
    new_columns = Index(['a', 'c', 'e', 'b', 'd'])

    fblock = _Block(floats, float_cols)
    oblock = _Block(objects, object_cols)

    blocks = [fblock, oblock]

    interleaved = _interleave(blocks, columns)
