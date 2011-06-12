import unittest

import numpy as np

from pandas import Index
from pandas.core.internals import *

from pandas.util.testing import (assert_almost_equal, randn)

class TestBlock(unittest.TestCase):

    def test_merge(self):
        pass

    def test_copy(self):
        pass

    def test_assign_columns(self):
        pass

    def test_reindex_index(self):
        pass

    def test_reindex_columns_from(self):
        n = 10
        floats = np.repeat(np.atleast_2d(np.arange(3.)), n, axis=0)
        block = make_block(floats, [0, 2, 4],
                           ['a', 'b', 'c', 'd', 'e'])

        new_cols = Index(['e', 'b', 'c', 'f'])

        reindexed = block.reindex_columns(new_cols)
        assert_almost_equal(reindexed.ref_locs, [0, 2])
        self.assertEquals(reindexed.values.shape[1], 2)
        self.assert_((reindexed.values[:, 0] == 2).all())
        self.assert_((reindexed.values[:, 0] == 1).all())

    def test_insert(self):
        pass

    def test_delete(self):
        pass

    def test_get(self):
        pass

    def test_set(self):
        pass

    def test_fillna(self):
        pass

    def test_repr(self):
        pass


class TestBlockManager(unittest.TestCase):

    def test_interleave(self):
        pass

    def test_consolidate(self):
        pass

    def test_xs(self):
        pass


if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--pdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

