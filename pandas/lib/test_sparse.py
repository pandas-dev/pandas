from unittest import TestCase

import nose
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal

from sparse import DenseIndex, BlockIndex, SparseVector
import sparse

class TestDenseIndex(TestCase):

    def setUp(self):
        pass

    def test_intersect(self):
        pass

class TestBlockIndex(TestCase):

    def setUp(self):
        pass

    def test_check_integrity(self):
        locs = []
        lengths = []

        # 0-length OK
        index = BlockIndex(0, locs, lengths)

        # also OK even though empty
        index = BlockIndex(1, locs, lengths)

        # block extend beyond end
        self.assertRaises(Exception, BlockIndex, 10, [5], [10])

        # block overlap
        self.assertRaises(Exception, BlockIndex, 10, [2, 5], [5, 3])

    def test_intersect(self):

        length = 20

        def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
            xindex = BlockIndex(length, xloc, xlen)
            yindex = BlockIndex(length, yloc, ylen)
            result = xindex.intersect(yindex)

            self.assert_(isinstance(result, BlockIndex))

            assert_equal(result.blocs, eloc)
            assert_equal(result.blengths, elen)

        # plain old test case

        xlocs = [0, 7, 15]
        xlengths = [3, 5, 5]
        ylocs = [2, 9, 14]
        ylengths = [2, 3, 5]
        exp_locs = [2, 9, 15]
        exp_lengths = [1, 3, 4]

        _check_case(xlocs, xlengths, ylocs, ylengths,
                    exp_locs, exp_lengths)

        # delete blocks
        xlocs = [0, 5]; xlengths = [4, 4]
        ylocs = [1]; ylengths = [4]
        exp_locs = [1]
        exp_lengths = [3]

        _check_case(xlocs, xlengths, ylocs, ylengths,
                    exp_locs, exp_lengths)

        # split blocks
        xlocs = [0]
        xlengths = [10]
        ylocs = [0, 5]
        ylengths = [3, 7]
        exp_locs = [0, 5]
        exp_lengths = [3, 5]

        _check_case(xlocs, xlengths, ylocs, ylengths,
                    exp_locs, exp_lengths)

        # skip block
        xlocs = [10]
        xlengths = [5]
        ylocs = [0, 12]
        ylengths = [5, 3]
        exp_locs = [12]
        exp_lengths = [3]
        _check_case(xlocs, xlengths, ylocs, ylengths,
                    exp_locs, exp_lengths)

        # no intersection
        xlocs = [0, 10]
        xlengths = [4, 6]
        ylocs = [5, 17]
        ylengths = [4, 2]
        exp_locs = []
        exp_lengths = []

        _check_case(xlocs, xlengths, ylocs, ylengths,
                    exp_locs, exp_lengths)

        # one or both is empty
        _check_case([0], [5], [], [], [], [])
        _check_case([], [], [], [], [], [])

    def test_to_dense(self):
        locs = [0, 10]
        lengths = [4, 6]
        exp_inds = [0, 1, 2, 3, 10, 11, 12, 13, 14, 15]

        block = BlockIndex(20, locs, lengths)
        dense = block.to_dense()

        assert_equal(dense.indices, exp_inds)

class TestSparseVector(TestCase):
    pass


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

