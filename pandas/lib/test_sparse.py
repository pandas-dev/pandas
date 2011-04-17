from unittest import TestCase

import nose
import numpy as np
import operator
from numpy.testing import assert_almost_equal, assert_equal

from sparse import DenseIndex, BlockIndex, SparseVector
import sparse

class TestDenseIndex(TestCase):

    def setUp(self):
        pass

    def test_intersect(self):
        pass



TEST_LENGTH = 20

plain_case = dict(xloc = [0, 7, 15],
                  xlen = [3, 5, 5],
                  yloc = [2, 9, 14],
                  ylen = [2, 3, 5],
                  eloc = [2, 9, 15],
                  elen = [1, 3, 4])
delete_blocks = dict(xloc = [0, 5],
                     xlen = [4, 4],
                     yloc = [1],
                     ylen = [4],
                     eloc = [1],
                     elen = [3])
split_blocks = dict(xloc = [0],
                    xlen = [10],
                    yloc = [0, 5],
                    ylen = [3, 7],
                    eloc = [0, 5],
                    elen = [3, 5])
skip_block = dict(xloc = [10],
                  xlen = [5],
                  yloc = [0, 12],
                  ylen = [5, 3],
                  eloc = [12],
                  elen = [3])

no_intersect = dict(xloc = [0, 10],
                    xlen = [4, 6],
                    yloc = [5, 17],
                    ylen = [4, 2],
                    eloc = [],
                    elen = [])

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

        def _check_case_dict(case):
            _check_case(case['xloc'], case['xlen'], case['yloc'], case['ylen'],
                        case['eloc'], case['elen'])

        def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
            xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
            yindex = BlockIndex(TEST_LENGTH, yloc, ylen)
            result = xindex.intersect(yindex)

            self.assert_(isinstance(result, BlockIndex))

            assert_equal(result.blocs, eloc)
            assert_equal(result.blengths, elen)

        _check_case_dict(plain_case)
        _check_case_dict(delete_blocks)
        _check_case_dict(split_blocks)
        _check_case_dict(skip_block)
        _check_case_dict(no_intersect)

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

    def _arith_op_tests(self, op):
        pass

# too cute? oh but how I abhor code duplication

check_ops = ['add', 'sub', 'mul', 'div']
def make_optestf(op):
    def f(self):
        self._arith_op_tests(getattr(operator, op))
    f.__name__ = 'test_%s' % op
    return f

for op in check_ops:
    f = make_optestf(op)
    setattr(TestSparseVector, f.__name__, f)
    del f

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

