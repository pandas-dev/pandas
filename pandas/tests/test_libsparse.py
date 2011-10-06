from unittest import TestCase

from pandas import Series

import nose
from numpy import nan
import numpy as np
import operator
from numpy.testing import assert_almost_equal, assert_equal

from pandas.core.sparse import SparseSeries
from pandas import DataFrame

from pandas._sparse import IntIndex, BlockIndex
import pandas._sparse as splib

TEST_LENGTH = 20

plain_case = dict(xloc = [0, 7, 15],
                  xlen = [3, 5, 5],
                  yloc = [2, 9, 14],
                  ylen = [2, 3, 5],
                  intersect_loc = [2, 9, 15],
                  intersect_len = [1, 3, 4])
delete_blocks = dict(xloc = [0, 5],
                     xlen = [4, 4],
                     yloc = [1],
                     ylen = [4],
                     intersect_loc = [1],
                     intersect_len = [3])
split_blocks = dict(xloc = [0],
                    xlen = [10],
                    yloc = [0, 5],
                    ylen = [3, 7],
                    intersect_loc = [0, 5],
                    intersect_len = [3, 5])
skip_block = dict(xloc = [10],
                  xlen = [5],
                  yloc = [0, 12],
                  ylen = [5, 3],
                  intersect_loc = [12],
                  intersect_len = [3])

no_intersect = dict(xloc = [0, 10],
                    xlen = [4, 6],
                    yloc = [5, 17],
                    ylen = [4, 2],
                    intersect_loc = [],
                    intersect_len = [])

def check_cases(_check_case):
    def _check_case_dict(case):
        _check_case(case['xloc'], case['xlen'], case['yloc'], case['ylen'],
                    case['intersect_loc'], case['intersect_len'])

    _check_case_dict(plain_case)
    _check_case_dict(delete_blocks)
    _check_case_dict(split_blocks)
    _check_case_dict(skip_block)
    _check_case_dict(no_intersect)

    # one or both is empty
    _check_case([0], [5], [], [], [], [])
    _check_case([], [], [], [], [], [])

def test_index_make_union():
    def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
        xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
        yindex = BlockIndex(TEST_LENGTH, yloc, ylen)
        bresult = xindex.make_union(yindex)
        assert(isinstance(bresult, BlockIndex))
        assert_equal(bresult.blocs, eloc)
        assert_equal(bresult.blengths, elen)

        ixindex = xindex.to_int_index()
        iyindex = yindex.to_int_index()
        iresult = ixindex.make_union(iyindex)
        assert(isinstance(iresult, IntIndex))
        assert_equal(iresult.indices, bresult.to_int_index().indices)

    """
    x: ----
    y:     ----
    r: --------
    """
    xloc = [0]; xlen = [5]
    yloc = [5]; ylen = [4]
    eloc = [0]; elen = [9]
    _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    """
    x: -----     -----
    y:   -----          --
    """
    xloc = [0, 10]; xlen = [5, 5]
    yloc = [2, 17]; ylen = [5, 2]
    eloc = [0, 10, 17]; elen = [7, 5, 2]
    _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    """
    x: ------
    y:    -------
    r: ----------
    """
    xloc = [1]; xlen = [5]
    yloc = [3]; ylen = [5]
    eloc = [1]; elen = [7]
    _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    """
    x: ------  -----
    y:    -------
    r: -------------
    """
    xloc = [2, 10]; xlen = [4, 4]
    yloc = [4]; ylen = [8]
    eloc = [2]; elen = [12]
    _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    """
    x: ---  -----
    y: -------
    r: -------------
    """
    xloc = [0, 5]; xlen = [3, 5]
    yloc = [0]; ylen = [7]
    eloc = [0]; elen = [10]
    _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    """
    x: ------  -----
    y:    -------  ---
    r: -------------
    """
    xloc = [2, 10]; xlen = [4, 4]
    yloc = [4, 13]; ylen = [8, 4]
    eloc = [2]; elen = [15]
    _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    """
    x: ----------------------
    y:   ----  ----   ---
    r: ----------------------
    """
    xloc = [2]; xlen = [15]
    yloc = [4, 9, 14]; ylen = [3, 2, 2]
    eloc = [2]; elen = [15]
    _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    """
    x: ----       ---
    y:       ---       ---
    """
    xloc = [0, 10]; xlen = [3, 3]
    yloc = [5, 15]; ylen = [2, 2]
    eloc = [0, 5, 10, 15]; elen = [3, 2, 3, 2]
    _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    # TODO: different-length index objects

def test_lookup():

    def _check(index):
        assert(index.lookup(0) == -1)
        assert(index.lookup(5) == 0)
        assert(index.lookup(7) == 2)
        assert(index.lookup(8) == -1)
        assert(index.lookup(9) == -1)
        assert(index.lookup(10) == -1)
        assert(index.lookup(11) == -1)
        assert(index.lookup(12) == 3)
        assert(index.lookup(17) == 8)
        assert(index.lookup(18) == -1)

    bindex = BlockIndex(20, [5, 12], [3, 6])
    iindex = bindex.to_int_index()

    _check(bindex)
    _check(iindex)

    # corner cases

def test_intersect():
    def _check_correct(a, b, expected):
        result = a.intersect(b)
        assert(result.equals(expected))

    def _check_length_exc(a, longer):
        nose.tools.assert_raises(Exception, a.intersect, longer)

    def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
        xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
        yindex = BlockIndex(TEST_LENGTH, yloc, ylen)
        expected = BlockIndex(TEST_LENGTH, eloc, elen)
        longer_index = BlockIndex(TEST_LENGTH + 1, yloc, ylen)

        _check_correct(xindex, yindex, expected)
        _check_correct(xindex.to_int_index(),
                       yindex.to_int_index(),
                       expected.to_int_index())

        _check_length_exc(xindex, longer_index)
        _check_length_exc(xindex.to_int_index(),
                          longer_index.to_int_index())

    check_cases(_check_case)

class TestBlockIndex(TestCase):

    def test_equals(self):
        index = BlockIndex(10, [0, 4], [2, 5])

        self.assert_(index.equals(index))
        self.assert_(not index.equals(BlockIndex(10, [0, 4], [2, 6])))

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

    def test_to_int_index(self):
        locs = [0, 10]
        lengths = [4, 6]
        exp_inds = [0, 1, 2, 3, 10, 11, 12, 13, 14, 15]

        block = BlockIndex(20, locs, lengths)
        dense = block.to_int_index()

        assert_equal(dense.indices, exp_inds)

    def test_to_block_index(self):
        index = BlockIndex(10, [0, 5], [4, 5])
        self.assert_(index.to_block_index() is index)

class TestIntIndex(TestCase):

    def test_equals(self):
        index = IntIndex(10, [0, 1, 2, 3, 4])
        self.assert_(index.equals(index))
        self.assert_(not index.equals(IntIndex(10, [0, 1, 2, 3])))

    def test_to_block_index(self):
        def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
            xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
            yindex = BlockIndex(TEST_LENGTH, yloc, ylen)

            # see if survive the round trip
            xbindex = xindex.to_int_index().to_block_index()
            ybindex = yindex.to_int_index().to_block_index()
            self.assert_(isinstance(xbindex, BlockIndex))
            self.assert_(xbindex.equals(xindex))
            self.assert_(ybindex.equals(yindex))
        check_cases(_check_case)

    def test_to_int_index(self):
        index = IntIndex(10, [2, 3, 4, 5, 6])
        self.assert_(index.to_int_index() is index)

class TestSparseOperators(TestCase):

    def _nan_op_tests(self, sparse_op, python_op):
        def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
            xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
            yindex = BlockIndex(TEST_LENGTH, yloc, ylen)

            xdindex = xindex.to_int_index()
            ydindex = yindex.to_int_index()

            x = np.arange(xindex.npoints) * 10. + 1
            y = np.arange(yindex.npoints) * 100. + 1

            result_block_vals, rb_index = sparse_op(x, xindex, y, yindex)
            result_int_vals, ri_index = sparse_op(x, xdindex, y, ydindex)

            self.assert_(rb_index.to_int_index().equals(ri_index))
            assert_equal(result_block_vals, result_int_vals)

            # check versus Series...
            xseries = Series(x, xdindex.indices)
            yseries = Series(y, ydindex.indices)
            series_result = python_op(xseries, yseries).valid()
            assert_equal(result_block_vals, series_result.values)
            assert_equal(result_int_vals, series_result.values)

        check_cases(_check_case)

    def _op_tests(self, sparse_op, python_op):
        def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
            xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
            yindex = BlockIndex(TEST_LENGTH, yloc, ylen)

            xdindex = xindex.to_int_index()
            ydindex = yindex.to_int_index()

            x = np.arange(xindex.npoints) * 10. + 1
            y = np.arange(yindex.npoints) * 100. + 1

            xfill = 0
            yfill = 2

            result_block_vals, rb_index = sparse_op(x, xindex, xfill, y, yindex, yfill)
            result_int_vals, ri_index = sparse_op(x, xdindex, xfill,
                                                  y, ydindex, yfill)

            self.assert_(rb_index.to_int_index().equals(ri_index))
            assert_equal(result_block_vals, result_int_vals)

            # check versus Series...
            xseries = Series(x, xdindex.indices)
            xseries = xseries.reindex(np.arange(TEST_LENGTH)).fillna(xfill)

            yseries = Series(y, ydindex.indices)
            yseries = yseries.reindex(np.arange(TEST_LENGTH)).fillna(yfill)

            series_result = python_op(xseries, yseries)
            series_result = series_result.reindex(ri_index.indices)

            assert_equal(result_block_vals, series_result.values)
            assert_equal(result_int_vals, series_result.values)

        check_cases(_check_case)

# too cute? oh but how I abhor code duplication

check_ops = ['add', 'sub', 'mul', 'truediv', 'floordiv']
def make_nanoptestf(op):
    def f(self):
        sparse_op = getattr(splib, 'sparse_nan%s' % op)
        python_op = getattr(operator, op)
        self._nan_op_tests(sparse_op, python_op)
    f.__name__ = 'test_nan%s' % op
    return f

def make_optestf(op):
    def f(self):
        sparse_op = getattr(splib, 'sparse_%s' % op)
        python_op = getattr(operator, op)
        self._op_tests(sparse_op, python_op)
    f.__name__ = 'test_%s' % op
    return f

for op in check_ops:
    f = make_nanoptestf(op)
    g = make_optestf(op)
    setattr(TestSparseOperators, f.__name__, f)
    setattr(TestSparseOperators, g.__name__, g)
    del f
    del g

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

