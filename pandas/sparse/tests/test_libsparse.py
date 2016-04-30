from pandas import Series

import nose  # noqa
import numpy as np
import operator
from numpy.testing import assert_equal
import pandas.util.testing as tm

from pandas import compat

from pandas.sparse.array import IntIndex, BlockIndex, _make_index
import pandas._sparse as splib

TEST_LENGTH = 20

plain_case = dict(xloc=[0, 7, 15], xlen=[3, 5, 5], yloc=[2, 9, 14],
                  ylen=[2, 3, 5], intersect_loc=[2, 9, 15],
                  intersect_len=[1, 3, 4])
delete_blocks = dict(xloc=[0, 5], xlen=[4, 4], yloc=[1], ylen=[4],
                     intersect_loc=[1], intersect_len=[3])
split_blocks = dict(xloc=[0], xlen=[10], yloc=[0, 5], ylen=[3, 7],
                    intersect_loc=[0, 5], intersect_len=[3, 5])
skip_block = dict(xloc=[10], xlen=[5], yloc=[0, 12], ylen=[5, 3],
                  intersect_loc=[12], intersect_len=[3])

no_intersect = dict(xloc=[0, 10], xlen=[4, 6], yloc=[5, 17], ylen=[4, 2],
                    intersect_loc=[], intersect_len=[])


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


class TestSparseIndexUnion(tm.TestCase):

    def test_index_make_union(self):
        def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
            xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
            yindex = BlockIndex(TEST_LENGTH, yloc, ylen)
            bresult = xindex.make_union(yindex)
            assert (isinstance(bresult, BlockIndex))
            assert_equal(bresult.blocs, eloc)
            assert_equal(bresult.blengths, elen)

            ixindex = xindex.to_int_index()
            iyindex = yindex.to_int_index()
            iresult = ixindex.make_union(iyindex)
            assert (isinstance(iresult, IntIndex))
            assert_equal(iresult.indices, bresult.to_int_index().indices)

        """
        x: ----
        y:     ----
        r: --------
        """
        xloc = [0]
        xlen = [5]
        yloc = [5]
        ylen = [4]
        eloc = [0]
        elen = [9]
        _check_case(xloc, xlen, yloc, ylen, eloc, elen)
        """
        x: -----     -----
        y:   -----          --
        """
        xloc = [0, 10]
        xlen = [5, 5]
        yloc = [2, 17]
        ylen = [5, 2]
        eloc = [0, 10, 17]
        elen = [7, 5, 2]
        _check_case(xloc, xlen, yloc, ylen, eloc, elen)
        """
        x: ------
        y:    -------
        r: ----------
        """
        xloc = [1]
        xlen = [5]
        yloc = [3]
        ylen = [5]
        eloc = [1]
        elen = [7]
        _check_case(xloc, xlen, yloc, ylen, eloc, elen)
        """
        x: ------  -----
        y:    -------
        r: -------------
        """
        xloc = [2, 10]
        xlen = [4, 4]
        yloc = [4]
        ylen = [8]
        eloc = [2]
        elen = [12]
        _check_case(xloc, xlen, yloc, ylen, eloc, elen)
        """
        x: ---  -----
        y: -------
        r: -------------
        """
        xloc = [0, 5]
        xlen = [3, 5]
        yloc = [0]
        ylen = [7]
        eloc = [0]
        elen = [10]
        _check_case(xloc, xlen, yloc, ylen, eloc, elen)
        """
        x: ------  -----
        y:    -------  ---
        r: -------------
        """
        xloc = [2, 10]
        xlen = [4, 4]
        yloc = [4, 13]
        ylen = [8, 4]
        eloc = [2]
        elen = [15]
        _check_case(xloc, xlen, yloc, ylen, eloc, elen)
        """
        x: ----------------------
        y:   ----  ----   ---
        r: ----------------------
        """
        xloc = [2]
        xlen = [15]
        yloc = [4, 9, 14]
        ylen = [3, 2, 2]
        eloc = [2]
        elen = [15]
        _check_case(xloc, xlen, yloc, ylen, eloc, elen)
        """
        x: ----       ---
        y:       ---       ---
        """
        xloc = [0, 10]
        xlen = [3, 3]
        yloc = [5, 15]
        ylen = [2, 2]
        eloc = [0, 5, 10, 15]
        elen = [3, 2, 3, 2]
        _check_case(xloc, xlen, yloc, ylen, eloc, elen)

    def test_intindex_make_union(self):
        a = IntIndex(5, np.array([0, 3, 4], dtype=np.int32))
        b = IntIndex(5, np.array([0, 2], dtype=np.int32))
        res = a.make_union(b)
        exp = IntIndex(5, np.array([0, 2, 3, 4], np.int32))
        self.assertTrue(res.equals(exp))

        a = IntIndex(5, np.array([], dtype=np.int32))
        b = IntIndex(5, np.array([0, 2], dtype=np.int32))
        res = a.make_union(b)
        exp = IntIndex(5, np.array([0, 2], np.int32))
        self.assertTrue(res.equals(exp))

        a = IntIndex(5, np.array([], dtype=np.int32))
        b = IntIndex(5, np.array([], dtype=np.int32))
        res = a.make_union(b)
        exp = IntIndex(5, np.array([], np.int32))
        self.assertTrue(res.equals(exp))

        a = IntIndex(5, np.array([0, 1, 2, 3, 4], dtype=np.int32))
        b = IntIndex(5, np.array([0, 1, 2, 3, 4], dtype=np.int32))
        res = a.make_union(b)
        exp = IntIndex(5, np.array([0, 1, 2, 3, 4], np.int32))
        self.assertTrue(res.equals(exp))

        a = IntIndex(5, np.array([0, 1], dtype=np.int32))
        b = IntIndex(4, np.array([0, 1], dtype=np.int32))
        with tm.assertRaises(ValueError):
            a.make_union(b)


class TestSparseIndexCommon(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_lookup(self):
        for kind in ['integer', 'block']:
            idx = _make_index(4, np.array([2, 3], dtype=np.int32), kind=kind)
            self.assertEqual(idx.lookup(-1), -1)
            self.assertEqual(idx.lookup(0), -1)
            self.assertEqual(idx.lookup(1), -1)
            self.assertEqual(idx.lookup(2), 0)
            self.assertEqual(idx.lookup(3), 1)
            self.assertEqual(idx.lookup(4), -1)

            idx = _make_index(4, np.array([], dtype=np.int32), kind=kind)

            for i in range(-1, 5):
                self.assertEqual(idx.lookup(i), -1)

            idx = _make_index(4, np.array([0, 1, 2, 3], dtype=np.int32),
                              kind=kind)
            self.assertEqual(idx.lookup(-1), -1)
            self.assertEqual(idx.lookup(0), 0)
            self.assertEqual(idx.lookup(1), 1)
            self.assertEqual(idx.lookup(2), 2)
            self.assertEqual(idx.lookup(3), 3)
            self.assertEqual(idx.lookup(4), -1)

            idx = _make_index(4, np.array([0, 2, 3], dtype=np.int32),
                              kind=kind)
            self.assertEqual(idx.lookup(-1), -1)
            self.assertEqual(idx.lookup(0), 0)
            self.assertEqual(idx.lookup(1), -1)
            self.assertEqual(idx.lookup(2), 1)
            self.assertEqual(idx.lookup(3), 2)
            self.assertEqual(idx.lookup(4), -1)

    def test_lookup_array(self):
        for kind in ['integer', 'block']:
            idx = _make_index(4, np.array([2, 3], dtype=np.int32), kind=kind)

            res = idx.lookup_array(np.array([-1, 0, 2], dtype=np.int32))
            exp = np.array([-1, -1, 0], dtype=np.int32)
            self.assert_numpy_array_equal(res, exp)

            res = idx.lookup_array(np.array([4, 2, 1, 3], dtype=np.int32))
            exp = np.array([-1, 0, -1, 1], dtype=np.int32)
            self.assert_numpy_array_equal(res, exp)

            idx = _make_index(4, np.array([], dtype=np.int32), kind=kind)
            res = idx.lookup_array(np.array([-1, 0, 2, 4], dtype=np.int32))
            exp = np.array([-1, -1, -1, -1], dtype=np.int32)

            idx = _make_index(4, np.array([0, 1, 2, 3], dtype=np.int32),
                              kind=kind)
            res = idx.lookup_array(np.array([-1, 0, 2], dtype=np.int32))
            exp = np.array([-1, 0, 2], dtype=np.int32)
            self.assert_numpy_array_equal(res, exp)

            res = idx.lookup_array(np.array([4, 2, 1, 3], dtype=np.int32))
            exp = np.array([-1, 2, 1, 3], dtype=np.int32)
            self.assert_numpy_array_equal(res, exp)

            idx = _make_index(4, np.array([0, 2, 3], dtype=np.int32),
                              kind=kind)
            res = idx.lookup_array(np.array([2, 1, 3, 0], dtype=np.int32))
            exp = np.array([1, -1, 2, 0], dtype=np.int32)
            self.assert_numpy_array_equal(res, exp)

            res = idx.lookup_array(np.array([1, 4, 2, 5], dtype=np.int32))
            exp = np.array([-1, -1, 1, -1], dtype=np.int32)
            self.assert_numpy_array_equal(res, exp)

    def test_lookup_basics(self):
        def _check(index):
            assert (index.lookup(0) == -1)
            assert (index.lookup(5) == 0)
            assert (index.lookup(7) == 2)
            assert (index.lookup(8) == -1)
            assert (index.lookup(9) == -1)
            assert (index.lookup(10) == -1)
            assert (index.lookup(11) == -1)
            assert (index.lookup(12) == 3)
            assert (index.lookup(17) == 8)
            assert (index.lookup(18) == -1)

        bindex = BlockIndex(20, [5, 12], [3, 6])
        iindex = bindex.to_int_index()

        _check(bindex)
        _check(iindex)

        # corner cases


def test_intersect():
    def _check_correct(a, b, expected):
        result = a.intersect(b)
        assert (result.equals(expected))

    def _check_length_exc(a, longer):
        nose.tools.assert_raises(Exception, a.intersect, longer)

    def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
        xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
        yindex = BlockIndex(TEST_LENGTH, yloc, ylen)
        expected = BlockIndex(TEST_LENGTH, eloc, elen)
        longer_index = BlockIndex(TEST_LENGTH + 1, yloc, ylen)

        _check_correct(xindex, yindex, expected)
        _check_correct(xindex.to_int_index(), yindex.to_int_index(),
                       expected.to_int_index())

        _check_length_exc(xindex, longer_index)
        _check_length_exc(xindex.to_int_index(), longer_index.to_int_index())

    if compat.is_platform_windows():
        raise nose.SkipTest("segfaults on win-64 when all tests are run")
    check_cases(_check_case)


class TestBlockIndex(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_block_internal(self):
        idx = _make_index(4, np.array([2, 3], dtype=np.int32), kind='block')
        self.assertIsInstance(idx, BlockIndex)
        self.assertEqual(idx.npoints, 2)
        tm.assert_numpy_array_equal(idx.blocs,
                                    np.array([2], dtype=np.int32))
        tm.assert_numpy_array_equal(idx.blengths,
                                    np.array([2], dtype=np.int32))

        idx = _make_index(4, np.array([], dtype=np.int32), kind='block')
        self.assertIsInstance(idx, BlockIndex)
        self.assertEqual(idx.npoints, 0)
        tm.assert_numpy_array_equal(idx.blocs,
                                    np.array([], dtype=np.int32))
        tm.assert_numpy_array_equal(idx.blengths,
                                    np.array([], dtype=np.int32))

        idx = _make_index(4, np.array([0, 1, 2, 3], dtype=np.int32),
                          kind='block')
        self.assertIsInstance(idx, BlockIndex)
        self.assertEqual(idx.npoints, 4)
        tm.assert_numpy_array_equal(idx.blocs,
                                    np.array([0], dtype=np.int32))
        tm.assert_numpy_array_equal(idx.blengths,
                                    np.array([4], dtype=np.int32))

        idx = _make_index(4, np.array([0, 2, 3], dtype=np.int32), kind='block')
        self.assertIsInstance(idx, BlockIndex)
        self.assertEqual(idx.npoints, 3)
        tm.assert_numpy_array_equal(idx.blocs,
                                    np.array([0, 2], dtype=np.int32))
        tm.assert_numpy_array_equal(idx.blengths,
                                    np.array([1, 2], dtype=np.int32))

    def test_equals(self):
        index = BlockIndex(10, [0, 4], [2, 5])

        self.assertTrue(index.equals(index))
        self.assertFalse(index.equals(BlockIndex(10, [0, 4], [2, 6])))

    def test_check_integrity(self):
        locs = []
        lengths = []

        # 0-length OK
        # TODO: index variables are not used...is that right?
        index = BlockIndex(0, locs, lengths)  # noqa

        # also OK even though empty
        index = BlockIndex(1, locs, lengths)  # noqa

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
        self.assertIs(index.to_block_index(), index)


class TestIntIndex(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_int_internal(self):
        idx = _make_index(4, np.array([2, 3], dtype=np.int32), kind='integer')
        self.assertIsInstance(idx, IntIndex)
        self.assertEqual(idx.npoints, 2)
        tm.assert_numpy_array_equal(idx.indices,
                                    np.array([2, 3], dtype=np.int32))

        idx = _make_index(4, np.array([], dtype=np.int32), kind='integer')
        self.assertIsInstance(idx, IntIndex)
        self.assertEqual(idx.npoints, 0)
        tm.assert_numpy_array_equal(idx.indices,
                                    np.array([], dtype=np.int32))

        idx = _make_index(4, np.array([0, 1, 2, 3], dtype=np.int32),
                          kind='integer')
        self.assertIsInstance(idx, IntIndex)
        self.assertEqual(idx.npoints, 4)
        tm.assert_numpy_array_equal(idx.indices,
                                    np.array([0, 1, 2, 3], dtype=np.int32))

    def test_equals(self):
        index = IntIndex(10, [0, 1, 2, 3, 4])
        self.assertTrue(index.equals(index))
        self.assertFalse(index.equals(IntIndex(10, [0, 1, 2, 3])))

    def test_to_block_index(self):
        def _check_case(xloc, xlen, yloc, ylen, eloc, elen):
            xindex = BlockIndex(TEST_LENGTH, xloc, xlen)
            yindex = BlockIndex(TEST_LENGTH, yloc, ylen)

            # see if survive the round trip
            xbindex = xindex.to_int_index().to_block_index()
            ybindex = yindex.to_int_index().to_block_index()
            tm.assertIsInstance(xbindex, BlockIndex)
            self.assertTrue(xbindex.equals(xindex))
            self.assertTrue(ybindex.equals(yindex))

        check_cases(_check_case)

    def test_to_int_index(self):
        index = IntIndex(10, [2, 3, 4, 5, 6])
        self.assertIs(index.to_int_index(), index)


class TestSparseOperators(tm.TestCase):

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

            result_block_vals, rb_index = sparse_op(x, xindex, xfill, y,
                                                    yindex, yfill)
            result_int_vals, ri_index = sparse_op(x, xdindex, xfill, y,
                                                  ydindex, yfill)

            self.assertTrue(rb_index.to_int_index().equals(ri_index))
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


def make_optestf(op):
    def f(self):
        sparse_op = getattr(splib, 'sparse_%s' % op)
        python_op = getattr(operator, op)
        self._op_tests(sparse_op, python_op)

    f.__name__ = 'test_%s' % op
    return f


for op in check_ops:
    g = make_optestf(op)
    setattr(TestSparseOperators, g.__name__, g)
    del g


if __name__ == '__main__':
    import nose  # noqa
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
