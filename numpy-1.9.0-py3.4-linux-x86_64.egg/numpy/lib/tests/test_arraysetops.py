"""Test functions for 1D array set operations.

"""
from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.testing import (
    run_module_suite, TestCase, assert_array_equal
    )
from numpy.lib.arraysetops import (
    ediff1d, intersect1d, setxor1d, union1d, setdiff1d, unique, in1d
    )


class TestSetOps(TestCase):

    def test_unique(self):

        def check_all(a, b, i1, i2, c, dt):
            base_msg = 'check {0} failed for type {1}'

            msg = base_msg.format('values', dt)
            v = unique(a)
            assert_array_equal(v, b, msg)

            msg = base_msg.format('return_index', dt)
            v, j = unique(a, 1, 0, 0)
            assert_array_equal(v, b, msg)
            assert_array_equal(j, i1, msg)

            msg = base_msg.format('return_inverse', dt)
            v, j = unique(a, 0, 1, 0)
            assert_array_equal(v, b, msg)
            assert_array_equal(j, i2, msg)

            msg = base_msg.format('return_counts', dt)
            v, j = unique(a, 0, 0, 1)
            assert_array_equal(v, b, msg)
            assert_array_equal(j, c, msg)

            msg = base_msg.format('return_index and return_inverse', dt)
            v, j1, j2 = unique(a, 1, 1, 0)
            assert_array_equal(v, b, msg)
            assert_array_equal(j1, i1, msg)
            assert_array_equal(j2, i2, msg)

            msg = base_msg.format('return_index and return_counts', dt)
            v, j1, j2 = unique(a, 1, 0, 1)
            assert_array_equal(v, b, msg)
            assert_array_equal(j1, i1, msg)
            assert_array_equal(j2, c, msg)

            msg = base_msg.format('return_inverse and return_counts', dt)
            v, j1, j2 = unique(a, 0, 1, 1)
            assert_array_equal(v, b, msg)
            assert_array_equal(j1, i2, msg)
            assert_array_equal(j2, c, msg)

            msg = base_msg.format(('return_index, return_inverse '
                                   'and return_counts'), dt)
            v, j1, j2, j3 = unique(a, 1, 1, 1)
            assert_array_equal(v, b, msg)
            assert_array_equal(j1, i1, msg)
            assert_array_equal(j2, i2, msg)
            assert_array_equal(j3, c, msg)

        a = [5, 7, 1, 2, 1, 5, 7]*10
        b = [1, 2, 5, 7]
        i1 = [2, 3, 0, 1]
        i2 = [2, 3, 0, 1, 0, 2, 3]*10
        c = np.multiply([2, 1, 2, 2], 10)

        # test for numeric arrays
        types = []
        types.extend(np.typecodes['AllInteger'])
        types.extend(np.typecodes['AllFloat'])
        types.append('datetime64[D]')
        types.append('timedelta64[D]')
        for dt in types:
            aa = np.array(a, dt)
            bb = np.array(b, dt)
            check_all(aa, bb, i1, i2, c, dt)

        # test for object arrays
        dt = 'O'
        aa = np.empty(len(a), dt)
        aa[:] = a
        bb = np.empty(len(b), dt)
        bb[:] = b
        check_all(aa, bb, i1, i2, c, dt)

        # test for structured arrays
        dt = [('', 'i'), ('', 'i')]
        aa = np.array(list(zip(a, a)), dt)
        bb = np.array(list(zip(b, b)), dt)
        check_all(aa, bb, i1, i2, c, dt)

        # test for ticket #2799
        aa = [1. + 0.j, 1 - 1.j, 1]
        assert_array_equal(np.unique(aa), [1. - 1.j, 1. + 0.j])

        # test for ticket #4785
        a = [(1, 2), (1, 2), (2, 3)]
        unq = [1, 2, 3]
        inv = [0, 1, 0, 1, 1, 2]
        a1 = unique(a)
        assert_array_equal(a1, unq)
        a2, a2_inv = unique(a, return_inverse=True)
        assert_array_equal(a2, unq)
        assert_array_equal(a2_inv, inv)

    def test_intersect1d(self):
        # unique inputs
        a = np.array([5, 7, 1, 2])
        b = np.array([2, 4, 3, 1, 5])

        ec = np.array([1, 2, 5])
        c = intersect1d(a, b, assume_unique=True)
        assert_array_equal(c, ec)

        # non-unique inputs
        a = np.array([5, 5, 7, 1, 2])
        b = np.array([2, 1, 4, 3, 3, 1, 5])

        ed = np.array([1, 2, 5])
        c = intersect1d(a, b)
        assert_array_equal(c, ed)

        assert_array_equal([], intersect1d([], []))

    def test_setxor1d(self):
        a = np.array([5, 7, 1, 2])
        b = np.array([2, 4, 3, 1, 5])

        ec = np.array([3, 4, 7])
        c = setxor1d(a, b)
        assert_array_equal(c, ec)

        a = np.array([1, 2, 3])
        b = np.array([6, 5, 4])

        ec = np.array([1, 2, 3, 4, 5, 6])
        c = setxor1d(a, b)
        assert_array_equal(c, ec)

        a = np.array([1, 8, 2, 3])
        b = np.array([6, 5, 4, 8])

        ec = np.array([1, 2, 3, 4, 5, 6])
        c = setxor1d(a, b)
        assert_array_equal(c, ec)

        assert_array_equal([], setxor1d([], []))

    def test_ediff1d(self):
        zero_elem = np.array([])
        one_elem = np.array([1])
        two_elem = np.array([1, 2])

        assert_array_equal([], ediff1d(zero_elem))
        assert_array_equal([0], ediff1d(zero_elem, to_begin=0))
        assert_array_equal([0], ediff1d(zero_elem, to_end=0))
        assert_array_equal([-1, 0], ediff1d(zero_elem, to_begin=-1, to_end=0))
        assert_array_equal([], ediff1d(one_elem))
        assert_array_equal([1], ediff1d(two_elem))

    def test_in1d(self):
        # we use two different sizes for the b array here to test the
        # two different paths in in1d().
        for mult in (1, 10):
            # One check without np.array, to make sure lists are handled correct
            a = [5, 7, 1, 2]
            b = [2, 4, 3, 1, 5] * mult
            ec = np.array([True, False, True, True])
            c = in1d(a, b, assume_unique=True)
            assert_array_equal(c, ec)

            a[0] = 8
            ec = np.array([False, False, True, True])
            c = in1d(a, b, assume_unique=True)
            assert_array_equal(c, ec)

            a[0], a[3] = 4, 8
            ec = np.array([True, False, True, False])
            c = in1d(a, b, assume_unique=True)
            assert_array_equal(c, ec)

            a = np.array([5, 4, 5, 3, 4, 4, 3, 4, 3, 5, 2, 1, 5, 5])
            b = [2, 3, 4] * mult
            ec = [False, True, False, True, True, True, True, True, True, False,
                  True, False, False, False]
            c = in1d(a, b)
            assert_array_equal(c, ec)

            b = b + [5, 5, 4] * mult
            ec = [True, True, True, True, True, True, True, True, True, True,
                  True, False, True, True]
            c = in1d(a, b)
            assert_array_equal(c, ec)

            a = np.array([5, 7, 1, 2])
            b = np.array([2, 4, 3, 1, 5] * mult)
            ec = np.array([True, False, True, True])
            c = in1d(a, b)
            assert_array_equal(c, ec)

            a = np.array([5, 7, 1, 1, 2])
            b = np.array([2, 4, 3, 3, 1, 5] * mult)
            ec = np.array([True, False, True, True, True])
            c = in1d(a, b)
            assert_array_equal(c, ec)

            a = np.array([5, 5])
            b = np.array([2, 2] * mult)
            ec = np.array([False, False])
            c = in1d(a, b)
            assert_array_equal(c, ec)

        a = np.array([5])
        b = np.array([2])
        ec = np.array([False])
        c = in1d(a, b)
        assert_array_equal(c, ec)

        assert_array_equal(in1d([], []), [])

    def test_in1d_char_array(self):
        a = np.array(['a', 'b', 'c', 'd', 'e', 'c', 'e', 'b'])
        b = np.array(['a', 'c'])

        ec = np.array([True, False, True, False, False, True, False, False])
        c = in1d(a, b)

        assert_array_equal(c, ec)

    def test_in1d_invert(self):
        "Test in1d's invert parameter"
        # We use two different sizes for the b array here to test the
        # two different paths in in1d().
        for mult in (1, 10):
            a = np.array([5, 4, 5, 3, 4, 4, 3, 4, 3, 5, 2, 1, 5, 5])
            b = [2, 3, 4] * mult
            assert_array_equal(np.invert(in1d(a, b)), in1d(a, b, invert=True))

    def test_in1d_ravel(self):
        # Test that in1d ravels its input arrays. This is not documented
        # behavior however. The test is to ensure consistentency.
        a = np.arange(6).reshape(2, 3)
        b = np.arange(3, 9).reshape(3, 2)
        long_b = np.arange(3, 63).reshape(30, 2)
        ec = np.array([False, False, False, True, True, True])

        assert_array_equal(in1d(a, b, assume_unique=True), ec)
        assert_array_equal(in1d(a, b, assume_unique=False), ec)
        assert_array_equal(in1d(a, long_b, assume_unique=True), ec)
        assert_array_equal(in1d(a, long_b, assume_unique=False), ec)

    def test_union1d(self):
        a = np.array([5, 4, 7, 1, 2])
        b = np.array([2, 4, 3, 3, 2, 1, 5])

        ec = np.array([1, 2, 3, 4, 5, 7])
        c = union1d(a, b)
        assert_array_equal(c, ec)

        assert_array_equal([], union1d([], []))

    def test_setdiff1d(self):
        a = np.array([6, 5, 4, 7, 1, 2, 7, 4])
        b = np.array([2, 4, 3, 3, 2, 1, 5])

        ec = np.array([6, 7])
        c = setdiff1d(a, b)
        assert_array_equal(c, ec)

        a = np.arange(21)
        b = np.arange(19)
        ec = np.array([19, 20])
        c = setdiff1d(a, b)
        assert_array_equal(c, ec)

        assert_array_equal([], setdiff1d([], []))

    def test_setdiff1d_char_array(self):
        a = np.array(['a', 'b', 'c'])
        b = np.array(['a', 'b', 's'])
        assert_array_equal(setdiff1d(a, b), np.array(['c']))

    def test_manyways(self):
        a = np.array([5, 7, 1, 2, 8])
        b = np.array([9, 8, 2, 4, 3, 1, 5])

        c1 = setxor1d(a, b)
        aux1 = intersect1d(a, b)
        aux2 = union1d(a, b)
        c2 = setdiff1d(aux2, aux1)
        assert_array_equal(c1, c2)


if __name__ == "__main__":
    run_module_suite()
