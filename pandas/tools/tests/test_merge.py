import nose
import unittest

import numpy as np

from pandas.tools.merge import merge
import pandas._sandbox as sbx

a_ = np.array

class TestMerge(unittest.TestCase):

    def setUp(self):
        pass

    def test_cython_left_outer_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype='i4')
        right = a_([1, 1, 0, 4, 2, 2, 1], dtype='i4')
        max_group = 5

        ls, rs = sbx.left_outer_join(left, right, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        exp_li = a_([0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
                     6, 6, 7, 7, 8, 8, 9, 10])
        exp_ri = a_([0, 0, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                     4, 5, 4, 5, 4, 5, -1, -1])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        self.assert_(np.array_equal(ls, exp_ls))
        self.assert_(np.array_equal(rs, exp_rs))

    def test_cython_right_outer_join(self):
        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype='i4')
        right = a_([1, 1, 0, 4, 2, 2, 1], dtype='i4')
        max_group = 5

        rs, ls  = sbx.left_outer_join(right, left, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        #            0        1        1        1
        exp_li = a_([0, 1, 2, 3, 4, 5, 3, 4, 5, 3, 4, 5,
        #            2        2        4
                     6, 7, 8, 6, 7, 8, -1])
        exp_ri = a_([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
                     4, 4, 4, 5, 5, 5, 6])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        self.assert_(np.array_equal(ls, exp_ls))
        self.assert_(np.array_equal(rs, exp_rs))

    def test_cython_inner_join(self):
        raise nose.SkipTest

        left = a_([0, 1, 2, 1, 2, 0, 0, 1, 2, 3, 3], dtype='i4')
        right = a_([1, 1, 0, 4, 2, 2, 1, 4], dtype='i4')
        max_group = 5

        ls, rs = sbx.inner_join(left, right, max_group)

        exp_ls = left.argsort(kind='mergesort')
        exp_rs = right.argsort(kind='mergesort')

        exp_li = a_([0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
                     6, 6, 7, 7, 8, 8])
        exp_ri = a_([0, 0, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                     4, 5, 4, 5, 4, 5])

        exp_ls = exp_ls.take(exp_li)
        exp_ls[exp_li == -1] = -1

        exp_rs = exp_rs.take(exp_ri)
        exp_rs[exp_ri == -1] = -1

        self.assert_(np.array_equal(ls, exp_ls))
        self.assert_(np.array_equal(rs, exp_rs))

    def test_cython_full_outer_join(self):
        pass

    def test_left_join(self):
        pass

    def test_merge_common(self):
        pass

    def test_merge_index(self):
        pass

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)



