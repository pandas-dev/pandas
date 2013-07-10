import os
import nose
import unittest

import numpy as np
from numpy.testing import assert_equal

from pandas.tools.util import cartesian_product

class TestCartesianProduct(unittest.TestCase):

    def test_simple(self):
        x, y = list('ABC'), [1, 22]
        result = cartesian_product([x, y])
        expected = [np.array(['A', 'A', 'B', 'B', 'C', 'C']),
                    np.array([ 1, 22,  1, 22,  1, 22])]
        assert_equal(result, expected)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
