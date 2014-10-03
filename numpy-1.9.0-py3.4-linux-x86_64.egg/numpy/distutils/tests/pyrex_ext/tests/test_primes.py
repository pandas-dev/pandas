from __future__ import division, absolute_import, print_function

import sys
from numpy.testing import *
from pyrex_ext.primes import primes

class TestPrimes(TestCase):
    def test_simple(self, level=1):
        l = primes(10)
        assert_equal(l, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29])


if __name__ == "__main__":
    run_module_suite()
