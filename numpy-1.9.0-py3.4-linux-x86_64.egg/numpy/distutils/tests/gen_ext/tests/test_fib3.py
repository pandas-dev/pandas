from __future__ import division, absolute_import, print_function

import sys
from numpy.testing import *
from gen_ext import fib3

class TestFib3(TestCase):
    def test_fib(self):
        assert_array_equal(fib3.fib(6), [0, 1, 1, 2, 3, 5])

if __name__ == "__main__":
    run_module_suite()
