from __future__ import division, absolute_import, print_function

import sys
from numpy.testing import *
from swig_ext import example

class TestExample(TestCase):
    def test_fact(self):
        assert_equal(example.fact(10), 3628800)

    def test_cvar(self):
        assert_equal(example.cvar.My_variable, 3.0)
        example.cvar.My_variable = 5
        assert_equal(example.cvar.My_variable, 5.0)


if __name__ == "__main__":
    run_module_suite()
