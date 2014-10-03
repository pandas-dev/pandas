from __future__ import division, absolute_import, print_function

from numpy.testing import assert_equal, TestCase
from numpy.core import ones
from numpy import matrix

class TestDot(TestCase):
    def test_matscalar(self):
        b1 = matrix(ones((3, 3), dtype=complex))
        assert_equal(b1*1.0, b1)
