from __future__ import division, absolute_import, print_function

import os
import math

from numpy.testing import *
from numpy import array

import util
import textwrap

def _path(*a):
    return os.path.join(*((os.path.dirname(__file__),) + a))

class TestMixed(util.F2PyTest):
    sources = [_path('src', 'mixed', 'foo.f'),
               _path('src', 'mixed', 'foo_fixed.f90'),
               _path('src', 'mixed', 'foo_free.f90')]

    @dec.slow
    def test_all(self):
        assert_( self.module.bar11() == 11)
        assert_( self.module.foo_fixed.bar12() == 12)
        assert_( self.module.foo_free.bar13() == 13)

    @dec.slow
    def test_docstring(self):
        expected = """
        a = bar11()

        Wrapper for ``bar11``.

        Returns
        -------
        a : int
        """
        assert_equal(self.module.bar11.__doc__, textwrap.dedent(expected).lstrip())

if __name__ == "__main__":
    import nose
    nose.runmodule()
