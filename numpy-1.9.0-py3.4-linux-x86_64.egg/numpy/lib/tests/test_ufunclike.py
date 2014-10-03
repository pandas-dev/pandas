from __future__ import division, absolute_import, print_function

import numpy.core as nx
import numpy.lib.ufunclike as ufl
from numpy.testing import (
    run_module_suite, TestCase, assert_, assert_equal, assert_array_equal
    )


class TestUfunclike(TestCase):

    def test_isposinf(self):
        a = nx.array([nx.inf, -nx.inf, nx.nan, 0.0, 3.0, -3.0])
        out = nx.zeros(a.shape, bool)
        tgt = nx.array([True, False, False, False, False, False])

        res = ufl.isposinf(a)
        assert_equal(res, tgt)
        res = ufl.isposinf(a, out)
        assert_equal(res, tgt)
        assert_equal(out, tgt)

    def test_isneginf(self):
        a = nx.array([nx.inf, -nx.inf, nx.nan, 0.0, 3.0, -3.0])
        out = nx.zeros(a.shape, bool)
        tgt = nx.array([False, True, False, False, False, False])

        res = ufl.isneginf(a)
        assert_equal(res, tgt)
        res = ufl.isneginf(a, out)
        assert_equal(res, tgt)
        assert_equal(out, tgt)

    def test_fix(self):
        a = nx.array([[1.0, 1.1, 1.5, 1.8], [-1.0, -1.1, -1.5, -1.8]])
        out = nx.zeros(a.shape, float)
        tgt = nx.array([[1., 1., 1., 1.], [-1., -1., -1., -1.]])

        res = ufl.fix(a)
        assert_equal(res, tgt)
        res = ufl.fix(a, out)
        assert_equal(res, tgt)
        assert_equal(out, tgt)
        assert_equal(ufl.fix(3.14), 3)

    def test_fix_with_subclass(self):
        class MyArray(nx.ndarray):
            def __new__(cls, data, metadata=None):
                res = nx.array(data, copy=True).view(cls)
                res.metadata = metadata
                return res

            def __array_wrap__(self, obj, context=None):
                obj.metadata = self.metadata
                return obj

        a = nx.array([1.1, -1.1])
        m = MyArray(a, metadata='foo')
        f = ufl.fix(m)
        assert_array_equal(f, nx.array([1, -1]))
        assert_(isinstance(f, MyArray))
        assert_equal(f.metadata, 'foo')

if __name__ == "__main__":
    run_module_suite()
