from __future__ import division, absolute_import, print_function

from numpy.testing import *
from numpy import logspace, linspace, dtype, array

class TestLogspace(TestCase):

    def test_basic(self):
        y = logspace(0, 6)
        assert_(len(y) == 50)
        y = logspace(0, 6, num=100)
        assert_(y[-1] == 10 ** 6)
        y = logspace(0, 6, endpoint=0)
        assert_(y[-1] < 10 ** 6)
        y = logspace(0, 6, num=7)
        assert_array_equal(y, [1, 10, 100, 1e3, 1e4, 1e5, 1e6])

    def test_dtype(self):
        y = logspace(0, 6, dtype='float32')
        assert_equal(y.dtype, dtype('float32'))
        y = logspace(0, 6, dtype='float64')
        assert_equal(y.dtype, dtype('float64'))
        y = logspace(0, 6, dtype='int32')
        assert_equal(y.dtype, dtype('int32'))


class TestLinspace(TestCase):

    def test_basic(self):
        y = linspace(0, 10)
        assert_(len(y) == 50)
        y = linspace(2, 10, num=100)
        assert_(y[-1] == 10)
        y = linspace(2, 10, endpoint=0)
        assert_(y[-1] < 10)

    def test_corner(self):
        y = list(linspace(0, 1, 1))
        assert_(y == [0.0], y)
        y = list(linspace(0, 1, 2.5))
        assert_(y == [0.0, 1.0])

    def test_type(self):
        t1 = linspace(0, 1, 0).dtype
        t2 = linspace(0, 1, 1).dtype
        t3 = linspace(0, 1, 2).dtype
        assert_equal(t1, t2)
        assert_equal(t2, t3)

    def test_dtype(self):
        y = linspace(0, 6, dtype='float32')
        assert_equal(y.dtype, dtype('float32'))
        y = linspace(0, 6, dtype='float64')
        assert_equal(y.dtype, dtype('float64'))
        y = linspace(0, 6, dtype='int32')
        assert_equal(y.dtype, dtype('int32'))

    def test_array_scalar(self):
        lim1 = array([-120, 100], dtype="int8")
        lim2 = array([120, -100], dtype="int8")
        lim3 = array([1200, 1000], dtype="uint16")
        t1 = linspace(lim1[0], lim1[1], 5)
        t2 = linspace(lim2[0], lim2[1], 5)
        t3 = linspace(lim3[0], lim3[1], 5)
        t4 = linspace(-120.0, 100.0, 5)
        t5 = linspace(120.0, -100.0, 5)
        t6 = linspace(1200.0, 1000.0, 5)
        assert_equal(t1, t4)
        assert_equal(t2, t5)
        assert_equal(t3, t6)

    def test_complex(self):
        lim1 = linspace(1 + 2j, 3 + 4j, 5)
        t1 = array([ 1.0+2.j ,  1.5+2.5j,  2.0+3.j ,  2.5+3.5j,  3.0+4.j])
        lim2 = linspace(1j, 10, 5)
        t2 = array([  0.0+1.j  ,   2.5+0.75j,   5.0+0.5j ,   7.5+0.25j,  10.0+0.j])
        assert_equal(lim1, t1)
        assert_equal(lim2, t2)
        
    def test_physical_quantities(self):
        class PhysicalQuantity(float):
            def __new__(cls, value):
                return float.__new__(cls, value)
                
            def __add__(self, x):
                assert_(isinstance(x, PhysicalQuantity))
                return PhysicalQuantity(float(x) + float(self))
            __radd__ = __add__
                
            def __sub__(self, x):
                assert_(isinstance(x, PhysicalQuantity))
                return PhysicalQuantity(float(self) - float(x))
                
            def __rsub__(self, x):
                assert_(isinstance(x, PhysicalQuantity))
                return PhysicalQuantity(float(x) - float(self))
                
            def __mul__(self, x):
                return PhysicalQuantity(float(x) * float(self))
            __rmul__ = __mul__
                
            def __div__(self, x):
                return PhysicalQuantity(float(self) / float(x))
                
            def __rdiv__(self, x):
                return PhysicalQuantity(float(x) / float(self))

        
        a = PhysicalQuantity(0.0)
        b = PhysicalQuantity(1.0)
        assert_equal(linspace(a, b), linspace(0.0, 1.0))