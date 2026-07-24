from numba import jit
import unittest
import numpy as np
import copy
from numba.tests.support import MemoryLeakMixin


class TestChainedAssign(MemoryLeakMixin, unittest.TestCase):
    def test_chain1(self):
        from numba.tests.chained_assign_usecases import chain1
        args = [
            [np.arange(2)],
            [np.arange(4, dtype=np.double)],
        ]
        self._test_template(chain1, args)

    def test_chain2(self):
        from numba.tests.chained_assign_usecases import chain2
        args = [
            [3],
            [3.0],
        ]
        self._test_template(chain2, args)

    def test_unpack1(self):
        from numba.tests.chained_assign_usecases import unpack1
        args = [
            [1, 3.0],
            [1.0, 3],
        ]
        self._test_template(unpack1, args)

    def test_unpack2(self):
        from numba.tests.chained_assign_usecases import unpack2
        args = [
            [np.array([2]), np.array([4.0])],
            [np.array([2.0]), np.array([4])],
        ]
        self._test_template(unpack2, args)

    def test_chain3(self):
        from numba.tests.chained_assign_usecases import chain3
        args = [
            [np.array([0]), np.array([1.5])],
            [np.array([0.5]), np.array([1])],
        ]
        self._test_template(chain3, args)

    def test_unpack3(self):
        from numba.tests.chained_assign_usecases import unpack3
        args = [
            [np.array([1])],
            [np.array([1.0])],
        ]
        self._test_template(unpack3, args)

    def test_unpack4(self):
        from numba.tests.chained_assign_usecases import unpack4
        args = [
            [np.array([1])],
            [np.array([1.0])],
        ]
        self._test_template(unpack4, args)

    def test_unpack5(self):
        from numba.tests.chained_assign_usecases import unpack5
        args = [
            [np.array([2])],
            [np.array([2.0])],
        ]
        self._test_template(unpack5, args)

    def test_unpack6(self):
        from numba.tests.chained_assign_usecases import unpack6
        args1 = 3.0, 2
        args2 = 3.0, 2.0
        self._test_template(unpack6, [args1, args2])

    def _test_template(self, pyfunc, argcases):
        cfunc = jit(pyfunc)
        for args in argcases:
            a1 = copy.deepcopy(args)
            a2 = copy.deepcopy(args)
            np.testing.assert_allclose(pyfunc(*a1), cfunc(*a2))


if __name__ == '__main__':
    unittest.main()

