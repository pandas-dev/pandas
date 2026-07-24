import numpy as np

import numba
from numba.tests.support import TestCase


class Issue455(object):
    """
    Test code from issue 455.
    """

    def __init__(self):
        self.f = []

    def create_f(self):
        code = """
        def f(x):
            n = x.shape[0]
            for i in range(n):
                x[i] = 1.
        """
        d = {}
        exec(code.strip(), d)
        self.f.append(numba.jit("void(f8[:])", nopython=True)(d['f']))

    def call_f(self):
        a = np.zeros(10)
        for f in self.f:
            f(a)
        return a


class TestDynFunc(TestCase):

    def test_issue_455(self):
        inst = Issue455()
        inst.create_f()
        a = inst.call_f()
        self.assertPreciseEqual(a, np.ones_like(a))


if __name__ == '__main__':
    unittest.main()
