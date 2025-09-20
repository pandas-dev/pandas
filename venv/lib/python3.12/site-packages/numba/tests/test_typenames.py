import numpy as np

from numba.core import types
import unittest


class TestTypeNames(unittest.TestCase):
    def test_numpy_integers(self):
        expect = getattr(types, "int%d" % (np.dtype("int").itemsize * 8))
        self.assertEqual(types.int_, expect)

        expect = getattr(types, "uint%d" % (np.dtype("uint").itemsize * 8))
        self.assertEqual(types.uint, expect)


if __name__ == '__main__':
    unittest.main()
