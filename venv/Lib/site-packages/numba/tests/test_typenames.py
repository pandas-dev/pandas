from numba.core import types
import unittest

import numpy as np


class TestTypeNames(unittest.TestCase):
    def test_numpy_integers(self):
        self.assertEqual(types.int_.bitwidth, np.dtype(np.intp).itemsize * 8)
        self.assertEqual(types.uint.bitwidth, np.dtype(np.uintp).itemsize * 8)


if __name__ == '__main__':
    unittest.main()
