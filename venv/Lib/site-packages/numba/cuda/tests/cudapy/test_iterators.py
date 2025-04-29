from numba import cuda
from numba.cuda.testing import unittest, CUDATestCase

import numpy as np


class TestIterators(CUDATestCase):

    def test_enumerate(self):
        @cuda.jit
        def enumerator(x, error):
            count = 0

            for i, v in enumerate(x):
                if count != i:
                    error[0] = 1
                if v != x[i]:
                    error[0] = 2

                count += 1

            if count != len(x):
                error[0] = 3

        x = np.asarray((10, 9, 8, 7, 6))
        error = np.zeros(1, dtype=np.int32)

        enumerator[1, 1](x, error)
        self.assertEqual(error[0], 0)

    def _test_twoarg_function(self, f):
        x = np.asarray((10, 9, 8, 7, 6))
        y = np.asarray((1, 2, 3, 4, 5))
        error = np.zeros(1, dtype=np.int32)

        f[1, 1](x, y, error)
        self.assertEqual(error[0], 0)

    def test_zip(self):
        @cuda.jit
        def zipper(x, y, error):
            i = 0

            for xv, yv in zip(x, y):
                if xv != x[i]:
                    error[0] = 1
                if yv != y[i]:
                    error[0] = 2

                i += 1

            if i != len(x):
                error[0] = 3

        self._test_twoarg_function(zipper)

    def test_enumerate_zip(self):
        @cuda.jit
        def enumerator_zipper(x, y, error):
            count = 0

            for i, (xv, yv) in enumerate(zip(x, y)):
                if i != count:
                    error[0] = 1
                if xv != x[i]:
                    error[0] = 2
                if yv != y[i]:
                    error[0] = 3

                count += 1

            if count != len(x):
                error[0] = 4

        self._test_twoarg_function(enumerator_zipper)

    def test_zip_enumerate(self):
        @cuda.jit
        def zipper_enumerator(x, y, error):
            count = 0

            for (i, xv), yv in zip(enumerate(x), y):
                if i != count:
                    error[0] = 1
                if xv != x[i]:
                    error[0] = 2
                if yv != y[i]:
                    error[0] = 3

                count += 1

            if count != len(x):
                error[0] = 4

        self._test_twoarg_function(zipper_enumerator)


if __name__ == '__main__':
    unittest.main()
