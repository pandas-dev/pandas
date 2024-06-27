# Contents in this file are referenced from the sphinx-generated docs.
# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.tests.support import captured_stdout, skip_parfors_unsupported
from numba import set_parallel_chunksize
from numba.tests.support import TestCase


@skip_parfors_unsupported
class ChunksizeExamplesTest(TestCase):

    _numba_parallel_test_ = False

    def setUp(self):
        set_parallel_chunksize(0)

    def tearDown(self):
        set_parallel_chunksize(0)

    def test_unbalanced_example(self):
        with captured_stdout():
            # magictoken.ex_unbalanced.begin
            from numba import (njit,
                               prange,
                               )
            import numpy as np

            @njit(parallel=True)
            def func1():
                n = 100
                vals = np.empty(n)
                # The work in each iteration of the following prange
                # loop is proportional to its index.
                for i in prange(n):
                    cur = i + 1
                    for j in range(i):
                        if cur % 2 == 0:
                            cur //= 2
                        else:
                            cur = cur * 3 + 1
                    vals[i] = cur
                return vals

            result = func1()
            # magictoken.ex_unbalanced.end
            self.assertPreciseEqual(result, func1.py_func())

    def test_chunksize_manual(self):
        with captured_stdout():
            # magictoken.ex_chunksize_manual.begin
            from numba import (njit,
                               prange,
                               set_parallel_chunksize,
                               get_parallel_chunksize,
                               )

            @njit(parallel=True)
            def func1(n):
                acc = 0
                print(get_parallel_chunksize()) # Will print 4.
                for i in prange(n):
                    print(get_parallel_chunksize()) # Will print 0.
                    acc += i
                print(get_parallel_chunksize()) # Will print 4.
                return acc

            @njit(parallel=True)
            def func2(n):
                acc = 0
                # This version gets the previous chunksize explicitly.
                old_chunksize = get_parallel_chunksize()
                set_parallel_chunksize(8)
                for i in prange(n):
                    acc += i
                set_parallel_chunksize(old_chunksize)
                return acc

            # This version saves the previous chunksize as returned
            # by set_parallel_chunksize.
            old_chunksize = set_parallel_chunksize(4)
            result1 = func1(12)
            result2 = func2(12)
            result3 = func1(12)
            set_parallel_chunksize(old_chunksize)
            # magictoken.ex_chunksize_manual.end
            self.assertPreciseEqual(result1, func1.py_func(12))
            self.assertPreciseEqual(result2, func2.py_func(12))
            self.assertPreciseEqual(result3, func1.py_func(12))

    def test_chunksize_with(self):
        with captured_stdout():
            # magictoken.ex_chunksize_with.begin
            from numba import njit, prange, parallel_chunksize

            @njit(parallel=True)
            def func1(n):
                acc = 0
                for i in prange(n):
                    acc += i
                return acc

            @njit(parallel=True)
            def func2(n):
                acc = 0
                with parallel_chunksize(8):
                    for i in prange(n):
                        acc += i
                return acc

            with parallel_chunksize(4):
                result1 = func1(12)
                result2 = func2(12)
                result3 = func1(12)
            # magictoken.ex_chunksize_with.end
            self.assertPreciseEqual(result1, func1.py_func(12))
            self.assertPreciseEqual(result2, func2.py_func(12))
            self.assertPreciseEqual(result3, func1.py_func(12))


if __name__ == '__main__':
    unittest.main()
