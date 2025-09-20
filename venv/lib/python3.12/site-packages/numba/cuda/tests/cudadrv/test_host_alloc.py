import numpy as np
from numba.cuda.cudadrv import driver
from numba import cuda
from numba.cuda.testing import unittest, ContextResettingTestCase


class TestHostAlloc(ContextResettingTestCase):
    def test_host_alloc_driver(self):
        n = 32
        mem = cuda.current_context().memhostalloc(n, mapped=True)

        dtype = np.dtype(np.uint8)
        ary = np.ndarray(shape=n // dtype.itemsize, dtype=dtype,
                         buffer=mem)

        magic = 0xab
        driver.device_memset(mem, magic, n)

        self.assertTrue(np.all(ary == magic))

        ary.fill(n)

        recv = np.empty_like(ary)

        driver.device_to_host(recv, mem, ary.size)

        self.assertTrue(np.all(ary == recv))
        self.assertTrue(np.all(recv == n))

    def test_host_alloc_pinned(self):
        ary = cuda.pinned_array(10, dtype=np.uint32)
        ary.fill(123)
        self.assertTrue(all(ary == 123))
        devary = cuda.to_device(ary)
        driver.device_memset(devary, 0, driver.device_memory_size(devary))
        self.assertTrue(all(ary == 123))
        devary.copy_to_host(ary)
        self.assertTrue(all(ary == 0))

    def test_host_alloc_mapped(self):
        ary = cuda.mapped_array(10, dtype=np.uint32)
        ary.fill(123)
        self.assertTrue(all(ary == 123))
        driver.device_memset(ary, 0, driver.device_memory_size(ary))
        self.assertTrue(all(ary == 0))
        self.assertTrue(sum(ary != 0) == 0)

    def test_host_operators(self):
        for ary in [cuda.mapped_array(10, dtype=np.uint32),
                    cuda.pinned_array(10, dtype=np.uint32)]:
            ary[:] = range(10)
            self.assertTrue(sum(ary + 1) == 55)
            self.assertTrue(sum((ary + 1) * 2 - 1) == 100)
            self.assertTrue(sum(ary < 5) == 5)
            self.assertTrue(sum(ary <= 5) == 6)
            self.assertTrue(sum(ary > 6) == 3)
            self.assertTrue(sum(ary >= 6) == 4)
            self.assertTrue(sum(ary ** 2) == 285)
            self.assertTrue(sum(ary // 2) == 20)
            self.assertTrue(sum(ary / 2.0) == 22.5)
            self.assertTrue(sum(ary % 2) == 5)


if __name__ == '__main__':
    unittest.main()
