import math

import numpy as np

from numba import cuda
from numba.cuda.testing import unittest
from numba.cuda.testing import skip_on_cudasim, CUDATestCase

from numba.cuda.random import \
    xoroshiro128p_uniform_float32, xoroshiro128p_normal_float32, \
    xoroshiro128p_uniform_float64, xoroshiro128p_normal_float64


# Distributions
UNIFORM = 1
NORMAL = 2


@cuda.jit
def rng_kernel_float32(states, out, count, distribution):
    thread_id = cuda.grid(1)

    for i in range(count):
        idx = thread_id * count + i

        if distribution == UNIFORM:
            out[idx] = xoroshiro128p_uniform_float32(states, thread_id)
        elif distribution == NORMAL:
            out[idx] = xoroshiro128p_normal_float32(states, thread_id)


@cuda.jit
def rng_kernel_float64(states, out, count, distribution):
    thread_id = cuda.grid(1)

    for i in range(count):
        idx = thread_id * count + i

        if distribution == UNIFORM:
            out[idx] = xoroshiro128p_uniform_float64(states, thread_id)
        elif distribution == NORMAL:
            out[idx] = xoroshiro128p_normal_float64(states, thread_id)


class TestCudaRandomXoroshiro128p(CUDATestCase):
    def test_create(self):
        states = cuda.random.create_xoroshiro128p_states(10, seed=1)
        s = states.copy_to_host()
        self.assertEqual(len(np.unique(s)), 10)

    def test_create_subsequence_start(self):
        states = cuda.random.create_xoroshiro128p_states(10, seed=1)
        s1 = states.copy_to_host()

        states = cuda.random.create_xoroshiro128p_states(10, seed=1,
                                                         subsequence_start=3)
        s2 = states.copy_to_host()

        # Starting seeds should match up with offset of 3
        np.testing.assert_array_equal(s1[3:], s2[:-3])

    def test_create_stream(self):
        stream = cuda.stream()
        states = cuda.random.create_xoroshiro128p_states(10, seed=1,
                                                         stream=stream)
        s = states.copy_to_host()
        self.assertEqual(len(np.unique(s)), 10)

    def check_uniform(self, kernel_func, dtype):
        states = cuda.random.create_xoroshiro128p_states(32 * 2, seed=1)
        out = np.zeros(2 * 32 * 32, dtype=np.float32)

        kernel_func[2, 32](states, out, 32, UNIFORM)
        self.assertAlmostEqual(out.min(), 0.0, delta=1e-3)
        self.assertAlmostEqual(out.max(), 1.0, delta=1e-3)
        self.assertAlmostEqual(out.mean(), 0.5, delta=1.5e-2)
        self.assertAlmostEqual(out.std(), 1.0 / (2 * math.sqrt(3)), delta=6e-3)

    def test_uniform_float32(self):
        self.check_uniform(rng_kernel_float32, np.float32)

    @skip_on_cudasim('skip test for speed under cudasim')
    def test_uniform_float64(self):
        self.check_uniform(rng_kernel_float64, np.float64)

    def check_normal(self, kernel_func, dtype):
        states = cuda.random.create_xoroshiro128p_states(32 * 2, seed=1)
        out = np.zeros(2 * 32 * 32, dtype=dtype)

        kernel_func[2, 32](states, out, 32, NORMAL)

        self.assertAlmostEqual(out.mean(), 0.0, delta=4e-3)
        self.assertAlmostEqual(out.std(), 1.0, delta=2e-3)

    def test_normal_float32(self):
        self.check_normal(rng_kernel_float32, np.float32)

    @skip_on_cudasim('skip test for speed under cudasim')
    def test_normal_float64(self):
        self.check_normal(rng_kernel_float64, np.float64)


if __name__ == '__main__':
    unittest.main()
