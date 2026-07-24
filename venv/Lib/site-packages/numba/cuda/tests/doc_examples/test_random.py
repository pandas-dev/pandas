# Contents in this file are referenced from the sphinx-generated docs.
# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.cuda.testing import CUDATestCase, skip_on_cudasim


@skip_on_cudasim("cudasim doesn't support cuda import at non-top-level")
class TestRandom(CUDATestCase):
    def test_ex_3d_grid(self):
        # magictoken.ex_3d_grid.begin
        from numba import cuda
        from numba.cuda.random import (create_xoroshiro128p_states,
                                       xoroshiro128p_uniform_float32)
        import numpy as np

        @cuda.jit
        def random_3d(arr, rng_states):
            # Per-dimension thread indices and strides
            startx, starty, startz = cuda.grid(3)
            stridex, stridey, stridez = cuda.gridsize(3)

            # Linearized thread index
            tid = (startz * stridey * stridex) + (starty * stridex) + startx

            # Use strided loops over the array to assign a random value to each entry
            for i in range(startz, arr.shape[0], stridez):
                for j in range(starty, arr.shape[1], stridey):
                    for k in range(startx, arr.shape[2], stridex):
                        arr[i, j, k] = xoroshiro128p_uniform_float32(rng_states, tid)

        # Array dimensions
        X, Y, Z = 701, 900, 719

        # Block and grid dimensions
        bx, by, bz = 8, 8, 8
        gx, gy, gz = 16, 16, 16

        # Total number of threads
        nthreads = bx * by * bz * gx * gy * gz

        # Initialize a state for each thread
        rng_states = create_xoroshiro128p_states(nthreads, seed=1)

        # Generate random numbers
        arr = cuda.device_array((X, Y, Z), dtype=np.float32)
        random_3d[(gx, gy, gz), (bx, by, bz)](arr, rng_states)
        # magictoken.ex_3d_grid.end

        # Some basic tests of the randomly-generated numbers
        host_arr = arr.copy_to_host()
        self.assertGreater(np.mean(host_arr), 0.49)
        self.assertLess(np.mean(host_arr), 0.51)
        self.assertTrue(np.all(host_arr <= 1.0))
        self.assertTrue(np.all(host_arr >= 0.0))


if __name__ == '__main__':
    unittest.main()
