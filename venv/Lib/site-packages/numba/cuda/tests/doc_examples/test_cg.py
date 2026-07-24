# Contents in this file are referenced from the sphinx-generated docs.
# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.cuda.testing import (CUDATestCase, skip_on_cudasim,
                                skip_if_cudadevrt_missing, skip_unless_cc_60,
                                skip_if_mvc_enabled)


@skip_if_cudadevrt_missing
@skip_unless_cc_60
@skip_if_mvc_enabled('CG not supported with MVC')
@skip_on_cudasim("cudasim doesn't support cuda import at non-top-level")
class TestCooperativeGroups(CUDATestCase):
    def test_ex_grid_sync(self):
        # magictoken.ex_grid_sync_kernel.begin
        from numba import cuda, int32
        import numpy as np

        sig = (int32[:,::1],)

        @cuda.jit(sig)
        def sequential_rows(M):
            col = cuda.grid(1)
            g = cuda.cg.this_grid()

            rows = M.shape[0]
            cols = M.shape[1]

            for row in range(1, rows):
                opposite = cols - col - 1
                # Each row's elements are one greater than the previous row
                M[row, col] = M[row - 1, opposite] + 1
                # Wait until all threads have written their column element,
                # and that the write is visible to all other threads
                g.sync()
        # magictoken.ex_grid_sync_kernel.end

        # magictoken.ex_grid_sync_data.begin
        # Empty input data
        A = np.zeros((1024, 1024), dtype=np.int32)
        # A somewhat arbitrary choice (one warp), but generally smaller block sizes
        # allow more blocks to be launched (noting that other limitations on
        # occupancy apply such as shared memory size)
        blockdim = 32
        griddim = A.shape[1] // blockdim
        # magictoken.ex_grid_sync_data.end

        # Skip this test if the grid size used in the example is too large for
        # a cooperative launch on the current GPU
        mb = sequential_rows.overloads[sig].max_cooperative_grid_blocks(blockdim)
        if mb < griddim:
            self.skipTest('Device does not support a large enough coop grid')

        # magictoken.ex_grid_sync_launch.begin
        # Kernel launch - this is implicitly a cooperative launch
        sequential_rows[griddim, blockdim](A)

        # What do the results look like?
        # print(A)
        #
        # [[   0    0    0 ...    0    0    0]
        #  [   1    1    1 ...    1    1    1]
        #  [   2    2    2 ...    2    2    2]
        #  ...
        #  [1021 1021 1021 ... 1021 1021 1021]
        #  [1022 1022 1022 ... 1022 1022 1022]
        #  [1023 1023 1023 ... 1023 1023 1023]]
        # magictoken.ex_grid_sync_launch.end

        # Sanity check - are the results what we expect?
        reference = np.tile(np.arange(1024), (1024, 1)).T
        np.testing.assert_equal(A, reference)


if __name__ == '__main__':
    unittest.main()
