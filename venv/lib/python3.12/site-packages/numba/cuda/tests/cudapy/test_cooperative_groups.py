from __future__ import print_function

import numpy as np

from numba import config, cuda, int32
from numba.cuda.testing import (unittest, CUDATestCase, skip_on_cudasim,
                                skip_unless_cc_60, skip_if_cudadevrt_missing,
                                skip_if_mvc_enabled)


@cuda.jit
def this_grid(A):
    cuda.cg.this_grid()
    A[0] = 1.0


@cuda.jit
def sync_group(A):
    g = cuda.cg.this_grid()
    g.sync()
    A[0] = 1.0


@cuda.jit
def no_sync(A):
    A[0] = cuda.grid(1)


def sequential_rows(M):
    # The grid writes rows one at a time. Each thread reads an element from
    # the previous row written by its "opposite" thread.
    #
    # A failure to sync the grid at each row would result in an incorrect
    # result as some threads could run ahead of threads in other blocks, or
    # fail to see the update to the previous row from their opposite thread.

    col = cuda.grid(1)
    g = cuda.cg.this_grid()

    rows = M.shape[0]
    cols = M.shape[1]

    for row in range(1, rows):
        opposite = cols - col - 1
        M[row, col] = M[row - 1, opposite] + 1
        g.sync()


@skip_if_cudadevrt_missing
@skip_if_mvc_enabled('CG not supported with MVC')
class TestCudaCooperativeGroups(CUDATestCase):
    @skip_unless_cc_60
    def test_this_grid(self):
        A = np.full(1, fill_value=np.nan)
        this_grid[1, 1](A)

        # Ensure the kernel executed beyond the call to cuda.this_grid()
        self.assertFalse(np.isnan(A[0]), 'Value was not set')

    @skip_unless_cc_60
    @skip_on_cudasim("Simulator doesn't differentiate between normal and "
                     "cooperative kernels")
    def test_this_grid_is_cooperative(self):
        A = np.full(1, fill_value=np.nan)
        this_grid[1, 1](A)

        # this_grid should have been determined to be cooperative
        for key, overload in this_grid.overloads.items():
            self.assertTrue(overload.cooperative)

    @skip_unless_cc_60
    def test_sync_group(self):
        A = np.full(1, fill_value=np.nan)
        sync_group[1, 1](A)

        # Ensure the kernel executed beyond the call to cuda.sync_group()
        self.assertFalse(np.isnan(A[0]), 'Value was not set')

    @skip_unless_cc_60
    @skip_on_cudasim("Simulator doesn't differentiate between normal and "
                     "cooperative kernels")
    def test_sync_group_is_cooperative(self):
        A = np.full(1, fill_value=np.nan)
        sync_group[1, 1](A)
        # sync_group should have been determined to be cooperative
        for key, overload in sync_group.overloads.items():
            self.assertTrue(overload.cooperative)

    @skip_on_cudasim("Simulator does not implement linking")
    def test_false_cooperative_doesnt_link_cudadevrt(self):
        """
        We should only mark a kernel as cooperative and link cudadevrt if the
        kernel uses grid sync. Here we ensure that one that doesn't use grid
        synsync isn't marked as such.
        """
        A = np.full(1, fill_value=np.nan)
        no_sync[1, 1](A)

        for key, overload in no_sync.overloads.items():
            self.assertFalse(overload.cooperative)
            for link in overload._codelibrary._linking_files:
                self.assertNotIn('cudadevrt', link)

    @skip_unless_cc_60
    def test_sync_at_matrix_row(self):
        if config.ENABLE_CUDASIM:
            # Use a small matrix to compute using a single block in a
            # reasonable amount of time
            shape = (32, 32)
        else:
            shape = (1024, 1024)
        A = np.zeros(shape, dtype=np.int32)
        blockdim = 32
        griddim = A.shape[1] // blockdim

        sig = (int32[:,::1],)
        c_sequential_rows = cuda.jit(sig)(sequential_rows)

        overload = c_sequential_rows.overloads[sig]
        mb = overload.max_cooperative_grid_blocks(blockdim)
        if griddim > mb:
            unittest.skip("GPU cannot support enough cooperative grid blocks")

        c_sequential_rows[griddim, blockdim](A)

        reference = np.tile(np.arange(shape[0]), (shape[1], 1)).T
        np.testing.assert_equal(A, reference)

    @skip_unless_cc_60
    def test_max_cooperative_grid_blocks(self):
        # The maximum number of blocks will vary based on the device so we
        # can't test for an expected value, but we can check that the function
        # doesn't error, and that varying the number of dimensions of the block
        # whilst keeping the total number of threads constant doesn't change
        # the maximum to validate some of the logic.
        sig = (int32[:,::1],)
        c_sequential_rows = cuda.jit(sig)(sequential_rows)
        overload = c_sequential_rows.overloads[sig]
        blocks1d = overload.max_cooperative_grid_blocks(256)
        blocks2d = overload.max_cooperative_grid_blocks((16, 16))
        blocks3d = overload.max_cooperative_grid_blocks((16, 4, 4))
        self.assertEqual(blocks1d, blocks2d)
        self.assertEqual(blocks1d, blocks3d)


if __name__ == '__main__':
    unittest.main()
