import unittest

from numba.cuda.testing import (CUDATestCase, skip_if_cudadevrt_missing,
                                skip_on_cudasim, skip_unless_cc_60,
                                skip_if_mvc_enabled)
from numba.tests.support import captured_stdout


@skip_if_cudadevrt_missing
@skip_unless_cc_60
@skip_if_mvc_enabled('CG not supported with MVC')
@skip_on_cudasim("cudasim doesn't support cuda import at non-top-level")
class TestLaplace(CUDATestCase):
    """
    Test simple vector addition
    """

    def setUp(self):
        # Prevent output from this test showing up when running the test suite
        self._captured_stdout = captured_stdout()
        self._captured_stdout.__enter__()
        super().setUp()

    def tearDown(self):
        # No exception type, value, or traceback
        self._captured_stdout.__exit__(None, None, None)
        super().tearDown()

    def test_ex_laplace(self):

        # set True to regenerate the figures that
        # accompany this example
        plot = False

        # ex_laplace.import.begin
        import numpy as np
        from numba import cuda
        # ex_laplace.import.end

        # ex_laplace.allocate.begin
        # Use an odd problem size.
        # This is so there can be an element truly in the "middle" for symmetry.
        size = 1001
        data = np.zeros(size)

        # Middle element is made very hot
        data[500] = 10000
        buf_0 = cuda.to_device(data)

        # This extra array is used for synchronization purposes
        buf_1 = cuda.device_array_like(buf_0)

        niter = 10000
        # ex_laplace.allocate.end

        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(16 * 0.66, 9 * 0.66))
            plt.plot(
                np.arange(len(buf_0)),
                buf_0.copy_to_host(),
                lw=3,
                marker="*",
                color='black'
            )

            plt.title('Initial State', fontsize=24)
            plt.xlabel('Position', fontsize=24)
            plt.ylabel('Temperature', fontsize=24)

            ax.set_xticks(ax.get_xticks(), fontsize=16)
            ax.set_yticks(ax.get_yticks(), fontsize=16)
            plt.xlim(0, len(data))
            plt.ylim(0, 10001)
            plt.savefig('laplace_initial.svg')

        # ex_laplace.kernel.begin
        @cuda.jit
        def solve_heat_equation(buf_0, buf_1, timesteps, k):
            i = cuda.grid(1)

            # Don't continue if our index is outside the domain
            if i >= len(buf_0):
                return

            # Prepare to do a grid-wide synchronization later
            grid = cuda.cg.this_grid()

            for step in range(timesteps):
                # Select the buffer from the previous timestep
                if (step % 2) == 0:
                    data = buf_0
                    next_data = buf_1
                else:
                    data = buf_1
                    next_data = buf_0

                # Get the current temperature associated with this point
                curr_temp = data[i]

                # Apply formula from finite difference equation
                if i == 0:
                    # Left wall is held at T = 0
                    next_temp = curr_temp + k * (data[i + 1] - (2 * curr_temp))
                elif i == len(data) - 1:
                    # Right wall is held at T = 0
                    next_temp = curr_temp + k * (data[i - 1] - (2 * curr_temp))
                else:
                    # Interior points are a weighted average of their neighbors
                    next_temp = curr_temp + k * (
                        data[i - 1] - (2 * curr_temp) + data[i + 1]
                    )

                # Write new value to the next buffer
                next_data[i] = next_temp

                # Wait for every thread to write before moving on
                grid.sync()
        # ex_laplace.kernel.end

        # ex_laplace.launch.begin
        solve_heat_equation.forall(len(data))(
            buf_0, buf_1, niter, 0.25
        )
        # ex_laplace.launch.end

        results = buf_1.copy_to_host()
        if plot:
            fig, ax = plt.subplots(figsize=(16 * 0.66, 9 * 0.66))
            plt.plot(
                np.arange(len(results)),
                results, lw=3,
                marker="*",
                color='black'
            )
            plt.title(f"T = {niter}", fontsize=24)
            plt.xlabel('Position', fontsize=24)
            plt.ylabel('Temperature', fontsize=24)

            ax.set_xticks(ax.get_xticks(), fontsize=16)
            ax.set_yticks(ax.get_yticks(), fontsize=16)

            plt.ylim(0, max(results))
            plt.xlim(0, len(results))
            plt.savefig('laplace_final.svg')

        # Integral over the domain should be equal to its initial value.
        # Note that this should match the initial value of data[500] above, but
        # we don't assign it to a variable because that would make the example
        # code look a bit oddly verbose.
        np.testing.assert_allclose(results.sum(), 10000)


if __name__ == "__main__":
    unittest.main()
