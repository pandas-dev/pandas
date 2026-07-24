import unittest

from numba.cuda.testing import CUDATestCase, skip_on_cudasim
from numba.tests.support import captured_stdout


@skip_on_cudasim("cudasim doesn't support cuda import at non-top-level")
class TestMonteCarlo(CUDATestCase):
    """
    Test monte-carlo integration
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

    def test_ex_montecarlo(self):
        # ex_montecarlo.import.begin
        import numba
        import numpy as np
        from numba import cuda
        from numba.cuda.random import (
            create_xoroshiro128p_states,
            xoroshiro128p_uniform_float32,
        )
        # ex_montecarlo.import.end

        # ex_montecarlo.define.begin
        # number of samples, higher will lead to a more accurate answer
        nsamps = 1000000
        # ex_montecarlo.define.end

        # ex_montecarlo.kernel.begin
        @cuda.jit
        def mc_integrator_kernel(out, rng_states, lower_lim, upper_lim):
            """
            kernel to draw random samples and evaluate the function to
            be integrated at those sample values
            """
            size = len(out)

            gid = cuda.grid(1)
            if gid < size:
                # draw a sample between 0 and 1 on this thread
                samp = xoroshiro128p_uniform_float32(rng_states, gid)

                # normalize this sample to the limit range
                samp = samp * (upper_lim - lower_lim) + lower_lim

                # evaluate the function to be
                # integrated at the normalized
                # value of the sample
                y = func(samp)
                out[gid] = y
        # ex_montecarlo.kernel.end

        # ex_montecarlo.callfunc.begin
        @cuda.reduce
        def sum_reduce(a, b):
            return a + b

        def mc_integrate(lower_lim, upper_lim, nsamps):
            """
            approximate the definite integral of `func` from
            `lower_lim` to `upper_lim`
            """
            out = cuda.to_device(np.zeros(nsamps, dtype="float32"))
            rng_states = create_xoroshiro128p_states(nsamps, seed=42)

            # jit the function for use in CUDA kernels

            mc_integrator_kernel.forall(nsamps)(
                out, rng_states, lower_lim, upper_lim
            )
            # normalization factor to convert
            # to the average: (b - a)/(N - 1)
            factor = (upper_lim - lower_lim) / (nsamps - 1)

            return sum_reduce(out) * factor
        # ex_montecarlo.callfunc.end

        # ex_montecarlo.launch.begin
        # define a function to integrate
        @numba.jit
        def func(x):
            return 1.0 / x

        mc_integrate(1, 2, nsamps)  # array(0.6929643, dtype=float32)
        mc_integrate(2, 3, nsamps)  # array(0.4054021, dtype=float32)
        # ex_montecarlo.launch.end

        # values computed independently using maple
        np.testing.assert_allclose(
            mc_integrate(1, 2, nsamps), 0.69315, atol=0.001
        )
        np.testing.assert_allclose(
            mc_integrate(2, 3, nsamps), 0.4055, atol=0.001
        )


if __name__ == "__main__":
    unittest.main()
