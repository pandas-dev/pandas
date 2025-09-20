import platform
import numpy as np
from numba import types
import unittest
from numba import njit
from numba.core import config
from numba.tests.support import TestCase

_DEBUG = False
if _DEBUG:
    from llvmlite import binding as llvm
    # Prints debug info from the LLVMs vectorizer
    llvm.set_option("", "--debug-only=loop-vectorize")


_skylake_env = {
    "NUMBA_CPU_NAME": "skylake-avx512",
    "NUMBA_CPU_FEATURES": "",
}


@unittest.skipIf(platform.machine() != 'x86_64', 'x86_64 only test')
class TestVectorization(TestCase):
    """
    Tests to assert that code which should vectorize does indeed vectorize
    """
    def gen_ir(self, func, args_tuple, fastmath=False):
        self.assertEqual(config.CPU_NAME, "skylake-avx512")
        self.assertEqual(config.CPU_FEATURES, "")

        jitted = njit(args_tuple, fastmath=fastmath)(func)
        return jitted.inspect_llvm(args_tuple)

    @TestCase.run_test_in_subprocess(envvars=_skylake_env)
    def test_nditer_loop(self):
        # see https://github.com/numba/numba/issues/5033
        def do_sum(x):
            acc = 0
            for v in np.nditer(x):
                acc += v.item()
            return acc

        llvm_ir = self.gen_ir(do_sum, (types.float64[::1],), fastmath=True)
        self.assertIn("vector.body", llvm_ir)
        self.assertIn("llvm.loop.isvectorized", llvm_ir)

    # SLP is off by default due to miscompilations, see #8705. Put this into a
    # subprocess to isolate any potential issues.
    @TestCase.run_test_in_subprocess(
        envvars={'NUMBA_SLP_VECTORIZE': '1', **_skylake_env},
    )
    def test_slp(self):
        # Sample translated from:
        # https://www.llvm.org/docs/Vectorizers.html#the-slp-vectorizer

        def foo(a1, a2, b1, b2, A):
            A[0] = a1 * (a1 + b1)
            A[1] = a2 * (a2 + b2)
            A[2] = a1 * (a1 + b1)
            A[3] = a2 * (a2 + b2)

        ty = types.float64
        llvm_ir = self.gen_ir(foo, ((ty,) * 4 + (ty[::1],)), fastmath=True)
        self.assertIn("2 x double", llvm_ir)

    @TestCase.run_test_in_subprocess(envvars=_skylake_env)
    def test_instcombine_effect(self):
        # Without instcombine running ahead of refprune, the IR has refops that
        # are trivially prunable (same BB) but the arguments are obfuscated
        # through aliases etc. The follow case triggers this situation as the
        # typed.List has a structproxy call for computing `len` and getting the
        # base pointer for use in iteration.

        def sum_sqrt_list(lst):
            acc = 0.0
            for item in lst:
                acc += np.sqrt(item)
            return acc

        llvm_ir = self.gen_ir(sum_sqrt_list, (types.ListType(types.float64),),
                              fastmath=True)
        self.assertIn("vector.body", llvm_ir)
        self.assertIn("llvm.loop.isvectorized", llvm_ir)


if __name__ == '__main__':
    unittest.main()
