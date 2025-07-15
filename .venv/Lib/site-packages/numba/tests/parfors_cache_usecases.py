import sys

import numpy as np

from numba import njit
from numba.tests.support import TestCase


@njit(parallel=True, cache=True)
def arrayexprs_case(arr):
    return arr / arr.sum()


@njit(parallel=True, cache=True)
def prange_case(arr):
    out = np.zeros_like(arr)
    c = 1 / arr.sum()
    for i in range(arr.size):
        out[i] = arr[i] * c
    return out


@njit(cache=True)
def caller_case(arr):
    return prange_case(arrayexprs_case(arr))


class _TestModule(TestCase):
    """
    Tests for functionality of this module's functions.
    Note this does not define any "test_*" method, instead check_module()
    should be called by hand.
    """
    def check_module(self, mod):
        total_cache_hits = 0
        for fn in [mod.arrayexprs_case, mod.prange_case, mod.caller_case]:
            arr = np.ones(20)
            np.testing.assert_allclose(
                fn(arr), fn.py_func(arr),
            )
            # Accumulate cache hits
            total_cache_hits += len(fn.stats.cache_hits)
        self.assertGreater(
            total_cache_hits, 0,
            msg="At least one dispatcher has used the cache",
        )

    def run_module(self, mod):
        # This just executes the module's functionality without asserting
        # anything about the cache, it's used in tests that ensure that
        # properties such as thread count aren't baked in to the cached object.
        for fn in [mod.arrayexprs_case, mod.prange_case, mod.caller_case]:
            arr = np.ones(20)
            np.testing.assert_allclose(
                fn(arr), fn.py_func(arr),
            )


def self_test():
    mod = sys.modules[__name__]
    _TestModule().check_module(mod)


def self_run():
    mod = sys.modules[__name__]
    _TestModule().run_module(mod)
