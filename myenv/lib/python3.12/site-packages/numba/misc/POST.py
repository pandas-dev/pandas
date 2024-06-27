""" Numba's POWER ON SELF TEST script. Used by CI to check:
0. That Numba imports ok!
1. That Numba can find an appropriate number of its own tests to run.
2. That Numba can manage to correctly compile and execute at least one thing.
"""
from numba.tests import test_runtests
from numba import njit


def _check_runtests():
    test_inst = test_runtests.TestCase()
    test_inst.test_default() # will raise an exception if there is a problem


def _check_cpu_compilation():
    @njit
    def foo(x):
        return x + 1

    result = foo(1)

    if result != 2:
        msg = ("Unexpected result from trial compilation. "
               f"Expected: 2, Got: {result}.")
        raise AssertionError(msg)


def check():
    _check_runtests()
    _check_cpu_compilation()


if __name__ == "__main__":
    check()
