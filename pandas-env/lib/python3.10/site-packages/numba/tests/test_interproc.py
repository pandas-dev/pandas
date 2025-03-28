import gc

from numba import jit, int32
import unittest


def foo(a, b):
    return a + b


def bar(a, b):
    return cfoo(a, b) + b

@jit
def inner(x, y):
    return x + y

@jit(nopython=True)
def outer(x, y):
    return inner(x, y)


class TestInterProc(unittest.TestCase):

    def test_bar_call_foo(self):
        global cfoo
        cfoo = jit((int32, int32), nopython=True)(foo)
        cbar = jit((int32, int32), nopython=True)(bar)
        self.assertEqual(cbar(1, 2), 1 + 2 + 2)

    def test_bar_call_foo_compiled_twice(self):
        # When a function is compiled twice, then called from another
        # compiled function, check that the right target is called.
        # (otherwise, LLVM would assert out or crash)
        global cfoo
        for i in range(2):
            cfoo = jit((int32, int32), nopython=True)(foo)
            gc.collect()
        cbar = jit((int32, int32), nopython=True)(bar)
        self.assertEqual(cbar(1, 2), 1 + 2 + 2)

    def test_callsite_compilation(self):
        self.assertEqual(outer(1, 2), 1 + 2)


if __name__ == '__main__':
    unittest.main()
