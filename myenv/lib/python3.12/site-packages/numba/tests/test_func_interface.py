import unittest
from numba import jit


class TestFuncInterface(unittest.TestCase):
    def test_jit_function_docstring(self):

        def add(x, y):
            '''Return sum of two numbers'''
            return x + y

        c_add = jit(add)
        self.assertEqual(c_add.__doc__, 'Return sum of two numbers')

    def test_jit_function_name(self):

        def add(x, y):
            return x + y

        c_add = jit(add)
        self.assertEqual(c_add.__name__, 'add')

    def test_jit_function_module(self):

        def add(x, y):
            return x + y

        c_add = jit(add)
        # Expected answer depends on how you run this test.
        # Compare to python function instead.
        self.assertEqual(c_add.__module__, add.__module__)

    def test_jit_function_code_object(self):
        def add(x, y):
            return x + y

        c_add = jit(add)
        self.assertEqual(c_add.__code__, add.__code__)
        self.assertEqual(c_add.func_code, add.__code__)


if __name__ == '__main__':
    unittest.main()
