from numba.tests.support import TestCase, numpy_support
from numba import njit, types
from numba.typed import List, Dict
import numpy as np


class TestConditionsAsPredicates(TestCase):

    def test_scalars(self):
        # checks that scalar types can be used as predicates
        dts = [np.int8, np.uint16, np.int64, np.float32, np.float64,
               np.complex128, int, float, complex, str, bool]
        for dt in dts:
            for c in 1, 0:
                x = dt(c)

                @njit
                def foo():
                    if x:
                        return 10
                    else:
                        return 20
                self.assertEqual(foo(), foo.py_func())
                self.assertEqual(foo(), 10 if c == 1 or dt is str else 20)

        # empty string
        @njit
        def foo(x):
            if x:
                return 10
            else:
                return 20
        s = ""
        self.assertEqual(foo(s), foo.py_func(s))

    def test_typed_list(self):
        @njit
        def foo(x):
            if x:
                return 10
            else:
                return 20

        # empty list
        z = List.empty_list(types.int64)
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 20)

        # non-empty list
        z.append(1)
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 10)

    def test_reflected_list(self):
        # non-empty
        @njit
        def foo(x):
            if x:
                return 10
            else:
                return 20

        z = [1]
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 10)

        # non-empty local
        @njit
        def foo():
            y = [1, 2]
            if y:
                return 10
            else:
                return 20

        self.assertEqual(foo(), foo.py_func())
        self.assertEqual(foo.py_func(), 10)

        # empty local
        @njit
        def foo():
            y = [1, 2]
            y.pop()
            y.pop()
            assert len(y) == 0
            if y:
                return 10
            else:
                return 20

        self.assertEqual(foo(), foo.py_func())
        self.assertEqual(foo.py_func(), 20)

    def test_reflected_set(self):
        # non-empty
        @njit
        def foo(x):
            if x:
                return 10
            else:
                return 20

        z = {1}
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 10)

        # non-empty local
        @njit
        def foo():
            y = {1, 2}
            if y:
                return 10
            else:
                return 20

        self.assertEqual(foo(), foo.py_func())
        self.assertEqual(foo.py_func(), 10)

        # empty local
        @njit
        def foo():
            y = {1, 2}
            y.pop()
            y.pop()
            assert len(y) == 0
            if y:
                return 10
            else:
                return 20

        self.assertEqual(foo(), foo.py_func())
        self.assertEqual(foo.py_func(), 20)

    def test_typed_dict(self):
        @njit
        def foo(x):
            if x:
                return 10
            else:
                return 20

        # empty
        z = Dict.empty(types.int64, types.int64)
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 20)

        # non-empty
        z[2] = 3
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 10)

    def test_arrays(self):
        @njit
        def foo(x):
            if x:
                return 10
            else:
                return 20

        # non-empty 0d, True
        z = np.array(1)
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 10)

        # non-empty 0d, False
        z = np.array(0)
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 20)
        # non-empty nd True
        z = np.array([[[1]]])
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 10)

        # non-empty nd False
        z = np.array([[[0]]])
        self.assertEqual(foo(z), foo.py_func(z))
        self.assertEqual(foo.py_func(z), 20)

        # various problems:

        # empty, NumPy warns or raises if NumPy >= 2.2
        z = np.empty(0)
        if numpy_support.numpy_version >= (2, 2):
            with self.assertRaises(ValueError) as raises:
                foo(z)
            msg = ("The truth value of an empty array is ambiguous."
                   " Use `array.size > 0` to check that an array is not empty.")
            self.assertIn(msg, str(raises.exception))
        else:
            self.assertEqual(foo(z), foo.py_func(z))
            self.assertEqual(foo.py_func(z), 20)

        # nd, NumPy raises
        z = np.array([1, 2])
        with self.assertRaises(ValueError) as raises:
            foo(z)

        msg = ("The truth value of an array with more than one element "
               "is ambiguous. Use a.any() or a.all()")
        self.assertIn(msg, str(raises.exception))
