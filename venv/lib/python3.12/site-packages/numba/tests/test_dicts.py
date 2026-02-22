import numpy as np
from numba import njit, jit
from numba.core.errors import TypingError
import unittest
from numba.tests.support import TestCase


def build_map():
    return {0: 1, 2: 3}

def build_map_from_local_vars():
    # There used to be a crash due to wrong IR generation for STORE_MAP
    x = TestCase
    return {0: x, x: 1}


class DictTestCase(TestCase):

    def check(self, pyfunc):
        cfunc = jit(forceobj=True)(pyfunc)
        self.assertPreciseEqual(pyfunc(), cfunc())

    def test_build_map(self):
        self.check(build_map)

    def test_build_map_from_local_vars(self):
        self.check(build_map_from_local_vars)


class TestCompiledDict(TestCase):
    """Testing `dict()` and `{}` usage that are redirected to
    `numba.typed.Dict`.
    """
    def test_use_dict(self):
        # Test dict()
        @njit
        def foo():
            d = dict()
            d[1] = 2
            return d

        d = foo()
        self.assertEqual(d, {1: 2})

    def test_use_dict_iterable_args(self):
        # Test dict(iterable)
        @njit
        def dict_iterable_1(a, b):
            d = dict(zip(a, b))
            return d

        @njit
        def dict_iterable_2():
            # from python docs
            return dict([('sape', 4139), ('guido', 4127), ('jack', 4098)])

        inps = (
            ([1, 2, 3], [4, 5, 6]),
            (np.arange(4), np.arange(4)),
            ([1, 2, 3], 'abc'),
            ([1, 2, 3, 4], 'abc'),
        )
        for a, b in inps:
            d = dict_iterable_1(a, b)
            self.assertEqual(d, dict(zip(a, b)))

        self.assertEqual(dict_iterable_2(), dict_iterable_2.py_func())

    def test_ctor_iterable_tuple(self):
        @njit
        def ctor():
            return dict(((1, 2), (1, 2)))

        expected = dict({1: 2})
        got = ctor()
        self.assertEqual(expected, got)

    def test_unsupported_dict_usage(self):
        # Test dict(dict())
        from numba.core.typing.dictdecl import _message_dict_support

        @njit
        def ctor1():
            d = dict()
            d[1] = 2
            return dict(d)

        @njit
        def ctor2():
            return dict(((1, 2), (3, 'a')))

        @njit
        def ctor3():
            return dict((('a', 'b', 'c'), ('d', 'e', 'f')))

        @njit
        def ctor4():
            return dict((({}, 1), ({}, 2)))

        _non_iter_args = "Non-iterable args used in dict(iterable)"
        _dict_upd_item_len = "dictionary update sequence element has length 3;"
        _unhashable_type = "Unhashable type"

        inputs = [
            (ctor1, TypingError, _message_dict_support),
            (ctor2, TypingError, _non_iter_args),
            (ctor3, TypingError, _dict_upd_item_len),
            (ctor4, TypingError, _unhashable_type),
        ]

        for func, exc, msg in inputs:
            with self.assertRaises(exc) as raises:
                func()

            self.assertIn(msg, str(raises.exception))

    def test_use_curlybraces(self):
        # Test {} with empty args
        @njit
        def foo():
            d = {}
            d[1] = 2
            return d

        d = foo()
        self.assertEqual(d, {1: 2})

    def test_use_curlybraces_with_init1(self):
        # Test {} with 1 item
        @njit
        def foo():
            return {1: 2}

        d = foo()
        self.assertEqual(d, {1: 2})

    def test_use_curlybraces_with_initmany(self):
        # Test {} with many items
        @njit
        def foo():
            return {1: 2.2, 3: 4.4, 5: 6.6}

        d = foo()
        self.assertEqual(d, {1: 2.2, 3: 4.4, 5: 6.6})

    def test_curlybraces_init_with_coercion(self):
        # Type coercion at dict init is tested
        @njit
        def foo():
            return {1: 2.2, 3: 4, 5: 6}

        self.assertEqual(foo(), foo.py_func())

    def test_use_curlybraces_with_manyvar(self):
        # Test using variable in {}
        @njit
        def foo(x, y):
            return {x: 1, y: x + y}

        x, y = 10, 20
        self.assertEqual(foo(x, y), foo.py_func(x, y))

    def test_mixed_curlybraces_and_dict(self):
        # Test mixed use of {} and dict()
        @njit
        def foo():
            k = dict()
            k[1] = {1: 3}
            k[2] = {4: 2}
            return k

        self.assertEqual(foo(), foo.py_func())

    def test_dict_use_with_none_value(self):
        # Test that NoneType cannot be used as value for Dict
        @njit
        def foo():
            k = {1: None}
            return k

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "Dict.value_type cannot be of type none",
            str(raises.exception),
        )


    def test_dict_use_with_optional_value(self):
        # Test that Optional cannot be used as value for Dict
        @njit
        def foo(choice):
            optional = 2.5 if choice else None
            k = {1: optional}
            return k

        with self.assertRaises(TypingError) as raises:
            foo(True)
        self.assertIn(
            "Dict.value_type cannot be of type OptionalType(float64)",
            str(raises.exception),
        )

    def test_dict_use_with_optional_key(self):
        # Test that Optional cannot be used as a key for Dict
        @njit
        def foo(choice):
            k = {2.5 if choice else None: 1}
            return k

        with self.assertRaises(TypingError) as raises:
            foo(True)
        self.assertIn(
            "Dict.key_type cannot be of type OptionalType(float64)",
            str(raises.exception),
        )

    def test_dict_use_with_none_key(self):
        # Test that NoneType cannot be used as a key for Dict
        @njit
        def foo():
            k = {None: 1}
            return k

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "Dict.key_type cannot be of type none",
            str(raises.exception),
        )

if __name__ == '__main__':
    unittest.main()
