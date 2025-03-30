"""
Tests for the as_numba_type() machinery.
"""
import typing as py_typing


import unittest

from numba.core import types
from numba.core.errors import TypingError
from numba.core.typing.typeof import typeof
from numba.core.typing.asnumbatype import as_numba_type, AsNumbaTypeRegistry
from numba.experimental.jitclass import jitclass
from numba.tests.support import TestCase


class TestAsNumbaType(TestCase):

    int_nb_type = typeof(0)
    float_nb_type = typeof(0.0)
    complex_nb_type = typeof(complex(0))
    str_nb_type = typeof("numba")
    bool_nb_type = typeof(True)
    none_nb_type = typeof(None)

    def test_simple_types(self):
        self.assertEqual(as_numba_type(int), self.int_nb_type)
        self.assertEqual(as_numba_type(float), self.float_nb_type)
        self.assertEqual(as_numba_type(complex), self.complex_nb_type)
        self.assertEqual(as_numba_type(str), self.str_nb_type)
        self.assertEqual(as_numba_type(bool), self.bool_nb_type)
        self.assertEqual(as_numba_type(type(None)), self.none_nb_type)

    def test_numba_types(self):
        numba_types = [
            types.intp,
            types.boolean,
            types.ListType(types.float64),
            types.DictType(
                types.intp, types.Tuple([types.float32, types.float32])
            ),
        ]

        for ty in numba_types:
            self.assertEqual(as_numba_type(ty), ty)

    def test_single_containers(self):
        self.assertEqual(
            as_numba_type(py_typing.List[float]),
            types.ListType(self.float_nb_type),
        )
        self.assertEqual(
            as_numba_type(py_typing.Dict[float, str]),
            types.DictType(self.float_nb_type, self.str_nb_type),
        )
        self.assertEqual(
            as_numba_type(py_typing.Set[complex]),
            types.Set(self.complex_nb_type),
        )
        self.assertEqual(
            as_numba_type(py_typing.Tuple[float, float]),
            types.Tuple([self.float_nb_type, self.float_nb_type]),
        )
        self.assertEqual(
            as_numba_type(py_typing.Tuple[float, complex]),
            types.Tuple([self.float_nb_type, self.complex_nb_type]),
        )

    def test_optional(self):
        self.assertEqual(
            as_numba_type(py_typing.Optional[float]),
            types.Optional(self.float_nb_type),
        )
        self.assertEqual(
            as_numba_type(py_typing.Union[str, None]),
            types.Optional(self.str_nb_type),
        )
        self.assertEqual(
            as_numba_type(py_typing.Union[None, bool]),
            types.Optional(self.bool_nb_type),
        )

        # Optional[x] is a special case of Union[x, None].  We raise a
        # TypingError if the right type is not NoneType.
        with self.assertRaises(TypingError) as raises:
            as_numba_type(py_typing.Union[int, float])
        self.assertIn("Cannot type Union that is not an Optional",
                      str(raises.exception))

    def test_nested_containers(self):
        IntList = py_typing.List[int]
        self.assertEqual(
            as_numba_type(py_typing.List[IntList]),
            types.ListType(types.ListType(self.int_nb_type)),
        )
        self.assertEqual(
            as_numba_type(py_typing.List[py_typing.Dict[float, bool]]),
            types.ListType(
                types.DictType(self.float_nb_type, self.bool_nb_type)
            ),
        )
        self.assertEqual(
            as_numba_type(
                py_typing.Set[py_typing.Tuple[py_typing.Optional[int], float]]),
            types.Set(types.Tuple(
                [types.Optional(self.int_nb_type), self.float_nb_type])),
        )

    def test_jitclass_registers(self):

        @jitclass
        class MyInt:
            x: int

            def __init__(self, value):
                self.x = value

        self.assertEqual(as_numba_type(MyInt), MyInt.class_type.instance_type)

    def test_type_alias(self):
        Pair = py_typing.Tuple[int, int]
        ListOfPairs = py_typing.List[Pair]

        pair_nb_type = types.Tuple((self.int_nb_type, self.int_nb_type))
        self.assertEqual(as_numba_type(Pair), pair_nb_type)
        self.assertEqual(
            as_numba_type(ListOfPairs), types.ListType(pair_nb_type)
        )

    def test_overwrite_type(self):
        as_numba_type = AsNumbaTypeRegistry()
        self.assertEqual(as_numba_type(float), self.float_nb_type)
        as_numba_type.register(float, types.float32)
        self.assertEqual(as_numba_type(float), types.float32)
        self.assertNotEqual(as_numba_type(float), self.float_nb_type)

    def test_any_throws(self):
        Any = py_typing.Any

        any_types = [
            py_typing.Optional[Any],
            py_typing.List[Any],
            py_typing.Set[Any],
            py_typing.Dict[float, Any],
            py_typing.Dict[Any, float],
            py_typing.Tuple[int, Any],
        ]

        for bad_py_type in any_types:
            with self.assertRaises(TypingError) as raises:
                as_numba_type(bad_py_type)
            self.assertIn(
                "Cannot infer numba type of python type",
                str(raises.exception),
            )

    def test_bad_union_throws(self):
        bad_unions = [
            py_typing.Union[str, int],
            py_typing.Union[int, type(None), py_typing.Tuple[bool, bool]],
        ]

        for bad_py_type in bad_unions:
            with self.assertRaises(TypingError) as raises:
                as_numba_type(bad_py_type)
            self.assertIn("Cannot type Union", str(raises.exception))


if __name__ == '__main__':
    unittest.main()
