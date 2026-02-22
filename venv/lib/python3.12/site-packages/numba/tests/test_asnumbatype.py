"""
Tests for the as_numba_type() machinery.
"""
import typing as py_typing


import unittest

from contextlib import ExitStack
from llvmlite import ir

from numba import jit
from numba.core import cgutils, types
from numba.core.datamodel.models import PrimitiveModel
from numba.core.errors import TypingError
from numba.core.extending import (register_model, type_callable, unbox,
                                  NativeValue)
from numba.core.types import Number
from numba.core.typing.typeof import typeof, typeof_impl
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
                "Cannot infer Numba type of Python type",
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

    def test_instance_check_usecase(self):
        # Demonstrates that registering the type class with as_numba_type
        # supports instance checks, at least for those subclasses supported by
        # the instance check (e.g. types.Number, etc.).
        #
        # To set up the test we need quite a lot of extension code to support
        # a new type independent of the existing types.

        # The Python class
        class bfloat16:
            def __init__(self, value):
                self._value = value

        # The Numba type class - we use a Number subclass both because it makes
        # sense for a new numeric type, and it's one of the types supported by
        # instance checks in Numba
        class _type_class_bfloat16(Number):
            def __init__(self):
                self.bitwidth = 16
                super().__init__(name="bfloat16")

        # The Numba type instance
        bfloat16_type = _type_class_bfloat16()

        # Register typing of the Python class for use as arguments and
        # constants
        @typeof_impl.register(bfloat16)
        def typeof_bfloat16(val, c):
            return bfloat16_type

        # A data model for the bfloat16 class. We don't need much actual
        # implementation so it doesn't matter too much what this is; a 16-bit
        # integer representation is sufficient.
        @register_model(_type_class_bfloat16)
        class _model_bfloat16(PrimitiveModel):
            def __init__(self, dmm, fe_type):
                be_type = ir.IntType(fe_type.bitwidth)
                super(_model_bfloat16, self).__init__(dmm, fe_type, be_type)

        # Ideally we pass in a value so we ensure that the instance check is
        # working with values dynamically passed in (preventing the whole check
        # being potentially optimized into a simple True or False). For this we
        # need an unboxing.
        @unbox(_type_class_bfloat16)
        def unbox_bfloat16(ty, obj, c):
            ll_type = c.context.get_argument_type(ty)
            val = cgutils.alloca_once(c.builder, ll_type)
            is_error_ptr = cgutils.alloca_once_value(c.builder,
                                                     cgutils.false_bit)

            with ExitStack() as stack:
                value_obj = c.pyapi.object_getattr_string(obj, "_value")

                with cgutils.early_exit_if_null(c.builder, stack, value_obj):
                    c.builder.store(cgutils.true_bit, is_error_ptr)

                value_native = c.unbox(types.uint16, value_obj)
                c.pyapi.decref(value_obj)

                with cgutils.early_exit_if(c.builder, stack,
                                           value_native.is_error):
                    c.builder.store(cgutils.true_bit, is_error_ptr)

                c.builder.store(value_native.value, val)

            return NativeValue(c.builder.load(val),
                               is_error=c.builder.load(is_error_ptr))

        # We never call bfloat16 to construct one inside a jitted function, but
        # we need this typing so that the type of the bfloat16 class can be
        # determined (it's the argument to the instance check).
        @type_callable(bfloat16)
        def type_bfloat16_ctor(context):
            # Note that the typer is never called in this test, because we
            # don't call bfloat16 - only the typing of it as a callable is
            # used.

            def typer(value):
                if isinstance(value, types.Integer):
                    return bfloat16_type

        # First we try the instance check without an as_numba_type
        # registration, to prove that it is necessary for instance checks to
        # work.
        @jit
        def instancecheck_no_ant_reg(x):
            return isinstance(x, bfloat16)

        # A "random" value to test with
        x_bf16 = bfloat16(0x4049)  # bfloat16(3.14)

        # Ensure the typing fails without the registration
        expected_message = r"Cannot infer Numba type of Python type.*bfloat16"
        with self.assertRaisesRegex(TypingError, expected_message):
            instancecheck_no_ant_reg(x_bf16)

        # Register the typing with as_numba_type so that we can expect instance
        # checks to work
        as_numba_type.register(bfloat16, bfloat16_type)

        # We define a new function to ensure all registrations are as-required
        @jit
        def instancecheck(x):
            return isinstance(x, bfloat16)

        # The instance check should be True for bfloat16 instances and False
        # otherwise.
        self.assertTrue(instancecheck(x_bf16))
        self.assertFalse(instancecheck(1))


if __name__ == '__main__':
    unittest.main()
