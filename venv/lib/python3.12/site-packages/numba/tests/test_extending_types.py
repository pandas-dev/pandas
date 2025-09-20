"""
Test extending types via the numba.extending.* API.
"""
import operator

from numba import njit, literally
from numba.core import types, cgutils
from numba.core.errors import TypingError, NumbaTypeError
from numba.core.extending import lower_builtin
from numba.core.extending import models, register_model
from numba.core.extending import make_attribute_wrapper
from numba.core.extending import type_callable
from numba.core.extending import overload
from numba.core.extending import typeof_impl

import unittest


def gen_mock_float():
    # Stub to overload, pretending to be `float`. The real `float` function is
    # not used as multiple registrations can collide.
    def mock_float(x):
        pass
    return mock_float


class TestExtTypDummy(unittest.TestCase):

    def setUp(self):
        class Dummy(object):
            def __init__(self, value):
                self.value = value

        class DummyType(types.Type):
            def __init__(self):
                super(DummyType, self).__init__(name='Dummy')

        dummy_type = DummyType()

        @register_model(DummyType)
        class DummyModel(models.StructModel):
            def __init__(self, dmm, fe_type):
                members = [
                    ('value', types.intp),
                ]
                models.StructModel.__init__(self, dmm, fe_type, members)

        make_attribute_wrapper(DummyType, 'value', 'value')

        @type_callable(Dummy)
        def type_dummy(context):
            def typer(value):
                return dummy_type
            return typer

        @lower_builtin(Dummy, types.intp)
        def impl_dummy(context, builder, sig, args):
            typ = sig.return_type
            [value] = args
            dummy = cgutils.create_struct_proxy(typ)(context, builder)
            dummy.value = value
            return dummy._getvalue()

        @typeof_impl.register(Dummy)
        def typeof_dummy(val, c):
            return DummyType()

        # Store attributes
        self.Dummy = Dummy
        self.DummyType = DummyType

    def _add_float_overload(self, mock_float_inst):
        @overload(mock_float_inst)
        def dummy_to_float(x):
            if isinstance(x, self.DummyType):
                def codegen(x):
                    return float(x.value)
                return codegen
            else:
                raise NumbaTypeError('cannot type float({})'.format(x))

    def test_overload_float(self):
        mock_float = gen_mock_float()
        self._add_float_overload(mock_float)
        Dummy = self.Dummy

        @njit
        def foo(x):
            return mock_float(Dummy(x))

        self.assertEqual(foo(123), float(123))

    def test_overload_float_error_msg(self):
        mock_float = gen_mock_float()
        self._add_float_overload(mock_float)

        @njit
        def foo(x):
            return mock_float(x)

        with self.assertRaises(TypingError) as raises:
            foo(1j)

        self.assertIn("cannot type float(complex128)",
                      str(raises.exception))

    def test_unboxing(self):
        """A test for the unboxing logic on unknown type
        """
        Dummy = self.Dummy

        @njit
        def foo(x):
            # pass a dummy object into another function
            bar(Dummy(x))

        # make sure a cpython wrapper is created
        @njit(no_cpython_wrapper=False)
        def bar(dummy_obj):
            pass

        foo(123)
        with self.assertRaises(TypeError) as raises:
            bar(Dummy(123))
        self.assertIn("can't unbox Dummy type", str(raises.exception))

    def test_boxing(self):
        """A test for the boxing logic on unknown type
        """
        Dummy = self.Dummy

        @njit
        def foo(x):
            return Dummy(x)

        with self.assertRaises(TypeError) as raises:
            foo(123)
        self.assertIn("cannot convert native Dummy to Python object",
                      str(raises.exception))

    def test_issue5565_literal_getitem(self):
        # the following test is adapted from
        # https://github.com/numba/numba/issues/5565
        Dummy, DummyType = self.Dummy, self.DummyType

        MAGIC_NUMBER = 12321

        @overload(operator.getitem)
        def dummy_getitem_ovld(self, idx):
            if not isinstance(self, DummyType):
                return None
            # suppose we can only support idx as literal argument
            if isinstance(idx, types.StringLiteral):
                def dummy_getitem_impl(self, idx):
                    return MAGIC_NUMBER
                return dummy_getitem_impl

            if isinstance(idx, types.UnicodeType):
                def dummy_getitem_impl(self, idx):
                    return literally(idx)
                return dummy_getitem_impl

            return None

        @njit
        def test_impl(x, y):
            return Dummy(x)[y]

        var = 'abc'
        self.assertEqual(test_impl(1, var), MAGIC_NUMBER)
