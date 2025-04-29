"""
Tests for numba.types.
"""


from collections import namedtuple
import gc
import os
import operator
import sys
import weakref

import numpy as np

from numba.core import types, typing, errors, sigutils
from numba.core.types.abstract import _typecache
from numba.core.types.functions import _header_lead
from numba.core.typing.templates import make_overload_template
from numba import jit, njit, typeof
from numba.core.extending import (overload, register_model, models, unbox,
                                  NativeValue, typeof_impl)
from numba.tests.support import TestCase, create_temp_module
from numba.tests.enum_usecases import Color, Shake, Shape
import unittest
from numba.np import numpy_support


try:
    import cPickle as pickle
except ImportError:
    import pickle


Point = namedtuple('Point', ('x', 'y'))

Rect = namedtuple('Rect', ('width', 'height'))


def gen(x):
    yield x + 1


class Dummy(object):
    pass


class TestTypes(TestCase):

    def test_equality(self):
        self.assertEqual(types.int32, types.int32)
        self.assertEqual(types.uint32, types.uint32)
        self.assertEqual(types.complex64, types.complex64)
        self.assertEqual(types.float32, types.float32)
        # Different signedness
        self.assertNotEqual(types.int32, types.uint32)
        # Different width
        self.assertNotEqual(types.int64, types.int32)
        self.assertNotEqual(types.float64, types.float32)
        self.assertNotEqual(types.complex64, types.complex128)
        # Different domain
        self.assertNotEqual(types.int64, types.float64)
        self.assertNotEqual(types.uint64, types.float64)
        self.assertNotEqual(types.complex64, types.float64)
        # Same arguments but different return types
        get_pointer = None
        sig_a = typing.signature(types.intp, types.intp)
        sig_b = typing.signature(types.voidptr, types.intp)
        a = types.ExternalFunctionPointer(sig=sig_a, get_pointer=get_pointer)
        b = types.ExternalFunctionPointer(sig=sig_b, get_pointer=get_pointer)
        self.assertNotEqual(a, b)
        # Different call convention
        a = types.ExternalFunctionPointer(sig=sig_a, get_pointer=get_pointer)
        b = types.ExternalFunctionPointer(sig=sig_a, get_pointer=get_pointer,
                                          cconv='stdcall')
        self.assertNotEqual(a, b)
        # Different get_pointer
        a = types.ExternalFunctionPointer(sig=sig_a, get_pointer=get_pointer)
        b = types.ExternalFunctionPointer(sig=sig_a, get_pointer=object())
        self.assertNotEqual(a, b)

        # Different template classes bearing the same name
        class DummyTemplate(object):
            key = "foo"
        a = types.BoundFunction(DummyTemplate, types.int32)

        class DummyTemplate(object):
            key = "bar"
        b = types.BoundFunction(DummyTemplate, types.int32)
        self.assertNotEqual(a, b)

        # Different dtypes
        self.assertNotEqual(types.DType(types.int32), types.DType(types.int64))

        # CPointer with same addrspace
        self.assertEqual(types.CPointer(types.float32),
                         types.CPointer(types.float32))

        # CPointer with different addrspace
        self.assertNotEqual(types.CPointer(types.float32, 0),
                            types.CPointer(types.float32, 1))

    def test_weaktype(self):
        d = Dummy()
        e = Dummy()
        a = types.Dispatcher(d)
        b = types.Dispatcher(d)
        c = types.Dispatcher(e)
        self.assertIs(a.dispatcher, d)
        self.assertIs(b.dispatcher, d)
        self.assertIs(c.dispatcher, e)
        # Equality of alive references
        self.assertTrue(a == b)
        self.assertFalse(a != b)
        self.assertTrue(a != c)
        self.assertFalse(a == c)
        z = types.int8
        self.assertFalse(a == z)
        self.assertFalse(b == z)
        self.assertFalse(c == z)
        self.assertTrue(a != z)
        self.assertTrue(b != z)
        self.assertTrue(c != z)
        # Hashing and mappings
        s = set([a, b, c])
        self.assertEqual(len(s), 2)
        self.assertIn(a, s)
        self.assertIn(b, s)
        self.assertIn(c, s)
        # Kill the references
        d = e = None
        gc.collect()
        with self.assertRaises(ReferenceError):
            a.dispatcher
        with self.assertRaises(ReferenceError):
            b.dispatcher
        with self.assertRaises(ReferenceError):
            c.dispatcher
        # Dead references are always unequal
        self.assertFalse(a == b)
        self.assertFalse(a == c)
        self.assertFalse(b == c)
        self.assertFalse(a == z)
        self.assertTrue(a != b)
        self.assertTrue(a != c)
        self.assertTrue(b != c)
        self.assertTrue(a != z)

    def test_interning(self):
        # Test interning and lifetime of dynamic types.
        a = types.Dummy('xyzzyx')
        code = a._code
        b = types.Dummy('xyzzyx')
        self.assertIs(b, a)
        wr = weakref.ref(a)
        del a
        gc.collect()
        c = types.Dummy('xyzzyx')
        self.assertIs(c, b)
        # The code is always the same
        self.assertEqual(c._code, code)
        del b, c
        gc.collect()
        self.assertIs(wr(), None)
        d = types.Dummy('xyzzyx')
        # The original code wasn't reused.
        self.assertNotEqual(d._code, code)

    def test_cache_trimming(self):
        # Test that the cache doesn't grow in size when types are
        # created and disposed of.
        cache = _typecache
        gc.collect()
        # Keep strong references to existing types, to avoid spurious failures
        existing_types = [wr() for wr in cache]  # noqa: F841
        cache_len = len(cache)
        a = types.Dummy('xyzzyx')
        b = types.Dummy('foox')
        self.assertEqual(len(cache), cache_len + 2)
        del a, b
        gc.collect()
        self.assertEqual(len(cache), cache_len)

    def test_array_notation(self):
        def check(arrty, scalar, ndim, layout):
            self.assertIs(arrty.dtype, scalar)
            self.assertEqual(arrty.ndim, ndim)
            self.assertEqual(arrty.layout, layout)

        def check_index_error(callable):
            with self.assertRaises(KeyError) as raises:
                callable()
            self.assertIn(
                "Can only index numba types with slices with no start or "
                "stop, got", str(raises.exception))

        scalar = types.int32
        check(scalar[:], scalar, 1, 'A')
        check(scalar[::1], scalar, 1, 'C')
        check(scalar[:, :], scalar, 2, 'A')
        check(scalar[:, ::1], scalar, 2, 'C')
        check(scalar[::1, :], scalar, 2, 'F')

        check_index_error(lambda: scalar[0])
        check_index_error(lambda: scalar[:, 4])
        check_index_error(lambda: scalar[::1, 1:])
        check_index_error(lambda: scalar[:2])
        check_index_error(lambda: list(scalar))

    def test_array_notation_for_dtype(self):
        def check(arrty, scalar, ndim, layout):
            self.assertIs(arrty.dtype, scalar)
            self.assertEqual(arrty.ndim, ndim)
            self.assertEqual(arrty.layout, layout)
        scalar = types.int32
        dtyped = types.DType(scalar)
        check(dtyped[:], scalar, 1, 'A')
        check(dtyped[::1], scalar, 1, 'C')
        check(dtyped[:, :], scalar, 2, 'A')
        check(dtyped[:, ::1], scalar, 2, 'C')
        check(dtyped[::1, :], scalar, 2, 'F')

    def test_call_notation(self):
        # Function call signature
        i = types.int32
        d = types.double
        self.assertEqual(i(), typing.signature(i))
        self.assertEqual(i(d), typing.signature(i, d))
        self.assertEqual(i(d, d), typing.signature(i, d, d))
        # Value cast
        self.assertPreciseEqual(i(42.5), 42)
        self.assertPreciseEqual(d(-5), -5.0)
        ty = types.NPDatetime('Y')
        self.assertPreciseEqual(ty('1900'), np.datetime64('1900', 'Y'))
        self.assertPreciseEqual(ty('NaT'), np.datetime64('NaT', 'Y'))
        ty = types.NPTimedelta('s')
        self.assertPreciseEqual(ty(5), np.timedelta64(5, 's'))
        self.assertPreciseEqual(ty('NaT'), np.timedelta64('NaT', 's'))
        ty = types.NPTimedelta('')
        self.assertPreciseEqual(ty(5), np.timedelta64(5))
        self.assertPreciseEqual(ty('NaT'), np.timedelta64('NaT'))

    def test_list_type_getitem(self):
        for listty in (types.int64, types.Array(types.float64, 1, 'C')):
            l_int = types.List(listty)
            self.assertTrue(isinstance(l_int, types.List))
            self.assertTrue(isinstance(l_int[0], type(listty)))

    def test_function_incompatible_templates(self):
        # issue 4345
        def func_stub():
            pass

        def func_stub2():
            pass

        def ol():
            pass

        template1 = make_overload_template(func_stub, ol, {}, True, 'never')
        template2 = make_overload_template(func_stub2, ol, {}, True, 'never')

        with self.assertRaises(ValueError) as raises:
            types.Function((template1, template2))
        self.assertIn("incompatible templates:", str(raises.exception))


class TestNumbers(TestCase):
    """
    Tests for number types.
    """

    def test_bitwidth(self):
        """
        All numeric types have bitwidth attribute
        """
        for ty in types.number_domain:
            self.assertTrue(hasattr(ty, "bitwidth"))

    def test_minval_maxval(self):
        self.assertEqual(types.int8.maxval, 127)
        self.assertEqual(types.int8.minval, -128)
        self.assertEqual(types.uint8.maxval, 255)
        self.assertEqual(types.uint8.minval, 0)
        self.assertEqual(types.int64.maxval, (1 << 63) - 1)
        self.assertEqual(types.int64.minval, -(1 << 63))
        self.assertEqual(types.uint64.maxval, (1 << 64) - 1)
        self.assertEqual(types.uint64.minval, 0)

    def test_from_bidwidth(self):
        f = types.Integer.from_bitwidth
        self.assertIs(f(32), types.int32)
        self.assertIs(f(8, signed=False), types.uint8)

    def test_ordering(self):
        def check_order(values):
            for i in range(len(values)):
                self.assertLessEqual(values[i], values[i])
                self.assertGreaterEqual(values[i], values[i])
                self.assertFalse(values[i] < values[i])
                self.assertFalse(values[i] > values[i])
                for j in range(i):
                    self.assertLess(values[j], values[i])
                    self.assertLessEqual(values[j], values[i])
                    self.assertGreater(values[i], values[j])
                    self.assertGreaterEqual(values[i], values[j])
                    self.assertFalse(values[i] < values[j])
                    self.assertFalse(values[i] <= values[j])
                    self.assertFalse(values[j] > values[i])
                    self.assertFalse(values[j] >= values[i])

        check_order([types.int8, types.int16, types.int32, types.int64])
        check_order([types.uint8, types.uint16, types.uint32, types.uint64])
        check_order([types.float32, types.float64])
        check_order([types.complex64, types.complex128])

        with self.assertRaises(TypeError):
            types.int8 <= types.uint32
        with self.assertRaises(TypeError):
            types.int8 <= types.float32
        with self.assertRaises(TypeError):
            types.float64 <= types.complex128


class TestNdIter(TestCase):

    def test_properties(self):
        def check(ty, dtypes, ndim, layout, indexers=None):
            self.assertEqual(ty.ndim, ndim)
            self.assertEqual(ty.layout, layout)
            self.assertEqual(ty.dtypes, dtypes)
            views = [types.Array(dtype, 0, "C") for dtype in dtypes]
            if len(views) > 1:
                self.assertEqual(
                    ty.yield_type,
                    types.BaseTuple.from_types(views))
            else:
                self.assertEqual(ty.yield_type, views[0])
            if indexers is not None:
                self.assertEqual(ty.indexers, indexers)

        f32 = types.float32
        c64 = types.complex64
        i16 = types.int16
        a = types.Array(f32, 1, "C")
        b = types.Array(f32, 2, "C")
        c = types.Array(c64, 2, "F")
        d = types.Array(i16, 2, "A")
        e = types.Array(i16, 0, "C")
        f = types.Array(f32, 1, "A")
        g = types.Array(f32, 0, "C")

        # 0-dim iterator
        ty = types.NumpyNdIterType((e,))
        check(ty, (i16,), 0, "C", [('0d', 0, 0, [0])])
        self.assertFalse(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((e, g))
        check(ty, (i16, f32), 0, "C", [('0d', 0, 0, [0, 1])])
        self.assertFalse(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((e, c64))
        check(ty, (i16, c64), 0, "C",
              [('0d', 0, 0, [0]), ('scalar', 0, 0, [1])])
        self.assertFalse(ty.need_shaped_indexing)

        # 1-dim iterator
        ty = types.NumpyNdIterType((a,))
        check(ty, (f32,), 1, "C",
              [('flat', 0, 1, [0])])
        self.assertFalse(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((a, a))
        check(ty, (f32, f32), 1, "C",
              [('flat', 0, 1, [0, 1])])
        self.assertFalse(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((a, e, e, c64))
        check(ty, (f32, i16, i16, c64), 1, "C",
              [('flat', 0, 1, [0]),  # a
               ('0d', 0, 0, [1, 2]),  # e, e
               ('scalar', 0, 0, [3]),  # c64
               ])
        self.assertFalse(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((a, f))
        check(ty, (f32, f32), 1, "C",
              [('flat', 0, 1, [0]), ('indexed', 0, 1, [1])])
        self.assertTrue(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((f,))
        check(ty, (f32,), 1, "C", [('indexed', 0, 1, [0])])
        self.assertTrue(ty.need_shaped_indexing)

        # 2-dim C-order iterator
        ty = types.NumpyNdIterType((b,))
        check(ty, (f32,), 2, "C", [('flat', 0, 2, [0])])
        self.assertFalse(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((b, c))
        check(
            ty, (f32, c64), 2, "C", [
                ('flat', 0, 2, [0]), ('indexed', 0, 2, [1])])
        self.assertTrue(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((d,))
        check(ty, (i16,), 2, "C", [('indexed', 0, 2, [0])])
        self.assertTrue(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((b, c, d, d, e))
        check(ty, (f32, c64, i16, i16, i16), 2, "C",
              [('flat', 0, 2, [0]),  # b
               ('indexed', 0, 2, [1, 2, 3]),  # c, d, d
               ('0d', 0, 0, [4]),  # e
               ])
        self.assertTrue(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((a, b, c, d, d, f))
        check(ty, (f32, f32, c64, i16, i16, f32), 2, "C",
              [('flat', 1, 2, [0]),  # a
               ('flat', 0, 2, [1]),  # b
               ('indexed', 0, 2, [2, 3, 4]),  # c, d, d
               ('indexed', 1, 2, [5]),  # f
               ])
        self.assertTrue(ty.need_shaped_indexing)

        # 2-dim F-order iterator
        ty = types.NumpyNdIterType((c,))
        check(ty, (c64,), 2, "F", [('flat', 0, 2, [0])])
        self.assertFalse(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((c, b, c, f))
        check(ty, (c64, f32, c64, f32), 2, "F",
              [('flat', 0, 2, [0, 2]),  # c, c
               ('indexed', 0, 2, [1]),  # b
               ('indexed', 0, 1, [3]),  # f
               ])
        self.assertTrue(ty.need_shaped_indexing)
        ty = types.NumpyNdIterType((b, c, c, d, d, a, e))
        check(ty, (f32, c64, c64, i16, i16, f32, i16), 2, "F",
              [('indexed', 0, 2, [0, 3, 4]),  # b, d, d
               ('flat', 0, 2, [1, 2]),  # c, c
               ('flat', 0, 1, [5]),  # a
               ('0d', 0, 0, [6]),  # e
               ])
        self.assertTrue(ty.need_shaped_indexing)


class TestPickling(TestCase):
    """
    Pickling and unpickling should preserve type identity (singleton-ness)
    and the _code attribute.  This is only a requirement for types that
    can be part of function signatures.
    """

    def predefined_types(self):
        """
        Yield all predefined type instances
        """
        for ty in types.__dict__.values():
            if isinstance(ty, types.Type):
                yield ty

    def check_pickling(self, orig):
        pickled = pickle.dumps(orig, protocol=-1)
        ty = pickle.loads(pickled)
        self.assertIs(ty, orig)
        self.assertGreaterEqual(ty._code, 0)

    def test_predefined_types(self):
        tys = list(self.predefined_types())
        self.assertIn(types.int16, tys)
        for ty in tys:
            self.check_pickling(ty)

    def test_atomic_types(self):
        for unit in ('M', 'ms'):
            ty = types.NPDatetime(unit)
            self.check_pickling(ty)
            ty = types.NPTimedelta(unit)
            self.check_pickling(ty)

    def test_arrays(self):
        for ndim in (0, 1, 2):
            for layout in ('A', 'C', 'F'):
                ty = types.Array(types.int16, ndim, layout)
                self.check_pickling(ty)

    def test_records(self):
        recordtype = np.dtype([('a', np.float64),
                               ('b', np.int32),
                               ('c', np.complex64),
                               ('d', (np.str_, 5))])
        ty = numpy_support.from_dtype(recordtype)
        self.check_pickling(ty)
        self.check_pickling(types.Array(ty, 1, 'A'))

    def test_optional(self):
        ty = types.Optional(types.int32)
        self.check_pickling(ty)

    def test_tuples(self):
        ty1 = types.UniTuple(types.int32, 3)
        self.check_pickling(ty1)
        ty2 = types.Tuple((types.int32, ty1))
        self.check_pickling(ty2)

    def test_namedtuples(self):
        ty1 = types.NamedUniTuple(types.intp, 2, Point)
        self.check_pickling(ty1)
        ty2 = types.NamedTuple((types.intp, types.float64), Point)
        self.check_pickling(ty2)

    def test_enums(self):
        ty1 = types.EnumMember(Color, types.int32)
        self.check_pickling(ty1)
        ty2 = types.EnumMember(Shake, types.int64)
        self.check_pickling(ty2)
        ty3 = types.IntEnumMember(Shape, types.int64)
        self.check_pickling(ty3)

    def test_lists(self):
        ty = types.List(types.int32)
        self.check_pickling(ty)

    def test_generator(self):
        cfunc = jit("(int32,)", nopython=True)(gen)
        sigs = list(cfunc.nopython_signatures)
        ty = sigs[0].return_type
        self.assertIsInstance(ty, types.Generator)
        self.check_pickling(ty)

    # call templates are not picklable
    @unittest.expectedFailure
    def test_external_function_pointers(self):
        from numba.core.typing import ctypes_utils
        from numba.tests.ctypes_usecases import c_sin, c_cos
        for fnptr in (c_sin, c_cos):
            ty = ctypes_utils.make_function_type(fnptr)
            self.assertIsInstance(ty, types.ExternalFunctionPointer)
            self.check_pickling(ty)


class TestSignatures(TestCase):

    def test_normalize_signature(self):
        f = sigutils.normalize_signature

        def check(sig, args, return_type):
            self.assertEqual(f(sig), (args, return_type))

        def check_error(sig, msg):
            with self.assertRaises(TypeError) as raises:
                f(sig)
            self.assertIn(msg, str(raises.exception))

        f32 = types.float32
        c64 = types.complex64
        i16 = types.int16
        a = types.Array(f32, 1, "C")

        check((c64,), (c64,), None)
        check((f32, i16), (f32, i16), None)
        check(a(i16), (i16,), a)
        check("int16(complex64)", (c64,), i16)
        check("(complex64, int16)", (c64, i16), None)
        check(typing.signature(i16, c64), (c64,), i16)

        msg = "invalid type in signature: expected a type instance"
        check_error((types.Integer,), msg)
        check_error((None,), msg)
        check_error([], "invalid signature")


class TestRecordDtype(unittest.TestCase):
    def test_record_type_equiv(self):
        rec_dt = np.dtype([('a', np.int32), ('b', np.float32)])
        rec_ty = typeof(rec_dt)
        art1 = rec_ty[::1]
        arr = np.zeros(5, dtype=rec_dt)
        art2 = typeof(arr)
        self.assertEqual(art2.dtype.dtype, rec_ty)
        self.assertEqual(art1, art2)

    def test_user_specified(self):
        rec_dt = np.dtype([('a', np.int32), ('b', np.float32)])
        rec_type = typeof(rec_dt)

        @jit((rec_type[:],), nopython=True)
        def foo(x):
            return x['a'], x['b']

        arr = np.zeros(1, dtype=rec_dt)
        arr[0]['a'] = 123
        arr[0]['b'] = 32.1

        a, b = foo(arr)

        self.assertEqual(a, arr[0]['a'])
        self.assertEqual(b, arr[0]['b'])


class TestDType(TestCase):
    def test_type_attr(self):
        # Test .type attribute of dtype
        def conv(arr, val):
            return arr.dtype.type(val)

        jit_conv = jit(nopython=True)(conv)

        def assert_matches(arr, val, exact):
            expect = conv(arr, val)
            got = jit_conv(arr, val)
            self.assertPreciseEqual(expect, exact)
            self.assertPreciseEqual(typeof(expect), typeof(got))
            self.assertPreciseEqual(expect, got)

        arr = np.zeros(5)
        assert_matches(arr.astype(np.intp), 1.2, 1)
        assert_matches(arr.astype(np.float64), 1.2, 1.2)
        assert_matches(arr.astype(np.complex128), 1.2, (1.2 + 0j))
        assert_matches(arr.astype(np.complex128), 1.2j, 1.2j)

    def test_kind(self):
        def tkind(A):
            return A.dtype.kind == 'f'
        jit_tkind = jit(nopython=True)(tkind)
        self.assertEqual(tkind(np.ones(3)), jit_tkind(np.ones(3)))
        self.assertEqual(tkind(np.ones(3, dtype=np.intp)),
                         jit_tkind(np.ones(3, dtype=np.intp)))

    def test_dtype_with_type(self):
        def impl():
            a = np.dtype(np.float64)
            return a.type(0)
        jit_impl = jit(nopython=True)(impl)
        self.assertEqual(impl(), jit_impl())

    def test_dtype_with_string(self):
        def impl():
            a = np.dtype('float64')
            return a.type(0)
        jit_impl = jit(nopython=True)(impl)
        self.assertEqual(impl(), jit_impl())


class TestIsInternalTypeMarker(TestCase):
    """Tests the use of the Type metaclass init correctly setting the flag on
    the `is_internal` attr of a concrete Type class
    """
    source_lines = """
from numba.core import types

class FooType(types.Type):
    def __init__(self):
        super(FooType, self).__init__(name='Foo')
"""

    def test_create_temp_module(self):
        sys_path_original = list(sys.path)
        sys_modules_original = dict(sys.modules)
        with create_temp_module(self.source_lines) as test_module:
            temp_module_dir = os.path.dirname(test_module.__file__)
            self.assertEqual(temp_module_dir, sys.path[0])
            self.assertEqual(sys.path[1:], sys_path_original)
            self.assertTrue(test_module.__name__ in sys.modules)
        # Test that modifications to sys.path / sys.modules are reverted
        self.assertEqual(sys.path, sys_path_original)
        self.assertEqual(sys.modules, sys_modules_original)

    def test_create_temp_module_with_exception(self):
        try:
            sys_path_original = list(sys.path)
            sys_modules_original = dict(sys.modules)
            with create_temp_module(self.source_lines):
                raise ValueError("Something went wrong!")
        except ValueError:
            # Test that modifications to sys.path / sys.modules are reverted
            self.assertEqual(sys.path, sys_path_original)
            self.assertEqual(sys.modules, sys_modules_original)

    def test_externally_defined_type_is_external(self):

        with create_temp_module(self.source_lines) as test_module:
            FooType = test_module.FooType
            self.assertFalse(FooType().is_internal)

            # set up an extension type
            class Foo(object):
                pass

            register_model(FooType)(models.OpaqueModel)

            @typeof_impl.register(Foo)
            def _typ_foo(val, c):
                return FooType()

            @unbox(FooType)
            def unbox_foo(typ, obj, c):
                return NativeValue(c.context.get_dummy_value())

            # function to overload
            def false_if_not_array(a):
                pass

            # Set up an overload which will accept all types irrespective of
            # whether they are from Numba's closed type system
            @overload(false_if_not_array)
            def ol_false_if_not_array(a):
                if isinstance(a, types.Array):
                    return lambda a : True
                else:
                    return lambda a : False

            @njit
            def call_false_if_not_array(a):
                return false_if_not_array(a)

            self.assertTrue(call_false_if_not_array(np.zeros(10)))
            self.assertFalse(call_false_if_not_array(10))

            # The extension type was accepted
            self.assertFalse(call_false_if_not_array(Foo()))

            # Now do the same sort of overload but put in a guard based on the
            # use of internal types

            def false_if_not_array_closed_system(a):
                pass

            @overload(false_if_not_array_closed_system)
            def ol_false_if_not_array_closed_system(a):
                if a.is_internal:  # guard
                    if isinstance(a, types.Array):
                        return lambda a : True
                    else:
                        return lambda a : False

            @njit
            def call_false_if_not_array_closed_system(a):
                return false_if_not_array_closed_system(a)

            self.assertTrue(call_false_if_not_array_closed_system(np.zeros(10)))
            self.assertFalse(call_false_if_not_array_closed_system(10))

            with self.assertRaises(errors.TypingError) as raises:
                call_false_if_not_array_closed_system(Foo())
            estr = str(raises.exception)
            self.assertIn(_header_lead, estr)
            self.assertIn("false_if_not_array_closed_system", estr)
            self.assertIn("(Foo)", estr)

    def test_mixin_against_real_example(self):
        # See issue #4970, this checks that unicode eq/ne now ignores extension
        # types.

        with create_temp_module(self.source_lines) as test_module:
            FooType = test_module.FooType
            self.assertFalse(FooType().is_internal)

            # set up an extension type
            class Foo(object):
                pass

            register_model(FooType)(models.OpaqueModel)

            @typeof_impl.register(Foo)
            def _typ_foo(val, c):
                return FooType()

            @unbox(FooType)
            def unbox_foo(typ, obj, c):
                return NativeValue(c.context.get_dummy_value())

            @overload(operator.eq)
            def foo_eq(a, b):
                if a == FooType():
                    return lambda a, b: "RAN CUSTOM EQ OVERLOAD"

            @overload(operator.ne)
            def foo_ne(a, b):
                if a == FooType():
                    return lambda a, b: "RAN CUSTOM NE OVERLOAD"

            @njit
            def f(a):
                return a == "A", a != "A"

            self.assertEqual(("RAN CUSTOM EQ OVERLOAD",
                              "RAN CUSTOM NE OVERLOAD"),
                             f(Foo()))


class TestIssues(TestCase):
    def test_omitted_type(self):
        # issue https://github.com/numba/numba/issues/5471
        def inner(a):
            pass

        @overload(inner)
        def inner_overload(a):
            if not isinstance(a, types.Literal):
                return
            return lambda a: a

        @njit
        def my_func(a='a'):
            return inner(a)

        @njit
        def f():
            return my_func()

        @njit
        def g():
            return my_func('b')

        self.assertEqual(f(), 'a')
        self.assertEqual(g(), 'b')

    def test_type_of_literal(self):
        # type(val) where val is a literal should not give a literal type.
        def inner(a):
            pass

        @overload(inner)
        def inner_overload(a):
            if not isinstance(a, types.Literal):
                return
            self.assertIsInstance(a, types.Literal)
            # NOTE: using 1.23 to ensure that the result is indeed an int.
            return lambda a: type(a)(a + 1.23)

        @njit
        def my_func(a=1):
            return inner(a)

        @njit
        def f():
            return my_func()

        @njit
        def g():
            return my_func(100)

        self.assertEqual(f(), 2)
        self.assertEqual(g(), 101)

    def test_issue_typeref_key(self):
        # issue https://github.com/numba/numba/issues/6336
        class NoUniqueNameType(types.Dummy):
            def __init__(self, param):
                super(NoUniqueNameType, self).__init__('NoUniqueNameType')
                self.param = param

            @property
            def key(self):
                return self.param

        no_unique_name_type_1 = NoUniqueNameType(1)
        no_unique_name_type_2 = NoUniqueNameType(2)

        for ty1 in (no_unique_name_type_1, no_unique_name_type_2):
            for ty2 in (no_unique_name_type_1, no_unique_name_type_2):
                self.assertIs(
                    types.TypeRef(ty1) == types.TypeRef(ty2), ty1 == ty2)

    def test_issue_list_type_key(self):
        # https://github.com/numba/numba/issues/6397
        class NoUniqueNameType(types.Dummy):
            def __init__(self, param):
                super(NoUniqueNameType, self).__init__('NoUniqueNameType')
                self.param = param

            @property
            def key(self):
                return self.param

        no_unique_name_type_1 = NoUniqueNameType(1)
        no_unique_name_type_2 = NoUniqueNameType(2)

        for ty1 in (no_unique_name_type_1, no_unique_name_type_2):
            for ty2 in (no_unique_name_type_1, no_unique_name_type_2):
                self.assertIs(
                    types.ListType(ty1) == types.ListType(ty2),  # noqa: E721
                    ty1 == ty2
                )


if __name__ == '__main__':
    unittest.main()
