import os, sys, subprocess
import dis
import itertools

import numpy as np

import numba
from numba import jit, njit
from numba.core import errors, ir, types, typing, typeinfer, utils
from numba.core.typeconv import Conversion
from numba.extending import overload_method

from numba.tests.support import TestCase, tag
from numba.tests.test_typeconv import CompatibilityTestMixin
from numba.core.untyped_passes import TranslateByteCode, IRProcessing
from numba.core.typed_passes import PartialTypeInference
from numba.core.compiler_machinery import FunctionPass, register_pass
import unittest


i8 = types.int8
i16 = types.int16
i32 = types.int32
i64 = types.int64
u8 = types.uint8
u16 = types.uint16
u32 = types.uint32
u64 = types.uint64
f32 = types.float32
f64 = types.float64
c64 = types.complex64
c128 = types.complex128

skip_unless_load_fast_and_clear = unittest.skipUnless(
    "LOAD_FAST_AND_CLEAR" in dis.opmap,
    "Requires LOAD_FAST_AND_CLEAR opcode",
)


class TestArgRetCasting(unittest.TestCase):
    def test_arg_ret_casting(self):
        def foo(x):
            return x

        args = (i32,)
        return_type = f32
        cfunc = njit(return_type(*args))(foo)
        cres = cfunc.overloads[args]
        self.assertTrue(isinstance(cfunc(123), float))
        self.assertEqual(cres.signature.args, args)
        self.assertEqual(cres.signature.return_type, return_type)

    def test_arg_ret_mismatch(self):
        def foo(x):
            return x

        args = (types.Array(i32, 1, 'C'),)
        return_type = f32
        try:
            njit(return_type(*args))(foo)
        except errors.TypingError as e:
            pass
        else:
            self.fail("Should complain about array casting to float32")

    def test_invalid_arg_type_forcing(self):
        def foo(iters):
            a = range(iters)
            return iters

        args = (u32,)
        return_type = u8
        cfunc = njit(return_type(*args))(foo)
        cres = cfunc.overloads[args]
        typemap = cres.type_annotation.typemap
        # Argument "iters" must be uint32
        self.assertEqual(typemap['iters'], u32)


class TestUnify(unittest.TestCase):
    """
    Tests for type unification with a typing context.
    """

    int_unify = {
        ('uint8', 'uint8'): 'uint8',
        ('int8', 'int8'): 'int8',
        ('uint16', 'uint16'): 'uint16',
        ('int16', 'int16'): 'int16',
        ('uint32', 'uint32'): 'uint32',
        ('int32', 'int32'): 'int32',
        ('uint64', 'uint64'): 'uint64',
        ('int64', 'int64'): 'int64',

        ('int8', 'uint8'): 'int16',
        ('int8', 'uint16'): 'int32',
        ('int8', 'uint32'): 'int64',

        ('uint8', 'int32'): 'int32',
        ('uint8', 'uint64'): 'uint64',

        ('int16', 'int8'): 'int16',
        ('int16', 'uint8'): 'int16',
        ('int16', 'uint16'): 'int32',
        ('int16', 'uint32'): 'int64',
        ('int16', 'int64'): 'int64',
        ('int16', 'uint64'): 'float64',

        ('uint16', 'uint8'): 'uint16',
        ('uint16', 'uint32'): 'uint32',
        ('uint16', 'int32'): 'int32',
        ('uint16', 'uint64'): 'uint64',

        ('int32', 'int8'): 'int32',
        ('int32', 'int16'): 'int32',
        ('int32', 'uint32'): 'int64',
        ('int32', 'int64'): 'int64',

        ('uint32', 'uint8'): 'uint32',
        ('uint32', 'int64'): 'int64',
        ('uint32', 'uint64'): 'uint64',

        ('int64', 'int8'): 'int64',
        ('int64', 'uint8'): 'int64',
        ('int64', 'uint16'): 'int64',

        ('uint64', 'int8'): 'float64',
        ('uint64', 'int32'): 'float64',
        ('uint64', 'int64'): 'float64',
    }

    def assert_unify(self, aty, bty, expected):
        ctx = typing.Context()
        template = "{0}, {1} -> {2} != {3}"
        for unify_func in ctx.unify_types, ctx.unify_pairs:
            unified = unify_func(aty, bty)
            self.assertEqual(unified, expected,
                             msg=template.format(aty, bty, unified, expected))
            unified = unify_func(bty, aty)
            self.assertEqual(unified, expected,
                             msg=template.format(bty, aty, unified, expected))

    def assert_unify_failure(self, aty, bty):
        self.assert_unify(aty, bty, None)

    def test_integer(self):
        ctx = typing.Context()
        for aty, bty in itertools.product(types.integer_domain,
                                          types.integer_domain):
            key = (str(aty), str(bty))
            try:
                expected = self.int_unify[key]
            except KeyError:
                expected = self.int_unify[key[::-1]]
            self.assert_unify(aty, bty, getattr(types, expected))

    def test_bool(self):
        aty = types.boolean
        for bty in types.integer_domain:
            self.assert_unify(aty, bty, bty)
        # Not sure about this one, but it respects transitivity
        for cty in types.real_domain:
            self.assert_unify(aty, cty, cty)

    def unify_number_pair_test(self, n):
        """
        Test all permutations of N-combinations of numeric types and ensure
        that the order of types in the sequence is irrelevant.
        """
        ctx = typing.Context()
        for tys in itertools.combinations(types.number_domain, n):
            res = [ctx.unify_types(*comb)
                   for comb in itertools.permutations(tys)]
            first_result = res[0]
            # Sanity check
            self.assertIsInstance(first_result, types.Number)
            # All results must be equal
            for other in res[1:]:
                self.assertEqual(first_result, other)

    def test_unify_number_pair(self):
        self.unify_number_pair_test(2)
        self.unify_number_pair_test(3)

    def test_none_to_optional(self):
        """
        Test unification of `none` and multiple number types to optional type
        """
        ctx = typing.Context()
        for tys in itertools.combinations(types.number_domain, 2):
            # First unify without none, to provide the control value
            tys = list(tys)
            expected = types.Optional(ctx.unify_types(*tys))
            results = [ctx.unify_types(*comb)
                       for comb in itertools.permutations(tys  + [types.none])]
            # All results must be equal
            for res in results:
                self.assertEqual(res, expected)

    def test_none(self):
        aty = types.none
        bty = types.none
        self.assert_unify(aty, bty, types.none)

    def test_optional(self):
        aty = types.Optional(i32)
        bty = types.none
        self.assert_unify(aty, bty, aty)
        aty = types.Optional(i32)
        bty = types.Optional(i64)
        self.assert_unify(aty, bty, bty)
        aty = types.Optional(i32)
        bty = i64
        self.assert_unify(aty, bty, types.Optional(i64))
        # Failure
        aty = types.Optional(i32)
        bty = types.Optional(types.slice3_type)
        self.assert_unify_failure(aty, bty)

    def test_tuple(self):
        aty = types.UniTuple(i32, 3)
        bty = types.UniTuple(i64, 3)
        self.assert_unify(aty, bty, types.UniTuple(i64, 3))
        # (Tuple, UniTuple) -> Tuple
        aty = types.UniTuple(i32, 2)
        bty = types.Tuple((i16, i64))
        self.assert_unify(aty, bty, types.Tuple((i32, i64)))
        aty = types.UniTuple(i64, 0)
        bty = types.Tuple(())
        self.assert_unify(aty, bty, bty)
        # (Tuple, Tuple) -> Tuple
        aty = types.Tuple((i8, i16, i32))
        bty = types.Tuple((i32, i16, i8))
        self.assert_unify(aty, bty, types.Tuple((i32, i16, i32)))
        aty = types.Tuple((i8, i32))
        bty = types.Tuple((i32, i8))
        self.assert_unify(aty, bty, types.Tuple((i32, i32)))
        aty = types.Tuple((i8, i16))
        bty = types.Tuple((i16, i8))
        self.assert_unify(aty, bty, types.Tuple((i16, i16)))
        # Different number kinds
        aty = types.UniTuple(f64, 3)
        bty = types.UniTuple(c64, 3)
        self.assert_unify(aty, bty, types.UniTuple(c128, 3))
        # Tuples of tuples
        aty = types.UniTuple(types.Tuple((u32, f32)), 2)
        bty = types.UniTuple(types.Tuple((i16, f32)), 2)
        self.assert_unify(aty, bty,
                          types.UniTuple(types.Tuple((i64, f32)), 2))
        # Failures
        aty = types.UniTuple(i32, 1)
        bty = types.UniTuple(types.slice3_type, 1)
        self.assert_unify_failure(aty, bty)
        aty = types.UniTuple(i32, 1)
        bty = types.UniTuple(i32, 2)
        self.assert_unify_failure(aty, bty)
        aty = types.Tuple((i8, types.slice3_type))
        bty = types.Tuple((i32, i8))
        self.assert_unify_failure(aty, bty)

    def test_optional_tuple(self):
        # Unify to optional tuple
        aty = types.none
        bty = types.UniTuple(i32, 2)
        self.assert_unify(aty, bty, types.Optional(types.UniTuple(i32, 2)))
        aty = types.Optional(types.UniTuple(i16, 2))
        bty = types.UniTuple(i32, 2)
        self.assert_unify(aty, bty, types.Optional(types.UniTuple(i32, 2)))
        # Unify to tuple of optionals
        aty = types.Tuple((types.none, i32))
        bty = types.Tuple((i16, types.none))
        self.assert_unify(aty, bty, types.Tuple((types.Optional(i16),
                                                 types.Optional(i32))))
        aty = types.Tuple((types.Optional(i32), i64))
        bty = types.Tuple((i16, types.Optional(i8)))
        self.assert_unify(aty, bty, types.Tuple((types.Optional(i32),
                                                 types.Optional(i64))))

    def test_arrays(self):
        aty = types.Array(i32, 3, "C")
        bty = types.Array(i32, 3, "A")
        self.assert_unify(aty, bty, bty)
        aty = types.Array(i32, 3, "C")
        bty = types.Array(i32, 3, "F")
        self.assert_unify(aty, bty, types.Array(i32, 3, "A"))
        aty = types.Array(i32, 3, "C")
        bty = types.Array(i32, 3, "C", readonly=True)
        self.assert_unify(aty, bty, bty)
        aty = types.Array(i32, 3, "A")
        bty = types.Array(i32, 3, "C", readonly=True)
        self.assert_unify(aty, bty,
                          types.Array(i32, 3, "A", readonly=True))
        # Failures
        aty = types.Array(i32, 2, "C")
        bty = types.Array(i32, 3, "C")
        self.assert_unify_failure(aty, bty)
        aty = types.Array(i32, 2, "C")
        bty = types.Array(u32, 2, "C")
        self.assert_unify_failure(aty, bty)

    def test_list(self):
        aty = types.List(types.undefined)
        bty = types.List(i32)
        self.assert_unify(aty, bty, bty)
        aty = types.List(i16)
        bty = types.List(i32)
        self.assert_unify(aty, bty, bty)
        aty = types.List(types.Tuple([i32, i16]))
        bty = types.List(types.Tuple([i16, i64]))
        cty = types.List(types.Tuple([i32, i64]))
        self.assert_unify(aty, bty, cty)
        # Different reflections
        aty = types.List(i16, reflected=True)
        bty = types.List(i32)
        cty = types.List(i32, reflected=True)
        self.assert_unify(aty, bty, cty)
        # Incompatible dtypes
        aty = types.List(i16)
        bty = types.List(types.Tuple([i16]))
        self.assert_unify_failure(aty, bty)

    def test_set(self):
        # Different reflections
        aty = types.Set(i16, reflected=True)
        bty = types.Set(i32)
        cty = types.Set(i32, reflected=True)
        self.assert_unify(aty, bty, cty)
        # Incompatible dtypes
        aty = types.Set(i16)
        bty = types.Set(types.Tuple([i16]))
        self.assert_unify_failure(aty, bty)

    def test_range(self):
        aty = types.range_state32_type
        bty = types.range_state64_type
        self.assert_unify(aty, bty, bty)


class TestTypeConversion(CompatibilityTestMixin, unittest.TestCase):
    """
    Test for conversion between types with a typing context.
    """

    def assert_can_convert(self, aty, bty, expected):
        ctx = typing.Context()
        got = ctx.can_convert(aty, bty)
        self.assertEqual(got, expected)

    def assert_cannot_convert(self, aty, bty):
        ctx = typing.Context()
        got = ctx.can_convert(aty, bty)
        self.assertIsNone(got)

    def test_convert_number_types(self):
        # Check that Context.can_convert() is compatible with the default
        # number conversion rules registered in the typeconv module
        # (which is used internally by the C _Dispatcher object).
        ctx = typing.Context()
        self.check_number_compatibility(ctx.can_convert)

    def test_tuple(self):
        # UniTuple -> UniTuple
        aty = types.UniTuple(i32, 3)
        bty = types.UniTuple(i64, 3)
        self.assert_can_convert(aty, aty, Conversion.exact)
        self.assert_can_convert(aty, bty, Conversion.promote)
        aty = types.UniTuple(i32, 3)
        bty = types.UniTuple(f64, 3)
        self.assert_can_convert(aty, bty, Conversion.safe)
        # Tuple -> Tuple
        aty = types.Tuple((i32, i32))
        bty = types.Tuple((i32, i64))
        self.assert_can_convert(aty, bty, Conversion.promote)
        # UniTuple <-> Tuple
        aty = types.UniTuple(i32, 2)
        bty = types.Tuple((i32, i64))
        self.assert_can_convert(aty, bty, Conversion.promote)
        self.assert_can_convert(bty, aty, Conversion.unsafe)
        # Empty tuples
        aty = types.UniTuple(i64, 0)
        bty = types.UniTuple(i32, 0)
        cty = types.Tuple(())
        self.assert_can_convert(aty, bty, Conversion.safe)
        self.assert_can_convert(bty, aty, Conversion.safe)
        self.assert_can_convert(aty, cty, Conversion.safe)
        self.assert_can_convert(cty, aty, Conversion.safe)
        # Failures
        aty = types.UniTuple(i64, 3)
        bty = types.UniTuple(types.none, 3)
        self.assert_cannot_convert(aty, bty)
        aty = types.UniTuple(i64, 2)
        bty = types.UniTuple(i64, 3)

    def test_arrays(self):
        # Different layouts
        aty = types.Array(i32, 3, "C")
        bty = types.Array(i32, 3, "A")
        self.assert_can_convert(aty, bty, Conversion.safe)
        aty = types.Array(i32, 2, "C")
        bty = types.Array(i32, 2, "F")
        self.assert_cannot_convert(aty, bty)
        # Different mutabilities
        aty = types.Array(i32, 3, "C")
        bty = types.Array(i32, 3, "C", readonly=True)
        self.assert_can_convert(aty, aty, Conversion.exact)
        self.assert_can_convert(bty, bty, Conversion.exact)
        self.assert_can_convert(aty, bty, Conversion.safe)
        self.assert_cannot_convert(bty, aty)
        # Various failures
        aty = types.Array(i32, 2, "C")
        bty = types.Array(i32, 3, "C")
        self.assert_cannot_convert(aty, bty)
        aty = types.Array(i32, 2, "C")
        bty = types.Array(i64, 2, "C")
        self.assert_cannot_convert(aty, bty)

    def test_optional(self):
        aty = types.int32
        bty = types.Optional(i32)
        self.assert_can_convert(types.none, bty, Conversion.promote)
        self.assert_can_convert(aty, bty, Conversion.promote)
        self.assert_cannot_convert(bty, types.none)
        self.assert_can_convert(bty, aty, Conversion.safe)  # XXX ???
        # Optional array
        aty = types.Array(i32, 2, "C")
        bty = types.Optional(aty)
        self.assert_can_convert(types.none, bty, Conversion.promote)
        self.assert_can_convert(aty, bty, Conversion.promote)
        self.assert_can_convert(bty, aty, Conversion.safe)
        aty = types.Array(i32, 2, "C")
        bty = types.Optional(aty.copy(layout="A"))
        self.assert_can_convert(aty, bty, Conversion.safe)  # C -> A
        self.assert_cannot_convert(bty, aty)                # A -> C
        aty = types.Array(i32, 2, "C")
        bty = types.Optional(aty.copy(layout="F"))
        self.assert_cannot_convert(aty, bty)
        self.assert_cannot_convert(bty, aty)


class TestResolveOverload(unittest.TestCase):
    """
    Tests for typing.Context.resolve_overload().
    """

    def assert_resolve_overload(self, cases, args, expected):
        ctx = typing.Context()
        got = ctx.resolve_overload("foo", cases, args, {})
        self.assertEqual(got, expected)

    def test_non_ambiguous_match(self):
        def check(args, expected):
            self.assert_resolve_overload(cases, args, expected)
            # Order shouldn't matter here
            self.assert_resolve_overload(cases[::-1], args, expected)

        cases = [i8(i8, i8), i32(i32, i32), f64(f64, f64)]
        # Exact match
        check((i8, i8), cases[0])
        check((i32, i32), cases[1])
        check((f64, f64), cases[2])
        # "Promote" conversion
        check((i8, i16), cases[1])
        check((i32, i8), cases[1])
        check((i32, i8), cases[1])
        check((f32, f32), cases[2])
        # "Safe" conversion
        check((u32, u32), cases[2])
        # "Unsafe" conversion
        check((i64, i64), cases[2])

    def test_ambiguous_match(self):
        # When the best match is ambiguous (there is a tie), the first
        # best case in original sequence order should be returned.
        def check(args, expected, expected_reverse):
            self.assert_resolve_overload(cases, args, expected)
            self.assert_resolve_overload(cases[::-1], args, expected_reverse)

        cases = [i16(i16, i16), i32(i32, i32), f64(f64, f64)]
        # Two "promote" conversions
        check((i8, i8), cases[0], cases[1])
        # Two "safe" conversions
        check((u16, u16), cases[1], cases[2])

        cases = [i32(i32, i32), f32(f32, f32)]
        # Two "unsafe" conversions
        check((u32, u32), cases[0], cases[1])

    def test_ambiguous_error(self):
        ctx = typing.Context()
        cases = [i16(i16, i16), i32(i32, i32)]
        with self.assertRaises(TypeError) as raises:
            ctx.resolve_overload("foo", cases, (i8, i8), {},
                                 allow_ambiguous=False)
        self.assertEqual(str(raises.exception).splitlines(),
                         ["Ambiguous overloading for foo (int8, int8):",
                          "(int16, int16) -> int16",
                          "(int32, int32) -> int32",
                          ])


class TestUnifyUseCases(unittest.TestCase):
    """
    Concrete cases where unification would fail.
    """

    @staticmethod
    def _actually_test_complex_unify():
        def pyfunc(a):
            res = 0.0
            for i in range(len(a)):
                res += a[i]
            return res

        argtys = (types.Array(c128, 1, 'C'),)
        cfunc = njit(argtys)(pyfunc)
        return (pyfunc, cfunc)

    def test_complex_unify_issue599(self):
        pyfunc, cfunc = self._actually_test_complex_unify()
        arg = np.array([1.0j])
        self.assertEqual(cfunc(arg), pyfunc(arg))

    def test_complex_unify_issue599_multihash(self):
        """
        Test issue #599 for multiple values of PYTHONHASHSEED.
        """
        env = os.environ.copy()
        for seedval in (1, 2, 1024):
            env['PYTHONHASHSEED'] = str(seedval)
            subproc = subprocess.Popen(
                [sys.executable, '-c',
                 'import numba.tests.test_typeinfer as test_mod\n' +
                 'test_mod.TestUnifyUseCases._actually_test_complex_unify()'],
                env=env)
            subproc.wait()
            self.assertEqual(subproc.returncode, 0, 'Child process failed.')

    def test_int_tuple_unify(self):
        """
        Test issue #493
        """
        def foo(an_int32, an_int64):
            a = an_int32, an_int32
            while True:  # infinite loop
                a = an_int32, an_int64
            return a

        args = (i32, i64)
        # Check if compilation is successful
        njit(args)(foo)


def issue_797(x0, y0, x1, y1, grid):
    nrows, ncols = grid.shape

    dx = abs(x1 - x0)
    dy = abs(y1 - y0)

    sx = 0
    if x0 < x1:
        sx = 1
    else:
        sx = -1
    sy = 0
    if y0 < y1:
        sy = 1
    else:
        sy = -1

    err = dx - dy

    while True:
        if x0 == x1 and y0 == y1:
            break

        if 0 <= x0 < nrows and 0 <= y0 < ncols:
            grid[x0, y0] += 1

        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x0 += sx
        if e2 < dx:
            err += dx
            y0 += sy


def issue_1080(a, b):
    if not a:
        return True
    return b


def list_unify_usecase1(n):
    res = 0
    x = []
    if n < 10:
        x.append(np.int32(n))
    else:
        for i in range(n):
            x.append(np.int64(i))
    x.append(5.0)

    # Note `i` and `j` may have different types (int64 vs. int32)
    for j in range(len(x)):
        res += j * x[j]
    for val in x:
        res += int(val) & len(x)
    while len(x) > 0:
        res += x.pop()
    return res

def list_unify_usecase2(n):
    res = []
    for i in range(n):
        if i & 1:
            res.append((i, 1.0))
        else:
            res.append((2.0, i))
    res.append((123j, 42))
    return res

def range_unify_usecase(v):
    if v:
        r = range(np.int32(3))
    else:
        r = range(np.int64(5))
    for x in r:
        return x

def issue_1394(a):
    if a:
        for i in range(a):
            a += i
        i = 1.2
    else:
        i = 3
    return a, i


class TestMiscIssues(TestCase):

    def test_issue_797(self):
        """https://github.com/numba/numba/issues/797#issuecomment-58592401

        Undeterministic triggering of tuple coercion error
        """
        foo = jit(nopython=True)(issue_797)
        g = np.zeros(shape=(10, 10), dtype=np.int32)
        foo(np.int32(0), np.int32(0), np.int32(1), np.int32(1), g)

    def test_issue_1080(self):
        """https://github.com/numba/numba/issues/1080

        Erroneous promotion of boolean args to int64
        """
        foo = jit(nopython=True)(issue_1080)
        foo(True, False)

    def test_list_unify1(self):
        """
        Exercise back-propagation of refined list type.
        """
        pyfunc = list_unify_usecase1
        cfunc = jit(nopython=True)(pyfunc)
        for n in [5, 100]:
            res = cfunc(n)
            self.assertPreciseEqual(res, pyfunc(n))

    def test_list_unify2(self):
        pyfunc = list_unify_usecase2
        cfunc = jit(nopython=True)(pyfunc)
        res = cfunc(3)
        # NOTE: the types will differ (Numba returns a homogeneous list with
        # converted values).
        self.assertEqual(res, pyfunc(3))

    def test_range_unify(self):
        pyfunc = range_unify_usecase
        cfunc = jit(nopython=True)(pyfunc)
        for v in (0, 1):
            res = cfunc(v)
            self.assertPreciseEqual(res, pyfunc(v))

    def test_issue_1394(self):
        pyfunc = issue_1394
        cfunc = jit(nopython=True)(pyfunc)
        for v in (0, 1, 2):
            res = cfunc(v)
            self.assertEqual(res, pyfunc(v))

    def test_issue_6293(self):
        """https://github.com/numba/numba/issues/6293

        Typer does not propagate return type to all return variables
        """
        @jit(nopython=True)
        def confuse_typer(x):
            if x == x:
                return int(x)
            else:
                return x

        confuse_typer.compile((types.float64,))
        cres = confuse_typer.overloads[(types.float64,)]
        typemap = cres.type_annotation.typemap
        return_vars = {}

        for block in cres.type_annotation.blocks.values():
            for inst in block.body:
                if isinstance(inst, ir.Return):
                    varname = inst.value.name
                    return_vars[varname] = typemap[varname]

        self.assertTrue(all(vt == types.float64 for vt in return_vars.values()))

    def test_issue_9162(self):
        @overload_method(types.Array, "aabbcc")
        def ol_aabbcc(self):

            def impl(self):
                return self.sum()

            return impl

        @jit
        def foo(ar):
            return ar.aabbcc()

        ar = np.ones(2)
        ret = foo(ar)

        overload = [value for value in foo.overloads.values()][0]
        typemap = overload.type_annotation.typemap
        calltypes = overload.type_annotation.calltypes
        for call_op in calltypes:
            name = call_op.list_vars()[0].name
            fc_ty = typemap[name]
            self.assertIsInstance(fc_ty, types.BoundFunction)
            tmplt = fc_ty.template
            info = tmplt.get_template_info(tmplt)
            py_file = info["filename"]
            self.assertIn("test_typeinfer.py", py_file)

    @skip_unless_load_fast_and_clear
    def test_load_fast_and_clear(self):
        @njit
        def foo(a):
            [x for x in (0,)]
            if a:
                # the test code cannot use a constant here due to constant
                # propagation issues.
                x = 3 + a
            x += 10
            return x

        self.assertEqual(foo(True), foo.py_func(True))
        # Interpreted version should raise an exception
        with self.assertRaises(UnboundLocalError):
            foo.py_func(False)
        # Compiled version returns 10 as x is zero initialized.
        self.assertEqual(foo(False), 10)

    @skip_unless_load_fast_and_clear
    def test_load_fast_and_clear_variant_2(self):
        @njit
        def foo():
            # The use of a literal False triggers different bytecode generation
            # necessary for this test. See test_load_fast_and_clear_variant_4.
            if False:
                x = 1
            [x for x in (1,)]
            # This return uses undefined variable
            return x

        with self.assertRaises(errors.TypingError) as raises:
            foo()
        self.assertIn("return value is undefined", str(raises.exception))

    @skip_unless_load_fast_and_clear
    def test_load_fast_and_clear_variant_3(self):
        @njit
        def foo():
            # The use of a literal False triggers different bytecode generation
            # necessary for this test. See test_load_fast_and_clear_variant_4.
            if False:
                x = 1
            [x for x in (1,)]
            # This print uses undefined variable
            print(1, 2, 3, x)

        with self.assertRaises(errors.TypingError) as raises:
            foo()
        self.assertIn("undefined variable used in call argument #4", str(raises.exception))

    @skip_unless_load_fast_and_clear
    def test_load_fast_and_clear_variant_4(self):
        @njit
        def foo(a):
            # This test variant is to show that non-literal boolean value here
            # produces a different behavior.
            if a:
                x = a
            [x for x in (1,)]
            return x
        self.assertEqual(foo(123), 123)
        self.assertEqual(foo(0), 0)


class TestFoldArguments(unittest.TestCase):
    def check_fold_arguments_list_inputs(self, func, args, kws):
        def make_tuple(*args):
            return args

        unused_handler = None

        pysig = utils.pysignature(func)
        names = list(pysig.parameters)

        with self.subTest(kind='dict'):
            folded_dict = typing.fold_arguments(
                pysig, args, kws, make_tuple, unused_handler, unused_handler,
            )
            # correct ordering
            for i, (j, k) in enumerate(zip(folded_dict, names)):
                (got_index, got_param, got_name) = j
                self.assertEqual(got_index, i)
                self.assertEqual(got_name, f'arg.{k}')

        kws = list(kws.items())
        with self.subTest(kind='list'):
            folded_list = typing.fold_arguments(
                pysig, args, kws, make_tuple, unused_handler, unused_handler,
            )
            self.assertEqual(folded_list, folded_dict)

    def test_fold_arguments_list_inputs(self):
        cases = [
            dict(
                func=lambda a, b, c, d: None,
                args=['arg.a', 'arg.b'],
                kws=dict(c='arg.c', d='arg.d')
            ),
            dict(
                func=lambda: None,
                args=[],
                kws=dict(),
            ),
            dict(
                func=lambda a: None,
                args=['arg.a'],
                kws={},
            ),
            dict(
                func=lambda a: None,
                args=[],
                kws=dict(a='arg.a'),
            ),
        ]
        for case in cases:
            with self.subTest(**case):
                self.check_fold_arguments_list_inputs(**case)


@register_pass(mutates_CFG=False, analysis_only=True)
class DummyCR(FunctionPass):
    """Dummy pass to add "cr" to compiler state to avoid errors in TyperCompiler since
    it doesn't have lowering.
    """

    _name = "dummy_cr"

    def __init__(self):
        FunctionPass.__init__(self)

    def run_pass(self, state):
        state.cr = 1  # arbitrary non-None value
        return True


class TyperCompiler(numba.core.compiler.CompilerBase):
    """A compiler pipeline that skips passes after typing (provides partial typing info
    but not lowering).
    """

    def define_pipelines(self):
        pm = numba.core.compiler_machinery.PassManager("custom_pipeline")
        pm.add_pass(TranslateByteCode, "analyzing bytecode")
        pm.add_pass(IRProcessing, "processing IR")
        pm.add_pass(PartialTypeInference, "do partial typing")
        pm.add_pass_after(DummyCR, PartialTypeInference)
        pm.finalize()
        return [pm]


def get_func_typing_errs(func, arg_types):
    """
    Get typing errors for function 'func'. It creates a pipeline that runs untyped
    passes as well as type inference.
    """
    typingctx = numba.core.registry.cpu_target.typing_context
    targetctx = numba.core.registry.cpu_target.target_context
    library = None
    return_type = None
    _locals = {}
    flags = numba.core.compiler.Flags()
    flags.nrt = True

    pipeline = TyperCompiler(
        typingctx, targetctx, library, arg_types, return_type, flags, _locals
    )
    pipeline.compile_extra(func)
    return pipeline.state.typing_errors


class TestPartialTypingErrors(unittest.TestCase):
    """
    Make sure partial typing stores type errors in compiler state properly
    """
    def test_partial_typing_error(self):
        # example with type unification error
        def impl(flag):
            if flag:
                a = 1
            else:
                a = str(1)
            return a

        typing_errs = get_func_typing_errs(impl, (types.bool_,))
        self.assertTrue(isinstance(typing_errs, list) and len(typing_errs) == 1)
        self.assertTrue(isinstance(typing_errs[0], errors.TypingError) and
                        "Cannot unify" in typing_errs[0].msg)


if __name__ == '__main__':
    unittest.main()
