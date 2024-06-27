import collections
import itertools

import numpy as np

from numba import njit, jit, typeof, literally
from numba.core import types, errors, utils
from numba.tests.support import TestCase, MemoryLeakMixin, tag
import unittest


Rect = collections.namedtuple('Rect', ('width', 'height'))

Point = collections.namedtuple('Point', ('x', 'y', 'z'))

Point2 = collections.namedtuple('Point2', ('x', 'y', 'z'))

Empty = collections.namedtuple('Empty', ())

def tuple_return_usecase(a, b):
    return a, b

def tuple_first(tup):
    a, b = tup
    return a

def tuple_second(tup):
    a, b = tup
    return b

def tuple_index(tup, idx):
    return tup[idx]

def tuple_index_static(tup):
    # Note the negative index
    return tup[-2]

def tuple_slice2(tup):
    return tup[1:-1]

def tuple_slice3(tup):
    return tup[1::2]

def len_usecase(tup):
    return len(tup)

def add_usecase(a, b):
    return a + b

def eq_usecase(a, b):
    return a == b

def ne_usecase(a, b):
    return a != b

def gt_usecase(a, b):
    return a > b

def ge_usecase(a, b):
    return a >= b

def lt_usecase(a, b):
    return a < b

def le_usecase(a, b):
    return a <= b

def in_usecase(a, b):
    return a in b

def bool_usecase(tup):
    return bool(tup), (3 if tup else 2)

def getattr_usecase(tup):
    return tup.z, tup.y, tup.x

def make_point(a, b, c):
    return Point(a, b, c)

def make_point_kws(a, b, c):
    return Point(z=c, y=b, x=a)

def make_point_nrt(n):
    r = Rect(list(range(n)), np.zeros(n + 1))
    # This also exercises attribute access
    p = Point(r, len(r.width), len(r.height))
    return p

def type_usecase(tup, *args):
    return type(tup)(*args)

def identity(tup):
    return tup

def index_method_usecase(tup, value):
    return tup.index(value)

def tuple_unpack_static_getitem_err():
    # see issue3895, `c` is imprecise
    a, b, c, d = [], [], [], 0.0
    a.append(1)
    b.append(1)
    return


class TestTupleLengthError(unittest.TestCase):

    def test_tuple_length_error(self):
        # issue 2195
        # raise an error on tuples greater than 1000 in length
        @njit
        def eattuple(tup):
            return len(tup)

        with self.assertRaises(errors.UnsupportedError) as raises:
            tup = tuple(range(1001))
            eattuple(tup)

        expected = "Tuple 'tup' length must be smaller than 1000"
        self.assertIn(expected, str(raises.exception))

class TestTupleTypeNotIterable(unittest.TestCase):
    '''
    issue 4369
    raise an error if 'type' is not iterable
    '''
    def test_namedtuple_types_exception(self):
        with self.assertRaises(errors.TypingError) as raises:
            types.NamedTuple(types.uint32, 'p')
        self.assertIn(
            "Argument 'types' is not iterable",
            str(raises.exception)
        )

    def test_tuple_types_exception(self):
        with self.assertRaises(errors.TypingError) as raises:
            types.Tuple((types.uint32))
        self.assertIn(
            "Argument 'types' is not iterable",
            str(raises.exception)
        )


class TestTupleReturn(TestCase):

    def test_array_tuple(self):
        aryty = types.Array(types.float64, 1, 'C')
        cfunc = njit((aryty, aryty))(tuple_return_usecase)
        a = b = np.arange(5, dtype='float64')
        ra, rb = cfunc(a, b)
        self.assertPreciseEqual(ra, a)
        self.assertPreciseEqual(rb, b)
        del a, b
        self.assertPreciseEqual(ra, rb)

    def test_scalar_tuple(self):
        scalarty = types.float32
        cfunc = njit((scalarty, scalarty))(tuple_return_usecase)
        a = b = 1
        ra, rb = cfunc(a, b)
        self.assertEqual(ra, a)
        self.assertEqual(rb, b)

    def test_hetero_tuple(self):
        alltypes = []
        allvalues = []

        alltypes.append((types.int32, types.int64))
        allvalues.append((1, 2))

        alltypes.append((types.float32, types.float64))
        allvalues.append((1.125, .25))

        alltypes.append((types.int32, types.float64))
        allvalues.append((1231, .5))

        for (ta, tb), (a, b) in zip(alltypes, allvalues):
            cfunc = njit((ta, tb))(tuple_return_usecase)
            ra, rb = cfunc(a, b)
            self.assertPreciseEqual((ra, rb), (a, b))


class TestTuplePassing(TestCase):

    def test_unituple(self):
        tuple_type = types.UniTuple(types.int32, 2)
        cf_first = njit((tuple_type,))(tuple_first)
        cf_second = njit((tuple_type,))(tuple_second)
        self.assertPreciseEqual(cf_first((4, 5)), 4)
        self.assertPreciseEqual(cf_second((4, 5)), 5)

    def test_hetero_tuple(self):
        tuple_type = types.Tuple((types.int64, types.float32))
        cf_first = njit((tuple_type,))(tuple_first)
        cf_second = njit((tuple_type,))(tuple_second)
        self.assertPreciseEqual(cf_first((2**61, 1.5)), 2**61)
        self.assertPreciseEqual(cf_second((2**61, 1.5)), 1.5)

    def test_size_mismatch(self):
        # Issue #1638: tuple size should be checked when unboxing
        tuple_type = types.UniTuple(types.int32, 2)
        cfunc = njit((tuple_type,))(tuple_first)
        entry_point = cfunc.overloads[cfunc.signatures[0]].entry_point
        with self.assertRaises(ValueError) as raises:
            entry_point((4, 5, 6))
        self.assertEqual(str(raises.exception),
                         ("size mismatch for tuple, "
                          "expected 2 element(s) but got 3"))


class TestOperations(TestCase):

    def test_len(self):
        pyfunc = len_usecase
        cfunc = njit((types.Tuple((types.int64, types.float32)),))(pyfunc)
        self.assertPreciseEqual(cfunc((4, 5)), 2)
        cfunc = njit((types.UniTuple(types.int64, 3),))(pyfunc)
        self.assertPreciseEqual(cfunc((4, 5, 6)), 3)

    def test_index_literal(self):
        # issue #6023, test non-static getitem with IntegerLiteral index
        def pyfunc(tup, idx):
            idx = literally(idx)
            return tup[idx]
        cfunc = njit(pyfunc)

        tup = (4, 3.1, 'sss')
        for i in range(len(tup)):
            self.assertPreciseEqual(cfunc(tup, i), tup[i])

    def test_index(self):
        pyfunc = tuple_index
        cfunc = njit((types.UniTuple(types.int64, 3), types.int64),)(pyfunc)
        tup = (4, 3, 6)
        for i in range(len(tup)):
            self.assertPreciseEqual(cfunc(tup, i), tup[i])

        # test negative indexing
        for i in range(len(tup) + 1):
            self.assertPreciseEqual(cfunc(tup, -i), tup[-i])

        # oob indexes, +ve then -ve
        with self.assertRaises(IndexError) as raises:
            cfunc(tup, len(tup))
        self.assertEqual("tuple index out of range", str(raises.exception))
        with self.assertRaises(IndexError) as raises:
            cfunc(tup, -(len(tup) + 1))
        self.assertEqual("tuple index out of range", str(raises.exception))

        # Test empty tuple, this is a bit unusual as `njit` will infer the empty
        # tuple arg as a types.Tuple and not match the compiled signature, this
        # is essentially because the test originally relied on
        # `compile_isolated`.
        args = (types.UniTuple(types.int64, 0), types.int64,)
        cr = njit(args)(pyfunc).overloads[args]
        with self.assertRaises(IndexError) as raises:
            cr.entry_point((), 0)
        self.assertEqual("tuple index out of range", str(raises.exception))

        # test uintp indexing (because, e.g., parfor generates unsigned prange)
        cfunc = njit((types.UniTuple(types.int64, 3), types.uintp,),)(pyfunc)
        for i in range(len(tup)):
            self.assertPreciseEqual(cfunc(tup, types.uintp(i)), tup[i])

        # With a compile-time static index (the code generation path is
        # different)
        pyfunc = tuple_index_static
        for typ in (types.UniTuple(types.int64, 4),
                    types.Tuple((types.int64, types.int32, types.int64, types.int32))):
            cfunc = njit((typ,))(pyfunc)
            tup = (4, 3, 42, 6)
            self.assertPreciseEqual(cfunc(tup), pyfunc(tup))

        typ = types.UniTuple(types.int64, 1)
        with self.assertTypingError():
            njit((typ,))(pyfunc)

        # test unpack, staticgetitem with imprecise type (issue #3895)
        pyfunc = tuple_unpack_static_getitem_err
        with self.assertTypingError() as raises:
            njit((),)(pyfunc)
        msg = ("Cannot infer the type of variable 'c', have imprecise type: "
               "list(undefined)<iv=None>.")
        self.assertIn(msg, str(raises.exception))

    def test_in(self):
        pyfunc = in_usecase
        cfunc = njit((types.int64, types.UniTuple(types.int64, 3),),)(pyfunc)
        tup = (4, 1, 5)
        for i in range(5):
            self.assertPreciseEqual(cfunc(i, tup), pyfunc(i, tup))

        # Test the empty case
        cfunc = njit((types.int64, types.Tuple([]),),)(pyfunc)
        self.assertPreciseEqual(cfunc(1, ()), pyfunc(1, ()))

    def check_slice(self, pyfunc):
        tup = (4, 5, 6, 7)
        cfunc = njit((types.UniTuple(types.int64, 4),),)(pyfunc)
        self.assertPreciseEqual(cfunc(tup), pyfunc(tup))
        args = types.Tuple((types.int64, types.int32, types.int64, types.int32))
        cfunc = njit((args,))(pyfunc)
        self.assertPreciseEqual(cfunc(tup), pyfunc(tup))

    def test_slice2(self):
        self.check_slice(tuple_slice2)

    def test_slice3(self):
        self.check_slice(tuple_slice3)

    def test_bool(self):
        pyfunc = bool_usecase
        cfunc = njit((types.Tuple((types.int64, types.int32)),),)(pyfunc)
        args = ((4, 5),)
        self.assertPreciseEqual(cfunc(*args), pyfunc(*args))
        cfunc = njit((types.UniTuple(types.int64, 3),),)(pyfunc)
        args = ((4, 5, 6),)
        self.assertPreciseEqual(cfunc(*args), pyfunc(*args))
        cfunc = njit((types.Tuple(()),),)(pyfunc)
        self.assertPreciseEqual(cfunc(()), pyfunc(()))

    def test_add(self):
        pyfunc = add_usecase
        samples = [(types.Tuple(()), ()),
                   (types.UniTuple(types.int32, 0), ()),
                   (types.UniTuple(types.int32, 1), (42,)),
                   (types.Tuple((types.int64, types.float32)), (3, 4.5)),
                   ]
        for (ta, a), (tb, b) in itertools.product(samples, samples):
            cfunc = njit((ta, tb),)(pyfunc)
            expected = pyfunc(a, b)
            got = cfunc(a, b)
            self.assertPreciseEqual(got, expected, msg=(ta, tb))

    def _test_compare(self, pyfunc):
        def eq(pyfunc, cfunc, args):
            self.assertIs(cfunc(*args), pyfunc(*args),
                          "mismatch for arguments %s" % (args,))

        # Same-sized tuples
        argtypes = [types.Tuple((types.int64, types.float32)),
                    types.UniTuple(types.int32, 2)]
        for ta, tb in itertools.product(argtypes, argtypes):
            cfunc = njit((ta, tb),)(pyfunc)
            for args in [((4, 5), (4, 5)),
                         ((4, 5), (4, 6)),
                         ((4, 6), (4, 5)),
                         ((4, 5), (5, 4))]:
                eq(pyfunc, cfunc, args)
        # Different-sized tuples
        argtypes = [types.Tuple((types.int64, types.float32)),
                    types.UniTuple(types.int32, 3)]
        cfunc = njit(tuple(argtypes),)(pyfunc)
        for args in [((4, 5), (4, 5, 6)),
                     ((4, 5), (4, 4, 6)),
                     ((4, 5), (4, 6, 7))]:
            eq(pyfunc, cfunc, args)

    def test_eq(self):
        self._test_compare(eq_usecase)

    def test_ne(self):
        self._test_compare(ne_usecase)

    def test_gt(self):
        self._test_compare(gt_usecase)

    def test_ge(self):
        self._test_compare(ge_usecase)

    def test_lt(self):
        self._test_compare(lt_usecase)

    def test_le(self):
        self._test_compare(le_usecase)


class TestNamedTuple(TestCase, MemoryLeakMixin):

    def test_unpack(self):
        def check(p):
            for pyfunc in tuple_first, tuple_second:
                cfunc = jit(nopython=True)(pyfunc)
                self.assertPreciseEqual(cfunc(p), pyfunc(p))

        # Homogeneous
        check(Rect(4, 5))
        # Heterogeneous
        check(Rect(4, 5.5))

    def test_len(self):
        def check(p):
            pyfunc = len_usecase
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        # Homogeneous
        check(Rect(4, 5))
        check(Point(4, 5, 6))
        # Heterogeneous
        check(Rect(4, 5.5))
        check(Point(4, 5.5, 6j))

    def test_index(self):
        pyfunc = tuple_index
        cfunc = jit(nopython=True)(pyfunc)

        p = Point(4, 5, 6)
        for i in range(len(p)):
            self.assertPreciseEqual(cfunc(p, i), pyfunc(p, i))

        # test uintp indexing (because, e.g., parfor generates unsigned prange)
        for i in range(len(p)):
            self.assertPreciseEqual(cfunc(p, types.uintp(i)), pyfunc(p, i))

    def test_bool(self):
        def check(p):
            pyfunc = bool_usecase
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        # Homogeneous
        check(Rect(4, 5))
        # Heterogeneous
        check(Rect(4, 5.5))
        check(Empty())

    def _test_compare(self, pyfunc):
        def eq(pyfunc, cfunc, args):
            self.assertIs(cfunc(*args), pyfunc(*args),
                          "mismatch for arguments %s" % (args,))

        cfunc = jit(nopython=True)(pyfunc)

        # Same-sized named tuples
        for a, b in [((4, 5), (4, 5)),
                     ((4, 5), (4, 6)),
                     ((4, 6), (4, 5)),
                     ((4, 5), (5, 4))]:
            eq(pyfunc, cfunc, (Rect(*a), Rect(*b)))

        # Different-sized named tuples
        for a, b in [((4, 5), (4, 5, 6)),
                     ((4, 5), (4, 4, 6)),
                     ((4, 5), (4, 6, 7))]:
            eq(pyfunc, cfunc, (Rect(*a), Point(*b)))

    def test_eq(self):
        self._test_compare(eq_usecase)

    def test_ne(self):
        self._test_compare(ne_usecase)

    def test_gt(self):
        self._test_compare(gt_usecase)

    def test_ge(self):
        self._test_compare(ge_usecase)

    def test_lt(self):
        self._test_compare(lt_usecase)

    def test_le(self):
        self._test_compare(le_usecase)

    def test_getattr(self):
        pyfunc = getattr_usecase
        cfunc = jit(nopython=True)(pyfunc)

        for args in (4, 5, 6), (4, 5.5, 6j):
            p = Point(*args)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

    def test_construct(self):
        def check(pyfunc):
            cfunc = jit(nopython=True)(pyfunc)
            for args in (4, 5, 6), (4, 5.5, 6j):
                expected = pyfunc(*args)
                got = cfunc(*args)
                self.assertIs(type(got), type(expected))
                self.assertPreciseEqual(got, expected)

        check(make_point)
        check(make_point_kws)

    def test_type(self):
        # Test the type() built-in on named tuples
        pyfunc = type_usecase
        cfunc = jit(nopython=True)(pyfunc)

        arg_tuples = [(4, 5, 6), (4, 5.5, 6j)]
        for tup_args, args in itertools.product(arg_tuples, arg_tuples):
            tup = Point(*tup_args)
            expected = pyfunc(tup, *args)
            got = cfunc(tup, *args)
            self.assertIs(type(got), type(expected))
            self.assertPreciseEqual(got, expected)

    def test_literal_unification(self):
        # Test for #3565.
        @jit(nopython=True)
        def Data1(value):
            return Rect(value, -321)

        @jit(nopython=True)
        def call(i, j):
            if j == 0:
                # In the error, `result` is typed to `Rect(int, LiteralInt)`
                # because of the `-321` literal.  This doesn't match the
                # `result` type in the other branch.
                result = Data1(i)
            else:
                # `result` is typed to be `Rect(int, int)`
                result = Rect(i, j)
            return result

        r = call(123, 1321)
        self.assertEqual(r, Rect(width=123, height=1321))
        r = call(123, 0)
        self.assertEqual(r, Rect(width=123, height=-321))

    def test_string_literal_in_ctor(self):
        # Test for issue #3813

        @jit(nopython=True)
        def foo():
            return Rect(10, 'somestring')

        r = foo()
        self.assertEqual(r, Rect(width=10, height='somestring'))

    def test_dispatcher_mistreat(self):
        # Test for issue #5215 that mistreat namedtuple as tuples
        @jit(nopython=True)
        def foo(x):
            return x

        in1 = (1, 2,  3)
        out1 = foo(in1)
        self.assertEqual(in1, out1)

        in2 = Point(1, 2, 3)
        out2 = foo(in2)
        self.assertEqual(in2, out2)

        # Check the signatures
        self.assertEqual(len(foo.nopython_signatures), 2)
        self.assertEqual(foo.nopython_signatures[0].args[0], typeof(in1))
        self.assertEqual(foo.nopython_signatures[1].args[0], typeof(in2))

        # Differently named
        in3 = Point2(1, 2, 3)
        out3 = foo(in3)
        self.assertEqual(in3, out3)
        self.assertEqual(len(foo.nopython_signatures), 3)
        self.assertEqual(foo.nopython_signatures[2].args[0], typeof(in3))


class TestTupleNRT(TestCase, MemoryLeakMixin):
    def test_tuple_add(self):
        def pyfunc(x):
            a = np.arange(3)
            return (a,) + (x,)

        cfunc = jit(nopython=True)(pyfunc)
        x = 123
        expect_a, expect_x = pyfunc(x)
        got_a, got_x = cfunc(x)
        np.testing.assert_equal(got_a, expect_a)
        self.assertEqual(got_x, expect_x)


class TestNamedTupleNRT(TestCase, MemoryLeakMixin):

    def test_return(self):
        # Check returning a namedtuple with a list inside it
        pyfunc = make_point_nrt
        cfunc = jit(nopython=True)(pyfunc)

        for arg in (3, 0):
            expected = pyfunc(arg)
            got = cfunc(arg)
            self.assertIs(type(got), type(expected))
            self.assertPreciseEqual(got, expected)


class TestConversions(TestCase):
    """
    Test implicit conversions between tuple types.
    """

    def check_conversion(self, fromty, toty, val):
        pyfunc = identity
        cfunc = njit(toty(fromty))(pyfunc)
        res = cfunc(val)
        self.assertEqual(res, val)

    def test_conversions(self):
        check = self.check_conversion
        fromty = types.UniTuple(types.int32, 2)
        check(fromty, types.UniTuple(types.float32, 2), (4, 5))
        check(fromty, types.Tuple((types.float32, types.int16)), (4, 5))
        aty = types.UniTuple(types.int32, 0)
        bty = types.Tuple(())
        check(aty, bty, ())
        check(bty, aty, ())

        with self.assertRaises(errors.TypingError) as raises:
            check(fromty, types.Tuple((types.float32,)), (4, 5))
        msg = "No conversion from UniTuple(int32 x 2) to UniTuple(float32 x 1)"
        self.assertIn(msg, str(raises.exception))


class TestMethods(TestCase):

    def test_index(self):
        pyfunc = index_method_usecase
        cfunc = jit(nopython=True)(pyfunc)
        self.assertEqual(cfunc((1, 2, 3), 2), 1)

        with self.assertRaises(ValueError) as raises:
            cfunc((1, 2, 3), 4)
        msg = 'tuple.index(x): x not in tuple'
        self.assertEqual(msg, str(raises.exception))


class TestTupleBuild(TestCase):

    def test_build_unpack(self):
        def check(p):
            pyfunc = lambda a: (1, *a)
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        # Homogeneous
        check((4, 5))
        # Heterogeneous
        check((4, 5.5))

    def test_build_unpack_assign_like(self):
        # see #6534
        def check(p):
            pyfunc = lambda a: (*a,)
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        # Homogeneous
        check((4, 5))
        # Heterogeneous
        check((4, 5.5))

    def test_build_unpack_fail_on_list_assign_like(self):
        # see #6534
        def check(p):
            pyfunc = lambda a: (*a,)
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        with self.assertRaises(errors.TypingError) as raises:
            check([4, 5])

        # Python 3.9 has a peephole rewrite due to large changes in tuple
        # unpacking. It results in a tuple + list situation from the above
        # so the error message reflects that. Catching this specific and
        # seemingly rare sequence in the peephole rewrite is prohibitively
        # hard. Should it be reported numerous times, revisit then.
        msg1 = "No implementation of function"
        self.assertIn(msg1, str(raises.exception))
        msg2 = "tuple(reflected list(" # ignore the rest of reflected list
                                        # part, it's repr is quite volatile.
        self.assertIn(msg2, str(raises.exception))

    def test_build_unpack_more(self):
        def check(p):
            pyfunc = lambda a: (1, *a, (1, 2), *a)
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        # Homogeneous
        check((4, 5))
        # Heterogeneous
        check((4, 5.5))

    def test_build_unpack_call(self):
        def check(p):
            @jit
            def inner(*args):
                return args
            pyfunc = lambda a: inner(1, *a)
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        # Homogeneous
        check((4, 5))
        # Heterogeneous
        check((4, 5.5))

    def test_build_unpack_call_more(self):
        def check(p):
            @jit
            def inner(*args):
                return args
            pyfunc = lambda a: inner(1, *a, *(1, 2), *a)
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        # Homogeneous
        check((4, 5))
        # Heterogeneous
        check((4, 5.5))

    def test_tuple_constructor(self):
        def check(pyfunc, arg):
            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(arg), pyfunc(arg))

        # empty
        check(lambda _: tuple(), ())
        # Homogeneous
        check(lambda a: tuple(a), (4, 5))
        # Heterogeneous
        check(lambda a: tuple(a), (4, 5.5))

    @unittest.skipIf(utils.PYVERSION < (3, 9), "needs Python 3.9+")
    def test_unpack_with_predicate_fails(self):
        # this fails as the list_to_tuple/list_extend peephole bytecode
        # rewriting needed for Python 3.9+ cannot yet traverse the CFG.
        @njit
        def foo():
            a = (1,)
            b = (3,2,  4)
            return (*(b if a[0] else (5, 6)),)

        with self.assertRaises(errors.UnsupportedError) as raises:
            foo()
        msg = "op_LIST_EXTEND at the start of a block"
        self.assertIn(msg, str(raises.exception))

    def test_build_unpack_with_calls_in_unpack(self):
        def check(p):
            def pyfunc(a):
                z = [1, 2]
                return (*a, z.append(3), z.extend(a), np.ones(3)), z

            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        check((4, 5))

    def test_build_unpack_complicated(self):
        def check(p):
            def pyfunc(a):
                z = [1, 2]
                return (*a, *(*a, a), *(a, (*(a, (1, 2), *(3,), *a),
                        (a, 1, (2, 3), *a, 1), (1,))),
                        *(z.append(4), z.extend(a))), z

            cfunc = jit(nopython=True)(pyfunc)
            self.assertPreciseEqual(cfunc(p), pyfunc(p))

        check((10, 20))


if __name__ == '__main__':
    unittest.main()
