import copy
import itertools
import operator
import unittest

import numpy as np

from numba import jit, njit
from numba.core import types, utils, errors
from numba.core.types.functions import _header_lead
from numba.tests.support import TestCase, tag, needs_blas
from numba.tests.matmul_usecase import (matmul_usecase, imatmul_usecase,
                                        DumbMatrix,)

Noflags = {'nopython': True}

force_pyobj_flags = {'forceobj': True}


def make_static_power(exp):
    def pow_usecase(x):
        return x ** exp
    return pow_usecase


class LiteralOperatorImpl(object):

    @staticmethod
    def add_usecase(x, y):
        return x + y

    @staticmethod
    def iadd_usecase(x, y):
        x += y
        return x

    @staticmethod
    def sub_usecase(x, y):
        return x - y

    @staticmethod
    def isub_usecase(x, y):
        x -= y
        return x

    @staticmethod
    def mul_usecase(x, y):
        return x * y

    @staticmethod
    def imul_usecase(x, y):
        x *= y
        return x

    @staticmethod
    def floordiv_usecase(x, y):
        return x // y

    @staticmethod
    def ifloordiv_usecase(x, y):
        x //= y
        return x

    @staticmethod
    def truediv_usecase(x, y):
        return x / y

    @staticmethod
    def itruediv_usecase(x, y):
        x /= y
        return x

    if matmul_usecase:
        matmul_usecase = staticmethod(matmul_usecase)
        imatmul_usecase = staticmethod(imatmul_usecase)

    @staticmethod
    def mod_usecase(x, y):
        return x % y

    @staticmethod
    def imod_usecase(x, y):
        x %= y
        return x

    @staticmethod
    def pow_usecase(x, y):
        return x ** y

    @staticmethod
    def ipow_usecase(x, y):
        x **= y
        return x

    @staticmethod
    def bitshift_left_usecase(x, y):
        return x << y

    @staticmethod
    def bitshift_ileft_usecase(x, y):
        x <<= y
        return x

    @staticmethod
    def bitshift_right_usecase(x, y):
        return x >> y

    @staticmethod
    def bitshift_iright_usecase(x, y):
        x >>= y
        return x

    @staticmethod
    def bitwise_and_usecase(x, y):
        return x & y

    @staticmethod
    def bitwise_iand_usecase(x, y):
        x &= y
        return x

    @staticmethod
    def bitwise_or_usecase(x, y):
        return x | y

    @staticmethod
    def bitwise_ior_usecase(x, y):
        x |= y
        return x

    @staticmethod
    def bitwise_xor_usecase(x, y):
        return x ^ y

    @staticmethod
    def bitwise_ixor_usecase(x, y):
        x ^= y
        return x

    @staticmethod
    def bitwise_not_usecase_binary(x, _unused):
        return ~x

    @staticmethod
    def bitwise_not_usecase(x):
        return ~x

    @staticmethod
    def not_usecase(x):
        return not(x)

    @staticmethod
    def negate_usecase(x):
        return -x

    @staticmethod
    def unary_positive_usecase(x):
        return +x

    @staticmethod
    def lt_usecase(x, y):
        return x < y

    @staticmethod
    def le_usecase(x, y):
        return x <= y

    @staticmethod
    def gt_usecase(x, y):
        return x > y

    @staticmethod
    def ge_usecase(x, y):
        return x >= y

    @staticmethod
    def eq_usecase(x, y):
        return x == y

    @staticmethod
    def ne_usecase(x, y):
        return x != y

    @staticmethod
    def in_usecase(x, y):
        return x in y

    @staticmethod
    def not_in_usecase(x, y):
        return x not in y

    @staticmethod
    def is_usecase(x, y):
        return x is y


class FunctionalOperatorImpl(object):

    @staticmethod
    def add_usecase(x, y):
        return operator.add(x, y)

    @staticmethod
    def iadd_usecase(x, y):
        return operator.iadd(x, y)

    @staticmethod
    def sub_usecase(x, y):
        return operator.sub(x, y)

    @staticmethod
    def isub_usecase(x, y):
        return operator.isub(x, y)

    @staticmethod
    def mul_usecase(x, y):
        return operator.mul(x, y)

    @staticmethod
    def imul_usecase(x, y):
        return operator.imul(x, y)

    @staticmethod
    def floordiv_usecase(x, y):
        return operator.floordiv(x, y)

    @staticmethod
    def ifloordiv_usecase(x, y):
        return operator.ifloordiv(x, y)

    @staticmethod
    def truediv_usecase(x, y):
        return operator.truediv(x, y)

    @staticmethod
    def itruediv_usecase(x, y):
        return operator.itruediv(x, y)

    @staticmethod
    def mod_usecase(x, y):
        return operator.mod(x, y)

    @staticmethod
    def imod_usecase(x, y):
        return operator.imod(x, y)

    @staticmethod
    def pow_usecase(x, y):
        return operator.pow(x, y)

    @staticmethod
    def ipow_usecase(x, y):
        return operator.ipow(x, y)

    @staticmethod
    def matmul_usecase(x, y):
        return operator.matmul(x, y)

    @staticmethod
    def imatmul_usecase(x, y):
        return operator.imatmul(x, y)

    @staticmethod
    def bitshift_left_usecase(x, y):
        return operator.lshift(x, y)

    @staticmethod
    def bitshift_ileft_usecase(x, y):
        return operator.ilshift(x, y)

    @staticmethod
    def bitshift_right_usecase(x, y):
        return operator.rshift(x, y)

    @staticmethod
    def bitshift_iright_usecase(x, y):
        return operator.irshift(x, y)

    @staticmethod
    def bitwise_and_usecase(x, y):
        return operator.and_(x, y)

    @staticmethod
    def bitwise_iand_usecase(x, y):
        return operator.iand(x, y)

    @staticmethod
    def bitwise_or_usecase(x, y):
        return operator.or_(x, y)

    @staticmethod
    def bitwise_ior_usecase(x, y):
        return operator.ior(x, y)

    @staticmethod
    def bitwise_xor_usecase(x, y):
        return operator.xor(x, y)

    @staticmethod
    def bitwise_ixor_usecase(x, y):
        return operator.ixor(x, y)

    @staticmethod
    def bitwise_not_usecase_binary(x, _unused):
        return operator.invert(x)

    @staticmethod
    def bitwise_not_usecase(x):
        return operator.invert(x)

    @staticmethod
    def not_usecase(x):
        return operator.not_(x)

    @staticmethod
    def negate_usecase(x):
        return operator.neg(x)

    @staticmethod
    def unary_positive_usecase(x):
        return operator.pos(x)

    @staticmethod
    def lt_usecase(x, y):
        return operator.lt(x, y)

    @staticmethod
    def le_usecase(x, y):
        return operator.le(x, y)

    @staticmethod
    def gt_usecase(x, y):
        return operator.gt(x, y)

    @staticmethod
    def ge_usecase(x, y):
        return operator.ge(x, y)

    @staticmethod
    def eq_usecase(x, y):
        return operator.eq(x, y)

    @staticmethod
    def ne_usecase(x, y):
        return operator.ne(x, y)

    @staticmethod
    def in_usecase(x, y):
        return operator.contains(y, x)

    @staticmethod
    def not_in_usecase(x, y):
        return not operator.contains(y, x)

    @staticmethod
    def is_usecase(x, y):
        return operator.is_(x, y)


class TestOperators(TestCase):
    """
    Test standard Python operators on scalars.

    NOTE: operators on array are generally tested in test_ufuncs.
    """

    op = LiteralOperatorImpl

    _bitwise_opnames = {
        'bitshift_left_usecase': operator.lshift,
        'bitshift_ileft_usecase': operator.ilshift,
        'bitshift_right_usecase': operator.rshift,
        'bitshift_iright_usecase': operator.irshift,
        'bitwise_and_usecase': operator.and_,
        'bitwise_iand_usecase': operator.iand,
        'bitwise_or_usecase': operator.or_,
        'bitwise_ior_usecase': operator.ior,
        'bitwise_xor_usecase': operator.xor,
        'bitwise_ixor_usecase': operator.ixor,
        'bitwise_not_usecase_binary': operator.invert,
    }

    def run_test_ints(self, pyfunc, x_operands, y_operands, types_list,
                      flags=force_pyobj_flags):
        for arg_types in types_list:
            cfunc = jit(arg_types, **flags)(pyfunc)
            for x, y in itertools.product(x_operands, y_operands):
                # For inplace ops, we check that the first operand
                # was correctly mutated.
                x_got = copy.copy(x)
                x_expected = copy.copy(x)
                got = cfunc(x_got, y)
                expected = pyfunc(x_expected, y)
                self.assertPreciseEqual(
                    got, expected,
                    msg="mismatch for (%r, %r) with types %s: %r != %r"
                        % (x, y, arg_types, got, expected))
                self.assertPreciseEqual(
                    x_got, x_expected,
                    msg="mismatch for (%r, %r) with types %s: %r != %r"
                        % (x, y, arg_types, x_got, x_expected))

    def run_test_floats(self, pyfunc, x_operands, y_operands, types_list,
                        flags=force_pyobj_flags):
        for arg_types in types_list:
            cfunc = jit(arg_types, **flags)(pyfunc)
            for x, y in itertools.product(x_operands, y_operands):
                # For inplace ops, we check that the first operand
                # was correctly mutated.
                x_got = copy.copy(x)
                x_expected = copy.copy(x)
                got = cfunc(x_got, y)
                expected = pyfunc(x_expected, y)
                np.testing.assert_allclose(got, expected, rtol=1e-5)
                np.testing.assert_allclose(x_got, x_expected, rtol=1e-5)

    def coerce_operand(self, op, numba_type):
        if hasattr(op, "dtype"):
            return numba_type.cast_python_value(op)
        elif numba_type in types.unsigned_domain:
            return abs(int(op.real))
        elif numba_type in types.integer_domain:
            return int(op.real)
        elif numba_type in types.real_domain:
            return float(op.real)
        else:
            return op

    def run_test_scalar_compare(self, pyfunc, flags=force_pyobj_flags,
                                ordered=True):
        ops = self.compare_scalar_operands
        types_list = self.compare_types
        if not ordered:
            types_list = types_list + self.compare_unordered_types
        for typ in types_list:
            cfunc = jit((typ, typ), **flags)(pyfunc)
            for x, y in itertools.product(ops, ops):
                x = self.coerce_operand(x, typ)
                y = self.coerce_operand(y, typ)
                expected = pyfunc(x, y)
                got = cfunc(x, y)
                # Scalar ops => scalar result
                self.assertIs(type(got), type(expected))
                self.assertEqual(got, expected,
                                 "mismatch with %r (%r, %r)"
                                 % (typ, x, y))


    #
    # Comparison operators
    #

    compare_scalar_operands = [-0.5, -1.0 + 1j, -1.0 + 2j, -0.5 + 1j, 1.5]
    compare_types = [types.int32, types.int64,
                     types.uint32, types.uint64,
                     types.float32, types.float64]
    compare_unordered_types = [types.complex64, types.complex128]

    def test_lt_scalar(self, flags=force_pyobj_flags):
        self.run_test_scalar_compare(self.op.lt_usecase, flags)

    def test_lt_scalar_npm(self):
        self.test_lt_scalar(flags=Noflags)

    def test_le_scalar(self, flags=force_pyobj_flags):
        self.run_test_scalar_compare(self.op.le_usecase, flags)

    def test_le_scalar_npm(self):
        self.test_le_scalar(flags=Noflags)

    def test_gt_scalar(self, flags=force_pyobj_flags):
        self.run_test_scalar_compare(self.op.gt_usecase, flags)

    def test_gt_scalar_npm(self):
        self.test_gt_scalar(flags=Noflags)

    def test_ge_scalar(self, flags=force_pyobj_flags):
        self.run_test_scalar_compare(self.op.ge_usecase, flags)

    def test_ge_scalar_npm(self):
        self.test_ge_scalar(flags=Noflags)

    def test_eq_scalar(self, flags=force_pyobj_flags):
        self.run_test_scalar_compare(self.op.eq_usecase, flags, ordered=False)

    def test_eq_scalar_npm(self):
        self.test_eq_scalar(flags=Noflags)

    def test_ne_scalar(self, flags=force_pyobj_flags):
        self.run_test_scalar_compare(self.op.ne_usecase, flags, ordered=False)

    def test_ne_scalar_npm(self):
        self.test_ne_scalar(flags=Noflags)

    def test_is_ellipsis(self):
        cfunc = njit((types.ellipsis, types.ellipsis))(self.op.is_usecase)
        self.assertTrue(cfunc(Ellipsis, Ellipsis))

    def test_is_void_ptr(self):
        # can't call this directly from python, as void cannot be unboxed
        cfunc_void = jit(
            (types.voidptr, types.voidptr), nopython=True
        )(self.op.is_usecase)

        # this wrapper performs the casts from int to voidptr for us
        @jit(nopython=True)
        def cfunc(x, y):
            return cfunc_void(x, y)

        self.assertTrue(cfunc(1, 1))
        self.assertFalse(cfunc(1, 2))

    #
    # Arithmetic operators
    #

    def run_binop_bools(self, pyfunc, flags=force_pyobj_flags):
        x_operands = [False, False, True, True]
        y_operands = [False, True, False, True]

        types_list = [(types.boolean, types.boolean)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

    def run_binop_ints(self, pyfunc, flags=force_pyobj_flags):
        x_operands = [-5, 0, 1, 2]
        y_operands = [-3, -1, 1, 3]

        types_list = [(types.int32, types.int32),
                      (types.int64, types.int64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = [2, 3]
        y_operands = [1, 2]

        types_list = [(types.byte, types.byte),
                      (types.uint32, types.uint32),
                      (types.uint64, types.uint64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

    def run_binop_floats(self, pyfunc, flags=force_pyobj_flags):
        x_operands = [-1.1, 0.0, 1.1]
        y_operands = [-1.5, 0.8, 2.1]

        types_list = [(types.float32, types.float32),
                      (types.float64, types.float64)]

        self.run_test_floats(pyfunc, x_operands, y_operands, types_list,
                             flags=flags)

    def run_binop_floats_floordiv(self, pyfunc, flags=force_pyobj_flags):
        self.run_binop_floats(pyfunc, flags=flags)

    def run_binop_complex(self, pyfunc, flags=force_pyobj_flags):
        x_operands = [-1.1 + 0.3j, 0.0 + 0.0j, 1.1j]
        y_operands = [-1.5 - 0.7j, 0.8j, 2.1 - 2.0j]

        types_list = [(types.complex64, types.complex64),
                      (types.complex128, types.complex128)]

        self.run_test_floats(pyfunc, x_operands, y_operands, types_list,
                             flags=flags)

    def generate_binop_tests(ns, usecases, tp_runners, npm_array=False):
        for usecase in usecases:
            for tp_name, runner_name in tp_runners.items():
                for nopython in (False, True):
                    test_name = "test_%s_%s" % (usecase, tp_name)
                    if nopython:
                        test_name += "_npm"
                    flags = Noflags if nopython else force_pyobj_flags
                    usecase_name = "%s_usecase" % usecase

                    def inner(self, runner_name=runner_name,
                              usecase_name=usecase_name, flags=flags):
                        runner = getattr(self, runner_name)
                        op_usecase = getattr(self.op, usecase_name)
                        runner(op_usecase, flags)

                    if nopython and 'array' in tp_name and not npm_array:
                        def test_meth(self):
                            with self.assertTypingError():
                                inner()
                    else:
                        test_meth = inner

                    test_meth.__name__ = test_name

                    if nopython:
                        test_meth = tag('important')(test_meth)

                    ns[test_name] = test_meth


    generate_binop_tests(locals(),
                         ('add', 'iadd', 'sub', 'isub', 'mul', 'imul'),
                         {'ints': 'run_binop_ints',
                          'floats': 'run_binop_floats',
                          'complex': 'run_binop_complex',
                          })

    generate_binop_tests(locals(),
                         ('truediv', 'itruediv'),
                         {'ints': 'run_binop_ints',
                          'floats': 'run_binop_floats',
                          'complex': 'run_binop_complex',
                          })

    # NOTE: floordiv and mod unsupported for complex numbers
    generate_binop_tests(locals(),
                         ('floordiv', 'ifloordiv', 'mod', 'imod'),
                         {'ints': 'run_binop_ints',
                          'floats': 'run_binop_floats_floordiv',
                          })

    def check_div_errors(self, usecase_name, msg, flags=force_pyobj_flags,
                         allow_complex=False):
        pyfunc = getattr(self.op, usecase_name)
        # Signed and unsigned division can take different code paths,
        # test them both.
        arg_types = [types.int32, types.uint32, types.float64]
        if allow_complex:
            arg_types.append(types.complex128)
        for tp in arg_types:
            cfunc = jit((tp, tp), **flags)(pyfunc)
            with self.assertRaises(ZeroDivisionError) as cm:
                cfunc(1, 0)
            # Test exception message if not in object mode
            if flags is not force_pyobj_flags:
                self.assertIn(msg, str(cm.exception))

    def test_truediv_errors(self, flags=force_pyobj_flags):
        self.check_div_errors("truediv_usecase", "division by zero", flags=flags,
                              allow_complex=True)

    def test_truediv_errors_npm(self):
        self.test_truediv_errors(flags=Noflags)

    def test_floordiv_errors(self, flags=force_pyobj_flags):
        self.check_div_errors("floordiv_usecase", "division by zero", flags=flags)

    def test_floordiv_errors_npm(self):
        self.test_floordiv_errors(flags=Noflags)

    def test_mod_errors(self, flags=force_pyobj_flags):
        self.check_div_errors("mod_usecase", "modulo by zero", flags=flags)

    def test_mod_errors_npm(self):
        self.test_mod_errors(flags=Noflags)

    def run_pow_ints(self, pyfunc, flags=force_pyobj_flags):
        x_operands = [-2, -1, 0, 1, 2]
        y_operands = [0, 1, 2]

        types_list = [(types.int32, types.int32),
                      (types.int64, types.int64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = [0, 1, 2]
        y_operands = [0, 1, 2]

        types_list = [(types.byte, types.byte),
                      (types.uint32, types.uint32),
                      (types.uint64, types.uint64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

    def run_pow_floats(self, pyfunc, flags=force_pyobj_flags):
        x_operands = [-222.222, -111.111, 111.111, 222.222]
        y_operands = [-2, -1, 0, 1, 2]

        types_list = [(types.float32, types.float32),
                      (types.float64, types.float64)]

        self.run_test_floats(pyfunc, x_operands, y_operands, types_list,
                             flags=flags)

        x_operands = [0.0]
        y_operands = [0, 1, 2]  # TODO native handling of 0 ** negative power

        types_list = [(types.float32, types.float32),
                      (types.float64, types.float64)]

        self.run_test_floats(pyfunc, x_operands, y_operands, types_list,
                             flags=flags)

    # XXX power operator is unsupported on complex numbers (see issue #488)
    generate_binop_tests(locals(),
                         ('pow', 'ipow'),
                         {'ints': 'run_pow_ints',
                          'floats': 'run_pow_floats',
                          })

    def test_add_complex(self, flags=force_pyobj_flags):
        pyfunc = self.op.add_usecase

        x_operands = [1+0j, 1j, -1-1j]
        y_operands = x_operands

        types_list = [(types.complex64, types.complex64),
                      (types.complex128, types.complex128),]

        self.run_test_floats(pyfunc, x_operands, y_operands, types_list,
                             flags=flags)

    def test_add_complex_npm(self):
        self.test_add_complex(flags=Noflags)

    def test_sub_complex(self, flags=force_pyobj_flags):
        pyfunc = self.op.sub_usecase

        x_operands = [1+0j, 1j, -1-1j]
        y_operands = [1, 2, 3]

        types_list = [(types.complex64, types.complex64),
                      (types.complex128, types.complex128),]

        self.run_test_floats(pyfunc, x_operands, y_operands, types_list,
                             flags=flags)

    def test_sub_complex_npm(self):
        self.test_sub_complex(flags=Noflags)

    def test_mul_complex(self, flags=force_pyobj_flags):
        pyfunc = self.op.mul_usecase

        x_operands = [1+0j, 1j, -1-1j]
        y_operands = [1, 2, 3]

        types_list = [(types.complex64, types.complex64),
                      (types.complex128, types.complex128),]

        self.run_test_floats(pyfunc, x_operands, y_operands, types_list,
                             flags=flags)

    def test_mul_complex_npm(self):
        self.test_mul_complex(flags=Noflags)

    def test_truediv_complex(self, flags=force_pyobj_flags):
        pyfunc = self.op.truediv_usecase

        x_operands = [1+0j, 1j, -1-1j]
        y_operands = [1, 2, 3]

        types_list = [(types.complex64, types.complex64),
                      (types.complex128, types.complex128),]

        self.run_test_floats(pyfunc, x_operands, y_operands, types_list,
                             flags=flags)

    def test_truediv_complex_npm(self):
        self.test_truediv_complex(flags=Noflags)

    def test_mod_complex(self, flags=force_pyobj_flags):
        pyfunc = self.op.mod_usecase
        cres = jit((types.complex64, types.complex64), **flags)(pyfunc)
        with self.assertRaises(TypeError) as raises:
            cres(4j, 2j)

        # error message depends on Python version.
        if utils.PYVERSION in ((3, 9),):
            msg = "can't mod complex numbers"
        elif utils.PYVERSION in ((3, 10), (3, 11), (3, 12)):
            msg = "unsupported operand type(s) for %"
        else:
            raise NotImplementedError(utils.PYVERSION)

        self.assertIn(msg, str(raises.exception))

    def test_mod_complex_npm(self):
        pyfunc = self.op.mod_usecase
        with self.assertTypingError():
            njit((types.complex64, types.complex64))(pyfunc)

    #
    # Matrix multiplication
    # (just check with simple values; computational tests are in test_linalg)
    #

    def check_matmul_objmode(self, pyfunc, inplace):
        # Use dummy objects, to work with any NumPy / SciPy version
        cfunc = jit((), **force_pyobj_flags)(pyfunc)
        a = DumbMatrix(3)
        b = DumbMatrix(4)
        got = cfunc(a, b)
        self.assertEqual(got.value, 12)
        if inplace:
            self.assertIs(got, a)
        else:
            self.assertIsNot(got, a)
            self.assertIsNot(got, b)

    def test_matmul(self):
        self.check_matmul_objmode(self.op.matmul_usecase, inplace=False)

    def test_imatmul(self):
        self.check_matmul_objmode(self.op.imatmul_usecase, inplace=True)

    @needs_blas
    def check_matmul_npm(self, pyfunc):
        arrty = types.Array(types.float32, 1, 'C')
        cfunc = njit((arrty, arrty))(pyfunc)
        a = np.float32([1, 2])
        b = np.float32([3, 4])
        got = cfunc(a, b)
        self.assertPreciseEqual(got, np.dot(a, b))
        # Never inplace
        self.assertIsNot(got, a)
        self.assertIsNot(got, b)

    def test_matmul_npm(self):
        self.check_matmul_npm(self.op.matmul_usecase)

    def test_imatmul_npm(self):
        with self.assertTypingError() as raises:
            self.check_matmul_npm(self.op.imatmul_usecase)

    #
    # Bitwise operators
    #

    def run_bitshift_left(self, pyfunc, flags=force_pyobj_flags):
        x_operands = [0, 1]
        y_operands = [0, 1, 2, 4, 8, 16, 31]

        types_list = [(types.uint32, types.uint32)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = [0, 1]
        y_operands = [0, 1, 2, 4, 8, 16, 32, 63]

        types_list = [(types.uint64, types.uint64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = [0, -1]
        y_operands = [0, 1, 2, 4, 8, 16, 31]

        types_list = [(types.int32, types.int32)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = [0, -1]
        y_operands = [0, 1, 2, 4, 8, 16, 32, 63]

        types_list = [(types.int64, types.int64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

    generate_binop_tests(locals(),
                         ('bitshift_left', 'bitshift_ileft'),
                         {'ints': 'run_bitshift_left',
                          })

    def run_bitshift_right(self, pyfunc, flags=force_pyobj_flags):
        x_operands = [0, 1, 2**32 - 1]
        y_operands = [0, 1, 2, 4, 8, 16, 31]

        types_list = [(types.uint32, types.uint32)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = [0, 1, 2**64 - 1]
        y_operands = [0, 1, 2, 4, 8, 16, 32, 63]

        types_list = [(types.uint64, types.uint64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = [0, 1, -(2**31)]
        y_operands = [0, 1, 2, 4, 8, 16, 31]

        types_list = [(types.int32, types.int32)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = [0, -1, -(2**31)]
        y_operands = [0, 1, 2, 4, 8, 16, 32, 63]

        types_list = [(types.int64, types.int64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

    generate_binop_tests(locals(),
                         ('bitshift_right', 'bitshift_iright'),
                         {'ints': 'run_bitshift_right',
                          })

    def run_logical(self, pyfunc, flags=force_pyobj_flags):
        x_operands = list(range(0, 8)) + [2**32 - 1]
        y_operands = list(range(0, 8)) + [2**32 - 1]

        types_list = [(types.uint32, types.uint32)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = list(range(0, 8)) + [2**64 - 1]
        y_operands = list(range(0, 8)) + [2**64 - 1]

        types_list = [(types.uint64, types.uint64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = list(range(-4, 4)) + [-(2**31), 2**31 - 1]
        y_operands = list(range(-4, 4)) + [-(2**31), 2**31 - 1]

        types_list = [(types.int32, types.int32)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = list(range(-4, 4)) + [-(2**63), 2**63 - 1]
        y_operands = list(range(-4, 4)) + [-(2**63), 2**63 - 1]

        types_list = [(types.int64, types.int64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

    generate_binop_tests(locals(),
                         ('bitwise_and', 'bitwise_iand',
                          'bitwise_or', 'bitwise_ior',
                          'bitwise_xor', 'bitwise_ixor'),
                         {'ints': 'run_logical',
                          'bools': 'run_binop_bools',
                          })

    #
    # Unary operators
    #

    def test_bitwise_not(self, flags=force_pyobj_flags):
        pyfunc = self.op.bitwise_not_usecase_binary

        x_operands = list(range(0, 8)) + [2**32 - 1]
        x_operands = [np.uint32(x) for x in x_operands]
        y_operands = [0]

        types_list = [(types.uint32, types.uint32)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = list(range(-4, 4)) + [-(2**31), 2**31 - 1]
        y_operands = [0]

        types_list = [(types.int32, types.int32)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = list(range(0, 8)) + [2**64 - 1]
        x_operands = [np.uint64(x) for x in x_operands]
        y_operands = [0]

        types_list = [(types.uint64, types.uint64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        x_operands = list(range(-4, 4)) + [-(2**63), 2**63 - 1]
        y_operands = [0]

        types_list = [(types.int64, types.int64)]

        self.run_test_ints(pyfunc, x_operands, y_operands, types_list,
                           flags=flags)

        # For booleans, we follow Numpy semantics (i.e. ~True == False,
        # not ~True == -2)
        values = [False, False, True, True]
        values = list(map(np.bool_, values))

        pyfunc = self.op.bitwise_not_usecase
        cfunc = jit((types.boolean,), **flags)(pyfunc)
        for val in values:
            self.assertPreciseEqual(pyfunc(val), cfunc(val))

    def test_bitwise_not_npm(self):
        self.test_bitwise_not(flags=Noflags)

    def test_bitwise_float(self):
        """
        Make sure that bitwise float operations are not allowed
        """
        def assert_reject_compile(pyfunc, argtypes, opname):
            msg = 'expecting TypingError when compiling {}'.format(pyfunc)
            with self.assertRaises(errors.TypingError, msg=msg) as raises:
                njit(argtypes)(pyfunc)
            # check error message
            fmt = _header_lead + ' {}'
            expecting = fmt.format(opname
                                   if isinstance(opname, str)
                                   else 'Function({})'.format(opname))
            self.assertIn(expecting, str(raises.exception))

        methods = [
            'bitshift_left_usecase',
            'bitshift_ileft_usecase',
            'bitshift_right_usecase',
            'bitshift_iright_usecase',
            'bitwise_and_usecase',
            'bitwise_iand_usecase',
            'bitwise_or_usecase',
            'bitwise_ior_usecase',
            'bitwise_xor_usecase',
            'bitwise_ixor_usecase',
            'bitwise_not_usecase_binary',
        ]

        for name in methods:
            pyfunc = getattr(self.op, name)
            assert_reject_compile(pyfunc, (types.float32, types.float32),
                                  opname=self._bitwise_opnames[name])

    def test_not(self):
        pyfunc = self.op.not_usecase

        values = [
            1,
            2,
            3,
            1.2,
            3.4j,
        ]

        cfunc = jit((), **force_pyobj_flags)(pyfunc)
        for val in values:
            self.assertEqual(pyfunc(val), cfunc(val))

    def test_not_npm(self):
        pyfunc = self.op.not_usecase
        # test native mode
        argtys = [
            types.int8,
            types.int32,
            types.int64,
            types.float32,
            types.complex128,
        ]
        values = [
            1,
            2,
            3,
            1.2,
            3.4j,
        ]
        for ty, val in zip(argtys, values):
            cfunc = njit((ty,))(pyfunc)
            self.assertEqual(cfunc.nopython_signatures[0].return_type,
                             types.boolean)
            self.assertEqual(pyfunc(val), cfunc(val))

    # XXX test_negate should check for negative and positive zeros and infinities

    def test_negate_npm(self):
        pyfunc = self.op.negate_usecase
        # test native mode
        argtys = [
            types.int8,
            types.int32,
            types.int64,
            types.float32,
            types.float64,
            types.complex128,
            types.boolean,
            types.boolean,
        ]
        values = [
            1,
            2,
            3,
            1.2,
            2.4,
            3.4j,
            True,
            False,
        ]
        for ty, val in zip(argtys, values):
            cfunc = njit((ty,))(pyfunc)
            self.assertAlmostEqual(pyfunc(val), cfunc(val))


    def test_negate(self):
        pyfunc = self.op.negate_usecase
        values = [
            1,
            2,
            3,
            1.2,
            3.4j,
            True,
            False,
        ]
        cfunc = jit((), **force_pyobj_flags)(pyfunc)
        for val in values:
            self.assertEqual(pyfunc(val), cfunc(val))

    def test_unary_positive_npm(self):
        pyfunc = self.op.unary_positive_usecase
        # test native mode
        argtys = [
            types.int8,
            types.int32,
            types.int64,
            types.float32,
            types.float64,
            types.complex128,
            types.boolean,
            types.boolean,
        ]
        values = [
            1,
            2,
            3,
            1.2,
            2.4,
            3.4j,
            True,
            False
        ]
        for ty, val in zip(argtys, values):
            cfunc = njit((ty,))(pyfunc)
            self.assertAlmostEqual(pyfunc(val), cfunc(val))

    def test_unary_positive(self):
        pyfunc = self.op.unary_positive_usecase
        values = [
            1,
            2,
            3,
            1.2,
            3.4j,
            True,
            False,
        ]
        cfunc = jit((), **force_pyobj_flags)(pyfunc)
        for val in values:
            self.assertEqual(pyfunc(val), cfunc(val))

    def _check_in(self, pyfunc, flags):
        dtype = types.int64
        cfunc = jit((dtype, types.UniTuple(dtype, 3)), **flags)(pyfunc)
        for i in (3, 4, 5, 6, 42):
            tup = (3, 42, 5)
            self.assertPreciseEqual(pyfunc(i, tup), cfunc(i, tup))

    def test_in(self, flags=force_pyobj_flags):
        self._check_in(self.op.in_usecase, flags)

    def test_in_npm(self):
        self.test_in(flags=Noflags)

    def test_not_in(self, flags=force_pyobj_flags):
        self._check_in(self.op.not_in_usecase, flags)

    def test_not_in_npm(self):
        self.test_not_in(flags=Noflags)


class TestOperatorModule(TestOperators):

    op = FunctionalOperatorImpl

    _bitwise_opnames = {
        'bitshift_left_usecase': operator.lshift,
        'bitshift_ileft_usecase': operator.ilshift,
        'bitshift_right_usecase': operator.rshift,
        'bitshift_iright_usecase': operator.irshift,
        'bitwise_and_usecase': operator.and_,
        'bitwise_iand_usecase': operator.iand,
        'bitwise_or_usecase': operator.or_,
        'bitwise_ior_usecase': operator.ior,
        'bitwise_xor_usecase': operator.xor,
        'bitwise_ixor_usecase': operator.ixor,
        'bitwise_not_usecase_binary': operator.invert,
    }


class TestMixedInts(TestCase):
    """
    Tests for operator calls with mixed integer types.
    """

    op = LiteralOperatorImpl

    int_samples = [0, 1, 3, 10, 42, 127, 10000, -1, -3, -10, -42, -127, -10000]

    int_types = [types.int8, types.uint8, types.int64, types.uint64]
    signed_types = [tp for tp in int_types if tp.signed]
    unsigned_types = [tp for tp in int_types if not tp.signed]
    type_pairs = list(itertools.product(int_types, int_types))
    signed_pairs = [(u, v) for u, v in type_pairs
                    if u.signed or v.signed]
    unsigned_pairs = [(u, v) for u, v in type_pairs
                      if not (u.signed or v.signed)]

    def int_in_dtype_range(self, val, tp):
        tp_info = np.iinfo(tp.key)
        return tp_info.min <= val <= tp_info.max

    def get_numpy_signed_upcast(self, *vals):
        bitwidth = max(v.dtype.itemsize * 8 for v in vals)
        bitwidth = max(bitwidth, types.intp.bitwidth)
        return getattr(np, "int%d" % bitwidth)

    def get_numpy_unsigned_upcast(self, *vals):
        bitwidth = max(v.dtype.itemsize * 8 for v in vals)
        bitwidth = max(bitwidth, types.intp.bitwidth)
        return getattr(np, "uint%d" % bitwidth)

    def get_typed_int(self, typ, val):
        return getattr(np, typ.name)(val)

    def get_control_signed(self, opname):
        op = getattr(operator, opname)
        def control_signed(a, b):
            tp = self.get_numpy_signed_upcast(a, b)
            return op(tp(a), tp(b))
        return control_signed

    def get_control_unsigned(self, opname):
        op = getattr(operator, opname)
        def control_unsigned(a, b):
            tp = self.get_numpy_unsigned_upcast(a, b)
            return op(tp(a), tp(b))
        return control_unsigned

    def run_binary(self, pyfunc, control_func, operands, types,
                   expected_type=int, force_type=lambda x: x,
                   **assertPreciseEqualArgs):
        for xt, yt in types:
            cfunc = njit((xt, yt))(pyfunc)
            for x, y in itertools.product(operands, operands):
                # Check if xt and yt are values with range of dtype x and y
                if not self.int_in_dtype_range(x, xt) or not self.int_in_dtype_range(y, yt):
                    continue
                # Get Numpy typed scalars for the given types and values
                x = self.get_typed_int(xt, x)
                y = self.get_typed_int(yt, y)
                expected = control_func(x, y)
                got = cfunc(x, y)
                self.assertIsInstance(got, expected_type)
                msg = ("mismatch for (%r, %r) with types %s"
                       % (x, y, (xt, yt)))
                got, expected = force_type(got), force_type(expected)
                self.assertPreciseEqual(got, expected, msg=msg,
                                        **assertPreciseEqualArgs)

    def run_unary(self, pyfunc, control_func, operands, types,
                  expected_type=int):
        for xt in types:
            cfunc = njit((xt,))(pyfunc)
            for x in operands:
                if not self.int_in_dtype_range(x, xt):
                    continue
                x = self.get_typed_int(xt, x)
                expected = control_func(x)
                got = cfunc(x)
                self.assertIsInstance(got, expected_type)
                self.assertPreciseEqual(
                    got, expected,
                    msg="mismatch for %r with type %s: %r != %r"
                        % (x, xt, got, expected))

    def run_arith_binop(self, pyfunc, opname, samples,
                        expected_type=int, force_type=lambda x: x,
                        **assertPreciseEqualArgs):
        self.run_binary(pyfunc, self.get_control_signed(opname),
                        samples, self.signed_pairs, expected_type,
                        force_type=force_type,
                        **assertPreciseEqualArgs)
        self.run_binary(pyfunc, self.get_control_unsigned(opname),
                        samples, self.unsigned_pairs, expected_type,
                        force_type=force_type,
                        **assertPreciseEqualArgs)

    def test_add(self):
        self.run_arith_binop(self.op.add_usecase, 'add', self.int_samples)

    def test_sub(self):
        self.run_arith_binop(self.op.sub_usecase, 'sub', self.int_samples)

    def test_mul(self):
        self.run_arith_binop(self.op.mul_usecase, 'mul', self.int_samples)

    def test_floordiv(self):
        samples = [x for x in self.int_samples if x != 0]
        self.run_arith_binop(self.op.floordiv_usecase, 'floordiv', samples)

    def test_mod(self):
        samples = [x for x in self.int_samples if x != 0]
        self.run_arith_binop(self.op.mod_usecase, 'mod', samples)

    def test_pow(self):
        extra_cast = {}
        if utils.PYVERSION == (3, 11):
            extra_cast["force_type"] = float
        pyfunc = self.op.pow_usecase
        # Only test with positive values, as otherwise trying to write the
        # control function in terms of Python or Numpy power turns out insane.
        samples = [x for x in self.int_samples if x >= 0]
        self.run_arith_binop(pyfunc, 'pow', samples, **extra_cast)

        # Now test all non-zero values, but only with signed types
        def control_signed(a, b):
            tp = self.get_numpy_signed_upcast(a, b)
            if b >= 0:
                return tp(a) ** tp(b)
            else:
                inv = tp(a) ** tp(-b)
                if inv == 0:
                    # Overflow
                    return 0
                return np.intp(1.0 / inv)
        samples = [x for x in self.int_samples if x != 0]
        signed_pairs = [(u, v) for u, v in self.type_pairs
                        if u.signed and v.signed]
        self.run_binary(pyfunc, control_signed,
                        samples, signed_pairs, **extra_cast)

    def test_truediv(self):

        def control(a, b):
            return float(a) / float(b)
        samples = [x for x in self.int_samples if x != 0]
        pyfunc = self.op.truediv_usecase

        # Note: there can be precision issues on x87
        # e.g. for `1 / 18446744073709541616`
        # -> 0x1.0000000000002p-64 vs. 0x1.0000000000003p-64.
        self.run_binary(pyfunc, control, samples, self.signed_pairs,
                        expected_type=float, prec='double')
        self.run_binary(pyfunc, control, samples, self.unsigned_pairs,
                        expected_type=float, prec='double')

    def test_and(self):
        self.run_arith_binop(self.op.bitwise_and_usecase, 'and_', self.int_samples)

    def test_or(self):
        self.run_arith_binop(self.op.bitwise_or_usecase, 'or_', self.int_samples)

    def test_xor(self):
        self.run_arith_binop(self.op.bitwise_xor_usecase, 'xor', self.int_samples)

    def run_shift_binop(self, pyfunc, opname):
        opfunc = getattr(operator, opname)
        def control_signed(a, b):
            tp = self.get_numpy_signed_upcast(a, b)
            return opfunc(tp(a), tp(b))
        def control_unsigned(a, b):
            tp = self.get_numpy_unsigned_upcast(a, b)
            return opfunc(tp(a), tp(b))

        samples = self.int_samples

        def check(xt, yt, control_func):
            cfunc = njit((xt, yt))(pyfunc)
            for x in samples:
                # Avoid shifting by more than the shiftand's bitwidth, as
                # we would hit undefined behaviour.
                maxshift = xt.bitwidth - 1
                for y in (0, 1, 3, 5, maxshift - 1, maxshift):
                    if not self.int_in_dtype_range(x, xt) or not self.int_in_dtype_range(y, yt):
                        continue
                    # Get Numpy typed scalars for the given types and values
                    x = self.get_typed_int(xt, x)
                    y = self.get_typed_int(yt, y)
                    expected = control_func(x, y)
                    got = cfunc(x, y)
                    msg = ("mismatch for (%r, %r) with types %s"
                           % (x, y, (xt, yt)))
                    self.assertPreciseEqual(got, expected, msg=msg)

        # For bitshifts, only the first operand's signedness matters
        # to choose the operation's signedness.
        signed_pairs = [(u, v) for u, v in self.type_pairs
                        if u.signed]
        unsigned_pairs = [(u, v) for u, v in self.type_pairs
                          if not u.signed]

        for xt, yt in signed_pairs:
            check(xt, yt, control_signed)
        for xt, yt in unsigned_pairs:
            check(xt, yt, control_unsigned)

    def test_lshift(self):
        self.run_shift_binop(self.op.bitshift_left_usecase, 'lshift')

    def test_rshift(self):
        self.run_shift_binop(self.op.bitshift_right_usecase, 'rshift')

    def test_unary_positive(self):
        def control(a):
            return a
        samples = self.int_samples
        pyfunc = self.op.unary_positive_usecase

        self.run_unary(pyfunc, control, samples, self.int_types)

    def test_unary_negative(self):
        def control_signed(a):
            tp = self.get_numpy_signed_upcast(a)
            return tp(-a)
        def control_unsigned(a):
            tp = self.get_numpy_unsigned_upcast(a)
            return tp(-a)
        samples = self.int_samples
        pyfunc = self.op.negate_usecase

        self.run_unary(pyfunc, control_signed, samples, self.signed_types)
        self.run_unary(pyfunc, control_unsigned, samples, self.unsigned_types)

    def test_invert(self):
        def control_signed(a):
            tp = self.get_numpy_signed_upcast(a)
            return tp(~a)
        def control_unsigned(a):
            tp = self.get_numpy_unsigned_upcast(a)
            return tp(~a)
        samples = self.int_samples
        pyfunc = self.op.bitwise_not_usecase

        self.run_unary(pyfunc, control_signed, samples, self.signed_types)
        self.run_unary(pyfunc, control_unsigned, samples, self.unsigned_types)


class TestMixedIntsOperatorModule(TestMixedInts):

    op = FunctionalOperatorImpl


class TestStaticPower(TestCase):
    """
    Test the ** operator with a static exponent, to exercise a
    dedicated optimization.
    """

    def _check_pow(self, exponents, values):
        for exp in exponents:
            # test against non-static version of the @jit-ed function
            regular_func = LiteralOperatorImpl.pow_usecase
            static_func = make_static_power(exp)

            static_cfunc = jit(nopython=True)(static_func)
            regular_cfunc = jit(nopython=True)(regular_func)
            for v in values:
                try:
                    expected = regular_cfunc(v, exp)
                except ZeroDivisionError:
                    with self.assertRaises(ZeroDivisionError):
                        static_cfunc(v)
                else:
                    got = static_cfunc(v)
                    self.assertPreciseEqual(expected, got, prec='double')

    def test_int_values(self):
        exponents = [1, 2, 3, 5, 17, 0, -1, -2, -3]
        vals = [0, 1, 3, -1, -4, np.int8(-3), np.uint16(4)]

        self._check_pow(exponents, vals)

    def test_real_values(self):
        exponents = [1, 2, 3, 5, 17, 0, -1, -2, -3, 0x111111, -0x111112]
        vals = [1.5, 3.25, -1.25, np.float32(-2.0), float('inf'), float('nan')]

        self._check_pow(exponents, vals)

class TestStringConstComparison(TestCase):
    """
    Test comparison of string constants
    """
    def test_eq(self):
        def test_impl1():
            s = 'test'
            return s == 'test'

        def test_impl2():
            s = 'test1'
            return s == 'test'

        cfunc1 = jit(nopython=True)(test_impl1)
        cfunc2 = jit(nopython=True)(test_impl2)
        self.assertEqual(test_impl1(), cfunc1())
        self.assertEqual(test_impl2(), cfunc2())

    def test_neq(self):
        def test_impl1():
            s = 'test'
            return s != 'test'

        def test_impl2():
            s = 'test1'
            return s != 'test'

        cfunc1 = jit(nopython=True)(test_impl1)
        cfunc2 = jit(nopython=True)(test_impl2)
        self.assertEqual(test_impl1(), cfunc1())
        self.assertEqual(test_impl2(), cfunc2())

class TestBooleanLiteralOperators(TestCase):
    """
    Test operators with Boolean constants
    """
    def test_eq(self):

        def test_impl1(b):
            return a_val == b

        def test_impl2(a):
            return a == b_val

        def test_impl3():
            r1 = True == True
            r2 = True == False
            r3 = False == True
            r4 = False == False
            return (r1, r2, r3, r4)

        for a_val, b in itertools.product([True, False], repeat=2):
            cfunc1 = jit(nopython=True)(test_impl1)
            self.assertEqual(test_impl1(b), cfunc1(b))

        for a, b_val in itertools.product([True, False], repeat=2):
            cfunc2 = jit(nopython=True)(test_impl2)
            self.assertEqual(test_impl2(a), cfunc2(a))

        cfunc3 = jit(nopython=True)(test_impl3)
        self.assertEqual(test_impl3(), cfunc3())

    def test_ne(self):

        def test_impl1(b):
            return a_val != b

        def test_impl2(a):
            return a != b_val

        def test_impl3():
            r1 = True != True
            r2 = True != False
            r3 = False != True
            r4 = False != False
            return (r1, r2, r3, r4)

        for a_val, b in itertools.product([True, False], repeat=2):
            cfunc1 = jit(nopython=True)(test_impl1)
            self.assertEqual(test_impl1(b), cfunc1(b))

        for a, b_val in itertools.product([True, False], repeat=2):
            cfunc2 = jit(nopython=True)(test_impl2)
            self.assertEqual(test_impl2(a), cfunc2(a))

        cfunc3 = jit(nopython=True)(test_impl3)
        self.assertEqual(test_impl3(), cfunc3())

    def test_is(self):

        def test_impl1(b):
            return a_val is b

        def test_impl2():
            r1 = True is True
            r2 = True is False
            r3 = False is True
            r4 = False is False
            return (r1, r2, r3, r4)

        for a_val, b in itertools.product([True, False], repeat=2):
            cfunc1 = jit(nopython=True)(test_impl1)
            self.assertEqual(test_impl1(b), cfunc1(b))

        cfunc2 = jit(nopython=True)(test_impl2)
        self.assertEqual(test_impl2(), cfunc2())

    def test_not(self):

        def test_impl():
            a, b = False, True
            return (not a, not b)

        cfunc = jit(nopython=True)(test_impl)
        self.assertEqual(test_impl(), cfunc())

    def test_bool(self):

        def test_impl():
            a, b = False, True
            return (bool(a), bool(b))

        cfunc = jit(nopython=True)(test_impl)
        self.assertEqual(test_impl(), cfunc())

    def test_bool_to_str(self):

        def test_impl():
            a, b = False, True
            return (str(a), str(b))

        cfunc = jit(nopython=True)(test_impl)
        self.assertEqual(test_impl(), cfunc())


if __name__ == '__main__':
    unittest.main()
