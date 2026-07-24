import functools
import itertools
import sys
import warnings
import threading
import operator

import numpy as np

import unittest
from numba import guvectorize, njit, typeof, vectorize
from numba.core import types
from numba.np.numpy_support import from_dtype
from numba.core.errors import LoweringError, TypingError
from numba.tests.support import TestCase, MemoryLeakMixin
from numba.core.typing.npydecl import supported_ufuncs
from numba.np import numpy_support
from numba.core.registry import cpu_target
from numba.core.base import BaseContext
from numba.np import ufunc_db

is32bits = tuple.__itemsize__ == 4
iswindows = sys.platform.startswith('win32')


def _unimplemented(func):
    """An 'expectedFailure' like decorator that only expects compilation errors
    caused by unimplemented functions that fail in no-python mode"""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except TypingError:
            raise unittest._ExpectedFailure(sys.exc_info())
        raise unittest._UnexpectedSuccess


def _make_ufunc_usecase(ufunc):
    ldict = {}
    arg_str = ','.join(['a{0}'.format(i) for i in range(ufunc.nargs)])
    func_str = 'def fn({0}):\n    np.{1}({0})'.format(arg_str, ufunc.__name__)
    exec(func_str, globals(), ldict)
    fn = ldict['fn']
    fn.__name__ = '{0}_usecase'.format(ufunc.__name__)
    return fn


def _make_unary_ufunc_op_usecase(ufunc_op):
    ldict = {}
    exec("def fn(x):\n    return {0}(x)".format(ufunc_op), globals(), ldict)
    fn = ldict["fn"]
    fn.__name__ = "usecase_{0}".format(hash(ufunc_op))
    return fn


def _make_binary_ufunc_op_usecase(ufunc_op):
    ldict = {}
    exec("def fn(x,y):\n    return x{0}y".format(ufunc_op), globals(), ldict)
    fn = ldict["fn"]
    fn.__name__ = "usecase_{0}".format(hash(ufunc_op))
    return fn


def _make_inplace_ufunc_op_usecase(ufunc_op):
    """Generates a function to be compiled that performs an inplace operation

    ufunc_op can be a string like '+=' or a function like operator.iadd
    """
    if isinstance(ufunc_op, str):
        ldict = {}
        exec("def fn(x,y):\n    x{0}y".format(ufunc_op), globals(), ldict)
        fn = ldict["fn"]
        fn.__name__ = "usecase_{0}".format(hash(ufunc_op))
    else:
        def inplace_op(x, y):
            ufunc_op(x, y)
        fn = inplace_op
    return fn


def _as_dtype_value(tyargs, args):
    """Convert python values into numpy scalar objects.
    """
    return [np.dtype(str(ty)).type(val) for ty, val in zip(tyargs, args)]


class BaseUFuncTest(MemoryLeakMixin):

    def setUp(self):
        super(BaseUFuncTest, self).setUp()
        self.inputs = [
            (np.uint32(0), types.uint32),
            (np.uint32(1), types.uint32),
            (np.int32(-1), types.int32),
            (np.int32(0), types.int32),
            (np.int32(1), types.int32),
            (np.uint64(0), types.uint64),
            (np.uint64(1), types.uint64),
            (np.int64(-1), types.int64),
            (np.int64(0), types.int64),
            (np.int64(1), types.int64),

            (np.float32(-0.5), types.float32),
            (np.float32(0.0), types.float32),
            (np.float32(0.5), types.float32),

            (np.float64(-0.5), types.float64),
            (np.float64(0.0), types.float64),
            (np.float64(0.5), types.float64),

            (np.array([0,1], dtype='u4'), types.Array(types.uint32, 1, 'C')),
            (np.array([0,1], dtype='u8'), types.Array(types.uint64, 1, 'C')),
            (np.array([-1,0,1], dtype='i4'), types.Array(types.int32, 1, 'C')),
            (np.array([-1,0,1], dtype='i8'), types.Array(types.int64, 1, 'C')),
            (np.array([-0.5, 0.0, 0.5], dtype='f4'),
             types.Array(types.float32, 1, 'C')),
            (np.array([-0.5, 0.0, 0.5], dtype='f8'),
             types.Array(types.float64, 1, 'C')),

            (np.array([0,1], dtype=np.int8), types.Array(types.int8, 1, 'C')),
            (np.array([0,1], dtype=np.int16), types.Array(types.int16, 1, 'C')),
            (np.array([0,1], dtype=np.uint8), types.Array(types.uint8, 1, 'C')),
            (np.array([0,1], dtype=np.uint16),
             types.Array(types.uint16, 1, 'C')),
        ]

    @functools.lru_cache(maxsize=None)
    def _compile(self, pyfunc, args, nrt=False):
        # NOTE: to test the implementation of Numpy ufuncs, we disable
        # rewriting of array expressions.
        return njit(args, _nrt=nrt, no_rewrites=True)(pyfunc)

    def _determine_output_type(self, input_type, int_output_type=None,
                               float_output_type=None):
        ty = input_type
        if isinstance(ty, types.Array):
            ndim = ty.ndim
            ty = ty.dtype
        else:
            ndim = 1

        if ty in types.signed_domain:
            if int_output_type:
                output_type = types.Array(int_output_type, ndim, 'C')
            else:
                output_type = types.Array(ty, ndim, 'C')
        elif ty in types.unsigned_domain:
            if int_output_type:
                output_type = types.Array(int_output_type, ndim, 'C')
            else:
                output_type = types.Array(ty, ndim, 'C')
        else:
            if float_output_type:
                output_type = types.Array(float_output_type, ndim, 'C')
            else:
                output_type = types.Array(ty, ndim, 'C')
        return output_type


class BasicUFuncTest(BaseUFuncTest):
    def _make_ufunc_usecase(self, ufunc):
        return _make_ufunc_usecase(ufunc)

    def basic_ufunc_test(self, ufunc, skip_inputs=None, additional_inputs=None,
                         int_output_type=None, float_output_type=None,
                         kinds='ifc', positive_only=False):

        # Necessary to avoid some Numpy warnings being silenced, despite
        # the simplefilter() call below.
        if skip_inputs is None:
            skip_inputs = []
        if additional_inputs is None:
            additional_inputs = []
        self.reset_module_warnings(__name__)

        pyfunc = self._make_ufunc_usecase(ufunc)

        inputs = list(self.inputs) + additional_inputs

        for input_tuple in inputs:
            input_operand = input_tuple[0]
            input_type = input_tuple[1]

            is_tuple = isinstance(input_operand, tuple)
            if is_tuple:
                args = input_operand
            else:
                args = (input_operand,) * ufunc.nin

            if input_type in skip_inputs:
                continue
            if positive_only and np.any(args[0] < 0):
                continue

            # Some ufuncs don't allow all kinds of arguments
            if (args[0].dtype.kind not in kinds):
                continue

            output_type = self._determine_output_type(
                input_type, int_output_type, float_output_type)

            input_types = (input_type,) * ufunc.nin
            output_types = (output_type,) * ufunc.nout
            argtys = input_types + output_types
            cfunc = self._compile(pyfunc, argtys)

            if isinstance(args[0], np.ndarray):
                results = [
                    np.zeros(args[0].shape,
                             dtype=out_ty.dtype.name)
                    for out_ty in output_types
                ]
                expected = [
                    np.zeros(args[0].shape, dtype=out_ty.dtype.name)
                    for out_ty in output_types
                ]
            else:
                results = [
                    np.zeros(1, dtype=out_ty.dtype.name)
                    for out_ty in output_types
                ]
                expected = [
                    np.zeros(1, dtype=out_ty.dtype.name)
                    for out_ty in output_types
                ]

            invalid_flag = False
            with warnings.catch_warnings(record=True) as warnlist:
                warnings.simplefilter('always')
                pyfunc(*args, *expected)

                warnmsg = "invalid value encountered"
                for thiswarn in warnlist:

                    if (issubclass(thiswarn.category, RuntimeWarning)
                            and str(thiswarn.message).startswith(warnmsg)):
                        invalid_flag = True

            cfunc(*args, *results)

            for expected_i, result_i in zip(expected, results):
                msg = '\n'.join(["ufunc '{0}' failed",
                                 "inputs ({1}):", "{2}",
                                 "got({3})", "{4}",
                                 "expected ({5}):", "{6}"
                                 ]).format(ufunc.__name__,
                                           input_type, input_operand,
                                           output_type, result_i,
                                           expected_i.dtype, expected_i)
                try:
                    np.testing.assert_array_almost_equal(
                        expected_i, result_i,
                        decimal=5,
                        err_msg=msg)
                except AssertionError:
                    if invalid_flag:
                        # Allow output to mismatch for invalid input
                        print("Output mismatch for invalid input",
                              input_tuple, result_i, expected_i)
                    else:
                        raise

    def signed_unsigned_cmp_test(self, comparison_ufunc):
        self.basic_ufunc_test(comparison_ufunc)

        if numpy_support.numpy_version < (1, 25):
            return

        # Test additional implementations that specifically handle signed /
        # unsigned comparisons added in NumPy 1.25:
        # https://github.com/numpy/numpy/pull/23713
        additional_inputs = (
            (np.int64(-1), np.uint64(0)),
            (np.int64(-1), np.uint64(1)),
            (np.int64(0), np.uint64(0)),
            (np.int64(0), np.uint64(1)),
            (np.int64(1), np.uint64(0)),
            (np.int64(1), np.uint64(1)),

            (np.uint64(0), np.int64(-1)),
            (np.uint64(0), np.int64(0)),
            (np.uint64(0), np.int64(1)),
            (np.uint64(1), np.int64(-1)),
            (np.uint64(1), np.int64(0)),
            (np.uint64(1), np.int64(1)),

            (np.array([-1, -1, 0, 0, 1, 1], dtype=np.int64),
             np.array([0, 1, 0, 1, 0, 1], dtype=np.uint64)),

            (np.array([0, 1, 0, 1, 0, 1], dtype=np.uint64),
             np.array([-1, -1, 0, 0, 1, 1], dtype=np.int64))
        )

        pyfunc = self._make_ufunc_usecase(comparison_ufunc)

        for a, b in additional_inputs:
            input_types = (typeof(a), typeof(b))
            output_type = types.Array(types.bool_, 1, 'C')
            argtys = input_types + (output_type,)
            cfunc = self._compile(pyfunc, argtys)

            if isinstance(a, np.ndarray):
                result = np.zeros(a.shape, dtype=np.bool_)
            else:
                result = np.zeros(1, dtype=np.bool_)

            expected = np.zeros_like(result)

            pyfunc(a, b, expected)
            cfunc(a, b, result)
            np.testing.assert_equal(expected, result)


class TestUFuncs(BasicUFuncTest, TestCase):
    def basic_int_ufunc_test(self, name=None):
        skip_inputs = [
            types.float32,
            types.float64,
            types.Array(types.float32, 1, 'C'),
            types.Array(types.float64, 1, 'C'),
        ]
        self.basic_ufunc_test(name, skip_inputs=skip_inputs)

    ############################################################################
    # Math operations

    def test_add_ufunc(self):
        self.basic_ufunc_test(np.add)

    def test_subtract_ufunc(self):
        self.basic_ufunc_test(np.subtract)

    def test_multiply_ufunc(self):
        self.basic_ufunc_test(np.multiply)

    def test_divide_ufunc(self):
        # Bear in mind that in python3 divide IS true_divide
        # so the out type for int types will be a double
        int_out_type = None
        int_out_type = types.float64

        self.basic_ufunc_test(np.divide,
                              int_output_type=int_out_type)

    def test_logaddexp_ufunc(self):
        self.basic_ufunc_test(np.logaddexp, kinds='f')

    def test_logaddexp2_ufunc(self):
        self.basic_ufunc_test(np.logaddexp2, kinds='f')

    def test_true_divide_ufunc(self):
        self.basic_ufunc_test(np.true_divide,
                              int_output_type=types.float64)

    def test_floor_divide_ufunc(self):
        self.basic_ufunc_test(np.floor_divide)

    def test_negative_ufunc(self):
        # NumPy ufunc has bug with uint32 as input and int64 as output,
        # so skip uint32 input.
        skip_inputs = [types.Array(types.uint32, 1, 'C'), types.uint32]
        self.basic_ufunc_test(np.negative, int_output_type=types.int64,
                              skip_inputs=skip_inputs)

    def test_positive_ufunc(self):
        self.basic_ufunc_test(np.positive)

    def test_power_ufunc(self):
        self.basic_ufunc_test(np.power, positive_only=True)

    def test_float_power_ufunc(self):
        self.basic_ufunc_test(np.float_power, kinds="fc")

    def test_gcd_ufunc(self):
        self.basic_ufunc_test(np.gcd, kinds="iu")

    def test_lcm_ufunc(self):
        self.basic_ufunc_test(np.lcm, kinds="iu")

    def test_remainder_ufunc(self):
        self.basic_ufunc_test(np.remainder)

    def test_mod_ufunc(self):
        additional_inputs = [
            ((np.uint64(np.iinfo(np.uint64).max), np.uint64(16)), types.uint64)
        ]
        self.basic_ufunc_test(np.mod, kinds='ifcu',
                              additional_inputs=additional_inputs)

    def test_fmod_ufunc(self):
        self.basic_ufunc_test(np.fmod)

    def test_abs_ufunc(self, ufunc=np.abs):
        additional_inputs = [
            (np.uint32(np.iinfo(np.uint32).max), types.uint32),
            (np.uint64(np.iinfo(np.uint64).max), types.uint64),
            (np.float32(np.finfo(np.float32).min), types.float32),
            (np.float64(np.finfo(np.float64).min), types.float64),
        ]
        self.basic_ufunc_test(ufunc,
                              additional_inputs=additional_inputs)

    def test_absolute_ufunc(self):
        self.test_abs_ufunc(ufunc=np.absolute)

    def test_fabs_ufunc(self):
        self.basic_ufunc_test(np.fabs, kinds='f')

    def test_rint_ufunc(self):
        self.basic_ufunc_test(np.rint, kinds='cf')

    def test_sign_ufunc(self):
        self.basic_ufunc_test(np.sign)

    def test_conj_ufunc(self):
        self.basic_ufunc_test(np.conj)

    def test_exp_ufunc(self):
        self.basic_ufunc_test(np.exp, kinds='cf')

    def test_exp2_ufunc(self):
        self.basic_ufunc_test(np.exp2, kinds='cf')

    def test_log_ufunc(self):
        self.basic_ufunc_test(np.log, kinds='cf')

    def test_log2_ufunc(self):
        self.basic_ufunc_test(np.log2, kinds='cf')

    def test_log10_ufunc(self):
        self.basic_ufunc_test(np.log10, kinds='cf')

    def test_expm1_ufunc(self):
        self.basic_ufunc_test(np.expm1, kinds='cf')

    def test_log1p_ufunc(self):
        self.basic_ufunc_test(np.log1p, kinds='cf')

    def test_sqrt_ufunc(self):
        self.basic_ufunc_test(np.sqrt, kinds='cf')

    def test_square_ufunc(self):
        self.basic_ufunc_test(np.square)

    def test_cbrt_ufunc(self):
        self.basic_ufunc_test(np.cbrt, kinds='f')

    def test_reciprocal_ufunc(self):
        # reciprocal for integers doesn't make much sense and is problematic
        # in the case of division by zero, as an inf will overflow float to
        # int conversions, which is undefined behavior.
        to_skip = [types.Array(types.uint32, 1, 'C'), types.uint32,
                   types.Array(types.int32, 1, 'C'), types.int32,
                   types.Array(types.uint64, 1, 'C'), types.uint64,
                   types.Array(types.int64, 1, 'C'), types.int64]
        self.basic_ufunc_test(np.reciprocal, skip_inputs=to_skip)

    def test_conjugate_ufunc(self):
        self.basic_ufunc_test(np.conjugate)

    ############################################################################
    # Trigonometric Functions

    def test_sin_ufunc(self):
        self.basic_ufunc_test(np.sin, kinds='cf')

    def test_cos_ufunc(self):
        self.basic_ufunc_test(np.cos, kinds='cf')

    def test_tan_ufunc(self):
        self.basic_ufunc_test(np.tan, kinds='cf')

    def test_arcsin_ufunc(self):
        self.basic_ufunc_test(np.arcsin, kinds='cf')

    def test_arccos_ufunc(self):
        self.basic_ufunc_test(np.arccos, kinds='cf')

    def test_arctan_ufunc(self):
        self.basic_ufunc_test(np.arctan, kinds='cf')

    def test_arctan2_ufunc(self):
        self.basic_ufunc_test(np.arctan2, kinds='cf')

    def test_hypot_ufunc(self):
        self.basic_ufunc_test(np.hypot, kinds='f')

    def test_sinh_ufunc(self):
        self.basic_ufunc_test(np.sinh, kinds='cf')

    def test_cosh_ufunc(self):
        self.basic_ufunc_test(np.cosh, kinds='cf')

    def test_tanh_ufunc(self):
        self.basic_ufunc_test(np.tanh, kinds='cf')

    def test_arcsinh_ufunc(self):
        self.basic_ufunc_test(np.arcsinh, kinds='cf')

    def test_arccosh_ufunc(self):
        self.basic_ufunc_test(np.arccosh, kinds='cf')

    def test_arctanh_ufunc(self):
        # arctanh is only valid is only finite in the range ]-1, 1[
        # This means that for any of the integer types it will produce
        # conversion from infinity/-infinity to integer. That's undefined
        # behavior in C, so the results may vary from implementation to
        # implementation. This means that the result from the compiler
        # used to compile NumPy may differ from the result generated by
        # llvm. Skipping the integer types in this test avoids failed
        # tests because of this.
        to_skip = [types.Array(types.uint32, 1, 'C'), types.uint32,
                   types.Array(types.int32, 1, 'C'), types.int32,
                   types.Array(types.uint64, 1, 'C'), types.uint64,
                   types.Array(types.int64, 1, 'C'), types.int64]

        self.basic_ufunc_test(np.arctanh, skip_inputs=to_skip, kinds='cf')

    def test_deg2rad_ufunc(self):
        self.basic_ufunc_test(np.deg2rad, kinds='f')

    def test_rad2deg_ufunc(self):
        self.basic_ufunc_test(np.rad2deg, kinds='f')

    def test_degrees_ufunc(self):
        self.basic_ufunc_test(np.degrees, kinds='f')

    def test_radians_ufunc(self):
        self.basic_ufunc_test(np.radians, kinds='f')

    ############################################################################
    # Bit-twiddling Functions

    def test_bitwise_and_ufunc(self):
        self.basic_int_ufunc_test(np.bitwise_and)

    def test_bitwise_or_ufunc(self):
        self.basic_int_ufunc_test(np.bitwise_or)

    def test_bitwise_xor_ufunc(self):
        self.basic_int_ufunc_test(np.bitwise_xor)

    def test_invert_ufunc(self):
        self.basic_int_ufunc_test(np.invert)

    def test_bitwise_not_ufunc(self):
        self.basic_int_ufunc_test(np.bitwise_not)

    # Note: there is no entry for left_shift and right_shift as this harness
    #       is not valid for them. This is so because left_shift and right
    #       shift implementation in NumPy has undefined behavior (in C-parlance)
    #       when the second argument is a negative (or bigger than the number
    #       of bits) value.
    #       Also, right_shift for negative first arguments also relies on
    #       implementation defined behavior, although numba warantees "sane"
    #       behavior (arithmetic shifts on signed integers, logic shifts on
    #       unsigned integers).

    ############################################################################
    # Comparison functions
    def test_greater_ufunc(self):
        self.signed_unsigned_cmp_test(np.greater)

    def test_greater_equal_ufunc(self):
        self.signed_unsigned_cmp_test(np.greater_equal)

    def test_less_ufunc(self):
        self.signed_unsigned_cmp_test(np.less)

    def test_less_equal_ufunc(self):
        self.signed_unsigned_cmp_test(np.less_equal)

    def test_not_equal_ufunc(self):
        self.signed_unsigned_cmp_test(np.not_equal)

    def test_equal_ufunc(self):
        self.signed_unsigned_cmp_test(np.equal)

    def test_logical_and_ufunc(self):
        self.basic_ufunc_test(np.logical_and)

    def test_logical_or_ufunc(self):
        self.basic_ufunc_test(np.logical_or)

    def test_logical_xor_ufunc(self):
        self.basic_ufunc_test(np.logical_xor)

    def test_logical_not_ufunc(self):
        self.basic_ufunc_test(np.logical_not)

    def test_maximum_ufunc(self):
        self.basic_ufunc_test(np.maximum)

    def test_minimum_ufunc(self):
        self.basic_ufunc_test(np.minimum)

    def test_fmax_ufunc(self):
        self.basic_ufunc_test(np.fmax)

    def test_fmin_ufunc(self):
        self.basic_ufunc_test(np.fmin)

    ############################################################################
    # Floating functions

    def bool_additional_inputs(self):
        return [
            (np.array([True, False], dtype=np.bool_),
             types.Array(types.bool_, 1, 'C')),
        ]

    def test_isfinite_ufunc(self):
        self.basic_ufunc_test(
            np.isfinite, kinds='ifcb',
            additional_inputs=self.bool_additional_inputs(),
        )

    def test_isinf_ufunc(self):
        self.basic_ufunc_test(
            np.isinf, kinds='ifcb',
            additional_inputs=self.bool_additional_inputs(),
        )

    def test_isnan_ufunc(self):
        self.basic_ufunc_test(
            np.isnan, kinds='ifcb',
            additional_inputs=self.bool_additional_inputs(),
        )

    def test_signbit_ufunc(self):
        self.basic_ufunc_test(np.signbit)

    def test_copysign_ufunc(self):
        self.basic_ufunc_test(np.copysign, kinds='f')

    def test_nextafter_ufunc(self):
        self.basic_ufunc_test(np.nextafter, kinds='f')

    @_unimplemented
    def test_modf_ufunc(self):
        self.basic_ufunc_test(np.modf, kinds='f')

    # Note: there is no entry for ldexp as this harness isn't valid for this
    #       ufunc. this is so because ldexp requires heterogeneous inputs.
    #       However, this ufunc is tested by the TestLoopTypes test classes.

    @_unimplemented
    def test_frexp_ufunc(self):
        self.basic_ufunc_test(np.frexp, kinds='f')

    def test_floor_ufunc(self):
        self.basic_ufunc_test(np.floor, kinds='f')

    def test_ceil_ufunc(self):
        self.basic_ufunc_test(np.ceil, kinds='f')

    def test_trunc_ufunc(self):
        self.basic_ufunc_test(np.trunc, kinds='f')

    def test_spacing_ufunc(self):
        # additional input to check inf behaviour as Numba uses a different alg
        # to NumPy
        additional = [(np.array([np.inf, -np.inf], dtype=np.float64),
                       types.Array(types.float64, 1, 'C')),]
        self.basic_ufunc_test(np.spacing, kinds='f',
                              additional_inputs=additional)

    ############################################################################
    # Other tests

    def binary_ufunc_mixed_types_test(self, ufunc):
        ufunc_name = ufunc.__name__
        ufunc = _make_ufunc_usecase(ufunc)
        inputs1 = [
            (1, types.uint64),
            (-1, types.int64),
            (0.5, types.float64),

            (np.array([0, 1], dtype='u8'), types.Array(types.uint64, 1, 'C')),
            (np.array([-1, 1], dtype='i8'), types.Array(types.int64, 1, 'C')),
            (np.array([-0.5, 0.5], dtype='f8'),
             types.Array(types.float64, 1, 'C'))]

        inputs2 = inputs1

        output_types = [types.Array(types.int64, 1, 'C'),
                        types.Array(types.float64, 1, 'C')]

        pyfunc = ufunc

        for vals in itertools.product(inputs1, inputs2, output_types):
            input1, input2, output_type = vals

            input1_operand = input1[0]
            input1_type = input1[1]

            input2_operand = input2[0]
            input2_type = input2[1]

            # Skip division by unsigned int because of NumPy bugs
            if ufunc_name == 'divide' and (
                    input2_type == types.Array(types.uint32, 1, 'C') or
                    input2_type == types.Array(types.uint64, 1, 'C')):
                continue

            # Skip some subtraction tests because of NumPy bugs
            if (ufunc_name == 'subtract'
                    and input1_type == types.Array(types.uint32, 1, 'C')
                    and input2_type == types.uint32
                    and types.Array(types.int64, 1, 'C')):
                continue
            if (ufunc_name == 'subtract'
                    and input1_type == types.Array(types.uint32, 1, 'C')
                    and input2_type == types.uint64
                    and types.Array(types.int64, 1, 'C')):
                continue

            if ((isinstance(input1_type, types.Array) or
                    isinstance(input2_type, types.Array)) and
                    not isinstance(output_type, types.Array)):
                continue

            args = (input1_type, input2_type, output_type)
            cfunc = self._compile(pyfunc, args)

            if isinstance(input1_operand, np.ndarray):
                result = np.zeros(input1_operand.size,
                                  dtype=output_type.dtype.name)
                expected = np.zeros(input1_operand.size,
                                    dtype=output_type.dtype.name)
            elif isinstance(input2_operand, np.ndarray):
                result = np.zeros(input2_operand.size,
                                  dtype=output_type.dtype.name)
                expected = np.zeros(input2_operand.size,
                                    dtype=output_type.dtype.name)
            else:
                result = np.zeros(1, dtype=output_type.dtype.name)
                expected = np.zeros(1, dtype=output_type.dtype.name)

            cfunc(input1_operand, input2_operand, result)
            pyfunc(input1_operand, input2_operand, expected)

            scalar_type = getattr(output_type, 'dtype', output_type)
            prec = ('single'
                    if scalar_type in (types.float32, types.complex64)
                    else 'double')
            self.assertPreciseEqual(expected, result, prec=prec)

    def test_broadcasting(self):

        # Test unary ufunc
        pyfunc = _make_ufunc_usecase(np.negative)

        input_operands = [
            np.arange(3, dtype='u8'),
            np.arange(3, dtype='u8').reshape(3,1),
            np.arange(3, dtype='u8').reshape(1,3),
            np.arange(3, dtype='u8').reshape(3,1),
            np.arange(3, dtype='u8').reshape(1,3),
            np.arange(3 * 3, dtype='u8').reshape(3,3)]

        output_operands = [
            np.zeros(3 * 3, dtype='i8').reshape(3,3),
            np.zeros(3 * 3, dtype='i8').reshape(3,3),
            np.zeros(3 * 3, dtype='i8').reshape(3,3),
            np.zeros(3 * 3 * 3, dtype='i8').reshape(3,3,3),
            np.zeros(3 * 3 * 3, dtype='i8').reshape(3,3,3),
            np.zeros(3 * 3 * 3, dtype='i8').reshape(3,3,3)]

        for x, result in zip(input_operands, output_operands):

            input_type = types.Array(types.uint64, x.ndim, 'C')
            output_type = types.Array(types.int64, result.ndim, 'C')
            args = (input_type, output_type)

            cfunc = self._compile(pyfunc, args)

            expected = np.zeros(result.shape, dtype=result.dtype)
            np.negative(x, expected)

            cfunc(x, result)

            self.assertPreciseEqual(result, expected)

        # Test binary ufunc
        pyfunc = _make_ufunc_usecase(np.add)

        input1_operands = [
            np.arange(3, dtype='u8'),
            np.arange(3 * 3, dtype='u8').reshape(3,3),
            np.arange(3 * 3 * 3, dtype='u8').reshape(3,3,3),
            np.arange(3, dtype='u8').reshape(3,1),
            np.arange(3, dtype='u8').reshape(1,3),
            np.arange(3, dtype='u8').reshape(3,1,1),
            np.arange(3 * 3, dtype='u8').reshape(3,3,1),
            np.arange(3 * 3, dtype='u8').reshape(3,1,3),
            np.arange(3 * 3, dtype='u8').reshape(1,3,3)]

        input2_operands = input1_operands

        for x, y in itertools.product(input1_operands, input2_operands):

            input1_type = types.Array(types.uint64, x.ndim, 'C')
            input2_type = types.Array(types.uint64, y.ndim, 'C')
            output_type = types.Array(types.uint64, max(x.ndim, y.ndim), 'C')
            args = (input1_type, input2_type, output_type)

            cfunc = self._compile(pyfunc, args)

            expected = np.add(x, y)
            result = np.zeros(expected.shape, dtype='u8')

            cfunc(x, y, result)
            self.assertPreciseEqual(result, expected)

    def test_implicit_output_npm(self):
        # Test for Issue #1078 (https://github.com/numba/numba/issues/1078) -
        # ensures that the output of a ufunc is an array.
        arr_ty = types.Array(types.uint64, 1, 'C')
        sig = (arr_ty, arr_ty)

        @njit((arr_ty, arr_ty))
        def myadd(a0, a1):
            return np.add(a0, a1)

        self.assertEqual(myadd.overloads[sig].signature.return_type, arr_ty)

    def test_broadcast_implicit_output_npm_nrt(self):
        def pyfunc(a0, a1):
            return np.add(a0, a1)

        input1_operands = [
            np.arange(3, dtype='u8'),
            np.arange(3 * 3, dtype='u8').reshape(3,3),
            np.arange(3 * 3 * 3, dtype='u8').reshape(3,3,3),
            np.arange(3, dtype='u8').reshape(3,1),
            np.arange(3, dtype='u8').reshape(1,3),
            np.arange(3, dtype='u8').reshape(3,1,1),
            np.arange(3 * 3, dtype='u8').reshape(3,3,1),
            np.arange(3 * 3, dtype='u8').reshape(3,1,3),
            np.arange(3 * 3, dtype='u8').reshape(1,3,3)]

        input2_operands = input1_operands

        for x, y in itertools.product(input1_operands, input2_operands):

            input1_type = types.Array(types.uint64, x.ndim, 'C')
            input2_type = types.Array(types.uint64, y.ndim, 'C')
            args = (input1_type, input2_type)

            cfunc = self._compile(pyfunc, args, nrt=True)

            expected = np.add(x, y)
            result = cfunc(x, y)
            np.testing.assert_array_equal(expected, result)

    def test_implicit_output_layout_binary(self):
        def pyfunc(a0, a1):
            return np.add(a0, a1)

        # C layout
        X = np.linspace(0, 1, 20).reshape(4, 5)
        # F layout
        Y = np.array(X, order='F')
        # A layout
        Z = X.reshape(5, 4).T[0]

        Xty = typeof(X)
        assert X.flags.c_contiguous and Xty.layout == 'C'
        Yty = typeof(Y)
        assert Y.flags.f_contiguous and Yty.layout == 'F'
        Zty = typeof(Z)
        assert Zty.layout == 'A'
        assert not Z.flags.c_contiguous
        assert not Z.flags.f_contiguous

        testcases = list(itertools.permutations([X, Y, Z], 2))
        testcases += [(X, X)]
        testcases += [(Y, Y)]
        testcases += [(Z, Z)]

        for arg0, arg1 in testcases:
            args = (typeof(arg0), typeof(arg1))
            cfunc = self._compile(pyfunc, args, nrt=True)
            expected = pyfunc(arg0, arg1)
            result = cfunc(arg0, arg1)

            self.assertEqual(expected.flags.c_contiguous,
                             result.flags.c_contiguous)
            self.assertEqual(expected.flags.f_contiguous,
                             result.flags.f_contiguous)
            np.testing.assert_array_equal(expected, result)

    def test_implicit_output_layout_unary(self):
        def pyfunc(a0):
            return np.sqrt(a0)

        # C layout
        X = np.linspace(0, 1, 20).reshape(4, 5)
        # F layout
        Y = np.array(X, order='F')
        # A layout
        Z = X.reshape(5, 4).T[0]

        Xty = typeof(X)
        assert X.flags.c_contiguous and Xty.layout == 'C'
        Yty = typeof(Y)
        assert Y.flags.f_contiguous and Yty.layout == 'F'
        Zty = typeof(Z)
        assert Zty.layout == 'A'
        assert not Z.flags.c_contiguous
        assert not Z.flags.f_contiguous

        for arg0 in [X, Y, Z]:
            args = (typeof(arg0),)
            cfunc = self._compile(pyfunc, args, nrt=True)
            expected = pyfunc(arg0)
            result = cfunc(arg0)

            self.assertEqual(expected.flags.c_contiguous,
                             result.flags.c_contiguous)
            self.assertEqual(expected.flags.f_contiguous,
                             result.flags.f_contiguous)
            np.testing.assert_array_equal(expected, result)


class TestArrayOperators(BaseUFuncTest, TestCase):

    def _check_results(self, expected, got):
        self.assertEqual(expected.dtype.kind, got.dtype.kind)
        np.testing.assert_array_almost_equal(expected, got)

    def unary_op_test(self, operator, nrt=True,
                      skip_inputs=None, additional_inputs=None,
                      int_output_type=None, float_output_type=None):
        if skip_inputs is None:
            skip_inputs = []
        if additional_inputs is None:
            additional_inputs = []
        operator_func = _make_unary_ufunc_op_usecase(operator)
        inputs = list(self.inputs)
        inputs.extend(additional_inputs)
        pyfunc = operator_func
        for input_tuple in inputs:
            input_operand, input_type = input_tuple

            if ((input_type in skip_inputs) or
                    (not isinstance(input_type, types.Array))):
                continue

            cfunc = self._compile(pyfunc, (input_type,), nrt=nrt)
            expected = pyfunc(input_operand)
            got = cfunc(input_operand)
            self._check_results(expected, got)

    def binary_op_test(self, operator, nrt=True,
                       skip_inputs=None, additional_inputs=None,
                       int_output_type=None, float_output_type=None,
                       positive_rhs=False):
        if skip_inputs is None:
            skip_inputs = []
        if additional_inputs is None:
            additional_inputs = []
        operator_func = _make_binary_ufunc_op_usecase(operator)
        inputs = list(self.inputs)
        inputs.extend(additional_inputs)
        pyfunc = operator_func
        # when generating arbitrary sequences, we use a fixed seed
        # for deterministic testing
        random_state = np.random.RandomState(1)
        for input_tuple in inputs:
            input_operand1, input_type = input_tuple
            input_dtype = numpy_support.as_dtype(
                getattr(input_type, "dtype", input_type))
            input_type1 = input_type

            if input_type in skip_inputs:
                continue

            if positive_rhs:
                zero = np.zeros(1, dtype=input_dtype)[0]
            # If we only use two scalars, the code generator will not
            # select the ufunctionalized operator, so we mix it up.
            if isinstance(input_type, types.Array):
                input_operand0 = input_operand1
                input_type0 = input_type
                if positive_rhs and np.any(input_operand1 < zero):
                    continue
            else:
                input_operand0 = (random_state.uniform(0, 100, 10)).astype(
                    input_dtype)
                input_type0 = typeof(input_operand0)
                if positive_rhs and input_operand1 < zero:
                    continue

            args = (input_type0, input_type1)
            cfunc = self._compile(pyfunc, args, nrt=nrt)
            expected = pyfunc(input_operand0, input_operand1)
            got = cfunc(input_operand0, input_operand1)
            self._check_results(expected, got)

    def bitwise_additional_inputs(self):
        # For bitwise operators, we want to check the results for boolean
        # arrays (see #1813).
        return [
            (True, types.boolean),
            (False, types.boolean),
            (np.array([True, False]), types.Array(types.boolean, 1, 'C')),
        ]

    def binary_int_op_test(self, *args, **kws):
        skip_inputs = kws.setdefault('skip_inputs', [])
        skip_inputs += [
            types.float32, types.float64,
            types.Array(types.float32, 1, 'C'),
            types.Array(types.float64, 1, 'C'),
        ]
        return self.binary_op_test(*args, **kws)

    def binary_bitwise_op_test(self, *args, **kws):
        additional_inputs = kws.setdefault('additional_inputs', [])
        additional_inputs += self.bitwise_additional_inputs()
        return self.binary_int_op_test(*args, **kws)

    def inplace_op_test(self, operator, lhs_values, rhs_values,
                        lhs_dtypes, rhs_dtypes, precise=True):
        operator_func = _make_inplace_ufunc_op_usecase(operator)
        pyfunc = operator_func

        if precise:
            assertion = self.assertPreciseEqual
        else:
            assertion = np.testing.assert_allclose

        # The left operand can only be an array, while the right operand
        # can be either an array or a scalar
        lhs_inputs = [np.array(lhs_values, dtype=dtype)
                      for dtype in lhs_dtypes]

        rhs_arrays = [np.array(rhs_values, dtype=dtype)
                      for dtype in rhs_dtypes]
        rhs_scalars = [dtype(v) for v in rhs_values for dtype in rhs_dtypes]
        rhs_inputs = rhs_arrays + rhs_scalars

        for lhs, rhs in itertools.product(lhs_inputs, rhs_inputs):
            lhs_type = typeof(lhs)
            rhs_type = typeof(rhs)
            args = (lhs_type, rhs_type)
            cfunc = self._compile(pyfunc, args)
            expected = lhs.copy()
            pyfunc(expected, rhs)
            got = lhs.copy()
            cfunc(got, rhs)
            assertion(got, expected)

    def inplace_float_op_test(self, operator, lhs_values, rhs_values,
                              precise=True):
        # Also accept integer inputs for the right operand (they should
        # be converted to float).
        return self.inplace_op_test(operator, lhs_values, rhs_values,
                                    (np.float32, np.float64),
                                    (np.float32, np.float64, np.int64),
                                    precise=precise)

    def inplace_int_op_test(self, operator, lhs_values, rhs_values):
        self.inplace_op_test(operator, lhs_values, rhs_values,
                             (np.int16, np.int32, np.int64),
                             (np.int16, np.uint32))

    def inplace_bitwise_op_test(self, operator, lhs_values, rhs_values):
        self.inplace_int_op_test(operator, lhs_values, rhs_values)
        self.inplace_op_test(operator, lhs_values, rhs_values,
                             (np.bool_,), (np.bool_, np.bool_))

    # ____________________________________________________________
    # Unary operators

    def test_unary_positive_array_op(self):
        self.unary_op_test('+')

    def test_unary_negative_array_op(self):
        self.unary_op_test('-')

    def test_unary_invert_array_op(self):
        self.unary_op_test('~',
                           skip_inputs=[types.float32, types.float64,
                                        types.Array(types.float32, 1, 'C'),
                                        types.Array(types.float64, 1, 'C')],
                           additional_inputs=self.bitwise_additional_inputs())

    # ____________________________________________________________
    # Inplace operators

    def test_inplace_add(self):
        self.inplace_float_op_test('+=', [-1, 1.5, 3], [-5, 0, 2.5])
        self.inplace_float_op_test(operator.iadd, [-1, 1.5, 3], [-5, 0, 2.5])

    def test_inplace_sub(self):
        self.inplace_float_op_test('-=', [-1, 1.5, 3], [-5, 0, 2.5])
        self.inplace_float_op_test(operator.isub, [-1, 1.5, 3], [-5, 0, 2.5])

    def test_inplace_mul(self):
        self.inplace_float_op_test('*=', [-1, 1.5, 3], [-5, 0, 2.5])
        self.inplace_float_op_test(operator.imul, [-1, 1.5, 3], [-5, 0, 2.5])

    def test_inplace_floordiv(self):
        self.inplace_float_op_test('//=', [-1, 1.5, 3], [-5, 1.25, 2.5])
        self.inplace_float_op_test(operator.ifloordiv, [-1, 1.5, 3],
                                   [-5, 1.25, 2.5])

    def test_inplace_div(self):
        self.inplace_float_op_test('/=', [-1, 1.5, 3], [-5, 0, 2.5])
        self.inplace_float_op_test(operator.itruediv, [-1, 1.5, 3],
                                   [-5, 1.25, 2.5])

    def test_inplace_remainder(self):
        self.inplace_float_op_test('%=', [-1, 1.5, 3], [-5, 2, 2.5])
        self.inplace_float_op_test(operator.imod, [-1, 1.5, 3], [-5, 2, 2.5])

    def test_inplace_pow(self):
        self.inplace_float_op_test('**=', [-1, 1.5, 3], [-5, 2, 2.5],
                                   precise=False)
        self.inplace_float_op_test(operator.ipow, [-1, 1.5, 3], [-5, 2, 2.5],
                                   precise=False)

    def test_inplace_and(self):
        self.inplace_bitwise_op_test('&=', [0, 1, 2, 3, 51],
                                     [0, 13, 16, 42, 255])
        self.inplace_bitwise_op_test(operator.iand, [0, 1, 2, 3, 51],
                                     [0, 13, 16, 42, 255])

    def test_inplace_or(self):
        self.inplace_bitwise_op_test('|=', [0, 1, 2, 3, 51],
                                     [0, 13, 16, 42, 255])
        self.inplace_bitwise_op_test(operator.ior, [0, 1, 2, 3, 51],
                                     [0, 13, 16, 42, 255])

    def test_inplace_xor(self):
        self.inplace_bitwise_op_test('^=', [0, 1, 2, 3, 51],
                                     [0, 13, 16, 42, 255])
        self.inplace_bitwise_op_test(operator.ixor, [0, 1, 2, 3, 51],
                                     [0, 13, 16, 42, 255])

    def test_inplace_lshift(self):
        self.inplace_int_op_test('<<=', [0, 5, -10, -51], [0, 1, 4, 14])
        self.inplace_int_op_test(operator.ilshift, [0, 5, -10, -51],
                                 [0, 1, 4, 14])

    def test_inplace_rshift(self):
        self.inplace_int_op_test('>>=', [0, 5, -10, -51], [0, 1, 4, 14])
        self.inplace_int_op_test(operator.irshift, [0, 5, -10, -51],
                                 [0, 1, 4, 14])

    def test_unary_positive_array_op_2(self):
        '''
        Verify that the unary positive operator copies values, and doesn't
        just alias to the input array (mirrors normal Numpy/Python
        interaction behavior).
        '''
        # Test originally from @gmarkall
        def f(a1):
            a2 = +a1
            a1[0] = 3
            a2[1] = 4
            return a2

        a1 = np.zeros(10)
        a2 = f(a1)
        self.assertTrue(a1[0] != a2[0] and a1[1] != a2[1])
        a3 = np.zeros(10)
        a4 = njit(f)(a3)
        self.assertTrue(a3[0] != a4[0] and a3[1] != a4[1])
        np.testing.assert_array_equal(a1, a3)
        np.testing.assert_array_equal(a2, a4)

    # ____________________________________________________________
    # Binary operators

    def test_add_array_op(self):
        self.binary_op_test('+')

    def test_subtract_array_op(self):
        self.binary_op_test('-')

    def test_multiply_array_op(self):
        self.binary_op_test('*')

    def test_divide_array_op(self):
        int_out_type = None
        int_out_type = types.float64
        self.binary_op_test('/', int_output_type=int_out_type)

    def test_floor_divide_array_op(self):
        # Avoid floating-point zeros as x // 0.0 can have varying results
        # depending on the algorithm (which changed across Numpy versions)
        self.inputs = [
            (np.uint32(1), types.uint32),
            (np.int32(-2), types.int32),
            (np.int32(0), types.int32),
            (np.uint64(4), types.uint64),
            (np.int64(-5), types.int64),
            (np.int64(0), types.int64),

            (np.float32(-0.5), types.float32),
            (np.float32(1.5), types.float32),

            (np.float64(-2.5), types.float64),
            (np.float64(3.5), types.float64),

            (np.array([1,2], dtype='u4'), types.Array(types.uint32, 1, 'C')),
            (np.array([3,4], dtype='u8'), types.Array(types.uint64, 1, 'C')),
            (np.array([-1,1,5], dtype='i4'), types.Array(types.int32, 1, 'C')),
            (np.array([-1,1,6], dtype='i8'), types.Array(types.int64, 1, 'C')),
            (np.array([-0.5, 1.5], dtype='f4'),
             types.Array(types.float32, 1, 'C')),
            (np.array([-2.5, 3.5], dtype='f8'),
             types.Array(types.float64, 1, 'C')),
        ]
        self.binary_op_test('//')

    def test_remainder_array_op(self):
        self.binary_op_test('%')

    def test_power_array_op(self):
        self.binary_op_test('**', positive_rhs=True)

    def test_left_shift_array_op(self):
        self.binary_int_op_test('<<', positive_rhs=True)

    def test_right_shift_array_op(self):
        self.binary_int_op_test('>>', positive_rhs=True)

    def test_bitwise_and_array_op(self):
        self.binary_bitwise_op_test('&')

    def test_bitwise_or_array_op(self):
        self.binary_bitwise_op_test('|')

    def test_bitwise_xor_array_op(self):
        self.binary_bitwise_op_test('^')

    def test_equal_array_op(self):
        self.binary_op_test('==')

    def test_greater_array_op(self):
        self.binary_op_test('>')

    def test_greater_equal_array_op(self):
        self.binary_op_test('>=')

    def test_less_array_op(self):
        self.binary_op_test('<')

    def test_less_equal_array_op(self):
        self.binary_op_test('<=')

    def test_not_equal_array_op(self):
        self.binary_op_test('!=')


class TestScalarUFuncs(TestCase):
    """check the machinery of ufuncs works when the result is an scalar.
    These are not exhaustive because:
    - the machinery to support this case is the same for all the functions of a
      given arity.
    - the result of the inner function itself is already tested in TestUFuncs
    """

    def run_ufunc(self, pyfunc, arg_types, arg_values):
        for tyargs, args in zip(arg_types, arg_values):
            cfunc = njit(tyargs)(pyfunc)
            got = cfunc(*args)
            expected = pyfunc(*_as_dtype_value(tyargs, args))

            msg = 'for args {0} typed {1}'.format(args, tyargs)

            # note: due to semantics of ufuncs, thing like adding a int32 to a
            # uint64 results in doubles (as neither int32 can be cast safely
            # to uint64 nor vice-versa, falling back to using the float version.
            # Modify in those cases the expected value (the numpy version does
            # not use typed integers as inputs so its result is an integer)
            special = set([
                (types.int32, types.uint64),
                (types.uint64, types.int32),
                (types.int64, types.uint64),
                (types.uint64, types.int64)
            ])
            if tyargs in special:
                expected = float(expected)
            else:
                # The numba version of scalar ufuncs return an actual value that
                # gets converted to a Python type, instead of using NumPy
                # scalars.  although in python 2 NumPy scalars are considered
                # and instance of the appropriate python type, in python 3 that
                # is no longer the case.  This is why the expected result is
                # casted to the appropriate Python type (which is actually the
                # expected behavior of the ufunc translation)
                if np.issubdtype(expected.dtype, np.inexact):
                    expected = float(expected)
                elif np.issubdtype(expected.dtype, np.integer):
                    expected = int(expected)
                elif np.issubdtype(expected.dtype, np.bool_):
                    expected = bool(expected)

            alltypes = tyargs + (cfunc.overloads[tyargs].signature.return_type,)

            # select the appropriate precision for comparison: note that an
            # argument typed at a lower precision can introduce precision
            # problems. For this reason the argument types must be taken into
            # account.
            if any([t == types.float32 for t in alltypes]):
                prec = 'single'
            elif any([t == types.float64 for t in alltypes]):
                prec = 'double'
            else:
                prec = 'exact'

            self.assertPreciseEqual(got, expected, msg=msg, prec=prec)

    def test_scalar_unary_ufunc(self):
        def _func(x):
            return np.sqrt(x)

        vals = [(2,), (2,), (1,), (2,), (.1,), (.2,)]
        tys = [(types.int32,), (types.uint32,),
               (types.int64,), (types.uint64,),
               (types.float32,), (types.float64,)]
        self.run_ufunc(_func, tys, vals)

    def test_scalar_binary_uniform_ufunc(self):
        def _func(x,y):
            return np.add(x,y)

        vals = [2, 2, 1, 2, .1, .2]
        tys = [types.int32, types.uint32,
               types.int64, types.uint64, types.float32, types.float64]
        self.run_ufunc(_func, zip(tys, tys), zip(vals, vals))

    def test_scalar_binary_mixed_ufunc(self):
        def _func(x,y):
            return np.add(x,y)

        vals = [2, 2, 1, 2, .1, .2]
        tys = [types.int32, types.uint32,
               types.int64, types.uint64,
               types.float32, types.float64]
        self.run_ufunc(_func, itertools.product(tys, tys),
                       itertools.product(vals, vals))


class TestUfuncIssues(TestCase):

    def test_issue_651(self):
        # Exercise the code path to make sure this does not fail
        @vectorize(["(float64,float64)"])
        def foo(x1, x2):
            return np.add(x1, x2) + np.add(x1, x2)

        a = np.arange(10, dtype='f8')
        b = np.arange(10, dtype='f8')
        self.assertPreciseEqual(foo(a, b), (a + b) + (a + b))

    def test_issue_2006(self):
        """
        <float32 ** int> should return float32, not float64.
        """
        def foo(x, y):
            return np.power(x, y)
        pyfunc = foo
        cfunc = njit(pyfunc)

        def check(x, y):
            got = cfunc(x, y)
            np.testing.assert_array_almost_equal(got, pyfunc(x, y))
            # Check the power operation conserved the input's dtype
            # (this is different from Numpy, whose behaviour depends on
            #  the *values* of the arguments -- see PyArray_CanCastArrayTo).
            self.assertEqual(got.dtype, x.dtype)

        xs = [np.float32([1, 2, 3]), np.complex64([1j, 2, 3 - 3j])]
        for x in xs:
            check(x, 3)
            check(x, np.uint64(3))
            check(x, np.int64([2, 2, 3]))


class _LoopTypesTester(TestCase):
    """Test code generation for the different loop types defined by ufunc.

    This test relies on class variables to configure the test. Subclasses
    of this class can just override some of these variables to check other
    ufuncs in a different compilation context. The variables supported are:

    _funcs: the ufuncs to test
    _skip_types: letter types that force skipping the loop when testing
                 if present in the NumPy ufunc signature.
    _supported_types: only test loops where all the types in the loop
                      signature are in this collection. If unset, all.

    Note that both, _skip_types and _supported_types must be met for a loop
    to be tested.

    The NumPy ufunc signature has a form like 'ff->f' (for a binary ufunc
    loop taking 2 floats and resulting in a float). In a NumPy ufunc object
    you can get a list of supported signatures by accessing the attribute
    'types'.
    """
    _skip_types = 'OegG'

    # Allowed deviation between Numpy and Numba results
    _ulps = {('arccos', 'F'): 2,
             ('arcsin', 'D'): 4,
             ('arcsin', 'F'): 4,
             ('log10', 'D'): 5,
             ('tanh', 'F'): 2,
             ('cbrt', 'd'): 2,
             ('logaddexp2', 'd'): 2,
             }

    def _arg_for_type(self, a_letter_type, index=0):
        """return a suitable array argument for testing the letter type"""
        # Note all possible arrays must have the same size, since they
        # may be used as inputs to the same func.
        if a_letter_type in 'bhilq':
            # an integral
            return np.array([1, 4, 0, -2], dtype=a_letter_type)
        if a_letter_type in 'BHILQ':
            return np.array([1, 2, 4, 0], dtype=a_letter_type)
        elif a_letter_type in '?':
            # a boolean
            return np.array([True, False, False, True], dtype=a_letter_type)
        elif a_letter_type[0] == 'm':
            # timedelta64
            if len(a_letter_type) == 1:
                a_letter_type = 'm8[D]'
            return np.array([2, -3, 'NaT', 0], dtype=a_letter_type)
        elif a_letter_type[0] == 'M':
            # datetime64
            if len(a_letter_type) == 1:
                a_letter_type = 'M8[D]'
            return np.array(['Nat', 1, 25, 0], dtype=a_letter_type)
        elif a_letter_type in 'fd':
            # floating point
            return np.array([1.5, -3.5, 0.0, float('nan')],
                            dtype=a_letter_type)
        elif a_letter_type in 'FD':
            # complex
            if sys.platform != 'win32':
                # Other platforms have better handling of negative zeros,
                # test them
                negzero = -(0.0 + 1.0j)
            else:
                negzero = 0.0 - 1.0j
            return np.array([negzero, 1.5 + 1.5j, 1j * float('nan'), 0j],
                            dtype=a_letter_type)
        else:
            raise RuntimeError("type %r not understood" % (a_letter_type,))

    def _check_loop(self, fn, ufunc, loop):
        # the letter types for the args
        letter_types = loop[:ufunc.nin] + loop[-ufunc.nout:]

        # ignore the loops containing an object argument. They will always
        # fail in no python mode. Usually the last loop in ufuncs is an all
        # object fallback
        supported_types = getattr(self, '_supported_types', [])
        if (supported_types and
                any(l not in supported_types for l in letter_types)):
            return
        skip_types = getattr(self, '_skip_types', [])
        if any(l in skip_types for l in letter_types):
            return
        # if the test case requires some types to be present, skip loops
        # not involving any of those types.
        required_types = getattr(self, '_required_types', [])
        if required_types and not any(l in letter_types
                                      for l in required_types):
            return

        self._check_ufunc_with_dtypes(fn, ufunc, letter_types)

    def _check_ufunc_with_dtypes(self, fn, ufunc, dtypes):
        # Arrays created with datetime and timedelta types (e.g. with np.array)
        # will have units, so in order to ensure that the dtypes of arguments
        # match the dtypes in the signature, we add units to unitless datetime
        # and timedelta types. This corresponds with the addition of units in
        # _arg_for_type() above.
        dtypes_with_units = []
        for t in dtypes:
            if t in ('m', 'M'):
                t = t + '8[D]'
            dtypes_with_units.append(t)

        arg_dty = [np.dtype(t) for t in dtypes_with_units]
        arg_nbty = tuple([types.Array(from_dtype(t), 1, 'C') for t in arg_dty])
        cfunc = njit(arg_nbty)(fn)

        # Ensure a good mix of input values
        c_args = [self._arg_for_type(t, index=index).repeat(2)
                  for index, t in enumerate(dtypes)]
        for arr in c_args:
            self.random.shuffle(arr)
        py_args = [a.copy() for a in c_args]

        cfunc(*c_args)
        fn(*py_args)

        # Check each array (including inputs, to ensure they weren't
        # mutated).
        for dtype, py_arg, c_arg in zip(arg_dty, py_args, c_args):
            py_arg, c_arg = self._fixup_results(dtype, py_arg, c_arg)
            typechar = c_arg.dtype.char
            ulps = self._ulps.get((ufunc.__name__, typechar), 1)
            prec = 'single' if typechar in 'fF' else 'exact'
            prec = 'double' if typechar in 'dD' else prec
            msg = '\n'.join(["ufunc '{0}' arrays differ ({1}):",
                             "args: {2}", "expected {3}", "got {4}"])
            msg = msg.format(ufunc.__name__, c_args, prec, py_arg, c_arg)
            self.assertPreciseEqual(py_arg, c_arg, prec=prec, msg=msg,
                                    ulps=ulps)

    def _fixup_results(self, dtype, py_arg, c_arg):
        return py_arg, c_arg

    @classmethod
    def _check_ufunc_loops(cls, ufunc):
        for loop in ufunc.types:
            cls._inject_test(ufunc, loop)

    @classmethod
    def _inject_test(cls, ufunc, loop):
        def test_template(self):
            fn = _make_ufunc_usecase(ufunc)
            self._check_loop(fn, ufunc, loop)
        setattr(cls, "test_{0}_{1}".format(ufunc.__name__,
                                           loop.replace('->', '_')),
                test_template)

    @classmethod
    def autogenerate(cls):
        for ufunc in cls._ufuncs:
            cls._check_ufunc_loops(ufunc)


class TestLoopTypesInt(_LoopTypesTester):
    _ufuncs = supported_ufuncs[:]
    # reciprocal and power need a special test due to issue #757
    _ufuncs.remove(np.power)
    _ufuncs.remove(np.reciprocal)
    _ufuncs.remove(np.left_shift) # has its own test class
    _ufuncs.remove(np.right_shift) # has its own test class
    # special test for bool subtract/negative
    _ufuncs.remove(np.subtract)
    _ufuncs.remove(np.negative)
    _required_types = '?bBhHiIlLqQ'
    _skip_types = 'fdFDmMO' + _LoopTypesTester._skip_types


TestLoopTypesInt.autogenerate()


class TestLoopTypesSubtractAndNegative(_LoopTypesTester):
    _ufuncs = [np.subtract, np.negative]
    _required_types = '?bBhHiIlLqQfdFD'
    _skip_types = 'mMO' + _LoopTypesTester._skip_types + '?'


TestLoopTypesSubtractAndNegative.autogenerate()


class TestLoopTypesReciprocal(_LoopTypesTester):
    _ufuncs = [np.reciprocal] # issue #757
    _required_types = 'bBhHiIlLqQfdFD'
    _skip_types = 'mMO' + _LoopTypesTester._skip_types

    def _arg_for_type(self, a_letter_type, index=0):
        res = super(self.__class__, self)._arg_for_type(a_letter_type,
                                                        index=index)
        if a_letter_type in 'bBhHiIlLqQ':
            # For integer reciprocal, avoid 0 as argument, as it triggers
            # undefined behavior that may differ in results from Numba
            # to the compiler used to compile NumPy.
            res[res == 0] = 42
        return res


TestLoopTypesReciprocal.autogenerate()


class TestLoopTypesPower(_LoopTypesTester):
    _ufuncs = [np.power] # issue #757
    _required_types = 'bBhHiIlLqQfdFD'
    _skip_types = 'mMO' + _LoopTypesTester._skip_types

    def _arg_for_type(self, a_letter_type, index=0):
        res = super(self.__class__, self)._arg_for_type(a_letter_type,
                                                        index=index)
        if a_letter_type in 'bBhHiIlLqQ' and index == 1:
            # For integer power, avoid a negative exponent, as it triggers
            # undefined behavior that may differ in results from Numba
            # to the compiler used to compile NumPy
            res[res < 0] = 3
        return res


TestLoopTypesPower.autogenerate()


class TestLoopTypesIntLeftShift(_LoopTypesTester):
    _ufuncs = [np.left_shift]
    _required_types = 'bBhHiIlLqQ'
    _skip_types = 'fdFDmMO' + _LoopTypesTester._skip_types

    def _arg_for_type(self, a_letter_type, index=0):
        res = super(self.__class__, self)._arg_for_type(a_letter_type,
                                                        index=index)
        # Shifting by a negative amount (argument with index 1) is undefined
        # behavior in C. It is also undefined behavior in numba. In the same
        # sense, it is also undefined behavior when the shift amount is larger
        # than the number of bits in the shifted integer.
        # To avoid problems in the test, the values are clamped (clipped) so
        # that 0 <= shift_amount < bitcount(shifted_integer)
        if index == 1:
            bit_count = res.dtype.itemsize * 8
            res = np.clip(res, 0, bit_count - 1)
        return res


TestLoopTypesIntLeftShift.autogenerate()


class TestLoopTypesIntRightShift(_LoopTypesTester):
    _ufuncs = [np.right_shift]
    _required_types = 'bBhHiIlLqQ'
    _skip_types = 'fdFDmMO' + _LoopTypesTester._skip_types

    def _arg_for_type(self, a_letter_type, index=0):
        res = super(self.__class__, self)._arg_for_type(a_letter_type,
                                                        index=index)
        # Shifting by a negative amount (argument with index 1) is undefined
        # behavior in C. It is also undefined behavior in numba. In the same
        # sense, it is also undefined behavior when the shift amount is larger
        # than the number of bits in the shifted integer.
        # To avoid problems in the test, the values are clamped (clipped) so
        # that 0 <= shift_amount < bitcount(shifted_integer)
        if index == 1:
            bit_count = res.dtype.itemsize * 8
            res = np.clip(res, 0, bit_count - 1)

        # Right shift has "implementation defined behavior" when the number
        # shifted is negative (in C). In numba, right shift for signed integers
        # is "arithmetic" while for unsigned integers is "logical".
        # This test compares against the NumPy implementation, that relies
        # on "implementation defined behavior", so the test could be a false
        # failure if the compiler used to compile NumPy doesn't follow the same
        # policy.
        # Hint: do not rely on right shifting negative numbers in NumPy.
        if index == 0:
            res = np.abs(res)
        return res


TestLoopTypesIntRightShift.autogenerate()


class TestLoopTypesFloorDivide(_LoopTypesTester):
    _ufuncs = [np.floor_divide, np.remainder, np.divmod]
    _required_types = 'bBhHiIlLqQfdFD'
    _skip_types = 'mMO' + _LoopTypesTester._skip_types

    def _fixup_results(self, dtype, py_arg, c_arg):
        if dtype.kind == 'f':
            # Discrepancies on floating-point floor division and remainder:
            # Numpy may return nan where Numba returns inf, e.g. 1. // 0.
            pred = (np.isinf(c_arg) & np.isnan(py_arg))
            # Numpy and Numba may differ in signed zeros, e.g. -0. // -1.
            pred |= (py_arg == 0.0) & (c_arg == 0.0)
            c_arg[pred] = py_arg[pred]
        return py_arg, c_arg


TestLoopTypesFloorDivide.autogenerate()


class TestLoopTypesFloat(_LoopTypesTester):
    _ufuncs = supported_ufuncs[:]
    if iswindows:
        _ufuncs.remove(np.signbit) # TODO: fix issue #758
    _ufuncs.remove(np.floor_divide) # has its own test class
    _ufuncs.remove(np.remainder) # has its own test class
    _ufuncs.remove(np.divmod) # has its own test class
    _ufuncs.remove(np.mod) # same as np.remainder
    _required_types = 'fd'
    _skip_types = 'FDmMO' + _LoopTypesTester._skip_types


TestLoopTypesFloat.autogenerate()


class TestLoopTypesComplex(_LoopTypesTester):
    _ufuncs = supported_ufuncs[:]

    # Test complex types
    # Every loop containing a complex argument must be tested
    _required_types = 'FD'
    _skip_types = 'mMO' + _LoopTypesTester._skip_types


TestLoopTypesComplex.autogenerate()


class TestLoopTypesDatetime(_LoopTypesTester):
    _ufuncs = supported_ufuncs[:]

    _ufuncs.remove(np.divmod)  # not implemented yet

    # NOTE: the full list of ufuncs supporting datetime64 and timedelta64
    # types in Numpy is:
    # ['absolute', 'add', 'divide', 'equal', 'floor_divide', 'fmax', 'fmin',
    #  'greater', 'greater_equal', 'less', 'less_equal', 'maximum',
    #  'minimum', 'multiply', 'negative', 'not_equal', 'sign', 'subtract',
    #  'true_divide']

    # Test datetime64 and timedelta64 types.
    _required_types = 'mM'

    # Test various units combinations (TestLoopTypes is only able to test
    # homogeneous units).

    def test_add(self):
        ufunc = np.add
        fn = _make_ufunc_usecase(ufunc)
        # heterogeneous inputs
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[s]', 'm8[m]', 'm8[s]'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[m]', 'm8[s]', 'm8[s]'])
        # heterogeneous inputs, scaled output
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[s]', 'm8[m]', 'm8[ms]'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[m]', 'm8[s]', 'm8[ms]'])
        # Cannot upscale result (Numpy would accept this)
        with self.assertRaises(LoweringError):
            self._check_ufunc_with_dtypes(fn, ufunc,
                                          ['m8[m]', 'm8[s]', 'm8[m]'])

    def test_subtract(self):
        ufunc = np.subtract
        fn = _make_ufunc_usecase(ufunc)
        # heterogeneous inputs
        self._check_ufunc_with_dtypes(fn, ufunc, ['M8[s]', 'M8[m]', 'm8[s]'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['M8[m]', 'M8[s]', 'm8[s]'])
        # heterogeneous inputs, scaled output
        self._check_ufunc_with_dtypes(fn, ufunc, ['M8[s]', 'M8[m]', 'm8[ms]'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['M8[m]', 'M8[s]', 'm8[ms]'])
        # Cannot upscale result (Numpy would accept this)
        with self.assertRaises(LoweringError):
            self._check_ufunc_with_dtypes(fn, ufunc,
                                          ['M8[m]', 'M8[s]', 'm8[m]'])

    def test_multiply(self):
        ufunc = np.multiply
        fn = _make_ufunc_usecase(ufunc)
        # scaled output
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[s]', 'q', 'm8[us]'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['q', 'm8[s]', 'm8[us]'])
        # Cannot upscale result (Numpy would accept this)
        with self.assertRaises(LoweringError):
            self._check_ufunc_with_dtypes(fn, ufunc, ['m8[s]', 'q', 'm8[m]'])

    def test_true_divide(self):
        ufunc = np.true_divide
        fn = _make_ufunc_usecase(ufunc)
        # heterogeneous inputs
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[m]', 'm8[s]', 'd'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[s]', 'm8[m]', 'd'])
        # scaled output
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[m]', 'q', 'm8[s]'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[m]', 'd', 'm8[s]'])
        # Cannot upscale result (Numpy would accept this)
        with self.assertRaises(LoweringError):
            self._check_ufunc_with_dtypes(fn, ufunc, ['m8[s]', 'q', 'm8[m]'])

    def test_floor_divide(self):
        ufunc = np.floor_divide
        fn = _make_ufunc_usecase(ufunc)
        # scaled output
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[m]', 'q', 'm8[s]'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[m]', 'd', 'm8[s]'])
        # Cannot upscale result (Numpy would accept this)
        with self.assertRaises(LoweringError):
            self._check_ufunc_with_dtypes(fn, ufunc, ['m8[s]', 'q', 'm8[m]'])

    def _check_comparison(self, ufunc):
        fn = _make_ufunc_usecase(ufunc)
        # timedelta
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[m]', 'm8[s]', '?'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['m8[s]', 'm8[m]', '?'])
        # datetime
        self._check_ufunc_with_dtypes(fn, ufunc, ['M8[m]', 'M8[s]', '?'])
        self._check_ufunc_with_dtypes(fn, ufunc, ['M8[s]', 'M8[m]', '?'])

    def test_comparisons(self):
        for ufunc in [np.equal, np.not_equal, np.less, np.less_equal,
                      np.greater, np.greater_equal]:
            self._check_comparison(ufunc)


TestLoopTypesDatetime.autogenerate()


class TestUFuncBadArgs(TestCase):
    def test_missing_args(self):
        def func(x):
            """error: np.add requires two args"""
            result = np.add(x)
            return result

        with self.assertRaises(TypingError):
            njit([types.float64(types.float64)])(func)

    def test_too_many_args(self):
        def func(x, out, out2):
            """error: too many args"""
            result = np.add(x, x, out, out2)
            return result

        array_type = types.Array(types.float64, 1, 'C')
        sig = array_type(array_type, array_type, array_type)

        with self.assertRaises(TypingError):
            njit(sig)(func)

    def test_no_scalar_result_by_reference(self):
        def func(x):
            """error: scalar as a return value is not supported"""
            y = 0
            np.add(x, x, y)

        with self.assertRaises(TypingError):
            njit([types.float64(types.float64)])(func)


class TestUFuncCompilationThreadSafety(TestCase):

    def test_lock(self):
        """
        Test that (lazy) compiling from several threads at once doesn't
        produce errors (see issue #2403).
        """
        errors = []

        @vectorize
        def foo(x):
            return x + 1

        def wrapper():
            try:
                a = np.ones((10,), dtype=np.float64)
                expected = np.ones((10,), dtype=np.float64) + 1.
                np.testing.assert_array_equal(foo(a), expected)
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=wrapper) for i in range(16)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertFalse(errors)


class TestUfuncOnContext(TestCase):
    def test_cpu_get_ufunc_info(self):
        # The CPU context defines get_ufunc_info that is the same as
        # ufunc_db.get_ufunc_info.
        targetctx = cpu_target.target_context
        # Check: get_ufunc_info returns a dict
        add_info = targetctx.get_ufunc_info(np.add)
        self.assertIsInstance(add_info, dict)
        # Check: it is the same as ufunc_db.get_ufunc_info
        expected = ufunc_db.get_ufunc_info(np.add)
        self.assertEqual(add_info, expected)
        # Check: KeyError raised on bad key
        badkey = object()
        with self.assertRaises(KeyError) as raises:
            ufunc_db.get_ufunc_info(badkey)
        self.assertEqual(raises.exception.args, (badkey,))

    def test_base_get_ufunc_info(self):
        # The BaseContext always raises NotImplementedError
        targetctx = BaseContext(cpu_target.typing_context, 'cpu')
        with self.assertRaises(NotImplementedError) as raises:
            targetctx.get_ufunc_info(np.add)
        self.assertRegex(
            str(raises.exception),
            r"<numba\..*\.BaseContext object at .*> does not support ufunc",
        )


class TestUfuncWriteInput(TestCase):
    def test_write_input_arg(self):
        @guvectorize(["void(float64[:], uint8[:])"], "(n)->(n)")
        def func(x, out):

            for i in range(x.size):
                # set every fourth element to 1
                if i % 4 == 0:
                    out[i] = 1

        x = np.random.rand(10, 5)
        out = np.zeros_like(x, dtype=np.int8)

        func(x, out)
        np.testing.assert_array_equal(
            np.array([True, False, False, False, True], dtype=np.bool_),
            out.any(axis=0))


if __name__ == '__main__':
    unittest.main()
