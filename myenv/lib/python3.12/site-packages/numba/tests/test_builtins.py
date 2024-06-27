import itertools
import functools
import sys
import operator
from collections import namedtuple

import numpy as np

import unittest
import warnings

from numba import jit, typeof, njit, typed
from numba.core import errors, types, config
from numba.tests.support import (TestCase, tag, ignore_internal_warnings,
                                 MemoryLeakMixin)
from numba.core.extending import overload_method, box


forceobj_flags = {'forceobj': True}

no_pyobj_flags = {'nopython': True, '_nrt': False}

nrt_no_pyobj_flags = {'nopython': True}


def abs_usecase(x):
    return abs(x)

def all_usecase(x, y):
    if x == None and y == None:
        return all([])
    elif x == None:
        return all([y])
    elif y == None:
        return all([x])
    else:
        return all([x, y])

def any_usecase(x, y):
    if x == None and y == None:
        return any([])
    elif x == None:
        return any([y])
    elif y == None:
        return any([x])
    else:
        return any([x, y])

def bool_usecase(x):
    return bool(x)

def complex_usecase(x, y):
    return complex(x, y)

def divmod_usecase(x, y):
    return divmod(x, y)

def enumerate_usecase():
    result = 0
    for i, j in enumerate((1., 2.5, 3.)):
        result += i * j
    return result

def enumerate_start_usecase():
    result = 0
    for i, j in enumerate((1., 2.5, 3.), 42):
        result += i * j
    return result

def enumerate_invalid_start_usecase():
    result = 0
    for i, j in enumerate((1., 2.5, 3.), 3.14159):
        result += i * j
    return result

def filter_usecase(x, filter_func):
    return filter(filter_func, x)

def float_usecase(x):
    return float(x)

def float_inf_usecase(x):
    d = {
        0: float('inf'),
        1: float('INF'),
        2: float('-inf'),
        3: float('-INF'),
        4: float('\r\nINF\r       '),
        5: float('       \r\n\t-INF'),
        6: float('1234.45'),
        7: float('\n-123.4\r'),
    }
    return d.get(x)

def format_usecase(x, y):
    return x.format(y)

def globals_usecase():
    return globals()

# NOTE: hash() is tested in test_hashing

def hex_usecase(x):
    return hex(x)

def str_usecase(x):
    return str(x)

def int_usecase(x, base):
    return int(x, base=base)

def iter_next_usecase(x):
    it = iter(x)
    return next(it), next(it)

def locals_usecase(x):
    y = 5
    return locals()['y']

def long_usecase(x, base):
    return long(x, base=base)

def map_usecase(x, map_func):
    return map(map_func, x)


def max_usecase1(x, y):
    return max(x, y)

def max_usecase2(x, y):
    return max([x, y])

def max_usecase3(x):
    return max(x)

def max_usecase4():
    return max(())


def min_usecase1(x, y):
    return min(x, y)

def min_usecase2(x, y):
    return min([x, y])

def min_usecase3(x):
    return min(x)

def min_usecase4():
    return min(())

def oct_usecase(x):
    return oct(x)

def reduce_usecase(reduce_func, x):
    return functools.reduce(reduce_func, x)

def round_usecase1(x):
    return round(x)

def round_usecase2(x, n):
    return round(x, n)

def sum_usecase(x):
    return sum(x)

def type_unary_usecase(a, b):
    return type(a)(b)

def truth_usecase(p):
    return operator.truth(p)

def unichr_usecase(x):
    return unichr(x)

def zip_usecase():
    result = 0
    for i, j in zip((1, 2, 3), (4.5, 6.7)):
        result += i * j
    return result

def zip_0_usecase():
    result = 0
    for i in zip():
        result += 1
    return result

def zip_1_usecase():
    result = 0
    for i, in zip((1, 2)):
        result += i
    return result


def zip_3_usecase():
    result = 0
    for i, j, k in zip((1, 2), (3, 4, 5), (6.7, 8.9)):
        result += i * j * k
    return result


def zip_first_exhausted():
    iterable = range(7)
    n = 3
    it = iter(iterable)
    # 1st iterator is shorter
    front = list(zip(range(n), it))
    # Make sure that we didn't skip one in `it`
    back = list(it)
    return front, back


def pow_op_usecase(x, y):
    return x ** y


def pow_usecase(x, y):
    return pow(x, y)


def sum_usecase(x):
    return sum(x)


def sum_kwarg_usecase(x, start=0):
    ret = sum(x, start)
    return sum(x, start=start), ret


def isinstance_usecase(a):
    if isinstance(a, (int, float)):
        if isinstance(a, int):
            return a + 1, 'int'
        if isinstance(a, float):
            return a + 2.0, 'float'
    elif isinstance(a, str):
        return a + ", world!", 'str'
    elif isinstance(a, complex):
        return a.imag, 'complex'
    elif isinstance(a, (tuple, list)):
        if isinstance(a, tuple):
            return 'tuple'
        else:
            return 'list'
    elif isinstance(a, set):
        return 'set'
    elif isinstance(a, bytes):
        return 'bytes'
    return 'no match'


def isinstance_dict():
    a = {1: 2, 3: 4}
    b = {'a': 10, 'b': np.zeros(3)}
    if isinstance(a, dict) and isinstance(b, dict):
        return 'dict'
    else:
        return 'not dict'


def isinstance_usecase_numba_types(a):
    if isinstance(a, typed.List):
        return 'typed list'
    elif isinstance(a, (types.int32, types.int64)):
        if isinstance(a, types.int32):
            return 'int32'
        else:
            return 'int64'
    elif isinstance(a, (types.float32, types.float64)):
        if isinstance(a, types.float32):
            return 'float32'
        elif isinstance(a, types.float64):
            return 'float64'
    elif isinstance(a, typed.Dict):
        return 'typed dict'
    else:
        return 'no match'


def isinstance_usecase_numba_types_2():
    # some types cannot be passed as argument to njit functions
    a = b'hello'
    b = range(1, 2)
    c = dict()
    c[2] = 3
    if isinstance(a, bytes) and \
            isinstance(b, range) and \
            isinstance(c, dict):
        return True
    return False


def invalid_isinstance_usecase(x):
    if isinstance(x, ('foo',)):
        return 'true branch'
    else:
        return 'false branch'


def isinstance_usecase_invalid_type(x):
    # this should be a valid call when x := float
    if isinstance(x, (float, 'not a type')):
        return True
    else:
        return False


def invalid_isinstance_usecase_phi_nopropagate(x):
    if x > 4:
        z = 10
    else:
        z = 'a'
    if isinstance(z, int):
        return True
    else:
        return False


def invalid_isinstance_usecase_phi_nopropagate2(a):
    # Numba issue #9125
    x = 0
    if isinstance(a, int):
        a = (a, a)

    for i in range(len(a)):
        x += i
    return x


def invalid_isinstance_optional_usecase(x):
    if x > 4:
        z = 10
    else:
        z = None
    if isinstance(z, int):
        return True
    else:
        return False

def invalid_isinstance_unsupported_type_usecase():
    ntpl = namedtuple('ntpl', ['a', 'b'])
    inst = ntpl(1, 2)
    def impl(x):
        return isinstance(inst, ntpl)
    return impl

class TestBuiltins(TestCase):

    def run_nullary_func(self, pyfunc, flags):
        cfunc = jit((), **flags)(pyfunc)
        expected = pyfunc()
        self.assertPreciseEqual(cfunc(), expected)

    def test_abs(self, flags=forceobj_flags):
        pyfunc = abs_usecase

        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in [-1, 0, 1]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

        cfunc = jit((types.float32,), **flags)(pyfunc)
        for x in [-1.1, 0.0, 1.1]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x), prec='single')

        complex_values = [-1.1 + 0.5j, 0.0 + 0j, 1.1 + 3j,
                          float('inf') + 1j * float('nan'),
                          float('nan') - 1j * float('inf')]
        cfunc = jit((types.complex64,), **flags)(pyfunc)
        for x in complex_values:
            self.assertPreciseEqual(cfunc(x), pyfunc(x), prec='single')
        cfunc = jit((types.complex128,), **flags)(pyfunc)
        for x in complex_values:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

        for unsigned_type in types.unsigned_domain:
            unsigned_values = [0, 10, 2, 2 ** unsigned_type.bitwidth - 1]
            cfunc = jit((unsigned_type,), **flags)(pyfunc)
            for x in unsigned_values:
                self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_abs_npm(self):
        self.test_abs(flags=no_pyobj_flags)

    def test_all(self, flags=forceobj_flags):
        pyfunc = all_usecase

        cfunc = jit((types.int32,types.int32), **flags)(pyfunc)
        x_operands = [-1, 0, 1, None]
        y_operands = [-1, 0, 1, None]
        for x, y in itertools.product(x_operands, y_operands):
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

    def test_all_npm(self):
        with self.assertTypingError():
            self.test_all(flags=no_pyobj_flags)

    def test_any(self, flags=forceobj_flags):
        pyfunc = any_usecase

        cfunc = jit((types.int32,types.int32), **flags)(pyfunc)
        x_operands = [-1, 0, 1, None]
        y_operands = [-1, 0, 1, None]
        for x, y in itertools.product(x_operands, y_operands):
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

    def test_any_npm(self):
        with self.assertTypingError():
            self.test_any(flags=no_pyobj_flags)

    def test_bool(self, flags=forceobj_flags):
        pyfunc = bool_usecase

        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in [-1, 0, 1]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))
        cfunc = jit((types.float64,), **flags)(pyfunc)
        for x in [0.0, -0.0, 1.5, float('inf'), float('nan')]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))
        cfunc = jit((types.complex128,), **flags)(pyfunc)
        for x in [complex(0, float('inf')), complex(0, float('nan'))]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_bool_npm(self):
        self.test_bool(flags=no_pyobj_flags)

    def test_bool_nonnumber(self, flags=forceobj_flags):
        pyfunc = bool_usecase

        cfunc = jit((types.string,), **flags)(pyfunc)
        for x in ['x', '']:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

        cfunc = jit((types.Dummy('list'),), **flags)(pyfunc)
        for x in [[1], []]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_bool_nonnumber_npm(self):
        with self.assertTypingError():
            self.test_bool_nonnumber(flags=no_pyobj_flags)

    def test_complex(self, flags=forceobj_flags):
        pyfunc = complex_usecase

        cfunc = jit((types.int32, types.int32), **flags)(pyfunc)

        x_operands = [-1, 0, 1]
        y_operands = [-1, 0, 1]
        for x, y in itertools.product(x_operands, y_operands):
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

    def test_complex_npm(self):
        self.test_complex(flags=no_pyobj_flags)

    def test_divmod_ints(self, flags=forceobj_flags):
        pyfunc = divmod_usecase

        cfunc = jit((types.int64, types.int64), **flags)(pyfunc)

        def truncate_result(x, bits=64):
            # Remove any extraneous bits (since Numba will return
            # a 64-bit result by definition)
            if x >= 0:
                x &= (1 << (bits - 1)) - 1
            return x

        denominators = [1, 3, 7, 15, -1, -3, -7, -15, 2**63 - 1, -2**63]
        numerators = denominators + [0]
        for x, y, in itertools.product(numerators, denominators):
            expected_quot, expected_rem = pyfunc(x, y)
            quot, rem = cfunc(x, y)
            f = truncate_result
            self.assertPreciseEqual((f(quot), f(rem)),
                                    (f(expected_quot), f(expected_rem)))

        for x in numerators:
            with self.assertRaises(ZeroDivisionError):
                cfunc(x, 0)

    def test_divmod_ints_npm(self):
        self.test_divmod_ints(flags=no_pyobj_flags)

    def test_divmod_floats(self, flags=forceobj_flags):
        pyfunc = divmod_usecase

        cfunc = jit((types.float64, types.float64), **flags)(pyfunc)

        denominators = [1., 3.5, 1e100, -2., -7.5, -1e101,
                        np.inf, -np.inf, np.nan]
        numerators = denominators + [-0.0, 0.0]
        for x, y, in itertools.product(numerators, denominators):
            expected_quot, expected_rem = pyfunc(x, y)
            quot, rem = cfunc(x, y)
            self.assertPreciseEqual((quot, rem), (expected_quot, expected_rem))

        for x in numerators:
            with self.assertRaises(ZeroDivisionError):
                cfunc(x, 0.0)

    def test_divmod_floats_npm(self):
        self.test_divmod_floats(flags=no_pyobj_flags)

    def test_enumerate(self, flags=forceobj_flags):
        self.run_nullary_func(enumerate_usecase, flags)

    def test_enumerate_npm(self):
        self.test_enumerate(flags=no_pyobj_flags)

    def test_enumerate_start(self, flags=forceobj_flags):
        self.run_nullary_func(enumerate_start_usecase, flags)

    def test_enumerate_start_npm(self):
        self.test_enumerate_start(flags=no_pyobj_flags)

    def test_enumerate_start_invalid_start_type(self):
        pyfunc = enumerate_invalid_start_usecase

        cfunc = jit((), **forceobj_flags)(pyfunc)

        with self.assertRaises(TypeError) as raises:
            cfunc()

        msg = "'float' object cannot be interpreted as an integer"
        self.assertIn(msg, str(raises.exception))

    def test_enumerate_start_invalid_start_type_npm(self):
        pyfunc = enumerate_invalid_start_usecase
        with self.assertRaises(errors.TypingError) as raises:
            jit((), **no_pyobj_flags)(pyfunc)
        msg = "Only integers supported as start value in enumerate"
        self.assertIn(msg, str(raises.exception))

    def test_filter(self, flags=forceobj_flags):
        pyfunc = filter_usecase
        argtys = (types.Dummy('list'), types.Dummy('function_ptr'))
        cfunc = jit(argtys, **flags)(pyfunc)

        filter_func = lambda x: x % 2
        x = [0, 1, 2, 3, 4]
        self.assertSequenceEqual(list(cfunc(x, filter_func)),
                                 list(pyfunc(x, filter_func)))

    def test_filter_npm(self):
        with self.assertTypingError():
            self.test_filter(flags=no_pyobj_flags)

    def test_float(self, flags=forceobj_flags):
        pyfunc = float_usecase

        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in [-1, 0, 1]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

        cfunc = jit((types.float32,), **flags)(pyfunc)
        for x in [-1.1, 0.0, 1.1]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x), prec='single')

        cfunc = jit((types.string,), **flags)(pyfunc)
        for x in ['-1.1', '0.0', '1.1', 'inf', '-inf', 'INF', '-INF']:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_float_npm(self):
        with self.assertTypingError():
            self.test_float(flags=no_pyobj_flags)

    def test_float_string_literal(self):
        pyfunc = float_inf_usecase
        cfunc = njit(pyfunc)
        for x in range(8):
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_format(self, flags=forceobj_flags):
        pyfunc = format_usecase

        cfunc = jit((types.string, types.int32,), **flags)(pyfunc)
        x = '{0}'
        for y in [-1, 0, 1]:
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

        cfunc = jit((types.string, types.float32,), **flags)(pyfunc)
        x = '{0}'
        for y in [-1.1, 0.0, 1.1]:
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

        cfunc = jit((types.string, types.string,), **flags)(pyfunc)
        x = '{0}'
        for y in ['a', 'b', 'c']:
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

    def test_format_npm(self):
        with self.assertTypingError():
            self.test_format(flags=no_pyobj_flags)

    def test_globals(self, flags=forceobj_flags):
        pyfunc = globals_usecase
        cfunc = jit((), **flags)(pyfunc)
        g = cfunc()
        self.assertIs(g, globals())

    def test_globals_npm(self):
        with self.assertTypingError():
            self.test_globals(flags=no_pyobj_flags)

    def test_globals_jit(self, flags=forceobj_flags):
        # Issue #416: weird behaviour of globals() in combination with
        # the @jit decorator.
        pyfunc = globals_usecase
        jitted = jit(**flags)(pyfunc)
        self.assertIs(jitted(), globals())
        self.assertIs(jitted(), globals())

    def test_globals_jit_npm(self):
        with self.assertTypingError():
            self.test_globals_jit(nopython=True)

    def test_hex(self, flags=forceobj_flags):
        pyfunc = hex_usecase

        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in [-1, 0, 1]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_hex_npm(self):
        with self.assertTypingError():
            self.test_hex(flags=no_pyobj_flags)

    def test_int_str(self):
        pyfunc = str_usecase

        small_inputs = [
            1234,
            1,
            0,
            10,
            1000,
        ]

        large_inputs = [
            123456789,
            2222222,
            1000000,
            ~0x0
        ]

        args = [*small_inputs, *large_inputs]

        typs = [
            types.int8,
            types.int16,
            types.int32,
            types.int64,
            types.uint,
            types.uint8,
            types.uint16,
            types.uint32,
            types.uint64,
        ]

        for typ in typs:
            cfunc = jit((typ,), **nrt_no_pyobj_flags)(pyfunc)
            for v in args:
                tp_info = np.iinfo(typ.key)
                if not (tp_info.min <= v <= tp_info.max):
                    continue
                self.assertPreciseEqual(cfunc(typ(v)), pyfunc(typ(v)))

                if typ.signed:
                    self.assertPreciseEqual(cfunc(typ(-v)), pyfunc(typ(-v)))

    def test_int(self, flags=forceobj_flags):
        pyfunc = int_usecase

        cfunc = jit((types.string, types.int32,), **flags)(pyfunc)

        x_operands = ['-1', '0', '1', '10']
        y_operands = [2, 8, 10, 16]
        for x, y in itertools.product(x_operands, y_operands):
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

    def test_int_npm(self):
        with self.assertTypingError():
            self.test_int(flags=no_pyobj_flags)

    def test_iter_next(self, flags=forceobj_flags):
        pyfunc = iter_next_usecase
        cfunc = jit((types.UniTuple(types.int32, 3),), **flags)(pyfunc)
        self.assertPreciseEqual(cfunc((1, 42, 5)), (1, 42))

        cfunc = jit((types.UniTuple(types.int32, 1),), **flags)(pyfunc)
        with self.assertRaises(StopIteration):
            cfunc((1,))

    def test_iter_next_npm(self):
        self.test_iter_next(flags=no_pyobj_flags)

    def test_locals(self, flags=forceobj_flags):
        pyfunc = locals_usecase
        with self.assertRaises(errors.ForbiddenConstruct):
            jit((types.int64,), **flags)(pyfunc)

    def test_locals_forceobj(self):
        self.test_locals(flags=forceobj_flags)

    def test_locals_npm(self):
        with self.assertTypingError():
            self.test_locals(flags=no_pyobj_flags)

    def test_map(self, flags=forceobj_flags):
        pyfunc = map_usecase
        argtys = (types.Dummy('list'), types.Dummy('function_ptr'))
        cfunc = jit(argtys, **flags)(pyfunc)

        map_func = lambda x: x * 2
        x = [0, 1, 2, 3, 4]
        self.assertSequenceEqual(list(cfunc(x, map_func)),
                                 list(pyfunc(x, map_func)))

    def test_map_npm(self):
        with self.assertTypingError():
            self.test_map(flags=no_pyobj_flags)

    #
    # min() and max()
    #

    def check_minmax_1(self, pyfunc, flags):
        cfunc = jit((types.int32, types.int32), **flags)(pyfunc)

        x_operands = [-1, 0, 1]
        y_operands = [-1, 0, 1]
        for x, y in itertools.product(x_operands, y_operands):
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

    def test_max_1(self, flags=forceobj_flags):
        """
        max(*args)
        """
        self.check_minmax_1(max_usecase1, flags)

    def test_min_1(self, flags=forceobj_flags):
        """
        min(*args)
        """
        self.check_minmax_1(min_usecase1, flags)

    def test_max_npm_1(self):
        self.test_max_1(flags=no_pyobj_flags)

    def test_min_npm_1(self):
        self.test_min_1(flags=no_pyobj_flags)

    def check_minmax_2(self, pyfunc, flags):
        cfunc = jit((types.int32, types.int32), **flags)(pyfunc)

        x_operands = [-1, 0, 1]
        y_operands = [-1, 0, 1]
        for x, y in itertools.product(x_operands, y_operands):
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

    def test_max_2(self, flags=forceobj_flags):
        """
        max(list)
        """
        self.check_minmax_2(max_usecase2, flags)

    def test_min_2(self, flags=forceobj_flags):
        """
        min(list)
        """
        self.check_minmax_2(min_usecase2, flags)

    def test_max_npm_2(self):
        with self.assertTypingError():
            self.test_max_2(flags=no_pyobj_flags)

    def test_min_npm_2(self):
        with self.assertTypingError():
            self.test_min_2(flags=no_pyobj_flags)

    def check_minmax_3(self, pyfunc, flags):
        def check(argty):
            cfunc = jit((argty,), **flags)(pyfunc)
            # Check that the algorithm matches Python's with a non-total order
            tup = (1.5, float('nan'), 2.5)
            for val in [tup, tup[::-1]]:
                self.assertPreciseEqual(cfunc(val), pyfunc(val))

        check(types.UniTuple(types.float64, 3))
        check(types.Tuple((types.float32, types.float64, types.float32)))

    def test_max_3(self, flags=forceobj_flags):
        """
        max(tuple)
        """
        self.check_minmax_3(max_usecase3, flags)

    def test_min_3(self, flags=forceobj_flags):
        """
        min(tuple)
        """
        self.check_minmax_3(min_usecase3, flags)

    def test_max_npm_3(self):
        self.test_max_3(flags=no_pyobj_flags)

    def test_min_npm_3(self):
        self.test_min_3(flags=no_pyobj_flags)

    def check_min_max_invalid_types(self, pyfunc, flags=forceobj_flags):
        cfunc = jit((types.int32, types.Dummy('list'),), **flags)(pyfunc)
        cfunc(1, [1])

    def test_max_1_invalid_types(self):
        with self.assertRaises(TypeError):
            self.check_min_max_invalid_types(max_usecase1)

    def test_max_1_invalid_types_npm(self):
        with self.assertTypingError():
            self.check_min_max_invalid_types(max_usecase1, flags=no_pyobj_flags)

    def test_min_1_invalid_types(self):
        with self.assertRaises(TypeError):
            self.check_min_max_invalid_types(min_usecase1)

    def test_min_1_invalid_types_npm(self):
        with self.assertTypingError():
            self.check_min_max_invalid_types(min_usecase1, flags=no_pyobj_flags)

    def check_minmax_bool1(self, pyfunc, flags):
        cfunc = jit((types.bool_, types.bool_), **flags)(pyfunc)

        operands = (False, True)
        for x, y in itertools.product(operands, operands):
            self.assertPreciseEqual(cfunc(x, y), pyfunc(x, y))

    def test_max_bool1(self, flags=forceobj_flags):
        # tests max(<booleans>)
        self.check_minmax_bool1(max_usecase1, flags)

    def test_min_bool1(self, flags=forceobj_flags):
        # tests min(<booleans>)
        self.check_minmax_bool1(min_usecase1, flags)

    # Test that max(1) and min(1) fail

    def check_min_max_unary_non_iterable(self, pyfunc, flags=forceobj_flags):
        cfunc = jit((types.int32,), **flags)(pyfunc)
        cfunc(1)

    def test_max_unary_non_iterable(self):
        with self.assertRaises(TypeError):
            self.check_min_max_unary_non_iterable(max_usecase3)

    def test_max_unary_non_iterable_npm(self):
        with self.assertTypingError():
            self.check_min_max_unary_non_iterable(max_usecase3)

    def test_min_unary_non_iterable(self):
        with self.assertRaises(TypeError):
            self.check_min_max_unary_non_iterable(min_usecase3)

    def test_min_unary_non_iterable_npm(self):
        with self.assertTypingError():
            self.check_min_max_unary_non_iterable(min_usecase3)

    # Test that max(()) and min(()) fail

    def check_min_max_empty_tuple(self, pyfunc, func_name):
        with self.assertTypingError() as raises:
            jit((), **no_pyobj_flags)(pyfunc)
        self.assertIn("%s() argument is an empty tuple" % func_name,
                      str(raises.exception))

    def test_max_empty_tuple(self):
        self.check_min_max_empty_tuple(max_usecase4, "max")

    def test_min_empty_tuple(self):
        self.check_min_max_empty_tuple(min_usecase4, "min")


    def test_oct(self, flags=forceobj_flags):
        pyfunc = oct_usecase

        cfunc = jit((types.int32,), **flags)(pyfunc)
        for x in [-8, -1, 0, 1, 8]:
            self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_oct_npm(self):
        with self.assertTypingError():
            self.test_oct(flags=no_pyobj_flags)

    def test_reduce(self, flags=forceobj_flags):
        pyfunc = reduce_usecase
        argtys = (types.Dummy('function_ptr'), types.Dummy('list'))
        cfunc = jit(argtys, **flags)(pyfunc)

        reduce_func = lambda x, y: x + y

        x = range(10)
        self.assertPreciseEqual(cfunc(reduce_func, x), pyfunc(reduce_func, x))

        x = [x + x/10.0 for x in range(10)]
        self.assertPreciseEqual(cfunc(reduce_func, x), pyfunc(reduce_func, x))

        x = [complex(x, x) for x in range(10)]
        self.assertPreciseEqual(cfunc(reduce_func, x), pyfunc(reduce_func, x))

    def test_reduce_npm(self):
        with self.assertTypingError():
            self.test_reduce(flags=no_pyobj_flags)

    def test_round1(self, flags=forceobj_flags):
        pyfunc = round_usecase1

        for tp in (types.float64, types.float32):
            cfunc = jit((tp,), **flags)(pyfunc)
            values = [-1.6, -1.5, -1.4, -0.5, 0.0, 0.1, 0.5, 0.6, 1.4, 1.5, 5.0]
            values += [-0.1, -0.0]
            for x in values:
                self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_round1_npm(self):
        self.test_round1(flags=no_pyobj_flags)

    def test_round2(self, flags=forceobj_flags):
        pyfunc = round_usecase2

        for tp in (types.float64, types.float32):
            prec = 'single' if tp is types.float32 else 'exact'
            cfunc = jit((tp, types.int32), **flags)(pyfunc)
            for x in [0.0, 0.1, 0.125, 0.25, 0.5, 0.75, 1.25,
                      1.5, 1.75, 2.25, 2.5, 2.75, 12.5, 15.0, 22.5]:
                for n in (-1, 0, 1, 2):
                    self.assertPreciseEqual(cfunc(x, n), pyfunc(x, n),
                                            prec=prec)
                    expected = pyfunc(-x, n)
                    self.assertPreciseEqual(cfunc(-x, n), pyfunc(-x, n),
                                            prec=prec)

    def test_round2_npm(self):
        self.test_round2(flags=no_pyobj_flags)

    def test_sum_objmode(self, flags=forceobj_flags):
        pyfunc = sum_usecase

        cfunc = jit((types.Dummy('list'),), **flags)(pyfunc)

        x = range(10)
        self.assertPreciseEqual(cfunc(x), pyfunc(x))

        x = [x + x/10.0 for x in range(10)]
        self.assertPreciseEqual(cfunc(x), pyfunc(x))

        x = [complex(x, x) for x in range(10)]
        self.assertPreciseEqual(cfunc(x), pyfunc(x))

    def test_sum(self):
        # In Python 3.8+ "start" can be specified as a kwarg, so test that too
        sum_default = njit(sum_usecase)
        sum_kwarg = njit(sum_kwarg_usecase)

        @njit
        def sum_range(sz, start=0):
            tmp = range(sz)
            ret = sum(tmp, start)
            return sum(tmp, start=start), ret

        ntpl = namedtuple('ntpl', ['a', 'b'])

        # check call with default kwarg, start=0
        def args():
            yield [*range(10)]
            yield [x + x/10.0 for x in range(10)]
            yield [x * 1j for x in range(10)]
            yield (1, 2, 3)
            yield (1, 2, 3j)
            # uints will likely end up as floats as `start` is signed, so just
            # test mixed signed ints
            yield (np.int64(32), np.int32(2), np.int8(3))
            tl = typed.List(range(5))
            yield tl
            yield np.ones(5)
            yield ntpl(100, 200)
            yield ntpl(100, 200j)

        for x in args():
            self.assertPreciseEqual(sum_default(x), sum_default.py_func(x))

        # Check the uint use case, as start is signed, NumPy will end up with
        # a float result whereas Numba will end up with an int (see integer
        # typing NBEP).
        x = (np.uint64(32), np.uint32(2), np.uint8(3))
        self.assertEqual(sum_default(x), sum_default.py_func(x))

        # check call with changing default kwarg, start
        def args_kws():
            yield [*range(10)], 12
            yield [x + x/10.0 for x in range(10)], 19j
            yield [x * 1j for x in range(10)], -2
            yield (1, 2, 3), 9
            yield (1, 2, 3j), -0
            # uints will likely end up as floats as `start` is signed, so just
            # test mixed signed ints
            yield (np.int64(32), np.int32(2), np.int8(3)), np.uint32(7)
            tl = typed.List(range(5))
            yield tl, 100
            yield np.ones((5, 5)), 10 * np.ones((5,))
            yield ntpl(100, 200), -50
            yield ntpl(100, 200j), 9

        for x, start in args_kws():
            self.assertPreciseEqual(sum_kwarg(x, start=start),
                                    sum_kwarg.py_func(x, start=start))

        # check call with range()
        for start in range(-3, 4):
            for sz in range(-3, 4):
                self.assertPreciseEqual(sum_range(sz, start=start),
                                        sum_range.py_func(sz, start=start))

    def test_sum_exceptions(self):
        sum_default = njit(sum_usecase)
        sum_kwarg = njit(sum_kwarg_usecase)

        # check start as string/bytes/bytearray is error
        msg = "sum() can't sum {}"

        with self.assertRaises(errors.TypingError) as raises:
            sum_kwarg((1, 2, 3), 'a')

        self.assertIn(msg.format('strings'), str(raises.exception))

        with self.assertRaises(errors.TypingError) as raises:
            sum_kwarg((1, 2, 3), b'123')

        self.assertIn(msg.format('bytes'), str(raises.exception))

        with self.assertRaises(errors.TypingError) as raises:
            sum_kwarg((1, 2, 3), bytearray(b'123'))

        self.assertIn(msg.format('bytearray'), str(raises.exception))

        # check invalid type has no impl
        with self.assertRaises(errors.TypingError) as raises:
            sum_default('abcd')

        self.assertIn('No implementation', str(raises.exception))

    def test_truth(self):
        pyfunc = truth_usecase
        cfunc = jit(nopython=True)(pyfunc)

        self.assertEqual(pyfunc(True), cfunc(True))
        self.assertEqual(pyfunc(False), cfunc(False))

    def test_type_unary(self):
        # Test type(val) and type(val)(other_val)
        pyfunc = type_unary_usecase
        cfunc = jit(nopython=True)(pyfunc)

        def check(*args):
            expected = pyfunc(*args)
            self.assertPreciseEqual(cfunc(*args), expected)

        check(1.5, 2)
        check(1, 2.5)
        check(1.5j, 2)
        check(True, 2)
        check(2.5j, False)

    def test_zip(self, flags=forceobj_flags):
        self.run_nullary_func(zip_usecase, flags)

    def test_zip_npm(self):
        self.test_zip(flags=no_pyobj_flags)

    def test_zip_1(self, flags=forceobj_flags):
        self.run_nullary_func(zip_1_usecase, flags)

    def test_zip_1_npm(self):
        self.test_zip_1(flags=no_pyobj_flags)

    def test_zip_3(self, flags=forceobj_flags):
        self.run_nullary_func(zip_3_usecase, flags)

    def test_zip_3_npm(self):
        self.test_zip_3(flags=no_pyobj_flags)

    def test_zip_0(self, flags=forceobj_flags):
        self.run_nullary_func(zip_0_usecase, flags)

    def test_zip_0_npm(self):
        self.test_zip_0(flags=no_pyobj_flags)

    def test_zip_first_exhausted(self, flags=forceobj_flags):
        """
        Test side effect to the input iterators when a left iterator has been
        exhausted before the ones on the right.
        """
        self.run_nullary_func(zip_first_exhausted, flags)

    def test_zip_first_exhausted_npm(self):
        self.test_zip_first_exhausted(flags=nrt_no_pyobj_flags)

    def test_pow_op_usecase(self):
        args = [
            (2, 3),
            (2.0, 3),
            (2, 3.0),
            (2j, 3.0j),
        ]

        for x, y in args:
            argtys = (typeof(x), typeof(y))
            cfunc = jit(argtys, **no_pyobj_flags)(pow_op_usecase)
            r = cfunc(x, y)
            self.assertPreciseEqual(r, pow_op_usecase(x, y))

    def test_pow_usecase(self):
        args = [
            (2, 3),
            (2.0, 3),
            (2, 3.0),
            (2j, 3.0j),
        ]

        for x, y in args:
            argtys = (typeof(x), typeof(y))
            cfunc = jit(argtys, **no_pyobj_flags)(pow_usecase)
            r = cfunc(x, y)
            self.assertPreciseEqual(r, pow_usecase(x, y))

    def _check_min_max(self, pyfunc):
        cfunc = njit()(pyfunc)
        expected = pyfunc()
        got = cfunc()
        self.assertPreciseEqual(expected, got)

    def test_min_max_iterable_input(self):

        @njit
        def frange(start, stop, step):
            i = start
            while i < stop:
                yield i
                i += step

        def sample_functions(op):
            yield lambda: op(range(10))
            yield lambda: op(range(4, 12))
            yield lambda: op(range(-4, -15, -1))
            yield lambda: op([6.6, 5.5, 7.7])
            yield lambda: op([(3, 4), (1, 2)])
            yield lambda: op(frange(1.1, 3.3, 0.1))
            yield lambda: op([np.nan, -np.inf, np.inf, np.nan])
            yield lambda: op([(3,), (1,), (2,)])

        for fn in sample_functions(op=min):
            self._check_min_max(fn)

        for fn in sample_functions(op=max):
            self._check_min_max(fn)


class TestOperatorMixedTypes(TestCase):

    def test_eq_ne(self):
        for opstr in ('eq', 'ne'):
            op = getattr(operator, opstr)

            @njit
            def func(a, b):
                return op(a, b)

            # all these things should evaluate to being equal or not, all should
            # survive typing.
            things = (1, 0, True, False, 1.0, 2.0, 1.1, 1j, None, "", "1")
            for x, y in itertools.product(things, things):
                self.assertPreciseEqual(func.py_func(x, y), func(x, y))

    def test_cmp(self):
        for opstr in ('gt', 'lt', 'ge', 'le', 'eq', 'ne'):
            op = getattr(operator, opstr)
            @njit
            def func(a, b):
                return op(a, b)

            # numerical things should all be comparable
            things = (1, 0, True, False, 1.0, 0.0, 1.1)
            for x, y in itertools.product(things, things):
                expected = func.py_func(x, y)
                got = func(x, y)
                message = ("%s %s %s does not match between Python and Numba"
                           % (x, opstr, y))
                self.assertEqual(expected, got, message)


class TestIsinstanceBuiltin(TestCase):
    def test_isinstance(self):
        pyfunc = isinstance_usecase
        cfunc = jit(nopython=True)(pyfunc)

        inputs = (
            3,              # int
            5.0,            # float
            "Hello",        # string
            b'world',       # bytes
            1j,             # complex
            [1, 2, 3],      # list
            (1, 3, 3, 3),   # UniTuple
            set([1, 2]),    # set
            (1, 'nba', 2),  # Heterogeneous Tuple
            # {'hello': 2},   # dict - doesn't work as input
            None,
        )

        for inpt in inputs:
            expected = pyfunc(inpt)
            got = cfunc(inpt)
            self.assertEqual(expected, got)

    def test_isinstance_dict(self):
        # Tests typed.Dict and LiteralStrKeyDict
        pyfunc = isinstance_dict
        cfunc = jit(nopython=True)(pyfunc)
        self.assertEqual(pyfunc(), cfunc())

    def test_isinstance_issue9125(self):
        pyfunc = invalid_isinstance_usecase_phi_nopropagate2
        cfunc = jit(nopython=True)(pyfunc)
        self.assertEqual(pyfunc(3), cfunc(3))

    def test_isinstance_numba_types(self):
        # This makes use of type aliasing between python scalars and NumPy
        # scalars, see also test_numba_types()
        pyfunc = isinstance_usecase_numba_types
        cfunc = jit(nopython=True)(pyfunc)

        inputs = (
            (types.int32(1), 'int32'),
            (types.int64(2), 'int64'),
            (types.float32(3.0), 'float32'),
            (types.float64(4.0), 'float64'),
            (types.complex64(5j), 'no match'),
            (typed.List([1, 2]), 'typed list'),
            (typed.Dict.empty(types.int64, types.int64), 'typed dict')
        )

        for inpt, expected in inputs:
            got = cfunc(inpt)
            self.assertEqual(expected, got)

    def test_isinstance_numba_types_2(self):
        pyfunc = isinstance_usecase_numba_types_2
        cfunc = jit(nopython=True)(pyfunc)
        self.assertEqual(pyfunc(), cfunc())

    def test_isinstance_invalid_type(self):
        pyfunc = isinstance_usecase_invalid_type
        cfunc = jit(nopython=True)(pyfunc)

        # valid type
        self.assertTrue(cfunc(3.4))

        # invalid type
        msg = 'Cannot infer numba type of python type'

        with self.assertRaises(errors.TypingError) as raises:
            cfunc(100)

        self.assertIn(msg, str(raises.exception))

    def test_isinstance_exceptions(self):
        fns = [
            (invalid_isinstance_usecase,
             'Cannot infer numba type of python type'),
            (invalid_isinstance_usecase_phi_nopropagate,
             ('isinstance() cannot determine the type of variable "z" due to a '
             'branch.')),
            (invalid_isinstance_optional_usecase,
             ('isinstance() cannot determine the type of variable "z" due to a '
             'branch.')),
            (invalid_isinstance_unsupported_type_usecase(),
             ('isinstance() does not support variables of type "ntpl(')),
        ]

        for fn, msg in fns:
            fn = njit(fn)

            with self.assertRaises(errors.TypingError) as raises:
                fn(100)

            self.assertIn(msg, str(raises.exception))

    def test_combinations(self):
        # Combinatorically test common classes and instances
        def gen_w_arg(clazz_type):
            def impl(x):
                return isinstance(x, clazz_type)
            return impl

        clazz_types = (int, float, complex, str, list, tuple, bytes, set, range,
                       np.int8, np.float32,)
        instances = (1, 2.3, 4j, '5', [6,], (7,), b'8', {9,}, None,
                     (10, 11, 12), (13, 'a', 14j), np.array([15, 16, 17]),
                     np.int8(18), np.float32(19),
                     typed.Dict.empty(types.unicode_type, types.float64),
                     typed.List.empty_list(types.complex128), np.ones(4))

        for ct in clazz_types:
            fn = njit(gen_w_arg(ct))
            for x in instances:
                expected = fn.py_func(x)
                got = fn(x)
                self.assertEqual(got, expected)

    def test_numba_types(self):
        # Check types which are Numba types, this would break without the jit
        # decorator in all cases except numba.typed containers.
        def gen_w_arg(clazz_type):
            def impl():
                return isinstance(1, clazz_type)
            return impl

        clazz_types = (types.Integer, types.Float, types.Array,)

        msg = "Numba type classes.*are not supported"
        for ct in clazz_types:
            fn = njit(gen_w_arg(ct))
            with self.assertRaises(errors.TypingError) as raises:
                fn()
            self.assertRegex(str(raises.exception), msg)

    def test_python_numpy_scalar_alias_problem(self):
        # There's a problem due to Python and NumPy scalars being aliased in the
        # type system. This is because e.g. int scalar values and NumPy np.intp
        # type alias to types.intp. This test merely records this fact.

        @njit
        def foo():
            return isinstance(np.intp(10), int)

        self.assertEqual(foo(), True)
        self.assertEqual(foo.py_func(), False)

        @njit
        def bar():
            return isinstance(1, np.intp)

        self.assertEqual(bar(), True)
        self.assertEqual(bar.py_func(), False)

    def test_branch_prune(self):
        # Check that isinstance branches are pruned allowing otherwise
        # impossible type specific specialisation.

        @njit
        def foo(x):
            if isinstance(x, str):
                return x + 'some_string'
            elif isinstance(x, complex):
                return np.imag(x)
            elif isinstance(x, tuple):
                return len(x)
            else:
                assert 0

        for x in ('string', 1 + 2j, ('a', 3, 4j)):
            expected = foo.py_func(x)
            got = foo(x)
            self.assertEqual(got, expected)


class TestGetattrBuiltin(MemoryLeakMixin, TestCase):
    # Tests the getattr() builtin

    def test_getattr_func_retty(self):

        @njit
        def foo(x):
            attr = getattr(x, '__hash__')
            return attr()

        for x in (1, 2.34, (5, 6, 7)):
            self.assertPreciseEqual(foo(x), foo.py_func(x))

    def test_getattr_value_retty(self):

        @njit
        def foo(x):
            return getattr(x, 'ndim')

        for x in range(3):
            tmp = np.empty((1, ) * x)
            self.assertPreciseEqual(foo(tmp), foo.py_func(tmp))

    def test_getattr_module_obj(self):
        # Consts on modules work ok

        @njit
        def foo():
            return getattr(np, 'pi')

        self.assertPreciseEqual(foo(), foo.py_func())

    def test_getattr_module_obj_not_implemented(self):
        # Functions on modules do not work at present

        @njit
        def foo():
            return getattr(np, 'cos')(1)

        with self.assertRaises(errors.TypingError) as raises:
            foo()

        msg = "Returning function objects is not implemented"
        self.assertIn(msg, str(raises.exception))

    def test_getattr_raises_attribute_error(self):

        invalid_attr = '__not_a_valid_attr__'

        @njit
        def foo(x):
            return getattr(x, invalid_attr)

        with self.assertRaises(AttributeError) as raises:
            foo(1.23)

        self.assertIn(f"'float64' has no attribute '{invalid_attr}'",
                      str(raises.exception))

    def test_getattr_with_default(self):
        # Checks returning a default works

        @njit
        def foo(x, default):
            return getattr(x, '__not_a_valid_attr__', default)

        for x, y in zip((1, 2.34, (5, 6, 7),), (None, 20, 'some_string')):
            self.assertPreciseEqual(foo(x, y), foo.py_func(x, y))

    def test_getattr_non_literal_str(self):

        @njit
        def foo(x, nonliteral_str):
            return getattr(x, nonliteral_str)

        with self.assertRaises(errors.TypingError) as raises:
            foo(1, '__hash__')

        msg = "argument 'name' must be a literal string"
        self.assertIn(msg, str(raises.exception))

    def test_getattr_no_optional_type_generated(self):

        @njit
        def default_hash():
            return 12345

        @njit
        def foo():
            hash_func = getattr(np.ones(1), "__not_a_valid_attr__",
                                default_hash)
            return hash_func() # Optionals have no call support

        self.assertPreciseEqual(foo(), foo.py_func())


class TestHasattrBuiltin(MemoryLeakMixin, TestCase):
    # Tests the hasattr() builtin

    def test_hasattr(self):

        @njit
        def foo(x):
            return hasattr(x, '__hash__'), hasattr(x, '__not_a_valid_attr__')

        ty = types.int64
        for x in (1, 2.34, (5, 6, 7), typed.Dict.empty(ty, ty),
                  typed.List.empty_list(ty), np.ones(4), 'ABC'):
            self.assertPreciseEqual(foo(x), foo.py_func(x))

    def test_hasattr_non_const_attr(self):
        # This tests that an error is raised in the case that a hasattr() call
        # is made on an attribute that cannot be resolved as a compile time
        # constant (there's a phi in the way!).

        @njit
        def foo(pred):
            if pred > 3:
                attr = "__hash__"
            else:
                attr = "__str__"

            hasattr(1, attr)

        with self.assertRaises(errors.NumbaTypeError) as raises:
            foo(6)

        msg = ('hasattr() cannot determine the type of variable '
               '"attr" due to a branch.')
        self.assertIn(msg, str(raises.exception))


class TestStrAndReprBuiltin(MemoryLeakMixin, TestCase):

    def test_str_default(self):

        @njit
        def foo():
            return str()

        self.assertEqual(foo(), foo.py_func())

    def test_str_object_kwarg(self):

        @njit
        def foo(x):
            return str(object=x)

        value = "a string"
        self.assertEqual(foo(value), foo.py_func(value))

    def test_str_calls_dunder_str(self):

        @njit
        def foo(x):
            return str(x)

        Dummy, DummyType = self.make_dummy_type()
        dummy = Dummy()
        string_repr = "this is the dummy object str"
        Dummy.__str__= lambda inst: string_repr

        @overload_method(DummyType, "__str__")
        def ol_dummy_string(dummy):
            def impl(dummy):
                return string_repr
            return impl

        @overload_method(DummyType, "__repr__")
        def ol_dummy_repr(dummy):
            def impl(dummy):
                return "SHOULD NOT BE CALLED"
            return impl

        self.assertEqual(foo(dummy), foo.py_func(dummy))

    def test_str_falls_back_to_repr(self):

        @njit
        def foo(x):
            return str(x)

        Dummy, DummyType = self.make_dummy_type()
        dummy = Dummy()
        string_repr = "this is the dummy object repr"
        Dummy.__repr__= lambda inst: string_repr

        @overload_method(DummyType, "__repr__")
        def ol_dummy_repr(dummy):
            def impl(dummy):
                return string_repr
            return impl

        self.assertEqual(foo(dummy), foo.py_func(dummy))

    def test_repr(self):
        @njit
        def foo(x):
            return repr(x), x

        for x in ("abc", False, 123):
            self.assertEqual(foo(x), foo.py_func(x))

    def test_repr_fallback(self):
        # checks str/repr fallback, there's no overloaded __str__ or __repr__
        # for the dummy type so it has to use generic '<object type:(type)>'
        # string for the `repr` call.

        Dummy, DummyType = self.make_dummy_type()
        dummy = Dummy()
        string_repr = f"<object type:{typeof(dummy)}>"
        Dummy.__repr__= lambda inst: string_repr

        @box(DummyType)
        def box_dummy(typ, obj, c):
            clazobj = c.pyapi.unserialize(c.pyapi.serialize_object(Dummy))
            return c.pyapi.call_function_objargs(clazobj, ())

        @njit
        def foo(x):
            return str(x)

        self.assertEqual(foo(dummy), foo.py_func(dummy))



if __name__ == '__main__':
    unittest.main()
