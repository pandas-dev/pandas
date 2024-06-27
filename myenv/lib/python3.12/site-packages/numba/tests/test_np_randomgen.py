import numba
import numpy as np
import sys
import itertools
import gc

from numba import types
from numba.tests.support import TestCase, MemoryLeakMixin
from numba.np.random.generator_methods import _get_proper_func
from numba.np.random.generator_core import next_uint32, next_uint64, next_double
from numpy.random import MT19937, Generator
from numba.core.errors import TypingError
from numba.tests.support import run_in_new_process_caching, SerialMixin


# TODO: Following testing tolerance adjustments should be reduced
# once NumPy Generator's fmadd issue described below is resolved:
# https://github.com/numba/numba/pull/8038#issuecomment-1165571368
# The progress is being tracked as one of the tasks in:
# https://github.com/numba/numba/issues/8519
adjusted_ulp_prec = 2048


class TestHelperFuncs(TestCase):
    def test_proper_func_provider(self):
        def test_32bit_func():
            return 32

        def test_64bit_func():
            return 64

        self.assertEqual(_get_proper_func(test_32bit_func, test_64bit_func,
                         np.float64)[0](), 64)
        self.assertEqual(_get_proper_func(test_32bit_func, test_64bit_func,
                         np.float32)[0](), 32)

        # With any other datatype it should return a TypingError
        with self.assertRaises(TypingError) as raises:
            _get_proper_func(test_32bit_func, test_64bit_func, np.int32)
        self.assertIn(
            'Argument dtype is not one of the expected type(s)',
            str(raises.exception)
        )
        with self.assertRaises(TypingError) as raises:
            _get_proper_func(test_32bit_func, test_64bit_func, types.float64)
        self.assertIn(
            'Argument dtype is not one of the expected type(s)',
            str(raises.exception)
        )

    def test_check_types(self):
        rng = np.random.default_rng(1)
        py_func = lambda x: x.normal(loc=(0,))
        numba_func = numba.njit(cache=True)(py_func)
        with self.assertRaises(TypingError) as raises:
            numba_func(rng)
        self.assertIn(
            'Argument loc is not one of the expected type(s): '
            + '[<class \'numba.core.types.scalars.Float\'>, '
            + '<class \'numba.core.types.scalars.Integer\'>, '
            + '<class \'int\'>, <class \'float\'>]',
            str(raises.exception)
        )

    def test_integers_arg_check(self):
        rng = np.random.default_rng(1)
        py_func = lambda x, low, high, dtype: \
            x.integers(low=low, high=high, dtype=dtype, endpoint=True)
        numba_func = numba.njit()(py_func)
        numba_func_low = numba.njit()(py_func)

        py_func = lambda x, low, high, dtype: \
            x.integers(low=low, high=high, dtype=dtype, endpoint=False)
        numba_func_endpoint_false = numba.njit()(py_func)

        cases = [
            # low, high, dtype
            (np.iinfo(np.uint8).min, np.iinfo(np.uint8).max, np.uint8),
            (np.iinfo(np.int8).min, np.iinfo(np.int8).max, np.int8),
            (np.iinfo(np.uint16).min, np.iinfo(np.uint16).max, np.uint16),
            (np.iinfo(np.int16).min, np.iinfo(np.int16).max, np.int16),
            (np.iinfo(np.uint32).min, np.iinfo(np.uint32).max, np.uint32),
            (np.iinfo(np.int32).min, np.iinfo(np.int32).max, np.int32),
        ]
        for low, high, dtype in cases:
            with self.subTest(low=low, high=high, dtype=dtype):
                with self.assertRaises(ValueError) as raises:
                    # min - 1
                    numba_func_low(rng, low - 1, high, dtype)
                self.assertIn(
                    'low is out of bounds',
                    str(raises.exception)
                )

                with self.assertRaises(ValueError) as raises:
                    # max + 1, endpoint=True
                    numba_func(rng, low, high + 1, dtype)
                self.assertIn(
                    'high is out of bounds',
                    str(raises.exception)
                )

                with self.assertRaises(ValueError) as raises:
                    # max + 2, endpoint=False
                    numba_func_endpoint_false(rng, low, high + 2, dtype)
                self.assertIn(
                    'high is out of bounds',
                    str(raises.exception)
                )

        low, high, dtype = (np.iinfo(np.uint64).min,
                            np.iinfo(np.uint64).max, np.uint64)
        with self.assertRaises(ValueError) as raises:
            # min - 1
            numba_func_low(rng, low - 1, high, dtype)
        self.assertIn(
            'low is out of bounds',
            str(raises.exception)
        )

        low, high, dtype = (np.iinfo(np.int64).min,
                            np.iinfo(np.int64).max, np.int64)
        with self.assertRaises(ValueError) as raises:
            # max + 1, endpoint=True
            numba_func(rng, low, high + 1, dtype)
        self.assertIn(
            'high is out of bounds',
            str(raises.exception)
        )

        with self.assertRaises(ValueError) as raises:
            # max + 2, endpoint=False
            numba_func_endpoint_false(rng, low, high + 2, dtype)
        self.assertIn(
            'high is out of bounds',
            str(raises.exception)
        )

        with self.assertRaises(ValueError) as raises:
            numba_func(rng, 105, 100, np.uint32)
        self.assertIn(
            'low is greater than high in given interval',
            str(raises.exception)
        )


def test_generator_caching():
    nb_rng = np.random.default_rng(1)
    np_rng = np.random.default_rng(1)
    py_func = lambda x: x.random(10)
    numba_func = numba.njit(cache=True)(py_func)
    assert np.allclose(np_rng.random(10), numba_func(nb_rng))


class TestRandomGenerators(MemoryLeakMixin, TestCase):
    def check_numpy_parity(self, distribution_func,
                           bitgen_type=None, seed=None,
                           test_size=None, test_dtype=None,
                           ulp_prec=5):

        distribution_func = numba.njit(distribution_func)
        if seed is None:
            seed = 1
        if bitgen_type is None:
            numba_rng_instance = np.random.default_rng(seed=seed)
            numpy_rng_instance = np.random.default_rng(seed=seed)
        else:
            numba_rng_instance = Generator(bitgen_type(seed))
            numpy_rng_instance = Generator(bitgen_type(seed))

        # Check parity for different size cases
        numba_res = distribution_func(numba_rng_instance,
                                      test_size, test_dtype)
        numpy_res = distribution_func.py_func(numpy_rng_instance,
                                              test_size, test_dtype)

        if (isinstance(numba_res, np.ndarray) and
            np.issubdtype(numba_res.dtype, np.floating)) \
                or isinstance(numba_res, float):
            # Float scalars and arrays
            np.testing.assert_array_max_ulp(numpy_res, numba_res,
                                            maxulp=ulp_prec, dtype=test_dtype)
        else:
            # Bool/int scalars and arrays
            np.testing.assert_equal(numba_res, numpy_res)

        # Check if the end state of both BitGenerators is same
        # after drawing the distributions
        numba_gen_state = numba_rng_instance.bit_generator.state['state']
        numpy_gen_state = numpy_rng_instance.bit_generator.state['state']

        for _state_key in numpy_gen_state:
            self.assertPreciseEqual(numba_gen_state[_state_key],
                                    numpy_gen_state[_state_key])

    def _test_bitgen_func_parity(self, func_name, bitgen_func, seed=1):
        numba_rng_instance = np.random.default_rng(seed=seed)
        numpy_rng_instance = np.random.default_rng(seed=seed)

        numpy_func = getattr(numpy_rng_instance.bit_generator.ctypes, func_name)
        numpy_res = numpy_func(numpy_rng_instance.bit_generator.ctypes.state)

        numba_func = numba.njit(lambda x: bitgen_func(x.bit_generator))
        numba_res = numba_func(numba_rng_instance)

        self.assertPreciseEqual(numba_res, numpy_res)

    def _check_invalid_types(self, dist_func, arg_list,
                             valid_args, invalid_args):
        rng = np.random.default_rng()
        for idx, _arg in enumerate(arg_list):
            curr_args = valid_args.copy()
            curr_args[idx] = invalid_args[idx]
            curr_args = [rng] + curr_args
            nb_dist_func = numba.njit(dist_func)
            with self.assertRaises(TypingError) as raises:
                nb_dist_func(*curr_args)
            self.assertIn(
                f'Argument {_arg} is not one of the expected type(s):',
                str(raises.exception)
            )

    def test_npgen_boxing_unboxing(self):
        rng_instance = np.random.default_rng()
        numba_func = numba.njit(lambda x: x)
        self.assertEqual(rng_instance, numba_func(rng_instance))
        self.assertEqual(id(rng_instance), id(numba_func(rng_instance)))

    def test_npgen_boxing_refcount(self):
        rng_instance = np.random.default_rng()
        no_box = numba.njit(lambda x:x.random())
        do_box = numba.njit(lambda x:x)

        y = do_box(rng_instance)
        gc.collect()
        ref_1 = sys.getrefcount(rng_instance)
        del y
        no_box(rng_instance)
        gc.collect()
        ref_2 = sys.getrefcount(rng_instance)

        self.assertEqual(ref_1, ref_2 + 1)

    def test_bitgen_funcs(self):
        func_names = ["next_uint32", "next_uint64", "next_double"]
        funcs = [next_uint32, next_uint64, next_double]

        for _func, _func_name in zip(funcs, func_names):
            with self.subTest(_func=_func, _func_name=_func_name):
                self._test_bitgen_func_parity(_func_name, _func)

    def test_integers(self):
        test_sizes = [None, (), (100,), (10, 20, 30)]
        test_dtypes = [np.int64, np.int32, np.int16, np.int8,
                       np.uint64, np.uint32, np.uint16, np.uint8]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.integers(0, 100)
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None, ulp_prec=0)

        dist_func = lambda x, size, dtype:\
            x.integers(5, 10, size=size, dtype=dtype)
        for _size in test_sizes:
            for _dtype in test_dtypes:
                for _bitgen in bitgen_types:
                    with self.subTest(_size=_size, _dtype=_dtype,
                                      _bitgen=_bitgen):
                        self.check_numpy_parity(dist_func, _bitgen,
                                                None, _size, _dtype, 0)

        # Checking dtype = bool seperately
        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:\
            x.integers(False, True, size=size, dtype=np.bool_)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size,
                                  _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, np.bool_, 0)

        # Test dtype casting for high and low
        dist_func = lambda x, size, dtype: \
            x.integers(np.uint8(0), np.int64(100))
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None)

        dist_func = lambda x, low, high, size, dtype, endpoint:\
            x.integers(low=low, high=high, size=size,
                       dtype=dtype, endpoint=endpoint)
        self._check_invalid_types(dist_func,
                                  ['low', 'high', 'size', 'dtype', 'endpoint'],
                                  [1, 5, (1,), np.int64, True],
                                  ['x', 'x', ('x',), np.float64, 'x'])

    # Testing .integers() dtype wise
    def test_integers_cases(self):
        cases = [
            # low, high, dtype
            (5, 6, np.uint64), # rng == 0 (rng stands for range)
            (5, 100, np.uint64), # rng <= 0xFFFFFFFF
            (0, 0xFFFFFFFFFF, np.uint64), # rng > 0xFFFFFFFF
            (0, 0xFFFFFFFFFFFFFFFF - 1, np.uint64),# rng == 0xFFFFFFFFFFFFFFFF-1
            (0, 0xFFFFFFFFFFFFFFFF, np.uint64), # rng == 0xFFFFFFFFFFFFFFFF

            (5, 6, np.int64), # rng == 0
            (5, 100, np.int64), # rng <= 0xFFFFFFFF
            (0, 0xFFFFFFFFFF, np.int64), # rng > 0xFFFFFFFF
            (0, 0xFFFFFFFFFFFFFFF - 1, np.int64), # rng == 0xFFFFFFFFFFFFFFF - 1
            (0, 0xFFFFFFFFFFFFFFF, np.int64), # rng == 0xFFFFFFFFFFFFFFF
            (-0xFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFF, np.int64), # min/max

            (5, 6, np.uint32), # rng == 0
            (5, 100, np.uint32), # rng < 0xFFFFFFFF
            (0, 0xFFFFFFFF - 1, np.uint32), # rng == 0xFFFFFFFF - 1
            (0, 0xFFFFFFFF, np.uint32), # rng == 0xFFFFFFFF

            (5, 6, np.int32), # rng == 0
            (5, 100, np.int32), # rng < 0xFFFFFFFF
            (0, 0xFFFFFFF - 1, np.int32), # rng == 0xFFFFFFF - 1
            (0, 0xFFFFFFF, np.int32), # rng == 0xFFFFFFF
            (-0xFFFFFFF, 0xFFFFFFF, np.int32),

            (5, 6, np.uint16), # rng == 0
            (5, 100, np.uint16), # rng < 0xFFFF
            (0, 0xFFFF - 1, np.uint16), # rng == 0xFFFF - 1
            (0, 0xFFFF, np.uint16), # rng == 0xFFFF

            (5, 6, np.int16), # rng == 0
            (5, 10, np.int16), # rng < 0xFFF
            (0, 0xFFF - 1, np.int16), # rng == 0xFFF - 1
            (0, 0xFFF, np.int16), # rng == 0xFFF
            (-0xFFF, 0xFFF, np.int16),

            (5, 6, np.uint8), # rng == 0
            (5, 10, np.uint8), # rng < 0xFF
            (0, 0xFF - 1, np.uint8), # rng == 0xFF - 1
            (0, 0xFF, np.uint8), # rng == 0xFF

            (5, 6, np.int8), # rng == 0
            (5, 10, np.int8), # rng < 0xF
            (0, 0xF - 1, np.int8), # rng == 0xF-1
            (0, 0xF, np.int8), # rng == 0xF
            (-0xF, 0xF, np.int8),
        ]
        size = (2, 3)

        for low, high, dtype in cases:
            with self.subTest(low=low, high=high, dtype=dtype):
                dist_func = lambda x, size, dtype:\
                    x.integers(low, high, size=size, dtype=dtype)
                self.check_numpy_parity(dist_func, None,
                                        None, size, dtype, 0)

    def test_random(self):
        test_sizes = [None, (), (100,), (10, 20, 30)]
        test_dtypes = [np.float32, np.float64]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.random()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None)

        dist_func = lambda x, size, dtype:x.random(size=size, dtype=dtype)

        for _size in test_sizes:
            for _dtype in test_dtypes:
                for _bitgen in bitgen_types:
                    with self.subTest(_size=_size, _dtype=_dtype,
                                      _bitgen=_bitgen):
                        self.check_numpy_parity(dist_func, _bitgen,
                                                None, _size, _dtype)
        dist_func = lambda x, size, dtype:\
            x.random(size=size, dtype=dtype)
        self._check_invalid_types(dist_func, ['size', 'dtype'],
                                  [(1,), np.float64], [('x',), 0.])

    def test_standard_normal(self):
        test_sizes = [None, (), (100,), (10, 20, 30)]
        test_dtypes = [np.float32, np.float64]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.standard_normal()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None)

        dist_func = lambda x, size, dtype:\
            x.standard_normal(size=size, dtype=dtype)

        for _size in test_sizes:
            for _dtype in test_dtypes:
                for _bitgen in bitgen_types:
                    with self.subTest(_size=_size, _dtype=_dtype,
                                      _bitgen=_bitgen):
                        self.check_numpy_parity(dist_func, _bitgen,
                                                None, _size, _dtype)
        dist_func = lambda x, size, dtype:\
            x.standard_normal(size=size, dtype=dtype)
        self._check_invalid_types(dist_func, ['size', 'dtype'],
                                  [(1,), np.float32], [('x',), 0])

    def test_standard_exponential(self):
        test_sizes = [None, (), (100,), (10, 20, 30)]
        test_dtypes = [np.float32, np.float64]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.standard_exponential()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None)

        dist_func = lambda x, size, dtype:\
            x.standard_exponential(size=size, dtype=dtype)

        for _size in test_sizes:
            for _dtype in test_dtypes:
                for _bitgen in bitgen_types:
                    with self.subTest(_size=_size, _dtype=_dtype,
                                      _bitgen=_bitgen):
                        self.check_numpy_parity(dist_func, _bitgen,
                                                None, _size, _dtype)

        dist_func = lambda x, method, size, dtype:\
            x.standard_exponential(method=method, size=size, dtype=dtype)
        self._check_invalid_types(dist_func, ['method', 'size', 'dtype'],
                                  ['zig', (1,), np.float32], [0, ('x',), 0])

    def test_standard_exponential_inv(self):
        test_sizes = [None, (), (100,), (10, 20, 30)]
        test_dtypes = [np.float32, np.float64]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:\
            x.standard_exponential(size=size, dtype=dtype,  method='inv')
        for _size in test_sizes:
            for _dtype in test_dtypes:
                for _bitgen in bitgen_types:
                    with self.subTest(_size=_size, _dtype=_dtype,
                                      _bitgen=_bitgen):
                        self.check_numpy_parity(dist_func, _bitgen,
                                                None, _size, _dtype)

    def test_standard_gamma(self):
        test_sizes = [None, (), (100,), (10, 20, 30)]
        test_dtypes = [np.float32, np.float64]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype: \
            x.standard_gamma(shape=5.0, size=size, dtype=dtype)
        for _size in test_sizes:
            for _dtype in test_dtypes:
                for _bitgen in bitgen_types:
                    with self.subTest(_size=_size, _dtype=_dtype,
                                      _bitgen=_bitgen):
                        self.check_numpy_parity(dist_func, _bitgen,
                                                None, _size, _dtype,
                                                adjusted_ulp_prec)
        dist_func = lambda x, shape, size, dtype:\
            x.standard_gamma(shape=shape, size=size, dtype=dtype)
        self._check_invalid_types(dist_func, ['shape', 'size', 'dtype'],
                                  [5.0, (1,), np.float32], ['x', ('x',), 0])

    def test_normal(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.normal()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None,
                                    ulp_prec=adjusted_ulp_prec)

        dist_func = lambda x, size, dtype:x.normal(loc=1.5, scale=3, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)
        dist_func = lambda x, loc, scale, size:\
            x.normal(loc=loc, scale=scale, size=size)
        self._check_invalid_types(dist_func, ['loc', 'scale', 'size'],
                                  [1.5, 3, (1,)], ['x', 'x', ('x',)])

    def test_uniform(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.uniform()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None,
                                    ulp_prec=adjusted_ulp_prec)

        dist_func = lambda x, size, dtype:x.uniform(low=1.5, high=3, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)
        dist_func = lambda x, low, high, size:\
            x.uniform(low=low, high=high, size=size)
        self._check_invalid_types(dist_func, ['low', 'high', 'size'],
                                  [1.5, 3, (1,)], ['x', 'x', ('x',)])

    def test_exponential(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.exponential()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None)

        dist_func = lambda x, size, dtype:x.exponential(scale=1.5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None)
        dist_func = lambda x, scale, size:\
            x.exponential(scale=scale, size=size)
        self._check_invalid_types(dist_func, ['scale', 'size'],
                                  [1.5, (1,)], ['x', ('x',)])

    def test_gamma(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.gamma(shape=5.0, scale=1.5,
                                                  size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)
        dist_func = lambda x, shape, scale, size:\
            x.gamma(shape=shape, scale=scale, size=size)
        self._check_invalid_types(dist_func, ['shape', 'scale', 'size'],
                                  [5.0, 1.5, (1,)], ['x', 'x', ('x',)])

    def test_beta(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.beta(a=1.5, b=2.5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, a, b, size:x.beta(a=a, b=b, size=size)
        self._check_invalid_types(dist_func, ['a', 'b', 'size'],
                                  [5.0, 1.5, (1,)], ['x', 'x', ('x',)])

    def test_f(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.f(dfnum=2, dfden=3, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, dfnum, dfden, size:\
            x.f(dfnum=dfnum, dfden=dfden, size=size)
        self._check_invalid_types(dist_func, ['dfnum', 'dfden', 'size'],
                                  [5, 1, (1,)], ['x', 'x', ('x',)])

    def test_chisquare(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.chisquare(df=2, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, df, size:\
            x.chisquare(df=df, size=size)
        self._check_invalid_types(dist_func, ['df', 'size'],
                                  [2, (1,)], ['x', ('x',)])

    def test_standard_cauchy(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.standard_cauchy()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None)

        dist_func = lambda x, size, dtype:x.standard_cauchy(size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None)

        dist_func = lambda x, size:x.standard_cauchy(size=size)
        self._check_invalid_types(dist_func, ['size'],
                                  [(1,)], [('x',)])

    def test_pareto(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.pareto(a=1.0, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None)

        dist_func = lambda x, a, size:x.pareto(a=a, size=size)
        self._check_invalid_types(dist_func, ['a', 'size'],
                                  [1, (1,)], ['x', ('x',)])

    def test_weibull(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.weibull(a=1.0, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None)

        dist_func = lambda x, a, size:x.weibull(a=a, size=size)
        self._check_invalid_types(dist_func, ['a', 'size'],
                                  [1, (1,)], ['x', ('x',)])

    def test_power(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.power(a=0.75, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None)

        dist_func = lambda x, a, size:x.power(a=a, size=size)
        self._check_invalid_types(dist_func, ['a', 'size'],
                                  [0.75, (1,)], ['x', ('x',)])

    def test_laplace(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.laplace()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None,
                                    ulp_prec=adjusted_ulp_prec)

        dist_func = lambda x, size, dtype:\
            x.laplace(loc=1.0, scale=1.5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, loc, scale, size:\
            x.laplace(loc=loc, scale=scale, size=size)
        self._check_invalid_types(dist_func, ['loc', 'scale', 'size'],
                                  [1.0, 1.5, (1,)], ['x', 'x', ('x',)])

    def test_logistic(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.logistic()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None,
                                    ulp_prec=adjusted_ulp_prec)

        dist_func = lambda x, size, dtype:\
            x.logistic(loc=1.0,scale=1.5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, loc, scale, size:\
            x.logistic(loc=loc, scale=scale, size=size)
        self._check_invalid_types(dist_func, ['loc', 'scale', 'size'],
                                  [1.0, 1.5, (1,)], ['x', 'x', ('x',)])

    def test_lognormal(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.lognormal()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None,
                                    ulp_prec=adjusted_ulp_prec)

        dist_func = lambda x, size, dtype:\
            x.lognormal(mean=5.0, sigma=1.5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, mean, sigma, size:\
            x.lognormal(mean=mean, sigma=sigma, size=size)
        self._check_invalid_types(dist_func, ['mean', 'sigma', 'size'],
                                  [1.0, 1.5, (1,)], ['x', 'x', ('x',)])

    def test_rayleigh(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        # Test with no arguments
        dist_func = lambda x, size, dtype:x.rayleigh()
        with self.subTest():
            self.check_numpy_parity(dist_func, test_size=None,
                                    test_dtype=None)

        dist_func = lambda x, size, dtype:x.rayleigh(scale=1.5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None)

        dist_func = lambda x, scale, size:x.rayleigh(scale=scale, size=size)
        self._check_invalid_types(dist_func, ['scale', 'size'],
                                  [1.5, (1,)], ['x', ('x',)])

    def test_standard_t(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.standard_t(df=2, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, df, size:x.standard_t(df=df, size=size)
        self._check_invalid_types(dist_func, ['df', 'size'],
                                  [2, (1,)], ['x', ('x',)])

    def test_wald(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.wald(mean=5.0, scale=1.5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, mean, scale, size:\
            x.wald(mean=mean, scale=scale, size=size)
        self._check_invalid_types(dist_func, ['mean', 'scale', 'size'],
                                  [1.0, 1.5, (1,)], ['x', 'x', ('x',)])

    def test_geometric(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.geometric(p=0.75, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, p, size:x.geometric(p=p, size=size)
        self._check_invalid_types(dist_func, ['p', 'size'],
                                  [0.75, (1,)], ['x', ('x',)])

    def test_zipf(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.zipf(a=1.5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None)

        dist_func = lambda x, a, size:x.zipf(a=a, size=size)
        self._check_invalid_types(dist_func, ['a', 'size'],
                                  [1, (1,)], ['x', ('x',)])

    def test_triangular(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:\
            x.triangular(left=0, mode=3, right=5, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None)

        dist_func = lambda x, left, mode, right, size:\
            x.triangular(left=left, mode=mode, right=right, size=size)
        self._check_invalid_types(dist_func, ['left', 'mode', 'right', 'size'],
                                  [0, 3, 5, (1,)], ['x', 'x', 'x', ('x',)])

    def test_poisson(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:x.poisson(lam=15, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, lam, size:x.poisson(lam=lam, size=size)
        self._check_invalid_types(dist_func, ['lam', 'size'],
                                  [15, (1,)], ['x', ('x',)])

    def test_negative_binomial(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:\
            x.negative_binomial(n=1, p=0.1, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, n, p, size:\
            x.negative_binomial(n=n, p=p, size=size)
        self._check_invalid_types(dist_func, ['n', 'p', 'size'],
                                  [1, 0.75, (1,)], ['x', 'x', ('x',)])

    # NumPy tests at:
    # https://github.com/numpy/numpy/blob/95e3e7f445407e4f355b23d6a9991d8774f0eb0c/numpy/random/tests/test_generator_mt19937.py#L936
    # Written in following format for semblance with existing Generator tests.
    def test_shuffle(self):
        test_sizes = [(10, 20, 30)]
        bitgen_types = [None, MT19937]
        axes = [0, 1, 2]

        for _size, _bitgen, _axis in itertools.product(test_sizes,
                                                       bitgen_types,
                                                       axes):
            with self.subTest(_size=_size, _bitgen=_bitgen, _axis=_axis):
                def dist_func(x, size, dtype):
                    arr = x.random(size=size)
                    x.shuffle(arr, axis=_axis)
                    return arr
                self.check_numpy_parity(dist_func, _bitgen,
                                        None, _size, None,
                                        0)

    def test_shuffle_empty(self):
        a = np.array([])
        b = np.array([])

        def dist_func(x, arr):
            x.shuffle(arr)
            return arr

        nb_func = numba.njit(dist_func)
        rng = lambda: np.random.default_rng(1)

        self.assertPreciseEqual(dist_func(rng(), a), nb_func(rng(), b))

    def test_shuffle_check(self):
        self.disable_leak_check()

        def dist_func(x, arr, axis):
            x.shuffle(arr, axis=axis)
            return arr

        self._check_invalid_types(dist_func, ['x', 'axis'],
                                  [np.array([3,4,5]), 0], ['x', 'x'])

        rng = np.random.default_rng(1)
        with self.assertRaises(IndexError) as raises:
            numba.njit(dist_func)(rng, np.array([3,4,5]), 2)
        self.assertIn(
            'Axis is out of bounds for the given array',
            str(raises.exception)
        )

    # NumPy tests at:
    # https://github.com/numpy/numpy/blob/95e3e7f445407e4f355b23d6a9991d8774f0eb0c/numpy/random/tests/test_generator_mt19937.py#L1030
    # Written in following format for semblance with existing Generator tests.
    def test_permutation(self):
        test_sizes = [(10, 20, 30)]
        bitgen_types = [None, MT19937]
        axes = [0, 1, 2, -1, -2]

        for _size, _bitgen, _axis in itertools.product(test_sizes,
                                                       bitgen_types,
                                                       axes):
            with self.subTest(_size=_size, _bitgen=_bitgen, _axis=_axis):
                def dist_func(x, size, dtype):
                    arr = x.random(size=size)
                    return x.permutation(arr, axis=1)
                self.check_numpy_parity(dist_func, _bitgen,
                                        None, _size, None,
                                        0)

        # Test that permutation is actually done on a copy of the array
        dist_func = numba.njit(lambda rng, arr: rng.permutation(arr))
        rng = np.random.default_rng()
        arr = rng.random(size=(10, 20))
        arr_cpy = arr.copy()
        dist_func(rng, arr)
        self.assertPreciseEqual(arr, arr_cpy)

    def test_permutation_exception(self):
        self.disable_leak_check()

        def dist_func(x, arr, axis):
            return x.permutation(arr, axis=axis)

        self._check_invalid_types(dist_func, ['x', 'axis'],
                                  [np.array([3,4,5]), 0], ['x', 'x'])

        rng = np.random.default_rng(1)
        with self.assertRaises(IndexError) as raises:
            numba.njit(dist_func)(rng, np.array([3,4,5]), 2)
        self.assertIn(
            'Axis is out of bounds for the given array',
            str(raises.exception)
        )
        with self.assertRaises(IndexError) as raises:
            numba.njit(dist_func)(rng, np.array([3,4,5]), -2)
        self.assertIn(
            'Axis is out of bounds for the given array',
            str(raises.exception)
        )

    def test_permutation_empty(self):
        a = np.array([])
        b = np.array([])

        def dist_func(x, arr):
            return x.permutation(arr)

        nb_func = numba.njit(dist_func)
        rng = lambda: np.random.default_rng(1)

        self.assertPreciseEqual(dist_func(rng(), a), nb_func(rng(), b))

    def test_noncentral_chisquare(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:\
            x.noncentral_chisquare(3.0, 20.0, size=size)
        for _size, _bitgen in itertools.product(test_sizes, bitgen_types):
            with self.subTest(_size=_size, _bitgen=_bitgen):
                self.check_numpy_parity(dist_func, _bitgen,
                                        None, _size, None)

        dist_func = lambda x, df, nonc, size:\
            x.noncentral_chisquare(df=df, nonc=nonc, size=size)
        valid_args = [3.0, 5.0, (1,)]
        self._check_invalid_types(dist_func, ['df', 'nonc', 'size'],
                                  valid_args, ['x', 'x', ('x',)])

        # Test argument bounds
        rng = np.random.default_rng()
        valid_args = [rng] + valid_args
        nb_dist_func = numba.njit(dist_func)
        with self.assertRaises(ValueError) as raises:
            curr_args = valid_args.copy()
            # Change df to an invalid value
            curr_args[1] = 0
            nb_dist_func(*curr_args)
        self.assertIn('df <= 0', str(raises.exception))
        with self.assertRaises(ValueError) as raises:
            curr_args = valid_args.copy()
            # Change nonc to an invalid value
            curr_args[2] = -1
            nb_dist_func(*curr_args)
        self.assertIn('nonc < 0', str(raises.exception))
        # Exceptions leak references
        self.disable_leak_check()

    def test_noncentral_f(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:\
            x.noncentral_f(3.0, 20.0, 3.0, size=size)
        for _size, _bitgen in itertools.product(test_sizes, bitgen_types):
            with self.subTest(_size=_size, _bitgen=_bitgen):
                self.check_numpy_parity(dist_func, _bitgen,
                                        None, _size, None,
                                        adjusted_ulp_prec)

        dist_func = lambda x, dfnum, dfden, nonc, size:\
            x.noncentral_f(dfnum=dfnum, dfden=dfden, nonc=nonc, size=size)
        valid_args = [3.0, 5.0, 3.0, (1,)]
        self._check_invalid_types(dist_func, ['dfnum', 'dfden', 'nonc', 'size'],
                                  valid_args, ['x', 'x', 'x', ('x',)])

        # Test argument bounds
        rng = np.random.default_rng()
        valid_args = [rng] + valid_args
        nb_dist_func = numba.njit(dist_func)
        with self.assertRaises(ValueError) as raises:
            curr_args = valid_args.copy()
            # Change dfnum to an invalid value
            curr_args[1] = 0
            nb_dist_func(*curr_args)
        self.assertIn('dfnum <= 0', str(raises.exception))
        with self.assertRaises(ValueError) as raises:
            curr_args = valid_args.copy()
            # Change dfden to an invalid value
            curr_args[2] = 0
            nb_dist_func(*curr_args)
        self.assertIn('dfden <= 0', str(raises.exception))
        with self.assertRaises(ValueError) as raises:
            curr_args = valid_args.copy()
            # Change nonc to an invalid value
            curr_args[3] = -1
            nb_dist_func(*curr_args)
        self.assertIn('nonc < 0', str(raises.exception))
        # Exceptions leak references
        self.disable_leak_check()

    def test_logseries(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:\
            x.logseries(0.3, size=size)
        for _size, _bitgen in itertools.product(test_sizes, bitgen_types):
            with self.subTest(_size=_size, _bitgen=_bitgen):
                self.check_numpy_parity(dist_func, _bitgen,
                                        None, _size, None)

        dist_func = lambda x, p, size:\
            x.logseries(p=p, size=size)
        valid_args = [0.3, (1,)]
        self._check_invalid_types(dist_func, ['p', 'size'],
                                  valid_args, ['x', ('x',)])

        # Test argument bounds
        rng = np.random.default_rng(1)
        valid_args = [rng] + valid_args
        nb_dist_func = numba.njit(dist_func)
        for _p in [-0.1, 1, np.nan]:
            with self.assertRaises(ValueError) as raises:
                curr_args = valid_args.copy()
                # Change p to an invalid negative, positive and nan value
                curr_args[1] = _p
                nb_dist_func(*curr_args)
            self.assertIn('p < 0, p >= 1 or p is NaN', str(raises.exception))
        # Exceptions leak references
        self.disable_leak_check()

    def test_binomial(self):
        # For this test dtype argument is never used, so we pass [None] as dtype
        # to make sure it runs only once with default system type.

        test_sizes = [None, (), (100,), (10, 20, 30)]
        bitgen_types = [None, MT19937]

        dist_func = lambda x, size, dtype:\
            x.binomial(n=1, p=0.1, size=size)
        for _size in test_sizes:
            for _bitgen in bitgen_types:
                with self.subTest(_size=_size, _bitgen=_bitgen):
                    self.check_numpy_parity(dist_func, _bitgen,
                                            None, _size, None,
                                            adjusted_ulp_prec)

        dist_func = lambda x, n, p, size:\
            x.binomial(n=n, p=p, size=size)
        self._check_invalid_types(dist_func, ['n', 'p', 'size'],
                                  [1, 0.75, (1,)], ['x', 'x', ('x',)])

    def test_binomial_cases(self):
        cases = [
            (1, 0.1), # p <= 0.5 && n * p <= 30
            (50, 0.9), # p > 0.5 && n * p <= 30
            (100, 0.4), # p <= 0.5 && n * p > 30
            (100, 0.9) # p > 0.5 && n * p > 30
        ]
        size = None

        for n, p in cases:
            with self.subTest(n=n, p=p):
                dist_func = lambda x, size, dtype:\
                    x.binomial(n, p, size=size)
                self.check_numpy_parity(dist_func, None,
                                        None, size, None, 0)


class TestGeneratorCaching(TestCase, SerialMixin):
    def test_randomgen_caching(self):
        nb_rng = np.random.default_rng(1)
        np_rng = np.random.default_rng(1)

        numba_func = numba.njit(lambda x: x.random(10), cache=True)
        self.assertPreciseEqual(np_rng.random(10), numba_func(nb_rng))
        # Run the function twice to make sure caching doesn't break anything.
        self.assertPreciseEqual(np_rng.random(10), numba_func(nb_rng))
        # Check that the function can be retrieved successfully from the cache.
        res = run_in_new_process_caching(test_generator_caching)
        self.assertEqual(res['exitcode'], 0)
