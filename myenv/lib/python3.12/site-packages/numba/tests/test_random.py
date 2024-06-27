import collections
import functools
import math
import multiprocessing
import os
import random
import subprocess
import sys
import threading
import itertools
from textwrap import dedent

import numpy as np

import unittest

import numba
from numba import jit, _helperlib, njit
from numba.core import types
from numba.tests.support import TestCase, compile_function, tag
from numba.core.errors import TypingError


# State size of the Mersenne Twister
N = 624


def get_py_state_ptr():
    return _helperlib.rnd_get_py_state_ptr()

def get_np_state_ptr():
    return _helperlib.rnd_get_np_state_ptr()


def numpy_randint1(a):
    return np.random.randint(a)

def numpy_randint2(a, b):
    return np.random.randint(a, b)

def random_randint(a, b):
    return random.randint(a, b)

def random_randrange1(a):
    return random.randrange(a)

def random_randrange2(a, b):
    return random.randrange(a, b)

def random_randrange3(a, b, c):
    return random.randrange(a, b, c)

def numpy_choice1(a):
    return np.random.choice(a)

def numpy_choice2(a, size):
    return np.random.choice(a, size=size)

def numpy_choice3(a, size, replace):
    return np.random.choice(a, size=size, replace=replace)

def numpy_multinomial2(n, pvals):
    return np.random.multinomial(n, pvals)

def numpy_multinomial3(n, pvals, size):
    return np.random.multinomial(n, pvals=pvals, size=size)

def numpy_dirichlet(alpha, size):
    return np.random.dirichlet(alpha, size=size)

def numpy_dirichlet_default(alpha):
    return np.random.dirichlet(alpha)

def numpy_noncentral_chisquare(df, nonc, size):
    return np.random.noncentral_chisquare(df, nonc, size=size)

def numpy_noncentral_chisquare_default(df, nonc):
    return np.random.noncentral_chisquare(df, nonc)

def numpy_check_rand(seed, a, b):
    np.random.seed(seed)
    expected = np.random.random((a, b))
    np.random.seed(seed)
    got = np.random.rand(a, b)
    return expected, got

def numpy_check_randn(seed, a, b):
    np.random.seed(seed)
    expected = np.random.standard_normal((a, b))
    np.random.seed(seed)
    got = np.random.randn(a, b)
    return expected, got

def jit_with_args(name, argstring):
    code = """def func(%(argstring)s):
        return %(name)s(%(argstring)s)
""" % locals()
    pyfunc = compile_function("func", code, globals())
    return jit(nopython=True)(pyfunc)

def jit_with_kwargs(name, kwarg_list):
    # Similar to jit_with_args, but uses keyword arguments
    call_args_with_kwargs = ','.join([f'{kw}={kw}' for kw in kwarg_list])
    signature = ','.join(kwarg_list)
    code = f"""def func({signature}):
        return {name}({call_args_with_kwargs})
"""
    pyfunc = compile_function("func", code, globals())
    return jit(nopython=True)(pyfunc)

def jit_nullary(name):
    return jit_with_args(name, "")

def jit_unary(name):
    return jit_with_args(name, "a")

def jit_binary(name):
    return jit_with_args(name, "a, b")

def jit_ternary(name):
    return jit_with_args(name, "a, b, c")


random_gauss = jit_binary("random.gauss")
random_random = jit_nullary("random.random")
random_seed = jit_unary("random.seed")

numpy_normal = jit_binary("np.random.normal")
numpy_random = jit_nullary("np.random.random")
numpy_seed = jit_unary("np.random.seed")


def _copy_py_state(r, ptr):
    """
    Copy state of Python random *r* to Numba state *ptr*.
    """
    mt = r.getstate()[1]
    ints, index = mt[:-1], mt[-1]
    _helperlib.rnd_set_state(ptr, (index, list(ints)))
    return ints, index

def _copy_np_state(r, ptr):
    """
    Copy state of Numpy random *r* to Numba state *ptr*.
    """
    ints, index = r.get_state()[1:3]
    _helperlib.rnd_set_state(ptr, (index, [int(x) for x in ints]))
    return ints, index

def sync_to_numpy(r):
    _ver, mt_st, _gauss_next = r.getstate()
    mt_pos = mt_st[-1]
    mt_ints = mt_st[:-1]
    assert len(mt_ints) == 624

    np_st = ('MT19937', np.array(mt_ints, dtype='uint32'), mt_pos)
    if _gauss_next is None:
        np_st += (0, 0.0)
    else:
        np_st += (1, _gauss_next)

    np.random.set_state(np_st)

# Pure Python equivalents of some of the Numpy distributions, using
# Python's basic generators.

def py_chisquare(r, df):
    return 2.0 * r.gammavariate(df / 2.0, 1.0)

def py_f(r, num, denom):
    return ((py_chisquare(r, num) * denom) /
            (py_chisquare(r, denom) * num))


class BaseTest(TestCase):

    def _follow_cpython(self, ptr, seed=2):
        r = random.Random(seed)
        _copy_py_state(r, ptr)
        return r

    def _follow_numpy(self, ptr, seed=2):
        r = np.random.RandomState(seed)
        _copy_np_state(r, ptr)
        return r


class TestInternals(BaseTest):
    """
    Test low-level internals of the implementation.
    """

    def _check_get_set_state(self, ptr):
        state = _helperlib.rnd_get_state(ptr)
        i, ints = state
        self.assertIsInstance(i, int)
        self.assertIsInstance(ints, list)
        self.assertEqual(len(ints), N)
        j = (i * 100007) % N
        ints = [i * 3 for i in range(N)]
        # Roundtrip
        _helperlib.rnd_set_state(ptr, (j, ints))
        self.assertEqual(_helperlib.rnd_get_state(ptr), (j, ints))

    def _check_shuffle(self, ptr):
        # We test shuffling against CPython
        r = random.Random()
        ints, index = _copy_py_state(r, ptr)
        # Force shuffling in CPython generator
        for i in range(index, N + 1, 2):
            r.random()
        _helperlib.rnd_shuffle(ptr)
        # Check new integer keys
        mt = r.getstate()[1]
        ints, index = mt[:-1], mt[-1]
        self.assertEqual(_helperlib.rnd_get_state(ptr)[1], list(ints))

    def _check_init(self, ptr):
        # We use the same integer seeding as Numpy
        # (CPython is different: it treats the integer as a byte array)
        r = np.random.RandomState()
        for i in [0, 1, 125, 2**32 - 5]:
            # Need to cast to a C-sized int (for Numpy <= 1.7)
            r.seed(np.uint32(i))
            st = r.get_state()
            ints = list(st[1])
            index = st[2]
            assert index == N  # sanity check
            _helperlib.rnd_seed(ptr, i)
            self.assertEqual(_helperlib.rnd_get_state(ptr), (index, ints))

    def _check_perturb(self, ptr):
        states = []
        for i in range(10):
            # Initialize with known state
            _helperlib.rnd_seed(ptr, 0)
            # Perturb with entropy
            _helperlib.rnd_seed(ptr, os.urandom(512))
            states.append(tuple(_helperlib.rnd_get_state(ptr)[1]))
        # No two identical states
        self.assertEqual(len(set(states)), len(states))

    def test_get_set_state(self):
        self._check_get_set_state(get_py_state_ptr())

    def test_shuffle(self):
        self._check_shuffle(get_py_state_ptr())

    def test_init(self):
        self._check_init(get_py_state_ptr())

    def test_perturb(self):
        self._check_perturb(get_py_state_ptr())


class TestRandom(BaseTest):

    # NOTE: there may be cascading imprecision issues (e.g. between x87-using
    # C code and SSE-using LLVM code), which is especially brutal for some
    # iterative algorithms with sensitive exit conditions.
    # Therefore we stick to hardcoded integers for seed values.

    def _check_random_seed(self, seedfunc, randomfunc):
        """
        Check seed()- and random()-like functions.
        """
        # Our seed() mimics NumPy's.
        r = np.random.RandomState()
        for i in [0, 1, 125, 2**32 - 1]:
            # Need to cast to a C-sized int (for Numpy <= 1.7)
            r.seed(np.uint32(i))
            seedfunc(i)
            # Be sure to trigger a reshuffle
            for j in range(N + 10):
                self.assertPreciseEqual(randomfunc(), r.uniform(0.0, 1.0))

    def test_random_random(self):
        self._check_random_seed(random_seed, random_random)

    def test_numpy_random(self):
        self._check_random_seed(numpy_seed, numpy_random)
        # Test aliases
        self._check_random_seed(numpy_seed, jit_nullary("np.random.random_sample"))
        self._check_random_seed(numpy_seed, jit_nullary("np.random.ranf"))
        self._check_random_seed(numpy_seed, jit_nullary("np.random.sample"))
        self._check_random_seed(numpy_seed, jit_nullary("np.random.rand"))

    def _check_random_sized(self, seedfunc, randomfunc):
        # Our seed() mimics NumPy's.
        r = np.random.RandomState()
        for i in [0, 1, 125, 2**32 - 1]:
            # Need to cast to a C-sized int (for Numpy <= 1.7)
            r.seed(np.uint32(i))
            seedfunc(i)
            for n in range(10):
                self.assertPreciseEqual(randomfunc(n), r.uniform(0.0, 1.0, n))

    def test_numpy_random_sized(self):
        self._check_random_sized(numpy_seed, jit_unary("np.random.random_sample"))
        self._check_random_sized(numpy_seed, jit_unary("np.random.ranf"))
        self._check_random_sized(numpy_seed, jit_unary("np.random.sample"))
        self._check_random_sized(numpy_seed, jit_unary("np.random.rand"))

    def test_independent_generators(self):
        # PRNGs for Numpy and Python are independent.
        N = 10
        random_seed(1)
        py_numbers = [random_random() for i in range(N)]
        numpy_seed(2)
        np_numbers = [numpy_random() for i in range(N)]
        random_seed(1)
        numpy_seed(2)
        pairs = [(random_random(), numpy_random()) for i in range(N)]
        self.assertPreciseEqual([p[0] for p in pairs], py_numbers)
        self.assertPreciseEqual([p[1] for p in pairs], np_numbers)

    def _check_getrandbits(self, func, ptr):
        """
        Check a getrandbits()-like function.
        """
        # Our implementation follows CPython's for bits <= 64.
        r = self._follow_cpython(ptr)
        for nbits in range(1, 65):
            expected = r.getrandbits(nbits)
            got = func(nbits)
            self.assertPreciseEqual(expected, got)
        self.assertRaises(OverflowError, func, 65)
        self.assertRaises(OverflowError, func, 9999999)
        self.assertRaises(OverflowError, func, -1)

    def test_random_getrandbits(self):
        self._check_getrandbits(jit_unary("random.getrandbits"), get_py_state_ptr())

    # Explanation for the large ulps value: on 32-bit platforms, our
    # LLVM-compiled functions use SSE but they are compared against
    # C functions which use x87.
    # On some distributions, the errors seem to accumulate dramatically.

    def _check_dist(self, func, pyfunc, argslist, niters=3,
                    prec='double', ulps=12, pydtype=None):
        assert len(argslist)
        for args in argslist:
            results = [func(*args) for i in range(niters)]
            pyresults = [(pyfunc(*args, dtype=pydtype) if pydtype else pyfunc(*args))
                         for i in range(niters)]
            self.assertPreciseEqual(results, pyresults, prec=prec, ulps=ulps,
                                    msg="for arguments %s" % (args,))

    def _check_dist_kwargs(self, func, pyfunc, kwargslist, niters=3,
                           prec='double', ulps=12, pydtype=None):
        assert len(kwargslist)
        for kwargs in kwargslist:
            results = [func(**kwargs) for i in range(niters)]
            pyresults = [(pyfunc(**kwargs, dtype=pydtype) if pydtype else pyfunc(**kwargs))
                         for i in range(niters)]
            self.assertPreciseEqual(results, pyresults, prec=prec, ulps=ulps,
                                    msg="for arguments %s" % (kwargs,))

    def _check_gauss(self, func2, func1, func0, ptr):
        """
        Check a gauss()-like function.
        """
        # Our implementation follows Numpy's.
        r = self._follow_numpy(ptr)
        if func2 is not None:
            self._check_dist(func2, r.normal,
                             [(1.0, 1.0), (2.0, 0.5), (-2.0, 0.5)],
                             niters=N // 2 + 10)
        if func1 is not None:
            self._check_dist(func1, r.normal, [(0.5,)])
        if func0 is not None:
            self._check_dist(func0, r.normal, [()])

    def test_random_gauss(self):
        self._check_gauss(jit_binary("random.gauss"), None, None, get_py_state_ptr())

    def test_random_normalvariate(self):
        # normalvariate() is really an alias to gauss() in Numba
        # (not in Python, though - they use different algorithms)
        self._check_gauss(jit_binary("random.normalvariate"), None, None,
                          get_py_state_ptr())

    def test_numpy_normal(self):
        self._check_gauss(jit_binary("np.random.normal"),
                          jit_unary("np.random.normal"),
                          jit_nullary("np.random.normal"),
                          get_np_state_ptr())

    def test_numpy_standard_normal(self):
        self._check_gauss(None, None, jit_nullary("np.random.standard_normal"),
                          get_np_state_ptr())

    def test_numpy_randn(self):
        self._check_gauss(None, None, jit_nullary("np.random.randn"),
                          get_np_state_ptr())

    def _check_lognormvariate(self, func2, func1, func0, ptr):
        """
        Check a lognormvariate()-like function.
        """
        # Our implementation follows Numpy's.
        r = self._follow_numpy(ptr)
        if func2 is not None:
            self._check_dist(func2, r.lognormal,
                             [(1.0, 1.0), (2.0, 0.5), (-2.0, 0.5)],
                             niters=N // 2 + 10)
        if func1 is not None:
            self._check_dist(func1, r.lognormal, [(0.5,)])
        if func0 is not None:
            self._check_dist(func0, r.lognormal, [()])

    def test_random_lognormvariate(self):
        self._check_lognormvariate(jit_binary("random.lognormvariate"),
                                   None, None, get_py_state_ptr())

    def test_numpy_lognormal(self):
        self._check_lognormvariate(jit_binary("np.random.lognormal"),
                                   jit_unary("np.random.lognormal"),
                                   jit_nullary("np.random.lognormal"),
                                   get_np_state_ptr())

    def _check_randrange(self, func1, func2, func3, ptr, max_width, is_numpy, tp=None):
        """
        Check a randrange()-like function.
        """
        # Sanity check
        ints = []
        for i in range(10):
            ints.append(func1(500000000))
            ints.append(func2(5, 500000000))
            if func3 is not None:
                ints.append(func3(5, 500000000, 3))
        if is_numpy:
            rr = self._follow_numpy(ptr).randint
        else:
            rr = self._follow_cpython(ptr).randrange
        widths = [w for w in [1, 5, 8, 5000, 2**40, 2**62 + 2**61] if w < max_width]
        pydtype = tp if is_numpy else None
        for width in widths:
            self._check_dist(func1, rr, [(width,)], niters=10,
                            pydtype=pydtype)
            self._check_dist(func2, rr, [(-2, 2 +width)], niters=10,
                            pydtype=pydtype)
            if func3 is not None:
                self.assertPreciseEqual(func3(-2, 2 + width, 6),
                                        rr(-2, 2 + width, 6))
                self.assertPreciseEqual(func3(2 + width, 2, -3),
                                        rr(2 + width, 2, -3))
        # Empty ranges
        self.assertRaises(ValueError, func1, 0)
        self.assertRaises(ValueError, func1, -5)
        self.assertRaises(ValueError, func2, 5, 5)
        self.assertRaises(ValueError, func2, 5, 2)
        if func3 is not None:
            self.assertRaises(ValueError, func3, 5, 7, -1)
            self.assertRaises(ValueError, func3, 7, 5, 1)

    def test_random_randrange(self):
        for tp, max_width in [(types.int64, 2**63), (types.int32, 2**31)]:
            cf1 = njit((tp,))(random_randrange1)
            cf2 = njit((tp, tp,),)(random_randrange2)
            cf3 = njit((tp, tp, tp,))(random_randrange3)
            self._check_randrange(cf1, cf2, cf3, get_py_state_ptr(),
                                  max_width, False)

    def test_numpy_randint(self):
        for tp, np_tp, max_width in [(types.int64, np.int64, 2**63),
                                     (types.int32, np.int32, 2**31)]:
            cf1 = njit((tp,))(numpy_randint1)
            cf2 = njit((tp, tp,))(numpy_randint2)
            self._check_randrange(cf1, cf2, None, get_np_state_ptr(), max_width,
                                  True, np_tp)

    def _check_randint(self, func, ptr, max_width):
        """
        Check a randint()-like function.
        """
        # Sanity check
        ints = []
        for i in range(10):
            ints.append(func(5, 500000000))
        self.assertEqual(len(ints), len(set(ints)), ints)

        r = self._follow_cpython(ptr)
        for args in [(1, 5), (13, 5000), (20, 2**62 + 2**61)]:
            if args[1] > max_width:
                continue
            self._check_dist(func, r.randint, [args], niters=10)
        # Empty ranges
        self.assertRaises(ValueError, func, 5, 4)
        self.assertRaises(ValueError, func, 5, 2)

    def test_random_randint(self):
        for tp, max_width in [(types.int64, 2**63), (types.int32, 2**31)]:
            cf = njit((tp, tp,))(random_randint)
            self._check_randint(cf, get_py_state_ptr(), max_width)

    def _check_uniform(self, func, ptr):
        """
        Check a uniform()-like function.
        """
        # Our implementation follows Python's.
        r = self._follow_cpython(ptr)
        self._check_dist(func, r.uniform,
                         [(1.5, 1e6), (-2.5, 1e3), (1.5, -2.5)])

    def _check_any_distrib_kwargs(self, func, ptr, distrib, paramlist):
        """
        Check any numpy distribution function. Does Numba use the same keyword
        argument names as Numpy?
        And given a fixed seed, do they both return the same samples?
        """
        # Our implementation follows Numpy's (not Python's)
        r = self._follow_numpy(ptr)
        distrib_method_of_numpy = getattr(r, distrib)
        self._check_dist_kwargs(func, distrib_method_of_numpy, paramlist)


    def test_random_uniform(self):
        self._check_uniform(jit_binary("random.uniform"), get_py_state_ptr())

    def test_numpy_uniform(self):
        self._check_uniform(jit_binary("np.random.uniform"), get_np_state_ptr())

    def test_numpy_uniform_kwargs(self):
        self._check_any_distrib_kwargs(
            jit_with_kwargs("np.random.uniform", ['low', 'high']),
            get_np_state_ptr(),
            'uniform',
            paramlist=[{'low': 1.5, 'high': 1e6},
                       {'low': -2.5, 'high': 1e3},
                       {'low': 1.5, 'high': -2.5}])

    def _check_triangular(self, func2, func3, ptr):
        """
        Check a triangular()-like function.
        """
        # Our implementation follows Python's.
        r = self._follow_cpython(ptr)
        if func2 is not None:
            self._check_dist(func2, r.triangular,
                             [(1.5, 3.5), (-2.5, 1.5), (1.5, 1.5)])
        self._check_dist(func3, r.triangular, [(1.5, 3.5, 2.2)])

    def test_random_triangular(self):
        self._check_triangular(jit_binary("random.triangular"),
                               jit_ternary("random.triangular"),
                               get_py_state_ptr())

    def test_numpy_triangular(self):
        triangular = jit_ternary("np.random.triangular")
        fixed_triangular = lambda l, r, m: triangular(l, m, r)
        self._check_triangular(None, fixed_triangular, get_np_state_ptr())

    def _check_gammavariate(self, func2, func1, ptr):
        """
        Check a gammavariate()-like function.
        """
        # Our implementation follows Python's.
        r = self._follow_cpython(ptr)
        if func2 is not None:
            self._check_dist(func2, r.gammavariate,
                             [(0.5, 2.5), (1.0, 1.5), (1.5, 3.5)])
        if func1 is not None:
            self.assertPreciseEqual(func1(1.5), r.gammavariate(1.5, 1.0))
        # Invalid inputs
        if func2 is not None:
            self.assertRaises(ValueError, func2, 0.0, 1.0)
            self.assertRaises(ValueError, func2, 1.0, 0.0)
            self.assertRaises(ValueError, func2, -0.5, 1.0)
            self.assertRaises(ValueError, func2, 1.0, -0.5)
        if func1 is not None:
            self.assertRaises(ValueError, func1, 0.0)
            self.assertRaises(ValueError, func1, -0.5)

    def test_random_gammavariate(self):
        self._check_gammavariate(jit_binary("random.gammavariate"), None,
                                 get_py_state_ptr())

    def test_numpy_gamma(self):
        self._check_gammavariate(jit_binary("np.random.gamma"),
                                 jit_unary("np.random.gamma"),
                                 get_np_state_ptr())
        self._check_gammavariate(None,
                                 jit_unary("np.random.standard_gamma"),
                                 get_np_state_ptr())

    def _check_betavariate(self, func, ptr):
        """
        Check a betavariate()-like function.
        """
        # Our implementation follows Python's.
        r = self._follow_cpython(ptr)
        self._check_dist(func, r.betavariate, [(0.5, 2.5)])
        # Invalid inputs
        self.assertRaises(ValueError, func, 0.0, 1.0)
        self.assertRaises(ValueError, func, 1.0, 0.0)
        self.assertRaises(ValueError, func, -0.5, 1.0)
        self.assertRaises(ValueError, func, 1.0, -0.5)

    def test_random_betavariate(self):
        self._check_betavariate(jit_binary("random.betavariate"), get_py_state_ptr())

    def test_numpy_beta(self):
        self._check_betavariate(jit_binary("np.random.beta"), get_np_state_ptr())

    def _check_vonmisesvariate(self, func, ptr):
        """
        Check a vonmisesvariate()-like function.
        """
        r = self._follow_cpython(ptr)
        self._check_dist(func, r.vonmisesvariate, [(0.5, 2.5)])

    def test_random_vonmisesvariate(self):
        self._check_vonmisesvariate(jit_binary("random.vonmisesvariate"),
                                    get_py_state_ptr())

    def test_numpy_vonmises(self):
        self._check_vonmisesvariate(jit_binary("np.random.vonmises"),
                                    get_np_state_ptr())

    def _check_expovariate(self, func, ptr):
        """
        Check a expovariate()-like function.  Note the second argument
        is inversed compared to np.random.exponential().
        """
        r = self._follow_numpy(ptr)
        for lambd in (0.2, 0.5, 1.5):
            for i in range(3):
                self.assertPreciseEqual(func(lambd), r.exponential(1 / lambd),
                                        prec='double')

    def test_random_expovariate(self):
        self._check_expovariate(jit_unary("random.expovariate"), get_py_state_ptr())

    def _check_exponential(self, func1, func0, ptr):
        """
        Check a exponential()-like function.
        """
        r = self._follow_numpy(ptr)
        if func1 is not None:
            self._check_dist(func1, r.exponential, [(0.5,), (1.0,), (1.5,)])
        if func0 is not None:
            self._check_dist(func0, r.exponential, [()])

    def test_numpy_exponential(self):
        self._check_exponential(jit_unary("np.random.exponential"),
                                jit_nullary("np.random.exponential"),
                                get_np_state_ptr())

    def test_numpy_standard_exponential(self):
        self._check_exponential(None,
                                jit_nullary("np.random.standard_exponential"),
                                get_np_state_ptr())

    def _check_paretovariate(self, func, ptr):
        """
        Check a paretovariate()-like function.
        """
        # Our implementation follows Python's.
        r = self._follow_cpython(ptr)
        self._check_dist(func, r.paretovariate, [(0.5,), (3.5,)])

    def test_random_paretovariate(self):
        self._check_paretovariate(jit_unary("random.paretovariate"), get_py_state_ptr())

    def test_numpy_pareto(self):
        pareto = jit_unary("np.random.pareto")
        fixed_pareto = lambda a: pareto(a) + 1.0
        self._check_paretovariate(fixed_pareto, get_np_state_ptr())

    def _check_weibullvariate(self, func2, func1, ptr):
        """
        Check a weibullvariate()-like function.
        """
        # Our implementation follows Python's.
        r = self._follow_cpython(ptr)
        if func2 is not None:
            self._check_dist(func2, r.weibullvariate, [(0.5, 2.5)])
        if func1 is not None:
            for i in range(3):
                self.assertPreciseEqual(func1(2.5),
                                        r.weibullvariate(1.0, 2.5))

    def test_random_weibullvariate(self):
        self._check_weibullvariate(jit_binary("random.weibullvariate"),
                                   None, get_py_state_ptr())

    def test_numpy_weibull(self):
        self._check_weibullvariate(None, jit_unary("np.random.weibull"),
                                   get_np_state_ptr())

    def test_numpy_binomial(self):
        # We follow Numpy's algorithm up to n*p == 30
        binomial = jit_binary("np.random.binomial")
        r = self._follow_numpy(get_np_state_ptr(), 0)
        self._check_dist(binomial, r.binomial, [(18, 0.25)])
        # Sanity check many values
        for n in (100, 1000, 10000):
            self.assertEqual(binomial(n, 0.0), 0)
            self.assertEqual(binomial(n, 1.0), n)
            for p in (0.0001, 0.1, 0.4, 0.49999, 0.5, 0.50001, 0.8, 0.9, 0.9999):
                r = binomial(n, p)
                if p > 0.5:
                    r = n - r
                    p = 1 - p
                self.assertGreaterEqual(r, 0)
                self.assertLessEqual(r, n)
                expected = p * n
                tol = 3 * n / math.sqrt(n)
                self.assertGreaterEqual(r, expected - tol, (p, n, r))
                self.assertLessEqual(r, expected + tol, (p, n, r))
        # Invalid values
        self.assertRaises(ValueError, binomial, -1, 0.5)
        self.assertRaises(ValueError, binomial, 10, -0.1)
        self.assertRaises(ValueError, binomial, 10, 1.1)

    def test_numpy_chisquare(self):
        chisquare = jit_unary("np.random.chisquare")
        r = self._follow_cpython(get_np_state_ptr())
        self._check_dist(chisquare,
                         functools.partial(py_chisquare, r),
                         [(1.5,), (2.5,)])

    def test_numpy_f(self):
        f = jit_binary("np.random.f")
        r = self._follow_cpython(get_np_state_ptr())
        self._check_dist(f, functools.partial(py_f, r),
                         [(0.5, 1.5), (1.5, 0.8)])

    def test_numpy_geometric(self):
        geom = jit_unary("np.random.geometric")
        # p out of domain
        self.assertRaises(ValueError, geom, -1.0)
        self.assertRaises(ValueError, geom, 0.0)
        self.assertRaises(ValueError, geom, 1.001)
        # Some basic checks
        N = 200
        r = [geom(1.0) for i in range(N)]
        self.assertPreciseEqual(r, [1] * N)
        r = [geom(0.9) for i in range(N)]
        n = r.count(1)
        self.assertGreaterEqual(n, N // 2)
        self.assertLess(n, N)
        self.assertFalse([i for i in r if i > 1000])  # unlikely
        r = [geom(0.4) for i in range(N)]
        self.assertTrue([i for i in r if i > 4])  # likely
        r = [geom(0.01) for i in range(N)]
        self.assertTrue([i for i in r if i > 50])  # likely
        r = [geom(1e-15) for i in range(N)]
        self.assertTrue([i for i in r if i > 2**32])  # likely

    def test_numpy_gumbel(self):
        gumbel = jit_binary("np.random.gumbel")
        r = self._follow_numpy(get_np_state_ptr())
        self._check_dist(gumbel, r.gumbel, [(0.0, 1.0), (-1.5, 3.5)])

    def test_numpy_gumbel_kwargs(self):
        self._check_any_distrib_kwargs(
            jit_with_kwargs("np.random.gumbel", ['loc', 'scale']),
            get_np_state_ptr(),
            distrib="gumbel",
            paramlist=[{'loc': 0.0, 'scale': 1.0},
                       {'loc': -1.5, 'scale': 3.5}])


    def test_numpy_hypergeometric(self):
        # Our implementation follows Numpy's up to nsamples = 10.
        hg = jit_ternary("np.random.hypergeometric")
        r = self._follow_numpy(get_np_state_ptr())
        self._check_dist(hg, r.hypergeometric,
                         [(1000, 5000, 10), (5000, 1000, 10)],
                         niters=30)
        # Sanity checks
        r = [hg(1000, 1000, 100) for i in range(100)]
        self.assertTrue(all(x >= 0 and x <= 100 for x in r), r)
        self.assertGreaterEqual(np.mean(r), 40.0)
        self.assertLessEqual(np.mean(r), 60.0)
        r = [hg(1000, 100000, 100) for i in range(100)]
        self.assertTrue(all(x >= 0 and x <= 100 for x in r), r)
        self.assertLessEqual(np.mean(r), 10.0)
        r = [hg(100000, 1000, 100) for i in range(100)]
        self.assertTrue(all(x >= 0 and x <= 100 for x in r), r)
        self.assertGreaterEqual(np.mean(r), 90.0)

    def test_numpy_laplace(self):
        r = self._follow_numpy(get_np_state_ptr())
        self._check_dist(jit_binary("np.random.laplace"), r.laplace,
                         [(0.0, 1.0), (-1.5, 3.5)])
        self._check_dist(jit_unary("np.random.laplace"), r.laplace,
                         [(0.0,), (-1.5,)])
        self._check_dist(jit_nullary("np.random.laplace"), r.laplace, [()])

    def test_numpy_logistic(self):
        r = self._follow_numpy(get_np_state_ptr())
        self._check_dist(jit_binary("np.random.logistic"), r.logistic,
                         [(0.0, 1.0), (-1.5, 3.5)])
        self._check_dist(jit_unary("np.random.logistic"), r.logistic,
                         [(0.0,), (-1.5,)])
        self._check_dist(jit_nullary("np.random.logistic"), r.logistic, [()])

    def test_numpy_logseries(self):
        r = self._follow_numpy(get_np_state_ptr())
        logseries = jit_unary("np.random.logseries")
        self._check_dist(logseries, r.logseries,
                         [(0.1,), (0.99,), (0.9999,)],
                         niters=50)
        # Numpy's logseries overflows on 32-bit builds, so instead
        # hardcode Numpy's (correct) output on 64-bit builds.
        r = self._follow_numpy(get_np_state_ptr(), seed=1)
        self.assertEqual([logseries(0.9999999999999) for i in range(10)],
                         [2022733531, 77296, 30, 52204, 9341294, 703057324,
                          413147702918, 1870715907, 16009330, 738])
        self.assertRaises(ValueError, logseries, 0.0)
        self.assertRaises(ValueError, logseries, -0.1)
        self.assertRaises(ValueError, logseries, 1.1)

    def test_numpy_poisson(self):
        r = self._follow_numpy(get_np_state_ptr())
        poisson = jit_unary("np.random.poisson")
        # Our implementation follows Numpy's.
        self._check_dist(poisson, r.poisson,
                         [(0.0,), (0.5,), (2.0,), (10.0,), (900.5,)],
                         niters=50)
        self.assertRaises(ValueError, poisson, -0.1)

    def test_numpy_negative_binomial(self):
        self._follow_numpy(get_np_state_ptr(), 0)
        negbin = jit_binary("np.random.negative_binomial")
        self.assertEqual([negbin(10, 0.9) for i in range(10)],
                         [2, 3, 1, 5, 2, 1, 0, 1, 0, 0])
        self.assertEqual([negbin(10, 0.1) for i in range(10)],
                         [55, 71, 56, 57, 56, 56, 34, 55, 101, 67])
        self.assertEqual([negbin(1000, 0.1) for i in range(10)],
                         [9203, 8640, 9081, 9292, 8938,
                          9165, 9149, 8774, 8886, 9117])
        m = np.mean([negbin(1000000000, 0.1)
                     for i in range(50)])
        self.assertGreater(m, 9e9 * 0.99)
        self.assertLess(m, 9e9 * 1.01)
        self.assertRaises(ValueError, negbin, 0, 0.5)
        self.assertRaises(ValueError, negbin, -1, 0.5)
        self.assertRaises(ValueError, negbin, 10, -0.1)
        self.assertRaises(ValueError, negbin, 10, 1.1)

    def test_numpy_power(self):
        r = self._follow_numpy(get_np_state_ptr())
        power = jit_unary("np.random.power")
        self._check_dist(power, r.power,
                         [(0.1,), (0.5,), (0.9,), (6.0,)])
        self.assertRaises(ValueError, power, 0.0)
        self.assertRaises(ValueError, power, -0.1)

    def test_numpy_rayleigh(self):
        r = self._follow_numpy(get_np_state_ptr())
        rayleigh1 = jit_unary("np.random.rayleigh")
        rayleigh0 = jit_nullary("np.random.rayleigh")
        self._check_dist(rayleigh1, r.rayleigh,
                         [(0.1,), (0.8,), (25.,), (1e3,)])
        self._check_dist(rayleigh0, r.rayleigh, [()])
        self.assertRaises(ValueError, rayleigh1, 0.0)
        self.assertRaises(ValueError, rayleigh1, -0.1)

    def test_numpy_standard_cauchy(self):
        r = self._follow_numpy(get_np_state_ptr())
        cauchy = jit_nullary("np.random.standard_cauchy")
        self._check_dist(cauchy, r.standard_cauchy, [()])

    def test_numpy_standard_t(self):
        # We use CPython's algorithm for the gamma dist and numpy's
        # for the normal dist.  Standard T calls both so we can't check
        # against either generator's output.
        r = self._follow_cpython(get_np_state_ptr())
        standard_t = jit_unary("np.random.standard_t")
        avg = np.mean([standard_t(5) for i in range(5000)])
        # Sanity check
        self.assertLess(abs(avg), 0.5)

    def test_numpy_wald(self):
        r = self._follow_numpy(get_np_state_ptr())
        wald = jit_binary("np.random.wald")
        self._check_dist(wald, r.wald, [(1.0, 1.0), (2.0, 5.0)])
        self.assertRaises(ValueError, wald, 0.0, 1.0)
        self.assertRaises(ValueError, wald, -0.1, 1.0)
        self.assertRaises(ValueError, wald, 1.0, 0.0)
        self.assertRaises(ValueError, wald, 1.0, -0.1)

    def test_numpy_wald_kwargs(self):
        numba_version = jit_with_kwargs("np.random.wald", ['mean', 'scale'])
        self._check_any_distrib_kwargs(numba_version,
                                       get_np_state_ptr(),
                                       distrib="wald",
                                       paramlist=[{'mean': 1.0, 'scale': 1.0},
                                                  {'mean': 2.0, 'scale': 5.0}])
        self.assertRaises(ValueError, numba_version, 0.0, 1.0)
        self.assertRaises(ValueError, numba_version, -0.1, 1.0)
        self.assertRaises(ValueError, numba_version, 1.0, 0.0)
        self.assertRaises(ValueError, numba_version, 1.0, -0.1)

    def test_numpy_zipf(self):
        r = self._follow_numpy(get_np_state_ptr())
        zipf = jit_unary("np.random.zipf")
        self._check_dist(zipf, r.zipf, [(1.5,), (2.5,)], niters=100)
        for val in (1.0, 0.5, 0.0, -0.1):
            self.assertRaises(ValueError, zipf, val)

    def _check_shuffle(self, func, ptr, is_numpy):
        """
        Check a shuffle()-like function for arrays.
        """
        arrs = [np.arange(20), np.arange(32).reshape((8, 4))]
        if is_numpy:
            r = self._follow_numpy(ptr)
        else:
            r = self._follow_cpython(ptr)
        for a in arrs:
            for i in range(3):
                got = a.copy()
                expected = a.copy()
                func(got)
                if is_numpy or len(a.shape) == 1:
                    r.shuffle(expected)
                    self.assertPreciseEqual(got, expected)
        # Test with an arbitrary buffer-providing object
        a = arrs[0]
        b = a.copy()
        func(memoryview(b))
        self.assertNotEqual(list(a), list(b))
        self.assertEqual(sorted(a), sorted(b))
        # Read-only object
        with self.assertTypingError():
            func(memoryview(b"xyz"))

    def test_random_shuffle(self):
        self._check_shuffle(jit_unary("random.shuffle"), get_py_state_ptr(), False)

    def test_numpy_shuffle(self):
        self._check_shuffle(jit_unary("np.random.shuffle"), get_np_state_ptr(), True)

    def _check_startup_randomness(self, func_name, func_args):
        """
        Check that the state is properly randomized at startup.
        """
        code = """if 1:
            from numba.tests import test_random
            func = getattr(test_random, %(func_name)r)
            print(func(*%(func_args)r))
            """ % (locals())
        numbers = set()
        for i in range(3):
            popen = subprocess.Popen([sys.executable, "-c", code],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = popen.communicate()
            if popen.returncode != 0:
                raise AssertionError("process failed with code %s: stderr follows\n%s\n"
                                     % (popen.returncode, err.decode()))
            numbers.add(float(out.strip()))
        self.assertEqual(len(numbers), 3, numbers)

    def test_random_random_startup(self):
        self._check_startup_randomness("random_random", ())

    def test_random_gauss_startup(self):
        self._check_startup_randomness("random_gauss", (1.0, 1.0))

    def test_numpy_random_startup(self):
        self._check_startup_randomness("numpy_random", ())

    def test_numpy_gauss_startup(self):
        self._check_startup_randomness("numpy_normal", (1.0, 1.0))

    def test_numpy_random_permutation(self):
        func = jit_unary("np.random.permutation")
        r = self._follow_numpy(get_np_state_ptr())
        for s in [5, 10, 15, 20]:
            a = np.arange(s)
            b = a.copy()
            # Test array version
            self.assertPreciseEqual(func(a), r.permutation(a))
            # Test int version
            self.assertPreciseEqual(func(s), r.permutation(s))
            # Permutation should not modify its argument
            self.assertPreciseEqual(a, b)
        # Check multi-dimensional arrays
        arrs = [np.arange(10).reshape(2, 5),
                np.arange(27).reshape(3, 3, 3),
                np.arange(36).reshape(2, 3, 3, 2)]
        for a in arrs:
            b = a.copy()
            self.assertPreciseEqual(func(a), r.permutation(a))
            self.assertPreciseEqual(a, b)


class TestRandomArrays(BaseTest):
    """
    Test array-producing variants of np.random.* functions.
    """

    def _compile_array_dist(self, funcname, nargs):
        qualname = "np.random.%s" % (funcname,)
        argstring = ', '.join('abcd'[:nargs])
        return jit_with_args(qualname, argstring)

    def _check_array_dist(self, funcname, scalar_args):
        """
        Check returning an array according to a given distribution.
        """
        cfunc = self._compile_array_dist(funcname, len(scalar_args) + 1)
        r = self._follow_numpy(get_np_state_ptr())
        pyfunc = getattr(r, funcname)
        for size in (8, (2, 3)):
            args = scalar_args + (size,)
            expected = pyfunc(*args)
            got = cfunc(*args)
            # Numpy may return int32s where we return int64s, adjust
            if (expected.dtype == np.dtype('int32')
                and got.dtype == np.dtype('int64')):
                expected = expected.astype(got.dtype)
            self.assertPreciseEqual(expected, got, prec='double', ulps=5)
        args = scalar_args + (None,)
        expected = pyfunc(*args)
        got = cfunc(*args)
        self.assertPreciseEqual(expected, got, prec='double', ulps=5)

    def _check_array_dist_gamma(self, funcname,  scalar_args, extra_pyfunc_args):
        """
        Check returning an array according to a given gamma distribution,
        where we use CPython's implementation rather than NumPy's.
        """
        cfunc = self._compile_array_dist(funcname, len(scalar_args) + 1)
        r = self._follow_cpython(get_np_state_ptr())
        pyfunc = getattr(r, "gammavariate")
        pyfunc_args = scalar_args + extra_pyfunc_args
        pyrandom = lambda *_args: pyfunc(*pyfunc_args)

        args = scalar_args + (None,)
        expected = pyrandom()
        got = cfunc(*args)
        self.assertPreciseEqual(expected, got, prec='double', ulps=5)
        for size in (8, (2, 3)):
            args = scalar_args + (size,)
            expected = np.empty(size)
            expected_flat = expected.flat
            for idx in range(expected.size):
                expected_flat[idx] = pyrandom()
            got = cfunc(*args)
            self.assertPreciseEqual(expected, got, prec='double', ulps=5)

    def _check_array_dist_self(self, funcname, scalar_args):
        """
        Check function returning an array against its scalar implementation.
        Because we use the CPython gamma distribution rather than the NumPy one,
        distributions which use the gamma distribution vary in ways that are
        difficult to compare. Instead, we compile both the array and scalar
        versions and check that the array is filled with the same values as
        we would expect from the scalar version.
        """
        @numba.njit
        def reset():
            np.random.seed(1234)

        array_func = self._compile_array_dist(funcname, len(scalar_args) + 1)

        qualname = "np.random.%s" % (funcname,)
        argstring = ', '.join('abcd'[:len(scalar_args)])
        scalar_func = jit_with_args(qualname, argstring)

        for size in (8, (2, 3)):
            args = scalar_args + (size,)
            reset()
            got = array_func(*args)
            reset()
            # We're just going to go with whatever type the array version
            # gives us and hope it's not Boolean or something useless.
            expected = np.empty(size, dtype=got.dtype)
            flat = expected.flat
            for idx in range(expected.size):
                flat[idx] = scalar_func(*scalar_args)
            self.assertPreciseEqual(expected, got, prec='double', ulps=5)

        reset()
        args = scalar_args + (None,)
        expected = scalar_func(*scalar_args)
        reset()
        got = array_func(*args)
        self.assertPreciseEqual(expected, got, prec='double', ulps=5)

    def test_numpy_randint(self):
        cfunc = self._compile_array_dist("randint", 3)
        low, high = 1000, 10000
        size = (30, 30)
        res = cfunc(low, high, size)
        self.assertIsInstance(res, np.ndarray)
        self.assertEqual(res.shape, size)
        self.assertIn(res.dtype, (np.dtype('int32'), np.dtype('int64')))
        self.assertTrue(np.all(res >= low))
        self.assertTrue(np.all(res < high))
        # Crude statistical tests
        mean = (low + high) / 2
        tol = (high - low) / 20
        self.assertGreaterEqual(res.mean(), mean - tol)
        self.assertLessEqual(res.mean(), mean + tol)

    def test_numpy_random_random(self):
        cfunc = self._compile_array_dist("random", 1)
        size = (30, 30)
        res = cfunc(size)
        self.assertIsInstance(res, np.ndarray)
        self.assertEqual(res.shape, size)
        self.assertEqual(res.dtype, np.dtype('float64'))
        # Results are within expected bounds
        self.assertTrue(np.all(res >= 0.0))
        self.assertTrue(np.all(res < 1.0))
        # Crude statistical tests
        self.assertTrue(np.any(res <= 0.1))
        self.assertTrue(np.any(res >= 0.9))
        mean = res.mean()
        self.assertGreaterEqual(mean, 0.45)
        self.assertLessEqual(mean, 0.55)

    # Sanity-check various distributions.  For convenience, we only check
    # those distributions that produce the exact same values as Numpy's.

    def test_numpy_beta(self):
        self._check_array_dist_self("beta", (0.5, 2.5))

    def test_numpy_binomial(self):
        self._check_array_dist("binomial", (20, 0.5))

    def test_numpy_chisquare(self):
        self._check_array_dist_self("chisquare", (1.5,))

    def test_numpy_exponential(self):
        self._check_array_dist("exponential", (1.5,))

    def test_numpy_f(self):
        self._check_array_dist_self("f", (0.5, 1.5))

    def test_numpy_gamma(self):
        self._check_array_dist_gamma("gamma", (2.0, 1.0), ())

    def test_numpy_geometric(self):
        self._check_array_dist("geometric", (1.0,))

    def test_numpy_gumbel(self):
        self._check_array_dist("gumbel", (1.5, 0.5))

    def test_numpy_hypergeometric(self):
        self._check_array_dist("hypergeometric", (1000, 5000, 10))

    def test_numpy_laplace(self):
        self._check_array_dist("laplace", (1.5, 0.5))

    def test_numpy_logistic(self):
        self._check_array_dist("logistic", (1.5, 0.5))

    def test_numpy_lognormal(self):
        self._check_array_dist("lognormal", (1.5, 2.0))

    def test_numpy_logseries(self):
        self._check_array_dist("logseries", (0.8,))

    def test_numpy_normal(self):
        self._check_array_dist("normal", (0.5, 2.0))

    def test_numpy_pareto(self):
        self._check_array_dist("pareto", (0.5,))

    def test_numpy_poisson(self):
        self._check_array_dist("poisson", (0.8,))

    def test_numpy_power(self):
        self._check_array_dist("power", (0.8,))

    def test_numpy_rand(self):
        cfunc = jit(nopython=True)(numpy_check_rand)
        expected, got = cfunc(42, 2, 3)
        self.assertEqual(got.shape, (2, 3))
        self.assertPreciseEqual(expected, got)

    def test_numpy_randn(self):
        cfunc = jit(nopython=True)(numpy_check_randn)
        expected, got = cfunc(42, 2, 3)
        self.assertEqual(got.shape, (2, 3))
        self.assertPreciseEqual(expected, got)

    def test_numpy_rayleigh(self):
        self._check_array_dist("rayleigh", (0.8,))

    def test_numpy_standard_cauchy(self):
        self._check_array_dist("standard_cauchy", ())

    def test_numpy_standard_exponential(self):
        self._check_array_dist("standard_exponential", ())

    def test_numpy_standard_gamma(self):
        self._check_array_dist_gamma("standard_gamma", (2.0,), (1.0,))

    def test_numpy_standard_normal(self):
        self._check_array_dist("standard_normal", ())

    def test_numpy_triangular(self):
        self._check_array_dist("triangular", (1.5, 2.2, 3.5))

    def test_numpy_uniform(self):
        self._check_array_dist("uniform", (0.1, 0.4))

    def test_numpy_wald(self):
        self._check_array_dist("wald", (0.1, 0.4))

    def test_numpy_vonmises(self):
        self._check_array_dist_self("vonmises", (0.5, 2.5))

    def test_numpy_zipf(self):
        self._check_array_dist("zipf", (2.5,))


class TestRandomChoice(BaseTest):
    """
    Test np.random.choice.
    """

    def _check_results(self, pop, res, replace=True):
        """
        Check basic expectations about a batch of samples.
        """
        spop = set(pop)
        sres = set(res)
        # All results are in the population
        self.assertLessEqual(sres, spop)
        # Sorted results are unlikely
        self.assertNotEqual(sorted(res), list(res))
        if replace:
            # Duplicates are likely
            self.assertLess(len(sres), len(res), res)
        else:
            # No duplicates
            self.assertEqual(len(sres), len(res), res)

    def _check_dist(self, pop, samples):
        """
        Check distribution of some samples.
        """
        # Sanity check that we have enough samples
        self.assertGreaterEqual(len(samples), len(pop) * 100)
        # Check equidistribution of samples
        expected_frequency = len(samples) / len(pop)
        c = collections.Counter(samples)
        for value in pop:
            n = c[value]
            self.assertGreaterEqual(n, expected_frequency * 0.5)
            self.assertLessEqual(n, expected_frequency * 2.0)

    def _accumulate_array_results(self, func, nresults):
        """
        Accumulate array results produced by *func* until they reach
        *nresults* elements.
        """
        res = []
        while len(res) < nresults:
            res += list(func().flat)
        return res[:nresults]

    def _check_choice_1(self, a, pop):
        """
        Check choice(a) against pop.
        """
        cfunc = jit(nopython=True)(numpy_choice1)
        n = len(pop)
        res = [cfunc(a) for i in range(n)]
        self._check_results(pop, res)
        dist = [cfunc(a) for i in range(n * 100)]
        self._check_dist(pop, dist)

    def test_choice_scalar_1(self):
        """
        Test choice(int)
        """
        n = 50
        pop = list(range(n))
        self._check_choice_1(n, pop)

    def test_choice_array_1(self):
        """
        Test choice(array)
        """
        pop = np.arange(50) * 2 + 100
        self._check_choice_1(pop, pop)

    def _check_array_results(self, func, pop, replace=True):
        """
        Check array results produced by *func* and their distribution.
        """
        n = len(pop)
        res = list(func().flat)
        self._check_results(pop, res, replace)
        dist = self._accumulate_array_results(func, n * 100)
        self._check_dist(pop, dist)

    def _check_choice_2(self, a, pop):
        """
        Check choice(a, size) against pop.
        """
        cfunc = jit(nopython=True)(numpy_choice2)
        n = len(pop)
        # Final sizes should be large enough, so as to stress
        # replacement
        sizes = [n - 10, (3, (n - 1) // 3), n * 10]

        for size in sizes:
            # Check result shape
            res = cfunc(a, size)
            expected_shape = size if isinstance(size, tuple) else (size,)
            self.assertEqual(res.shape, expected_shape)
            # Check results and their distribution
            self._check_array_results(lambda: cfunc(a, size), pop)

    def test_choice_scalar_2(self):
        """
        Test choice(int, size)
        """
        n = 50
        pop = np.arange(n)
        self._check_choice_2(n, pop)

    def test_choice_array_2(self):
        """
        Test choice(array, size)
        """
        pop = np.arange(50) * 2 + 100
        self._check_choice_2(pop, pop)

    def _check_choice_3(self, a, pop):
        """
        Check choice(a, size, replace) against pop.
        """
        cfunc = jit(nopython=True)(numpy_choice3)
        n = len(pop)
        # Final sizes should be close but slightly <= n, so as to stress
        # replacement (or not)
        sizes = [n - 10, (3, (n - 1) // 3)]
        replaces = [True, False]

        # Check result shapes
        for size in sizes:
            for replace in [True, False]:
                res = cfunc(a, size, replace)
                expected_shape = size if isinstance(size, tuple) else (size,)
                self.assertEqual(res.shape, expected_shape)

        # Check results for replace=True
        for size in sizes:
            self._check_array_results(lambda: cfunc(a, size, True), pop)
        # Check results for replace=False
        for size in sizes:
            self._check_array_results(lambda: cfunc(a, size, False), pop, False)

        # Can't ask for more samples than population size with replace=False
        for size in [n + 1, (3, n // 3 + 1)]:
            with self.assertRaises(ValueError):
                cfunc(a, size, False)

    def test_choice_scalar_3(self):
        """
        Test choice(int, size, replace)
        """
        n = 50
        pop = np.arange(n)
        self._check_choice_3(n, pop)

    def test_choice_array_3(self):
        """
        Test choice(array, size, replace)
        """
        pop = np.arange(50) * 2 + 100
        self._check_choice_3(pop, pop)

    def test_choice_follows_seed(self):
        # See issue #3888, np.random.choice must acknowledge the seed

        @jit(nopython=True)
        def numba_rands(n_to_return, choice_array):
            np.random.seed(1337)
            out = np.empty((n_to_return, 2), np.int32)
            for i in range(n_to_return):
                out[i] = np.random.choice(choice_array, 2, False)
            return out

        choice_array = np.random.randint(300, size=1000).astype(np.int32)
        tmp_np = choice_array.copy()
        expected = numba_rands.py_func(5, tmp_np)
        tmp_nb = choice_array.copy()
        got = numba_rands(5, tmp_nb)
        np.testing.assert_allclose(expected, got)
        # check no mutation
        np.testing.assert_allclose(choice_array, tmp_np)
        np.testing.assert_allclose(choice_array, tmp_nb)


class TestRandomMultinomial(BaseTest):
    """
    Test np.random.multinomial.
    """
    # A biased dice
    pvals = np.array([1, 1, 1, 2, 3, 1], dtype=np.float64)
    pvals /= pvals.sum()

    def _check_sample(self, n, pvals, sample):
        """
        Check distribution of some samples.
        """
        self.assertIsInstance(sample, np.ndarray)
        self.assertEqual(sample.shape, (len(pvals),))
        self.assertIn(sample.dtype, (np.dtype('int32'), np.dtype('int64')))
        # Statistical properties
        self.assertEqual(sample.sum(), n)
        for p, nexp in zip(pvals, sample):
            self.assertGreaterEqual(nexp, 0)
            self.assertLessEqual(nexp, n)
            pexp = float(nexp) / n
            self.assertGreaterEqual(pexp, p * 0.5)
            self.assertLessEqual(pexp, p * 2.0)

    def test_multinomial_2(self):
        """
        Test multinomial(n, pvals)
        """
        cfunc = jit(nopython=True)(numpy_multinomial2)
        n, pvals = 1000, self.pvals
        res = cfunc(n, pvals)
        self._check_sample(n, pvals, res)
        # pvals as list
        pvals = list(pvals)
        res = cfunc(n, pvals)
        self._check_sample(n, pvals, res)
        # A case with extreme probabilities
        n = 1000000
        pvals = np.array([1, 0, n // 100, 1], dtype=np.float64)
        pvals /= pvals.sum()
        res = cfunc(n, pvals)
        self._check_sample(n, pvals, res)

    def test_multinomial_3_int(self):
        """
        Test multinomial(n, pvals, size: int)
        """
        cfunc = jit(nopython=True)(numpy_multinomial3)
        n, pvals = 1000, self.pvals
        k = 10
        res = cfunc(n, pvals, k)
        self.assertEqual(res.shape[0], k)
        for sample in res:
            self._check_sample(n, pvals, sample)

    def test_multinomial_3_tuple(self):
        """
        Test multinomial(n, pvals, size: tuple)
        """
        cfunc = jit(nopython=True)(numpy_multinomial3)
        n, pvals = 1000, self.pvals
        k = (3, 4)
        res = cfunc(n, pvals, k)
        self.assertEqual(res.shape[:-1], k)
        for sample in res.reshape((-1, res.shape[-1])):
            self._check_sample(n, pvals, sample)


class TestRandomDirichlet(BaseTest):
    alpha = np.array([1, 1, 1, 2], dtype=np.float64)

    def _check_sample(self, alpha, size, sample):

        """Check output structure"""
        self.assertIsInstance(sample, np.ndarray)
        self.assertEqual(sample.dtype, np.float64)
        if size is None:
            self.assertEqual(sample.size, len(alpha))
        elif type(size) is int:
            self.assertEqual(sample.shape, (size, len(alpha)))
        else:
            self.assertEqual(sample.shape, size + (len(alpha),))

        """Check statistical properties"""
        for val in np.nditer(sample):
            self.assertGreaterEqual(val, 0)
            self.assertLessEqual(val, 1)
        if size is None:
            self.assertAlmostEqual(sample.sum(), 1, places=5)
        else:
            for totals in np.nditer(sample.sum(axis=-1)):
                self.assertAlmostEqual(totals, 1, places=5)

    def test_dirichlet_default(self):
        """
        Test dirichlet(alpha, size=None)
        """
        cfunc = jit(nopython=True)(numpy_dirichlet_default)
        alphas = (
            self.alpha,
            tuple(self.alpha),
            np.array([1, 1, 10000, 1], dtype=np.float64),
            np.array([1, 1, 1.5, 1], dtype=np.float64),
        )
        for alpha in alphas:
            res = cfunc(alpha)
            self._check_sample(alpha, None, res)

    def test_dirichlet(self):
        """
        Test dirichlet(alpha, size=None)
        """
        cfunc = jit(nopython=True)(numpy_dirichlet)
        sizes = (None, (10,), (10, 10))
        alphas = (
            self.alpha,
            tuple(self.alpha),
            np.array([1, 1, 10000, 1], dtype=np.float64),
            np.array([1, 1, 1.5, 1], dtype=np.float64),
        )

        for alpha, size in itertools.product(alphas, sizes):
            res = cfunc(alpha, size)
            self._check_sample(alpha, size, res)

    def test_dirichlet_exceptions(self):
        cfunc = jit(nopython=True)(numpy_dirichlet)
        alpha = tuple((0, 1, 1))
        with self.assertRaises(ValueError) as raises:
            cfunc(alpha, 1)
        self.assertIn("dirichlet: alpha must be > 0.0", str(raises.exception))
        
        alpha = self.alpha
        sizes = (True, 3j, 1.5, (1.5, 1), (3j, 1), (3j, 3j), (np.int8(3), np.int64(7)))
        for size in sizes:
            with self.assertRaises(TypingError) as raises:
                cfunc(alpha, size)
            self.assertIn(
                "np.random.dirichlet(): size should be int or "
                "tuple of ints or None, got",
                str(raises.exception),
            )

class TestRandomNoncentralChiSquare(BaseTest):

    def _check_sample(self, size, sample):

        # Check output structure
        if size is not None:
            self.assertIsInstance(sample, np.ndarray)
            self.assertEqual(sample.dtype, np.float64)
            
            if isinstance(size, int):
                self.assertEqual(sample.shape, (size,))
            else:
                self.assertEqual(sample.shape, size)
        else:
             self.assertIsInstance(sample, float)

        # Check statistical properties
        for val in np.nditer(sample):
            self.assertGreaterEqual(val, 0)

    def test_noncentral_chisquare_default(self):
        """
        Test noncentral_chisquare(df, nonc, size=None)
        """
        cfunc = jit(nopython=True)(numpy_noncentral_chisquare_default)
        inputs = (
            (0.5, 1), # test branch when df < 1
            (1, 5),
            (5, 1),
            (100000, 1),
            (1, 10000),
        )
        for df, nonc in inputs:
            res = cfunc(df, nonc)
            self._check_sample(None, res)
            res = cfunc(df, np.nan) # test branch when nonc is nan
            self.assertTrue(np.isnan(res))


    def test_noncentral_chisquare(self):
        """
        Test noncentral_chisquare(df, nonc, size)
        """
        cfunc = jit(nopython=True)(numpy_noncentral_chisquare)
        sizes = (None, 10, (10,), (10, 10))
        inputs = (
            (0.5, 1),
            (1, 5),
            (5, 1),
            (100000, 1),
            (1, 10000),
        )

        for (df, nonc), size in itertools.product(inputs, sizes):
            res = cfunc(df, nonc, size)
            self._check_sample(size, res)
            res = cfunc(df, np.nan, size) # test branch when nonc is nan
            self.assertTrue(np.isnan(res).all())

    def test_noncentral_chisquare_exceptions(self):
        cfunc = jit(nopython=True)(numpy_noncentral_chisquare)
        df, nonc = 0, 1
        with self.assertRaises(ValueError) as raises:
            cfunc(df, nonc, 1)
        self.assertIn("df <= 0", str(raises.exception))
        
        df, nonc = 1, -1
        with self.assertRaises(ValueError) as raises:
            cfunc(df, nonc, 1)
        self.assertIn("nonc < 0", str(raises.exception))        

        df, nonc = 1, 1
        sizes = (True, 3j, 1.5, (1.5, 1), (3j, 1), (3j, 3j), (np.int8(3), np.int64(7)))
        for size in sizes:
            with self.assertRaises(TypingError) as raises:
                cfunc(df, nonc, size)
            self.assertIn(
                "np.random.noncentral_chisquare(): size should be int or "
                "tuple of ints or None, got",
                str(raises.exception),
            )

@jit(nopython=True, nogil=True)
def py_extract_randomness(seed, out):
    if seed != 0:
        random.seed(seed)
    for i in range(out.size):
        out[i] = random.getrandbits(32)

_randint_limit = 1 << 32

@jit(nopython=True, nogil=True)
def np_extract_randomness(seed, out):
    if seed != 0:
        np.random.seed(seed)
    s = 0
    for i in range(out.size):
        out[i] = np.random.randint(_randint_limit)



class ConcurrencyBaseTest(TestCase):

    # Enough iterations for:
    # 1. Mersenne-Twister state shuffles to occur (once every 624)
    # 2. Race conditions to be plausible
    # 3. Nice statistical properties to emerge
    _extract_iterations = 100000

    def setUp(self):
        # Warm up, to avoid compiling in the threads
        args = (42, self._get_output(1))
        py_extract_randomness(*args)
        np_extract_randomness(*args)

    def _get_output(self, size):
        return np.zeros(size, dtype=np.uint32)

    def check_output(self, out):
        """
        Check statistical properties of output.
        """
        # Output should follow a uniform distribution in [0, 1<<32)
        expected_avg = 1 << 31
        expected_std = (1 << 32) / np.sqrt(12)
        rtol = 0.05  # given enough iterations
        np.testing.assert_allclose(out.mean(), expected_avg, rtol=rtol)
        np.testing.assert_allclose(out.std(), expected_std, rtol=rtol)

    def check_several_outputs(self, results, same_expected):
        # Outputs should have the expected statistical properties
        # (an uninitialized PRNG or a PRNG whose internal state was
        #  corrupted by a race condition could produce bogus randomness)
        for out in results:
            self.check_output(out)

        # Check all threads gave either the same sequence or
        # distinct sequences
        if same_expected:
            expected_distinct = 1
        else:
            expected_distinct = len(results)

        heads = {tuple(out[:5]) for out in results}
        tails = {tuple(out[-5:]) for out in results}
        sums = {out.sum() for out in results}
        self.assertEqual(len(heads), expected_distinct, heads)
        self.assertEqual(len(tails), expected_distinct, tails)
        self.assertEqual(len(sums), expected_distinct, sums)


class TestThreads(ConcurrencyBaseTest):
    """
    Check the PRNG behaves well with threads.
    """

    def extract_in_threads(self, nthreads, extract_randomness, seed):
        """
        Run *nthreads* threads extracting randomness with the given *seed*
        (no seeding if 0).
        """
        results = [self._get_output(self._extract_iterations)
                   for i in range(nthreads + 1)]

        def target(i):
            # The PRNG will be seeded in thread
            extract_randomness(seed=seed, out=results[i])

        threads = [threading.Thread(target=target, args=(i,))
                   for i in range(nthreads)]

        for th in threads:
            th.start()
        # Exercise main thread as well
        target(nthreads)
        for th in threads:
            th.join()

        return results

    def check_thread_safety(self, extract_randomness):
        """
        When initializing the PRNG the same way, each thread
        should produce the same sequence of random numbers,
        using independent states, regardless of parallel
        execution.
        """
        # Note the seed value doesn't matter, as long as it's
        # the same for all threads
        results = self.extract_in_threads(15, extract_randomness, seed=42)

        # All threads gave the same sequence
        self.check_several_outputs(results, same_expected=True)

    def check_implicit_initialization(self, extract_randomness):
        """
        The PRNG in new threads should be implicitly initialized with
        system entropy, if seed() wasn't called.
        """
        results = self.extract_in_threads(4, extract_randomness, seed=0)

        # All threads gave a different, valid random sequence
        self.check_several_outputs(results, same_expected=False)

    def test_py_thread_safety(self):
        self.check_thread_safety(py_extract_randomness)

    def test_np_thread_safety(self):
        self.check_thread_safety(np_extract_randomness)

    def test_py_implicit_initialization(self):
        self.check_implicit_initialization(py_extract_randomness)

    def test_np_implicit_initialization(self):
        self.check_implicit_initialization(np_extract_randomness)


@unittest.skipIf(os.name == 'nt', "Windows is not affected by fork() issues")
class TestProcesses(ConcurrencyBaseTest):
    """
    Check the PRNG behaves well in child processes.
    """

    # Avoid nested multiprocessing AssertionError
    # ("daemonic processes are not allowed to have children")
    _numba_parallel_test_ = False


    def extract_in_processes(self, nprocs, extract_randomness):
        """
        Run *nprocs* processes extracting randomness
        without explicit seeding.
        """
        q = multiprocessing.Queue()
        results = []

        def target_inner():
            out = self._get_output(self._extract_iterations)
            extract_randomness(seed=0, out=out)
            return out

        def target():
            try:
                out = target_inner()
                q.put(out)
            except Exception as e:
                # Ensure an exception in a child gets reported
                # in the parent.
                q.put(e)
                raise

        if hasattr(multiprocessing, 'get_context'):
            # The test works only in fork context.
            mpc = multiprocessing.get_context('fork')
        else:
            mpc = multiprocessing
        procs = [mpc.Process(target=target)
                 for i in range(nprocs)]
        for p in procs:
            p.start()
        # Need to dequeue before joining, otherwise the large size of the
        # enqueued objects will lead to deadlock.
        for i in range(nprocs):
            results.append(q.get(timeout=5))
        for p in procs:
            p.join()

        # Exercise parent process as well; this will detect if the
        # same state was reused for one of the children.
        results.append(target_inner())
        for res in results:
            if isinstance(res, Exception):
                self.fail("Exception in child: %s" % (res,))

        return results

    def check_implicit_initialization(self, extract_randomness):
        """
        The PRNG in new processes should be implicitly initialized
        with system entropy, to avoid reproducing the same sequences.
        """
        results = self.extract_in_processes(2, extract_randomness)

        # All processes gave a different, valid random sequence
        self.check_several_outputs(results, same_expected=False)

    def test_py_implicit_initialization(self):
        self.check_implicit_initialization(py_extract_randomness)

    def test_np_implicit_initialization(self):
        self.check_implicit_initialization(np_extract_randomness)


class TestNumPyRandomAPI(TestCase):

    def test_call_by_name(self):
        # Checks that the NumPy impls in Numba can be used via call-by-name
        # args, see issue numba#9053.
        #
        # Checking call-by-name has to be done somewhat manually as the NumPy
        # numpy.random.* functions do not have signatures, see numpy#8734.

        # Herein, it doesn't matter what the values are, the names and types
        # just have to make sense.
        data = {"np.random.beta": {'a': 1., 'b': 2., 'size': 3},
                "np.random.binomial": {'n': 1, 'p': 0.3, 'size': 3},
                "np.random.chisquare": {'df': 2., 'size': 3},
                "np.random.choice": {'a': 2, 'size': 3},
                "np.random.dirichlet": {'alpha': (2,), 'size': 3},
                "np.random.exponential": {'scale': 1., 'size': 3},
                "np.random.f": {'dfnum': 1., 'dfden': 2., 'size': 3},
                "np.random.gamma": {'shape': 2, 'scale': 2.0, 'size': 3},
                "np.random.geometric": {'p': 1., 'size': 3},
                "np.random.gumbel": {'loc': 0., 'scale': 1., 'size': 3},
                "np.random.hypergeometric": {'ngood': 1, 'nbad': 1,
                                             'nsample': 1, 'size': 3},
                "np.random.laplace": {'loc': 0., 'scale': 1., 'size': 3},
                "np.random.logistic": {'loc': 0., 'scale': 1., 'size': 3},
                "np.random.lognormal": {'mean': 0., 'sigma': 1., 'size': 3},
                "np.random.logseries": {'p': 0.5, 'size': 3},
                "np.random.multinomial": {'n': 1, 'pvals': (1,), 'size': 3},
                "np.random.negative_binomial": {'n': 1, 'p': 0.5},
                "np.random.noncentral_chisquare": {'df': 1., 'nonc': 1.,
                                                   'size': 3},
                "np.random.normal": {'loc': 0., 'scale': 1., 'size': 3},
                "np.random.pareto": {'a': 2., 'size': 3},
                # NOTE: The NumPy impl of permutation "takes no keyword
                # arguments".
                # "np.random.permutation": {'x': (1, 2, 3)},
                "np.random.poisson": {'lam': 1., 'size': 3},
                "np.random.power": {'a': 2., 'size': 3},
                # NOTE: The NumPy impl of rand essentially takes *args so kwargs
                # are unsupported.
                # "np.random.rand": {'d0': 1, 'd1': 2, ...}}
                "np.random.randint": {'low': 1, 'high': 2, 'size': 3},
                # NOTE: The NumPy impl of randn essentially takes *args so
                # kwargs are unsupported.
                # "np.random.randn":  {'d0': 1, 'd1': 2, ...}}
                "np.random.random": {'size': 3},
                "np.random.random_sample": {'size': 3},
                "np.random.ranf": {'size': 3},
                "np.random.rayleigh": {'scale': 1., 'size': 3},
                "np.random.sample": {'size': 3},
                "np.random.seed": {'seed': 4},
                # NOTE: The NumPy impl of shuffle "takes no keyword arguments".
                # "np.random.shuffle"
                "np.random.standard_cauchy": {'size': 3},
                "np.random.standard_exponential": {'size': 3},
                "np.random.standard_gamma": {'shape': 2., 'size': 3},
                "np.random.standard_normal": {'size': 3},
                "np.random.standard_t": {'df': 2., 'size': 3},
                "np.random.triangular": {'left': 1., 'mode': 2., 'right': 3.,
                                         'size': 3},
                "np.random.uniform": {'low': 1., 'high': 2., 'size': 3},
                "np.random.vonmises": {'mu': 1., 'kappa': 2., 'size': 3},
                "np.random.wald": {'mean': 1., 'scale': 2., 'size': 3},
                "np.random.weibull": {'a': 1., 'size': 3},
                "np.random.zipf": {'a': 2., 'size': 3},}

        for fn, args in data.items():
            argstr = ', '.join([f'{k}={v}' for k, v in args.items()])
            template = dedent(f"""
                def foo():
                    return {fn}({argstr})
                """)
            l = {}
            exec(template, {'np': np}, l)
            # The answer doesn't matter, these are tested in the tests above,
            # the purpose of this test is to ensure that the code compiles with
            # the args presented via name, i.e. the overloads are defined
            # correctly with respect to the public API of the function.
            func = l['foo']
            func()
            njit(func).compile(())


if __name__ == "__main__":
    unittest.main()
