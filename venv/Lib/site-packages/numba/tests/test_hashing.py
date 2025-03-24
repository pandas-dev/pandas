# -*- coding: utf-8 -*-
"""
Test hashing of various supported types.
"""

import unittest

import os
import sys
import subprocess
from collections import defaultdict
from textwrap import dedent

import numpy as np

from numba import jit, config, typed, typeof
from numba.core import types, utils
import unittest
from numba.tests.support import (TestCase, skip_unless_py10_or_later,
                                 run_in_subprocess)

from numba.cpython.unicode import compile_time_get_string_data
from numba.cpython import hashing
from numba.np.numpy_support import numpy_version


def hash_usecase(x):
    return hash(x)


class TestHashingSetup(TestCase):

    def test_warn_on_fnv(self):
        # FNV hash alg variant is not supported, check Numba warns
        work = """
        import sys
        import warnings
        from collections import namedtuple

        # hash_info is a StructSequence, mock as a named tuple
        fields = ["width", "modulus", "inf", "nan", "imag", "algorithm",
                  "hash_bits", "seed_bits", "cutoff"]

        hinfo = sys.hash_info
        FAKE_HASHINFO = namedtuple('FAKE_HASHINFO', fields)

        fd = dict()
        for f in fields:
            fd[f] = getattr(hinfo, f)

        fd['algorithm'] = 'fnv'

        fake_hashinfo = FAKE_HASHINFO(**fd)

        # replace the hashinfo with the fnv version
        sys.hash_info = fake_hashinfo
        with warnings.catch_warnings(record=True) as warns:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            from numba import njit
            @njit
            def foo():
                hash(1)
            foo()
            assert len(warns) > 0
            expect = "FNV hashing is not implemented in Numba. See PEP 456"
            for w in warns:
                if expect in str(w.message):
                    break
            else:
                raise RuntimeError("Expected warning not found")
        """
        subprocess.check_call([sys.executable, '-c', dedent(work)])


class TestHashAlgs(TestCase):
    # This tests Numba hashing replication against cPython "gold", i.e. the
    # actual hash values for given inputs, algs and PYTHONHASHSEEDs
    # Test adapted from:
    # https://github.com/python/cpython/blob/9dda9020abcf0d51d59b283a89c58c8e1fb0f574/Lib/test/test_hash.py#L197-L264
    # and
    # https://github.com/python/cpython/blob/9dda9020abcf0d51d59b283a89c58c8e1fb0f574/Lib/test/test_hash.py#L174-L189

    # 32bit little, 64bit little, 32bit big, 64bit big
    known_hashes = {
        'djba33x': [ # only used for small strings
            # seed 0, 'abc'
            [193485960, 193485960,  193485960, 193485960],
            # seed 42, 'abc'
            [-678966196, 573763426263223372, -820489388, -4282905804826039665],
            ],
        'siphash13': [
            # NOTE: PyUCS2 layout depends on endianness
            # seed 0, 'abc'
            [69611762, -4594863902769663758, 69611762, -4594863902769663758],
            # seed 42, 'abc'
            [-975800855, 3869580338025362921, -975800855, 3869580338025362921],
            # seed 42, 'abcdefghijk'
            [-595844228, 7764564197781545852, -595844228, 7764564197781545852],
            # seed 0, 'äú∑ℇ'
            [-1093288643, -2810468059467891395, -1041341092, 4925090034378237276],
            # seed 42, 'äú∑ℇ'
            [-585999602, -2845126246016066802, -817336969, -2219421378907968137],
        ],
        'siphash24': [
            # NOTE: PyUCS2 layout depends on endianness
            # seed 0, 'abc'
            [1198583518, 4596069200710135518, 1198583518, 4596069200710135518],
            # seed 42, 'abc'
            [273876886, -4501618152524544106, 273876886, -4501618152524544106],
            # seed 42, 'abcdefghijk'
            [-1745215313, 4436719588892876975, -1745215313, 4436719588892876975],
            # seed 0, 'äú∑ℇ'
            [493570806, 5749986484189612790, -1006381564, -5915111450199468540],
            # seed 42, 'äú∑ℇ'
            [-1677110816, -2947981342227738144, -1860207793, -4296699217652516017],
        ],
    }

    def get_expected_hash(self, position, length):
        if length < sys.hash_info.cutoff:
            algorithm = "djba33x"
        else:
            algorithm = sys.hash_info.algorithm
        IS_64BIT = not config.IS_32BITS
        if sys.byteorder == 'little':
            platform = 1 if IS_64BIT else 0
        else:
            assert(sys.byteorder == 'big')
            platform = 3 if IS_64BIT else 2
        return self.known_hashes[algorithm][position][platform]

    def get_hash_command(self, repr_):
        return 'print(hash(eval(%a)))' % repr_

    def get_hash(self, repr_, seed=None):
        env = os.environ.copy()
        if seed is not None:
            env['PYTHONHASHSEED'] = str(seed)
        else:
            env.pop('PYTHONHASHSEED', None)
        out, _ = run_in_subprocess(code=self.get_hash_command(repr_),
                                   env=env)
        stdout = out.decode().strip()
        return int(stdout)

    def test_against_cpython_gold(self):

        args = (('abc', 0, 0), ('abc', 42, 1), ('abcdefghijk', 42, 2),
                ('äú∑ℇ', 0, 3), ('äú∑ℇ', 42, 4),)

        for input_str, seed, position in args:
            with self.subTest(input_str=input_str, seed=seed):
                got = self.get_hash(repr(input_str), seed=seed)
                expected = self.get_expected_hash(position, len(input_str))
                self.assertEqual(got, expected)


class BaseTest(TestCase):

    def setUp(self):
        self.cfunc = jit(nopython=True)(hash_usecase)

    def check_hash_values(self, values):
        cfunc = self.cfunc
        for val in list(values):
            nb_hash = cfunc(val)
            self.assertIsInstance(nb_hash, int)
            try:
                self.assertEqual(nb_hash, hash(val))
            except AssertionError as e:
                print("val, nb_hash, hash(val)")
                print(val, nb_hash, hash(val))
                print("abs(val), hashing._PyHASH_MODULUS - 1")
                print(abs(val), hashing._PyHASH_MODULUS - 1)
                raise e

    def int_samples(self, typ=np.int64):
        for start in (0, -50, 60000, 1 << 32):
            info = np.iinfo(typ)
            if not info.min <= start <= info.max:
                continue
            n = 100
            yield range(start, start + n)
            yield range(start, start + 100 * n, 100)
            yield range(start, start + 128 * n, 128)
            yield [-1]

    def float_samples(self, typ):
        info = np.finfo(typ)

        for start in (0, 10, info.max ** 0.5, info.max / 1000.0):
            n = 100
            min_step = max(info.tiny, start * info.resolution)
            for step in (1.2, min_step ** 0.5, min_step):
                if step < min_step:
                    continue
                a = np.linspace(start, start + n * step, n)
                a = a.astype(typ)
                yield a
                yield -a
                yield a + a.mean()

        # Infs, nans, zeros, magic -1
        a = [0.0, 0.5, -0.0, -1.0, float('inf'), -float('inf'),]

        # Python 3.10 has a hash for nan based on the pointer to the PyObject
        # containing the nan, skip this input and use explicit test instead.

        yield typ(a)

    def complex_samples(self, typ, float_ty):
        for real in self.float_samples(float_ty):
            for imag in self.float_samples(float_ty):
                # Ensure equal sizes
                real = real[:len(imag)]
                imag = imag[:len(real)]
                a = real + typ(1j) * imag
                # Python 3.10 has a hash for nan based on the pointer to the
                # PyObject containing the nan, skip input that ends up as nan
                if not np.any(np.isnan(a)):
                    yield a


class TestNumberHashing(BaseTest):
    """
    Test hashing of number types.
    """

    def setUp(self):
        if numpy_version >= (2, 0):
            # Temporarily set promotions state to legacy,
            # to ensure overflow logic works
            self.initial_state = np._get_promotion_state()
            np._set_promotion_state("legacy")

        return super().setUp()

    def tearDown(self) -> None:
        if numpy_version >= (2, 0):
            # Reset numpy promotion state to initial state
            # since the setting is global
            np._set_promotion_state(self.initial_state)

        return super().tearDown()

    def check_floats(self, typ):
        for a in self.float_samples(typ):
            self.assertEqual(a.dtype, np.dtype(typ))
            self.check_hash_values(a)

    def check_complex(self, typ, float_ty):
        for a in self.complex_samples(typ, float_ty):
            self.assertEqual(a.dtype, np.dtype(typ))
            self.check_hash_values(a)

    def test_floats(self):
        self.check_floats(np.float32)
        self.check_floats(np.float64)

    def test_complex(self):
        self.check_complex(np.complex64, np.float32)
        self.check_complex(np.complex128, np.float64)

    def test_bool(self):
        self.check_hash_values([False, True])

    def test_ints(self):
        minmax = []

        for ty in [np.int8, np.uint8, np.int16, np.uint16,
                   np.int32, np.uint32, np.int64, np.uint64]:
            for a in self.int_samples(ty):
                self.check_hash_values(a)
            info = np.iinfo(ty)
            # check hash(-1) = -2
            # check hash(0) = 0
            self.check_hash_values([ty(-1)])
            self.check_hash_values([ty(0)])
            signed = 'uint' not in str(ty)
            # check bit shifting patterns from min through to max
            sz = ty().itemsize
            for x in [info.min, info.max]:
                shifts = 8 * sz
                # x is a python int, do shifts etc as a python int and init
                # numpy type from that to avoid numpy type rules
                y = x
                for i in range(shifts):
                    twiddle1 = 0xaaaaaaaaaaaaaaaa
                    twiddle2 = 0x5555555555555555
                    vals = [y]
                    for tw in [twiddle1, twiddle2]:
                        val = y & twiddle1
                        if val < sys.maxsize:
                            vals.append(val)
                    for v in vals:
                        self.check_hash_values([ty(v)])
                    if signed:  # try the same with flipped signs
                        # negated signed INT_MIN will overflow
                        for v in vals:
                            if v != info.min:
                                self.check_hash_values([ty(-v)])
                    if x == 0:  # unsigned min is 0, shift up
                        y = (y | 1) << 1
                    else:  # everything else shift down
                        y = y >> 1

        # these straddle the branch between returning the int as the hash and
        # doing the PyLong hash alg
        self.check_hash_values([np.int64(0x1ffffffffffffffe)])
        self.check_hash_values([np.int64(0x1fffffffffffffff)])
        self.check_hash_values([np.uint64(0x1ffffffffffffffe)])
        self.check_hash_values([np.uint64(0x1fffffffffffffff)])

        # check some values near sys int mins
        self.check_hash_values([np.int64(-0x7fffffffffffffff)])
        self.check_hash_values([np.int64(-0x7ffffffffffffff6)])
        self.check_hash_values([np.int64(-0x7fffffffffffff9c)])
        self.check_hash_values([np.int32(-0x7fffffff)])
        self.check_hash_values([np.int32(-0x7ffffff6)])
        self.check_hash_values([np.int32(-0x7fffff9c)])

    @skip_unless_py10_or_later
    def test_py310_nan_hash(self):
        # On Python 3.10+ nan's hash to a value which is based on the pointer to
        # the PyObject containing the nan. Numba cannot replicate as there's no
        # object, it instead produces equivalent behaviour, i.e. hashes to
        # something "unique".

        # Run 10 hashes, make sure that the "uniqueness" is sufficient that
        # there's more than one hash value. Not much more can be done!
        x = [float('nan') for i in range(10)]
        out = set([self.cfunc(z) for z in x])
        self.assertGreater(len(out), 1)


class TestTupleHashing(BaseTest):
    """
    Test hashing of tuples.
    """

    def setUp(self):
        if numpy_version >= (2, 0):
            # Temporarily set promotions state to legacy,
            # to ensure overflow logic works
            self.initial_state = np._get_promotion_state()
            np._set_promotion_state("legacy")

        return super().setUp()

    def tearDown(self) -> None:
        if numpy_version >= (2, 0):
            # Reset numpy promotion state to initial state
            # since the setting is global
            np._set_promotion_state(self.initial_state)

        return super().tearDown()

    def check_tuples(self, value_generator, split):
        for values in value_generator:
            tuples = [split(a) for a in values]
            self.check_hash_values(tuples)

    def test_homogeneous_tuples(self):
        typ = np.uint64

        def split2(i):
            """
            Split i's bits into 2 integers.
            """
            i = typ(i)
            return (i & typ(0x5555555555555555),
                    i & typ(0xaaaaaaaaaaaaaaaa),
                    )

        def split3(i):
            """
            Split i's bits into 3 integers.
            """
            i = typ(i)
            return (i & typ(0x2492492492492492),
                    i & typ(0x4924924924924924),
                    i & typ(0x9249249249249249),
                    )

        self.check_tuples(self.int_samples(), split2)
        self.check_tuples(self.int_samples(), split3)

        # Check exact. Sample values from:
        # https://github.com/python/cpython/blob/b738237d6792acba85b1f6e6c8993a812c7fd815/Lib/test/test_tuple.py#L80-L93
        # Untypable empty tuples are replaced with (7,).
        self.check_hash_values([(7,), (0,), (0, 0), (0.5,),
                                (0.5, (7,), (-2, 3, (4, 6)))])

    def test_heterogeneous_tuples(self):
        modulo = 2**63

        def split(i):
            a = i & 0x5555555555555555
            b = (i & 0xaaaaaaaa) ^ ((i >> 32) & 0xaaaaaaaa)
            return np.int64(a), np.float64(b * 0.0001)

        self.check_tuples(self.int_samples(), split)


class TestUnicodeHashing(BaseTest):

    def test_basic_unicode(self):
        kind1_string = "abcdefghijklmnopqrstuvwxyz"
        for i in range(len(kind1_string)):
            self.check_hash_values([kind1_string[:i]])

        sep = "眼"
        kind2_string = sep.join(list(kind1_string))
        for i in range(len(kind2_string)):
            self.check_hash_values([kind2_string[:i]])

        sep = "🐍⚡"
        kind4_string = sep.join(list(kind1_string))
        for i in range(len(kind4_string)):
            self.check_hash_values([kind4_string[:i]])

        empty_string = ""
        self.check_hash_values(empty_string)

    def test_hash_passthrough(self):
        # no `hash` call made, this just checks that `._hash` is correctly
        # passed through from an already existing string
        kind1_string = "abcdefghijklmnopqrstuvwxyz"

        @jit(nopython=True)
        def fn(x):
            return x._hash

        hash_value = compile_time_get_string_data(kind1_string)[-1]
        self.assertTrue(hash_value != -1)
        self.assertEqual(fn(kind1_string), hash_value)

    def test_hash_passthrough_call(self):
        # check `x._hash` and hash(x) are the same
        kind1_string = "abcdefghijklmnopqrstuvwxyz"

        @jit(nopython=True)
        def fn(x):
            return x._hash, hash(x)

        hash_value = compile_time_get_string_data(kind1_string)[-1]
        self.assertTrue(hash_value != -1)
        self.assertEqual(fn(kind1_string), (hash_value, hash_value))

    @unittest.skip("Needs hash computation at const unpickling time")
    def test_hash_literal(self):
        # a strconst always seem to have an associated hash value so the hash
        # member of the returned value should contain the correct hash
        @jit(nopython=True)
        def fn():
            x = "abcdefghijklmnopqrstuvwxyz"
            return x
        val = fn()
        tmp = hash("abcdefghijklmnopqrstuvwxyz")
        self.assertEqual(tmp, (compile_time_get_string_data(val)[-1]))

    def test_hash_on_str_creation(self):
        # In cPython some? new strings do not have a cached hash until hash() is
        # called
        def impl(do_hash):
            const1 = "aaaa"
            const2 = "眼眼眼眼"
            new = const1 + const2
            if do_hash:
                hash(new)
            return new

        jitted = jit(nopython=True)(impl)

        # do not compute the hash, cPython will have no cached hash, but Numba
        # will
        compute_hash = False
        expected = impl(compute_hash)
        got = jitted(compute_hash)
        a = (compile_time_get_string_data(expected))
        b = (compile_time_get_string_data(got))
        self.assertEqual(a[:-1], b[:-1])
        self.assertTrue(a[-1] != b[-1])

        # now with compute hash enabled, cPython will have a cached hash as will
        # Numba
        compute_hash = True
        expected = impl(compute_hash)
        got = jitted(compute_hash)
        a = (compile_time_get_string_data(expected))
        b = (compile_time_get_string_data(got))
        self.assertEqual(a, b)


class TestUnhashable(TestCase):
    # Tests that unhashable types behave correctly and raise a TypeError at
    # runtime.

    def test_hash_unhashable(self):
        unhashables = (typed.Dict().empty(types.int64, types.int64),
                       typed.List().empty_list(types.int64),
                       np.ones(4))
        cfunc = jit(nopython=True)(hash_usecase)
        for ty in unhashables:
            with self.assertRaises(TypeError) as raises:
                cfunc(ty)
            expected = f"unhashable type: '{str(typeof(ty))}'"
            self.assertIn(expected, str(raises.exception))

    def test_no_generic_hash(self):
        # In CPython, if there's no attr `__hash__` on an object, a hash of the
        # object's pointer is returned (see: _Py_HashPointer in the CPython
        # source). Numba has no access to such objects and can't create them
        # either, so it catches this case and raises an exception.

        @jit(nopython=True)
        def foo():
            hash(np.cos)

        with self.assertRaises(TypeError) as raises:
            foo()

        expected = ("No __hash__ is defined for object ")
        self.assertIn(expected, str(raises.exception))


if __name__ == "__main__":
    unittest.main()
