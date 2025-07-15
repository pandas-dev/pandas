"""
Test np.datetime64 and np.timedelta64 support.
"""

# NOTE: datetime64 and timedelta64 ufuncs are tested in test_ufuncs.


import contextlib
import itertools
import re
import unittest
import warnings

import numpy as np

from numba import jit, vectorize, njit
from numba.np.numpy_support import numpy_version
from numba.core import types, config
from numba.core.errors import TypingError
from numba.tests.support import TestCase, tag, skip_parfors_unsupported
from numba.np import npdatetime_helpers, numpy_support

TIMEDELTA_M = np.dtype('timedelta64[M]')
TIMEDELTA_Y = np.dtype('timedelta64[Y]')

def value_unit(val):
    ty = numpy_support.from_dtype(val.dtype)
    return ty.unit


date_units = ('Y', 'M')
time_units = ('W', 'D', 'h', 'm', 's', 'ms', 'us', 'ns', 'ps', 'fs', 'as')
# All except generic ("")
all_units = date_units + time_units


def add_usecase(x, y):
    return x + y

def sub_usecase(x, y):
    return x - y

def mul_usecase(x, y):
    return x * y

def div_usecase(x, y):
    return x / y

def floordiv_usecase(x, y):
    return x // y

def eq_usecase(x, y):
    return x == y

def ne_usecase(x, y):
    return x != y

def lt_usecase(x, y):
    return x < y

def le_usecase(x, y):
    return x <= y

def gt_usecase(x, y):
    return x > y

def ge_usecase(x, y):
    return x >= y

def pos_usecase(x):
    return +x

def neg_usecase(x):
    return -x

def abs_usecase(x):
    return abs(x)

def hash_usecase(x):
    return hash(x)

def min_usecase(x, y):
    return min(x, y)

def max_usecase(x, y):
    return max(x, y)

def int_cast_usecase(x):
    return int(x)

def make_add_constant(const):
    def add_constant(x):
        return x + const
    return add_constant


class TestModuleHelpers(TestCase):
    """
    Test the various helpers in numba.npdatetime_helpers.
    """

    def test_can_cast_timedelta(self):
        f = npdatetime_helpers.can_cast_timedelta_units
        for a, b in itertools.product(date_units, time_units):
            self.assertFalse(f(a, b), (a, b))
            self.assertFalse(f(b, a), (a, b))
        for unit in all_units:
            self.assertFalse(f(unit, ''))
            self.assertTrue(f('', unit))
        for unit in all_units + ('',):
            self.assertTrue(f(unit, unit))

        def check_units_group(group):
            for i, a in enumerate(group):
                for b in group[:i]:
                    # large into smaller is ok
                    self.assertTrue(f(b, a))
                    # small into larger is not
                    self.assertFalse(f(a, b))

        check_units_group(date_units)
        check_units_group(time_units)

    def test_timedelta_conversion(self):
        f = npdatetime_helpers.get_timedelta_conversion_factor
        for unit in all_units + ('',):
            self.assertEqual(f(unit, unit), 1)
        for unit in all_units:
            self.assertEqual(f('', unit), 1)
        for a, b in itertools.product(time_units, date_units):
            self.assertIs(f(a, b), None)
            self.assertIs(f(b, a), None)

        def check_units_group(group):
            for i, a in enumerate(group):
                for b in group[:i]:
                    self.assertGreater(f(b, a), 1, (b, a))
                    self.assertIs(f(a, b), None)

        check_units_group(date_units)
        check_units_group(time_units)

        # Check some hand-picked values
        self.assertEqual(f('Y', 'M'), 12)
        self.assertEqual(f('W', 'h'), 24 * 7)
        self.assertEqual(f('W', 'm'), 24 * 7 * 60)
        self.assertEqual(f('W', 'us'), 24 * 7 * 3600 * 1000 * 1000)

    def test_datetime_timedelta_scaling(self):
        f = npdatetime_helpers.get_datetime_timedelta_conversion
        def check_error(dt_unit, td_unit):
            with self.assertRaises(RuntimeError):
                f(dt_unit, td_unit)
        # Cannot combine a Y or M timedelta64 with a finer-grained datetime64
        for dt_unit, td_unit in itertools.product(time_units, date_units):
            check_error(dt_unit, td_unit)
        # Sanity check that all other unit pairs can be converted, we'll
        # check individual results below
        for dt_unit, td_unit in itertools.product(time_units, time_units):
            f(dt_unit, td_unit)
        for dt_unit, td_unit in itertools.product(date_units, time_units):
            f(dt_unit, td_unit)
        for dt_unit, td_unit in itertools.product(date_units, date_units):
            f(dt_unit, td_unit)
        # No-op conversions
        for unit in all_units:
            self.assertEqual(f(unit, unit), (unit, 1, 1))
            self.assertEqual(f(unit, ''), (unit, 1, 1))
            self.assertEqual(f('', unit), ('', 1, 1))
        self.assertEqual(f('', ''), ('', 1, 1))
        # "Regular" values
        self.assertEqual(f('Y', 'M'), ('M', 12, 1))
        self.assertEqual(f('M', 'Y'), ('M', 1, 12))
        self.assertEqual(f('W', 'D'), ('D', 7, 1))
        self.assertEqual(f('D', 'W'), ('D', 1, 7))
        self.assertEqual(f('W', 's'), ('s', 7 * 24 * 3600, 1))
        self.assertEqual(f('s', 'W'), ('s', 1, 7 * 24 * 3600))
        self.assertEqual(f('s', 'as'), ('as', 1000 ** 6, 1))
        self.assertEqual(f('as', 's'), ('as', 1, 1000 ** 6))
        # "Interesting" values
        self.assertEqual(f('Y', 'D'), ('D', 97 + 400 * 365, 400))
        self.assertEqual(f('Y', 'W'), ('W', 97 + 400 * 365, 400 * 7))
        self.assertEqual(f('M', 'D'), ('D', 97 + 400 * 365, 400 * 12))
        self.assertEqual(f('M', 'W'), ('W', 97 + 400 * 365, 400 * 12 * 7))
        self.assertEqual(f('Y', 's'), ('s', (97 + 400 * 365) *  24 * 3600, 400))
        self.assertEqual(f('M', 's'), ('s', (97 + 400 * 365) *  24 * 3600, 400 * 12))

    def test_combine_datetime_timedelta_units(self):
        f = npdatetime_helpers.combine_datetime_timedelta_units
        for unit in all_units:
            self.assertEqual(f(unit, unit), unit)
            self.assertEqual(f('', unit), unit)
            self.assertEqual(f(unit, ''), unit)
        self.assertEqual(f('', ''), '')
        for dt_unit, td_unit in itertools.product(time_units, date_units):
            self.assertIs(f(dt_unit, td_unit), None)
        for dt_unit, td_unit in itertools.product(date_units, time_units):
            self.assertEqual(f(dt_unit, td_unit), td_unit)

    def test_same_kind(self):
        f = npdatetime_helpers.same_kind
        for u in all_units:
            self.assertTrue(f(u, u))
        A = ('Y', 'M', 'W', 'D')
        B = ('h', 'm', 's', 'ms', 'us', 'ns', 'ps', 'fs', 'as')
        for a, b in itertools.product(A, A):
            self.assertTrue(f(a, b))
        for a, b in itertools.product(B, B):
            self.assertTrue(f(a, b))
        for a, b in itertools.product(A, B):
            self.assertFalse(f(a, b))
            self.assertFalse(f(b, a))


TD = np.timedelta64
DT = np.datetime64


class TestMiscCompiling(TestCase):

    def test_jit_explicit_signature(self):
        def _check_explicit_signature(sig):
            f = jit(sig, nopython=True)(add_usecase)
            # Just a sanity check
            args = DT(1, 'ms'), TD(2, 'us')
            expected = add_usecase(*args)
            self.assertPreciseEqual(f(*args), expected)

        # Test passing the signature in object form
        sig = types.NPDatetime('us')(types.NPDatetime('ms'), types.NPTimedelta('us'))
        _check_explicit_signature(sig)
        # Same with the signature in string form
        sig = "NPDatetime('us')(NPDatetime('ms'), NPTimedelta('us'))"
        _check_explicit_signature(sig)

    def test_vectorize_explicit_signature(self):
        def _check_explicit_signature(sig):
            f = vectorize([sig], nopython=True)(mul_usecase)
            # This isn't really right but we can't do better than this,
            # since Numpy's ufuncs don't store the metadata of return types.
            # Related to https://github.com/numpy/numpy/issues/5429
            self.assertPreciseEqual(f(TD(2), 3), TD(6))

        # Test passing the signature in object form (issue #917)
        sig = types.NPTimedelta('s')(types.NPTimedelta('s'), types.int64)
        _check_explicit_signature(sig)
        # Same with the signature in string form
        sig = "NPTimedelta('s')(NPTimedelta('s'), int64)"
        _check_explicit_signature(sig)

    def test_constant_datetime(self):
        def check(const):
            pyfunc = make_add_constant(const)
            f = jit(nopython=True)(pyfunc)
            x = TD(4, 'D')
            expected = pyfunc(x)
            self.assertPreciseEqual(f(x), expected)
        check(DT('2001-01-01'))
        check(DT('NaT', 'D'))

    def test_constant_timedelta(self):
        def check(const):
            pyfunc = make_add_constant(const)
            f = jit(nopython=True)(pyfunc)
            x = TD(4, 'D')
            expected = pyfunc(x)
            self.assertPreciseEqual(f(x), expected)
        check(TD(4, 'D'))
        check(TD(-4, 'D'))
        check(TD('NaT', 'D'))


class TestTimedeltaArithmetic(TestCase):

    jitargs = dict(forceobj=True)

    def jit(self, pyfunc):
        return jit(**self.jitargs)(pyfunc)

    def test_add(self):
        f = self.jit(add_usecase)
        def check(a, b, expected):
            self.assertPreciseEqual(f(a, b), expected)
            self.assertPreciseEqual(f(b, a), expected)

        check(TD(1), TD(2), TD(3))
        check(TD(1, 's'), TD(2, 's'), TD(3, 's'))
        # Implicit unit promotion
        check(TD(1, 's'), TD(2, 'us'), TD(1000002, 'us'))
        check(TD(1, 'W'), TD(2, 'D'), TD(9, 'D'))
        # NaTs
        check(TD('NaT'), TD(1), TD('NaT'))
        check(TD('NaT', 's'), TD(1, 'D'), TD('NaT', 's'))
        check(TD('NaT', 's'), TD(1, 'ms'), TD('NaT', 'ms'))
        # Cannot add days and months
        with self.assertRaises((TypeError, TypingError)):
            f(TD(1, 'M'), TD(1, 'D'))

    def test_sub(self):
        f = self.jit(sub_usecase)
        def check(a, b, expected):
            self.assertPreciseEqual(f(a, b), expected)
            self.assertPreciseEqual(f(b, a), -expected)

        check(TD(3), TD(2), TD(1))
        check(TD(3, 's'), TD(2, 's'), TD(1, 's'))
        # Implicit unit promotion
        check(TD(3, 's'), TD(2, 'us'), TD(2999998, 'us'))
        check(TD(1, 'W'), TD(2, 'D'), TD(5, 'D'))
        # NaTs
        check(TD('NaT'), TD(1), TD('NaT'))
        check(TD('NaT', 's'), TD(1, 'D'), TD('NaT', 's'))
        check(TD('NaT', 's'), TD(1, 'ms'), TD('NaT', 'ms'))
        # Cannot sub days to months
        with self.assertRaises((TypeError, TypingError)):
            f(TD(1, 'M'), TD(1, 'D'))

    def test_mul(self):
        f = self.jit(mul_usecase)
        def check(a, b, expected):
            self.assertPreciseEqual(f(a, b), expected)
            self.assertPreciseEqual(f(b, a), expected)

        # non-int64 int * timedelta64
        check(TD(3), np.uint32(2), TD(6))
        # int * timedelta64
        check(TD(3), 2, TD(6))
        check(TD(3, 'ps'), 2, TD(6, 'ps'))
        check(TD('NaT', 'ps'), 2, TD('NaT', 'ps'))
        # float * timedelta64
        check(TD(7), 1.5, TD(10))
        check(TD(-7), 1.5, TD(-10))
        check(TD(7, 'ps'), -1.5, TD(-10, 'ps'))
        check(TD(-7), -1.5, TD(10))
        check(TD('NaT', 'ps'), -1.5, TD('NaT', 'ps'))
        check(TD(7, 'ps'), float('nan'), TD('NaT', 'ps'))
        # wraparound on overflow
        check(TD(2**62, 'ps'), 16, TD(0, 'ps'))

    def test_div(self):
        div = self.jit(div_usecase)
        floordiv = self.jit(floordiv_usecase)
        def check(a, b, expected):
            self.assertPreciseEqual(div(a, b), expected)
            self.assertPreciseEqual(floordiv(a, b), expected)

        # timedelta64 / non-int64 int
        check(TD(-3, 'ps'), np.uint32(2), TD(-1, 'ps'))
        # timedelta64 / int
        check(TD(3), 2, TD(1))
        check(TD(-3, 'ps'), 2, TD(-1, 'ps'))
        check(TD('NaT', 'ps'), 2, TD('NaT', 'ps'))
        check(TD(3, 'ps'), 0, TD('NaT', 'ps'))
        check(TD('NaT', 'ps'), 0, TD('NaT', 'ps'))
        # timedelta64 / float
        check(TD(7), 0.5, TD(14))
        check(TD(-7, 'ps'), 1.5, TD(-4, 'ps'))
        check(TD('NaT', 'ps'), 2.5, TD('NaT', 'ps'))
        check(TD(3, 'ps'), 0.0, TD('NaT', 'ps'))
        check(TD('NaT', 'ps'), 0.0, TD('NaT', 'ps'))
        check(TD(3, 'ps'), float('nan'), TD('NaT', 'ps'))
        check(TD('NaT', 'ps'), float('nan'), TD('NaT', 'ps'))

    def test_homogeneous_div(self):
        div = self.jit(div_usecase)
        def check(a, b, expected):
            self.assertPreciseEqual(div(a, b), expected)

        # timedelta64 / timedelta64
        check(TD(7), TD(3), 7. / 3.)
        check(TD(7, 'us'), TD(3, 'ms'), 7. / 3000.)
        check(TD(7, 'ms'), TD(3, 'us'), 7000. / 3.)
        check(TD(7), TD(0), float('+inf'))
        check(TD(-7), TD(0), float('-inf'))
        check(TD(0), TD(0), float('nan'))
        # NaTs
        check(TD('nat'), TD(3), float('nan'))
        check(TD(3), TD('nat'), float('nan'))
        check(TD('nat'), TD(0), float('nan'))
        # Cannot div months with days
        with self.assertRaises((TypeError, TypingError)):
            div(TD(1, 'M'), TD(1, 'D'))

    def test_eq_ne(self):
        eq = self.jit(eq_usecase)
        ne = self.jit(ne_usecase)
        def check(a, b, expected):
            expected_val = expected
            not_expected_val = not expected

            # all NaT comparisons are False, including NaT==NaT,
            # conversely != is True
            if np.isnat(a) or np.isnat(a):
                expected_val = False
                not_expected_val = True

            self.assertPreciseEqual(eq(a, b), expected_val)
            self.assertPreciseEqual(eq(b, a), expected_val)
            self.assertPreciseEqual(ne(a, b), not_expected_val)
            self.assertPreciseEqual(ne(b, a), not_expected_val)

        check(TD(1), TD(2), False)
        check(TD(1), TD(1), True)
        check(TD(1, 's'), TD(2, 's'), False)
        check(TD(1, 's'), TD(1, 's'), True)
        check(TD(2000, 's'), TD(2, 's'), False)
        check(TD(2000, 'ms'), TD(2, 's'), True)
        check(TD(1, 'Y'), TD(12, 'M'), True)
        # NaTs
        check(TD('Nat'), TD('Nat'), True)
        check(TD('Nat', 'ms'), TD('Nat', 's'), True)
        check(TD('Nat'), TD(1), False)
        # Incompatible units => timedeltas compare unequal
        if numpy_version < (1, 25):
            check(TD(1, 'Y'), TD(365, 'D'), False)
            check(TD(1, 'Y'), TD(366, 'D'), False)
            # ... except when both are NaT!
            check(TD('NaT', 'W'), TD('NaT', 'D'), True)
        else:
            # incompatible units raise
            # The exception is different depending on Python mode
            with self.assertRaises((TypeError, TypingError)):
                eq(TD(1, 'Y'), TD(365, 'D'))
            with self.assertRaises((TypeError, TypingError)):
                ne(TD(1, 'Y'), TD(365, 'D'))

    def test_lt_ge(self):
        lt = self.jit(lt_usecase)
        ge = self.jit(ge_usecase)
        def check(a, b, expected):
            expected_val = expected
            not_expected_val = not expected

            # since np 1.16 all NaT magnitude comparisons including equality
            # are False (as NaT == NaT is now False)
            if np.isnat(a) or np.isnat(a):
                expected_val = False
                not_expected_val = False

            self.assertPreciseEqual(lt(a, b), expected_val)
            self.assertPreciseEqual(ge(a, b), not_expected_val)

        check(TD(1), TD(2), True)
        check(TD(1), TD(1), False)
        check(TD(2), TD(1), False)
        check(TD(1, 's'), TD(2, 's'), True)
        check(TD(1, 's'), TD(1, 's'), False)
        check(TD(2, 's'), TD(1, 's'), False)
        check(TD(1, 'm'), TD(61, 's'), True)
        check(TD(1, 'm'), TD(60, 's'), False)
        # NaTs
        check(TD('Nat'), TD('Nat'), False)
        check(TD('Nat', 'ms'), TD('Nat', 's'), False)
        check(TD('Nat'), TD(-(2**63)+1), True)
        # Incompatible units => exception raised
        with self.assertRaises((TypeError, TypingError)):
            lt(TD(1, 'Y'), TD(365, 'D'))
        with self.assertRaises((TypeError, TypingError)):
            ge(TD(1, 'Y'), TD(365, 'D'))
        # ... even when both are NaT
        with self.assertRaises((TypeError, TypingError)):
            lt(TD('NaT', 'Y'), TD('NaT', 'D'))
        with self.assertRaises((TypeError, TypingError)):
            ge(TD('NaT', 'Y'), TD('NaT', 'D'))

    def test_le_gt(self):
        le = self.jit(le_usecase)
        gt = self.jit(gt_usecase)
        def check(a, b, expected):
            expected_val = expected
            not_expected_val = not expected

            # since np 1.16 all NaT magnitude comparisons including equality
            # are False (as NaT == NaT is now False)
            if np.isnat(a) or np.isnat(a):
                expected_val = False
                not_expected_val = False
            self.assertPreciseEqual(le(a, b), expected_val)
            self.assertPreciseEqual(gt(a, b), not_expected_val)

        check(TD(1), TD(2), True)
        check(TD(1), TD(1), True)
        check(TD(2), TD(1), False)
        check(TD(1, 's'), TD(2, 's'), True)
        check(TD(1, 's'), TD(1, 's'), True)
        check(TD(2, 's'), TD(1, 's'), False)
        check(TD(1, 'm'), TD(61, 's'), True)
        check(TD(1, 'm'), TD(60, 's'), True)
        check(TD(1, 'm'), TD(59, 's'), False)
        # NaTs
        check(TD('Nat'), TD('Nat'), True)
        check(TD('Nat', 'ms'), TD('Nat', 's'), True)
        check(TD('Nat'), TD(-(2**63)+1), True)
        # Incompatible units => exception raised
        with self.assertRaises((TypeError, TypingError)):
            le(TD(1, 'Y'), TD(365, 'D'))
        with self.assertRaises((TypeError, TypingError)):
            gt(TD(1, 'Y'), TD(365, 'D'))
        # ... even when both are NaT
        with self.assertRaises((TypeError, TypingError)):
            le(TD('NaT', 'Y'), TD('NaT', 'D'))
        with self.assertRaises((TypeError, TypingError)):
            gt(TD('NaT', 'Y'), TD('NaT', 'D'))

    def test_pos(self):
        pos = self.jit(pos_usecase)
        def check(a):
            self.assertPreciseEqual(pos(a), +a)

        check(TD(3))
        check(TD(-4))
        check(TD(3, 'ms'))
        check(TD(-4, 'ms'))
        check(TD('NaT'))
        check(TD('NaT', 'ms'))

    def test_neg(self):
        neg = self.jit(neg_usecase)
        def check(a):
            self.assertPreciseEqual(neg(a), -a)

        check(TD(3))
        check(TD(-4))
        check(TD(3, 'ms'))
        check(TD(-4, 'ms'))
        check(TD('NaT'))
        check(TD('NaT', 'ms'))

    def test_abs(self):
        f = self.jit(abs_usecase)
        def check(a):
            self.assertPreciseEqual(f(a), abs(a))

        check(TD(3))
        check(TD(-4))
        check(TD(3, 'ms'))
        check(TD(-4, 'ms'))
        check(TD('NaT'))
        check(TD('NaT', 'ms'))

    def test_hash(self):
        f = self.jit(hash_usecase)
        def check(a):
            if numpy_version >= (2, 2):
                # Generic timedeltas (those without a unit)
                # are no longer hashable beyond NumPy 2.2
                # Non-generic timedeltas will have dtype name
                # as timedelta64[<unit>]
                if a.dtype.name == 'timedelta64':
                    return

                # If the function is not being compiled in objmode
                # then the hash should be equal to the hash of the
                # integer representation of the timedelta
                if self.jitargs.get('nopython', False):
                    self.assertPreciseEqual(f(a), a.astype(int))
                else:
                    self.assertPreciseEqual(f(a), hash(a))
            else:
                self.assertPreciseEqual(f(a), hash(a))

        TD_CASES = ((3,), (-4,), (3, 'ms'), (-4, 'ms'), (27, 'D'),
                    (2, 'D'), (2, 'W'), (2, 'Y'), (3, 'W'),
                    (365, 'D'), (10000, 'D'), (-10000, 'D'),
                    ('NaT',), ('NaT', 'ms'), ('NaT', 'D'), (-1,))
        DT_CASES = (('2014',), ('2016',), ('2000',), ('2014-02',),
                    ('2014-03',), ('2014-04',), ('2016-02',), ('2000-12-31',),
                    ('2014-01-16',), ('2014-01-05',), ('2014-01-07',),
                    ('2014-01-06',), ('2014-02-02',), ('2014-02-27',),
                    ('2014-02-16',), ('2014-03-01',), ('2000-01-01T01:02:03.002Z',),
                    ('2000-01-01T01:02:03Z',), ('NaT',))

        for case, typ in zip(TD_CASES + DT_CASES,
                             (TD,) * len(TD_CASES) + (DT,) * len(TD_CASES)):
            check(typ(*case))

        if numpy_version >= (2, 2):
            with self.assertRaises(ValueError) as raises:
                f(TD(3))
            self.assertIn("Can't hash generic timedelta64", str(raises.exception))

    def _test_min_max(self, usecase):
        f = self.jit(usecase)
        def check(a, b):
            self.assertPreciseEqual(f(a, b), usecase(a, b))

        for cases in (
            (TD(0), TD(1), TD(2), TD('NaT')),
            (TD(0, 's'), TD(1, 's'), TD(2, 's'), TD('NaT', 's')),
        ):
            for a, b in itertools.product(cases, cases):
                check(a, b)

    def test_min(self):
        self._test_min_max(min_usecase)

    def test_max(self):
        self._test_min_max(max_usecase)


class TestTimedeltaArithmeticNoPython(TestTimedeltaArithmetic):

    jitargs = dict(nopython=True)

    def test_int_cast(self):
        f = self.jit(int_cast_usecase)
        def check(a):
            self.assertPreciseEqual(f(a), int(a))

        for (delta, unit) in ((3, 'ns'), (-4, 'ns'), (30000, 'ns'),
                        (-40000000, 'ns'), (1, 'Y')):
            check(TD(delta, unit).astype('timedelta64[ns]'))

        for time in ('2014', '2016', '2000', '2014-02', '2014-03', '2014-04',
                     '2016-02', '2000-12-31', '2014-01-16', '2014-01-05',
                     '2014-01-07', '2014-01-06', '2014-02-02', '2014-02-27',
                     '2014-02-16', '2014-03-01', '2000-01-01T01:02:03.002Z',
                     '2000-01-01T01:02:03Z'):
            check(DT(time).astype('datetime64[ns]'))

        with self.assertRaises(TypingError, msg=('Only datetime64[ns] can be ' +
                                                 'converted, but got ' +
                                                 'datetime64[y]')):
            f(DT('2014'))


class TestDatetimeArithmetic(TestCase):

    jitargs = dict(forceobj=True)

    def jit(self, pyfunc):
        return jit(**self.jitargs)(pyfunc)

    @contextlib.contextmanager
    def silence_numpy_warnings(self):
        # Numpy can raise warnings when combining e.g. a generic timedelta64
        # with a non-generic datetime64.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',
                                    message='Implicitly casting between incompatible kinds',
                                    category=DeprecationWarning)
            yield

    def test_add_sub_timedelta(self):
        """
        Test `datetime64 + timedelta64` and `datetime64 - timedelta64`.
        """
        add = self.jit(add_usecase)
        sub = self.jit(sub_usecase)
        def check(a, b, expected):
            with self.silence_numpy_warnings():
                self.assertPreciseEqual(add(a, b), expected, (a, b))
                self.assertPreciseEqual(add(b, a), expected, (a, b))
                self.assertPreciseEqual(sub(a, -b), expected, (a, b))
                # Did we get it right?
                self.assertPreciseEqual(a + b, expected)

        # Y + ...
        check(DT('2014'), TD(2, 'Y'), DT('2016'))
        check(DT('2014'), TD(2, 'M'), DT('2014-03'))
        check(DT('2014'), TD(3, 'W'), DT('2014-01-16', 'W'))
        check(DT('2014'), TD(4, 'D'), DT('2014-01-05'))
        check(DT('2000'), TD(365, 'D'), DT('2000-12-31'))
        # M + ...
        check(DT('2014-02'), TD(2, 'Y'), DT('2016-02'))
        check(DT('2014-02'), TD(2, 'M'), DT('2014-04'))
        check(DT('2014-02'), TD(2, 'D'), DT('2014-02-03'))
        # W + ...
        check(DT('2014-01-07', 'W'), TD(2, 'W'), DT('2014-01-16', 'W'))
        # D + ...
        check(DT('2014-02-02'), TD(27, 'D'), DT('2014-03-01'))
        check(DT('2012-02-02'), TD(27, 'D'), DT('2012-02-29'))
        check(DT('2012-02-02'), TD(2, 'W'), DT('2012-02-16'))
        # s + ...
        check(DT('2000-01-01T01:02:03Z'), TD(2, 'h'), DT('2000-01-01T03:02:03Z'))
        check(DT('2000-01-01T01:02:03Z'), TD(2, 'ms'), DT('2000-01-01T01:02:03.002Z'))
        # More thorough checking with leap years and faraway years
        for dt_str in ('600', '601', '604', '801',
                       '1900', '1904', '2200', '2300', '2304',
                       '2400', '6001'):
            for dt_suffix in ('', '-01', '-12'):
                dt = DT(dt_str + dt_suffix)
                for td in [TD(2, 'D'), TD(2, 'W'),
                           TD(100, 'D'), TD(10000, 'D'),
                           TD(-100, 'D'), TD(-10000, 'D'),
                           TD(100, 'W'), TD(10000, 'W'),
                           TD(-100, 'W'), TD(-10000, 'W'),
                           TD(100, 'M'), TD(10000, 'M'),
                           TD(-100, 'M'), TD(-10000, 'M')]:
                    self.assertEqual(add(dt, td), dt + td, (dt, td))
                    self.assertEqual(add(td, dt), dt + td, (dt, td))
                    self.assertEqual(sub(dt, -td), dt + td, (dt, td))

        # NaTs
        check(DT('NaT'), TD(2), DT('NaT'))
        check(DT('NaT', 's'), TD(2, 'h'), DT('NaT', 's'))
        check(DT('NaT', 's'), TD(2, 'ms'), DT('NaT', 'ms'))
        check(DT('2014'), TD('NaT', 'W'), DT('NaT', 'W'))
        check(DT('2014-01-01'), TD('NaT', 'W'), DT('NaT', 'D'))
        check(DT('NaT', 's'), TD('NaT', 'ms'), DT('NaT', 'ms'))

        # Cannot add datetime days and timedelta months or years
        for f in (add, sub):
            with self.assertRaises((TypeError, TypingError)):
                f(DT(1, '2014-01-01'), TD(1, 'Y'))
            with self.assertRaises((TypeError, TypingError)):
                f(DT(1, '2014-01-01'), TD(1, 'M'))

    def datetime_samples(self):
        dt_years = ['600', '601', '604', '1968', '1969', '1973',
                   '2000', '2004', '2005', '2100', '2400', '2401']
        dt_suffixes = ['', '-01', '-12', '-02-28', '-12-31',
                       '-01-05T12:30:56Z', '-01-05T12:30:56.008Z']
        dts = [DT(a + b) for (a, b) in itertools.product(dt_years, dt_suffixes)]
        dts += [DT(s, 'W') for s in dt_years]
        return dts

    def test_datetime_difference(self):
        """
        Test `datetime64 - datetime64`.
        """
        sub = self.jit(sub_usecase)
        def check(a, b, expected=None):
            with self.silence_numpy_warnings():
                self.assertPreciseEqual(sub(a, b), a - b, (a, b))
                self.assertPreciseEqual(sub(b, a), b - a, (a, b))
                # Did we get it right?
                self.assertPreciseEqual(a - b, expected)

        check(DT('2014'), DT('2017'), TD(-3, 'Y'))
        check(DT('2014-02'), DT('2017-01'), TD(-35, 'M'))
        check(DT('2014-02-28'), DT('2015-03-01'), TD(-366, 'D'))
        # NaTs
        check(DT('NaT', 'M'), DT('2000'), TD('NaT', 'M'))
        check(DT('NaT', 'M'), DT('2000-01-01'), TD('NaT', 'D'))
        check(DT('NaT'), DT('NaT'), TD('NaT'))
        # Test many more values
        with self.silence_numpy_warnings():
            dts = self.datetime_samples()
            for a, b in itertools.product(dts, dts):
                if (not npdatetime_helpers.same_kind(value_unit(a), value_unit(b))):
                    continue
                self.assertPreciseEqual(sub(a, b), a - b, (a, b))

    def test_comparisons(self):
        # Test all datetime comparisons all at once
        eq = self.jit(eq_usecase)
        ne = self.jit(ne_usecase)
        lt = self.jit(lt_usecase)
        le = self.jit(le_usecase)
        gt = self.jit(gt_usecase)
        ge = self.jit(ge_usecase)

        def check_eq(a, b, expected):
            expected_val = expected
            not_expected_val = not expected

            # since np 1.16 all NaT comparisons bar != are False, including
            # NaT==NaT
            if np.isnat(a) or np.isnat(b):
                expected_val = False
                not_expected_val = True
                self.assertFalse(le(a, b), (a, b))
                self.assertFalse(ge(a, b), (a, b))
                self.assertFalse(le(b, a), (a, b))
                self.assertFalse(ge(b, a), (a, b))
                self.assertFalse(lt(a, b), (a, b))
                self.assertFalse(gt(a, b), (a, b))
                self.assertFalse(lt(b, a), (a, b))
                self.assertFalse(gt(b, a), (a, b))

            with self.silence_numpy_warnings():
                self.assertPreciseEqual(eq(a, b), expected_val, (a, b, expected))
                self.assertPreciseEqual(eq(b, a), expected_val, (a, b, expected))
                self.assertPreciseEqual(ne(a, b), not_expected_val, (a, b, expected))
                self.assertPreciseEqual(ne(b, a), not_expected_val, (a, b, expected))
                if expected_val:
                    # If equal, then equal-ordered comparisons are true
                    self.assertTrue(le(a, b), (a, b))
                    self.assertTrue(ge(a, b), (a, b))
                    self.assertTrue(le(b, a), (a, b))
                    self.assertTrue(ge(b, a), (a, b))
                    # and strictly ordered comparisons are false
                    self.assertFalse(lt(a, b), (a, b))
                    self.assertFalse(gt(a, b), (a, b))
                    self.assertFalse(lt(b, a), (a, b))
                    self.assertFalse(gt(b, a), (a, b))
                # Did we get it right?
                self.assertPreciseEqual(a == b, expected_val)

        def check_lt(a, b, expected):
            expected_val = expected
            not_expected_val = not expected

            # since np 1.16 all NaT magnitude comparisons including equality
            # are False (as NaT == NaT is now False)
            if np.isnat(a) or np.isnat(b):
                expected_val = False
                not_expected_val = False

            with self.silence_numpy_warnings():
                lt = self.jit(lt_usecase)
                self.assertPreciseEqual(lt(a, b), expected_val, (a, b, expected))
                self.assertPreciseEqual(gt(b, a), expected_val, (a, b, expected))
                self.assertPreciseEqual(ge(a, b), not_expected_val, (a, b, expected))
                self.assertPreciseEqual(le(b, a), not_expected_val, (a, b, expected))
                if expected_val:
                    # If true, then values are not equal
                    check_eq(a, b, False)
                # Did we get it right?
                self.assertPreciseEqual(a < b, expected_val)

        check_eq(DT('2014'), DT('2017'), False)
        check_eq(DT('2014'), DT('2014-01'), True)
        check_eq(DT('2014'), DT('2014-01-01'), True)
        check_eq(DT('2014'), DT('2014-01-01', 'W'), True)
        check_eq(DT('2014-01'), DT('2014-01-01', 'W'), True)
        # Yes, it's not transitive
        check_eq(DT('2014-01-01'), DT('2014-01-01', 'W'), False)
        check_eq(DT('2014-01-02'), DT('2014-01-06', 'W'), True)
        # with times
        check_eq(DT('2014-01-01T00:01:00Z', 's'),
              DT('2014-01-01T00:01Z', 'm'), True)
        check_eq(DT('2014-01-01T00:01:01Z', 's'),
              DT('2014-01-01T00:01Z', 'm'), False)
        # NaTs
        check_lt(DT('NaT', 'Y'), DT('2017'), True)
        check_eq(DT('NaT'), DT('NaT'), True)

        # Check comparison between various units
        dts = self.datetime_samples()
        for a in dts:
            # Take a number of smaller units
            a_unit = a.dtype.str.split('[')[1][:-1]
            i = all_units.index(a_unit)
            units = all_units[i:i+6]
            for unit in units:
                # Force conversion
                b = a.astype('M8[%s]' % unit)
                if (not npdatetime_helpers.same_kind(value_unit(a),
                                                     value_unit(b))):
                    continue
                check_eq(a, b, True)
                check_lt(a, b + np.timedelta64(1, unit), True)
                check_lt(b - np.timedelta64(1, unit), a, True)

    def _test_min_max(self, usecase):
        f = self.jit(usecase)
        def check(a, b):
            self.assertPreciseEqual(f(a, b), usecase(a, b))

        for cases in (
            (DT(0, 'ns'), DT(1, 'ns'), DT(2, 'ns'), DT('NaT', 'ns')),
            (DT(0, 's'), DT(1, 's'), DT(2, 's'), DT('NaT', 's')),
        ):
            for a, b in itertools.product(cases, cases):
                check(a, b)

    def test_min(self):
        self._test_min_max(min_usecase)

    def test_max(self):
        self._test_min_max(max_usecase)

class TestDatetimeArithmeticNoPython(TestDatetimeArithmetic):

    jitargs = dict(nopython=True)


class TestMetadataScalingFactor(TestCase):
    """
    Tests than non-1 scaling factors are not supported in datetime64
    and timedelta64 dtypes.
    """

    def test_datetime(self, jitargs={'forceobj':True}):
        eq = jit(**jitargs)(eq_usecase)
        self.assertTrue(eq(DT('2014', '10Y'), DT('2010')))

    def test_datetime_npm(self):
        with self.assertTypingError():
            self.test_datetime(jitargs={'nopython':True})

    def test_timedelta(self, jitargs={'forceobj':True}):
        eq = jit(**jitargs)(eq_usecase)
        self.assertTrue(eq(TD(2, '10Y'), TD(20, 'Y')))

    def test_timedelta_npm(self):
        with self.assertTypingError():
            self.test_timedelta(jitargs={'nopython':True})


class TestDatetimeDeltaOps(TestCase):
    def test_div(self):
        """
        Test the division of a timedelta by numeric types
        """
        def arr_div(a, b):
            return a / b

        py_func = arr_div
        cfunc = njit(arr_div)
        test_cases = [
            (np.ones(3, TIMEDELTA_M), np.ones(3, TIMEDELTA_M)),
            (np.ones(3, TIMEDELTA_M), np.ones(3, TIMEDELTA_Y)),
            (np.ones(3, TIMEDELTA_Y), np.ones(3, TIMEDELTA_M)),
            (np.ones(3, TIMEDELTA_Y), np.ones(3, TIMEDELTA_Y)),
            (np.ones(3, TIMEDELTA_M), 1),
            (np.ones(3, TIMEDELTA_M), np.ones(3, np.int64)),
            (np.ones(3, TIMEDELTA_M), np.ones(3, np.float64)),
        ]
        for a, b in test_cases:
            self.assertTrue(np.array_equal(py_func(a, b), cfunc(a, b)))


class TestDatetimeArrayOps(TestCase):

    def _test_td_add_or_sub(self, operation, parallel):
        """
        Test the addition/subtraction of a datetime array with a timedelta type
        """
        def impl(a, b):
            return operation(a, b)

        arr_one = np.array([
                        np.datetime64("2011-01-01"),
                        np.datetime64("1971-02-02"),
                        np.datetime64("2021-03-03"),
                        np.datetime64("2004-12-07"),
                    ], dtype="datetime64[ns]")
        arr_two = np.array([
                        np.datetime64("2011-01-01"),
                        np.datetime64("1971-02-02"),
                        np.datetime64("2021-03-03"),
                        np.datetime64("2004-12-07"),
                    ], dtype="datetime64[D]")
        py_func = impl
        cfunc = njit(parallel=parallel)(impl)
        test_cases = [
            (arr_one, np.timedelta64(1000)),
            (arr_two, np.timedelta64(1000)),
            (arr_one, np.timedelta64(-54557)),
            (arr_two, np.timedelta64(-54557)),
        ]
        # np.add is commutative so test the reversed order
        if operation is np.add:
            test_cases.extend([
                (np.timedelta64(1000), arr_one),
                (np.timedelta64(1000), arr_two),
                (np.timedelta64(-54557), arr_one),
                (np.timedelta64(-54557), arr_two),
            ])
        for a, b in test_cases:
            self.assertTrue(np.array_equal(py_func(a, b), cfunc(a, b)))

    def test_add_td(self):
        self._test_td_add_or_sub(np.add, False)

    @skip_parfors_unsupported
    def test_add_td_parallel(self):
        self._test_td_add_or_sub(np.add, True)

    def test_sub_td(self):
        self._test_td_add_or_sub(np.subtract, False)

    @skip_parfors_unsupported
    def test_sub_td_parallel(self):
        self._test_td_add_or_sub(np.subtract, True)

    def _test_add_sub_td_no_match(self, operation):
        """
        Tests that attempting to add/sub a datetime64 and timedelta64
        with types that cannot be cast raises a reasonable exception.
        """
        @njit
        def impl(a, b):
            return operation(a, b)

        fname = operation.__name__
        expected = re.escape((f"ufunc '{fname}' is not supported between "
                              "datetime64[ns] and timedelta64[M]"))
        with self.assertRaisesRegex((TypingError, TypeError), expected):
            impl(
                np.array([np.datetime64("2011-01-01"),],
                         dtype="datetime64[ns]"),
                np.timedelta64(1000,'M')
            )

    def test_add_td_no_match(self):
        self._test_add_sub_td_no_match(np.add)

    def test_sub_td_no_match(self):
        self._test_add_sub_td_no_match(np.subtract)

    def _get_testcases(self):
        test_cases = [
            np.array([
                DT(0, "ns"),
                DT(1, "ns"),
                DT(2, "ns"),
                DT(3, "ns"),
            ]),
            np.array([
                DT("2011-01-01", "ns"),
                DT("1971-02-02", "ns"),
                DT("1900-01-01", "ns"),
                DT("2021-03-03", "ns"),
                DT("2004-12-07", "ns"),
            ]),
            np.array([
                DT("2011-01-01", "D"),
                DT("1971-02-02", "D"),
                DT("1900-01-01", "D"),
                DT("2021-03-03", "D"),
                DT("2004-12-07", "D"),
            ]),
            np.array([
                DT("2011-01-01", "ns"),
                DT("1971-02-02", "ns"),
                DT("1900-01-01", "ns"),
                DT("2021-03-03", "ns"),
                DT("2004-12-07", "ns"),
                DT("NaT", "ns"),
            ]),
            np.array([
                DT("NaT", "ns"),
                DT("2011-01-01", "ns"),
                DT("1971-02-02", "ns"),
                DT("1900-01-01", "ns"),
                DT("2021-03-03", "ns"),
                DT("2004-12-07", "ns"),
            ]),
            np.array([
                DT("1971-02-02", "ns"),
                DT("NaT", "ns"),
            ]),
            np.array([
                DT("NaT", "ns"),
                DT("NaT", "ns"),
                DT("NaT", "ns"),
            ]),
            np.array([
                TD(1, "ns"),
                TD(2, "ns"),
                TD(3, "ns"),
                TD(4, "ns"),
            ]),
            np.array([
                TD(1, "D"),
                TD(2, "D"),
                TD(3, "D"),
                TD(4, "D"),
            ]),
            np.array([
                TD("NaT", "ns"),
                TD(1, "ns"),
                TD(2, "ns"),
                TD(3, "ns"),
                TD(4, "ns"),
            ]),
            np.array([
                TD(1, "ns"),
                TD(2, "ns"),
                TD(3, "ns"),
                TD(4, "ns"),
                TD("NaT", "ns"),
            ]),
            np.array([
                TD("NaT", "ns"),
            ]),
            np.array([
                TD("NaT", "ns"),
                TD("NaT", "ns"),
                TD("NaT", "ns"),
            ]),
        ]
        return test_cases

    def _test_min_max(self, operation, parallel, method):
        if method:
            if operation is np.min:
                def impl(arr):
                    return arr.min()
            else:
                def impl(arr):
                    return arr.max()
        else:
            def impl(arr):
                return operation(arr)

        py_func = impl
        cfunc = njit(parallel=parallel)(impl)

        test_cases = self._get_testcases()
        for arr in test_cases:
            py_res = py_func(arr)
            c_res = cfunc(arr)
            if np.isnat(py_res) or np.isnat(c_res):
                self.assertTrue(np.isnat(py_res))
                self.assertTrue(np.isnat(c_res))
            else:
                self.assertEqual(py_res, c_res)

    def test_min_func(self):
        self._test_min_max(min, False, False)

    def test_np_min_func(self):
        self._test_min_max(np.min, False, False)

    def test_min_method(self):
        self._test_min_max(np.min, False, True)

    def test_max_func(self):
        self._test_min_max(max, False, False)

    def test_np_max_func(self):
        self._test_min_max(np.max, False, False)

    def test_max_method(self):
        self._test_min_max(np.max, False, True)

    @skip_parfors_unsupported
    def test_min_func_parallel(self):
        self._test_min_max(np.min, True, False)

    @skip_parfors_unsupported
    def test_min_method_parallel(self):
        self._test_min_max(np.min, True, True)

    @skip_parfors_unsupported
    def test_max_func_parallel(self):
        self._test_min_max(np.max, True, False)

    @skip_parfors_unsupported
    def test_max_method_parallel(self):
        self._test_min_max(np.max, True, True)

    def test_searchsorted_datetime(self):
        from .test_np_functions import (
            searchsorted, searchsorted_left, searchsorted_right,
        )
        pyfunc_list = [searchsorted, searchsorted_left, searchsorted_right]
        cfunc_list = [jit(fn) for fn in pyfunc_list]

        def check(pyfunc, cfunc, a, v):
            expected = pyfunc(a, v)
            got = cfunc(a, v)
            self.assertPreciseEqual(expected, got)

        cases = self._get_testcases()
        for pyfunc, cfunc in zip(pyfunc_list, cfunc_list):
            for arr in cases:
                arr = np.sort(arr)
                for n in range(1, min(3, arr.size) + 1):
                    idx = np.random.randint(0, arr.size, n)
                    vs = arr[idx]
                    if n == 1:
                        [v] = vs
                        check(pyfunc, cfunc, arr, v)
                    check(pyfunc, cfunc, arr, vs)



class TestDatetimeTypeOps(TestCase):
    def test_isinstance_datetime(self):
        @njit
        def is_complex(a):
            return isinstance(a, complex)
        @njit
        def is_datetime(a):
            return isinstance(a, np.datetime64)
        @njit
        def is_timedelta(a):
            return isinstance(a, np.timedelta64)

        dt_a = np.datetime64(1, 'ns')
        dt_b = np.datetime64(2, 'ns')
        td_c = dt_b - dt_a

        def check(jit_func, x):
            with self.subTest(f'{jit_func.__name__}({type(x).__name__})'):
                got = jit_func(x)
                expect = jit_func.py_func(x)
                self.assertEqual(got, expect)

        fns = [
            is_complex,
            is_datetime,
            is_timedelta,
        ]
        args = [
            dt_a,
            dt_b,
            td_c,
        ]
        for fn, arg in itertools.product(fns, args):
            check(fn, arg)


if __name__ == '__main__':
    unittest.main()
