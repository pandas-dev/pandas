"""
Testing numba implementation of the numba dictionary.

The tests here only check that the numba typing and codegen are working
correctly.  Detailed testing of the underlying dictionary operations is done
in test_dictimpl.py.
"""

import sys
import warnings

import numpy as np

from numba import njit, literally
from numba import int32, int64, float32, float64
from numba import typeof
from numba.typed import Dict, dictobject, List
from numba.typed.typedobjectutils import _sentry_safe_cast
from numba.core.errors import TypingError
from numba.core import types
from numba.tests.support import (TestCase, MemoryLeakMixin, unittest,
                                 override_config, forbid_codegen)
from numba.experimental import jitclass
from numba.extending import overload


class TestDictObject(MemoryLeakMixin, TestCase):
    def test_dict_bool(self):
        """
        Exercise bool(dict)
        """
        @njit
        def foo(n):
            d = dictobject.new_dict(int32, float32)
            for i in range(n):
                d[i] = i + 1
            return bool(d)

        # Insert nothing
        self.assertEqual(foo(n=0), False)
        # Insert 1 entry
        self.assertEqual(foo(n=1), True)
        # Insert 2 entries
        self.assertEqual(foo(n=2), True)
        # Insert 100 entries
        self.assertEqual(foo(n=100), True)

    def test_dict_create(self):
        """
        Exercise dictionary creation, insertion and len
        """
        @njit
        def foo(n):
            d = dictobject.new_dict(int32, float32)
            for i in range(n):
                d[i] = i + 1
            return len(d)

        # Insert nothing
        self.assertEqual(foo(n=0), 0)
        # Insert 1 entry
        self.assertEqual(foo(n=1), 1)
        # Insert 2 entries
        self.assertEqual(foo(n=2), 2)
        # Insert 100 entries
        self.assertEqual(foo(n=100), 100)

    def test_dict_get(self):
        """
        Exercise dictionary creation, insertion and get
        """
        @njit
        def foo(n, targets):
            d = dictobject.new_dict(int32, float64)
            # insertion loop
            for i in range(n):
                d[i] = i
            # retrieval loop
            output = []
            for t in targets:
                output.append(d.get(t))
            return output

        self.assertEqual(foo(5, [0, 1, 9]), [0, 1, None])
        self.assertEqual(foo(10, [0, 1, 9]), [0, 1, 9])
        self.assertEqual(foo(10, [-1, 9, 1]), [None, 9, 1])

    def test_dict_get_with_default(self):
        """
        Exercise dict.get(k, d) where d is set
        """
        @njit
        def foo(n, target, default):
            d = dictobject.new_dict(int32, float64)
            # insertion loop
            for i in range(n):
                d[i] = i
            # retrieval loop
            return d.get(target, default)

        self.assertEqual(foo(5, 3, -1), 3)
        self.assertEqual(foo(5, 5, -1), -1)

    def test_dict_getitem(self):
        """
        Exercise dictionary __getitem__
        """
        @njit
        def foo(keys, vals, target):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v

            # lookup
            return d[target]

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]
        self.assertEqual(foo(keys, vals, 1), 0.1)
        self.assertEqual(foo(keys, vals, 2), 0.2)
        self.assertEqual(foo(keys, vals, 3), 0.3)
        # check no leak so far
        self.assert_no_memory_leak()
        # disable leak check for exception test
        self.disable_leak_check()
        with self.assertRaisesRegex(KeyError, "0"):
            foo(keys, vals, 0)
        with self.assertRaisesRegex(KeyError, "4"):
            foo(keys, vals, 4)

    def test_dict_popitem(self):
        """
        Exercise dictionary .popitem
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v

            # popitem
            return d.popitem()

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]
        for i in range(1, len(keys)):
            self.assertEqual(
                foo(keys[:i], vals[:i]),
                (keys[i - 1], vals[i - 1]),
            )

    def test_dict_popitem_many(self):
        """
        Exercise dictionary .popitem
        """

        @njit
        def core(d, npop):
            # popitem
            keysum, valsum = 0, 0
            for _ in range(npop):
                k, v = d.popitem()
                keysum += k
                valsum -= v
            return keysum, valsum

        @njit
        def foo(keys, vals, npop):
            d = dictobject.new_dict(int32, int32)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v

            return core(d, npop)

        keys = [1, 2, 3]
        vals = [10, 20, 30]

        for i in range(len(keys)):
            self.assertEqual(
                foo(keys, vals, npop=3),
                core.py_func(dict(zip(keys, vals)), npop=3),
            )

        # check no leak so far
        self.assert_no_memory_leak()
        # disable leak check for exception test
        self.disable_leak_check()

        with self.assertRaises(KeyError):
            foo(keys, vals, npop=4)

    def test_dict_pop(self):
        """
        Exercise dictionary .pop
        """
        @njit
        def foo(keys, vals, target):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v

            # popitem
            return d.pop(target, None), len(d)

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]

        self.assertEqual(foo(keys, vals, 1), (0.1, 2))
        self.assertEqual(foo(keys, vals, 2), (0.2, 2))
        self.assertEqual(foo(keys, vals, 3), (0.3, 2))
        self.assertEqual(foo(keys, vals, 0), (None, 3))

        # check no leak so far
        self.assert_no_memory_leak()
        # disable leak check for exception test
        self.disable_leak_check()

        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            # popitem
            return d.pop(0)

        with self.assertRaises(KeyError):
            foo()

    def test_dict_pop_many(self):
        """
        Exercise dictionary .pop
        """

        @njit
        def core(d, pops):
            total = 0
            for k in pops:
                total += k + d.pop(k, 0.123) + len(d)
                total *= 2
            return total

        @njit
        def foo(keys, vals, pops):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            # popitem
            return core(d, pops)

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]
        pops = [2, 3, 3, 1, 0, 2, 1, 0, -1]

        self.assertEqual(
            foo(keys, vals, pops),
            core.py_func(dict(zip(keys, vals)), pops),
        )

    def test_dict_delitem(self):
        @njit
        def foo(keys, vals, target):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            del d[target]
            return len(d), d.get(target)

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]
        self.assertEqual(foo(keys, vals, 1), (2, None))
        self.assertEqual(foo(keys, vals, 2), (2, None))
        self.assertEqual(foo(keys, vals, 3), (2, None))
        # check no leak so far
        self.assert_no_memory_leak()
        # disable leak check for exception test
        self.disable_leak_check()
        with self.assertRaises(KeyError):
            foo(keys, vals, 0)

    def test_dict_clear(self):
        """
        Exercise dict.clear
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            b4 = len(d)
            # clear
            d.clear()
            return b4, len(d)

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]
        self.assertEqual(foo(keys, vals), (3, 0))

    def test_dict_items(self):
        """
        Exercise dict.items
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            out = []
            for kv in d.items():
                out.append(kv)
            return out

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]

        self.assertEqual(
            foo(keys, vals),
            list(zip(keys, vals)),
        )

        # Test .items() on empty dict
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            out = []
            for kv in d.items():
                out.append(kv)
            return out

        self.assertEqual(foo(), [])

    def test_dict_keys(self):
        """
        Exercise dict.keys
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            out = []
            for k in d.keys():
                out.append(k)
            return out

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]

        self.assertEqual(
            foo(keys, vals),
            keys,
        )

    def test_dict_keys_len(self):
        """
        Exercise len(dict.keys())
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            return len(d.keys())

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]

        self.assertEqual(
            foo(keys, vals),
            len(keys),
        )

    def test_dict_values(self):
        """
        Exercise dict.values
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            out = []
            for v in d.values():
                out.append(v)
            return out

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]

        self.assertEqual(
            foo(keys, vals),
            vals,
        )

    def test_dict_values_len(self):
        """
        Exercise len(dict.values())
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            return len(d.values())

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]

        self.assertEqual(
            foo(keys, vals),
            len(vals),
        )

    def test_dict_items_len(self):
        """
        Exercise len(dict.items())
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            return len(d.items())

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]
        self.assertPreciseEqual(
            foo(keys, vals),
            len(vals),
        )

    def test_dict_iter(self):
        """
        Exercise iter(dict)
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            out = []
            for k in d:
                out.append(k)
            return out

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]

        self.assertEqual(
            foo(keys, vals),
            [1, 2, 3]
        )

    def test_dict_contains(self):
        """
        Exercise operator.contains
        """
        @njit
        def foo(keys, vals, checklist):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            out = []
            for k in checklist:
                out.append(k in d)
            return out

        keys = [1, 2, 3]
        vals = [0.1, 0.2, 0.3]

        self.assertEqual(
            foo(keys, vals, [2, 3, 4, 1, 0]),
            [True, True, False, True, False],
        )

    def test_dict_copy(self):
        """
        Exercise dict.copy
        """
        @njit
        def foo(keys, vals):
            d = dictobject.new_dict(int32, float64)
            # insertion
            for k, v in zip(keys, vals):
                d[k] = v
            return list(d.copy().items())

        keys = list(range(20))
        vals = [x + i / 100 for i, x in enumerate(keys)]
        out = foo(keys, vals)
        self.assertEqual(out, list(zip(keys, vals)))

    def test_dict_setdefault(self):
        """
        Exercise dict.setdefault
        """
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            d.setdefault(1, 1.2)  # used because key is not in
            a = d.get(1)
            d[1] = 2.3
            b = d.get(1)
            d[2] = 3.4
            d.setdefault(2, 4.5)  # not used because key is in
            c = d.get(2)
            return a, b, c

        self.assertEqual(foo(), (1.2, 2.3, 3.4))

    def test_dict_equality(self):
        """
        Exercise dict.__eq__ and .__ne__
        """
        @njit
        def foo(na, nb, fa, fb):
            da = dictobject.new_dict(int32, float64)
            db = dictobject.new_dict(int32, float64)
            for i in range(na):
                da[i] = i * fa
            for i in range(nb):
                db[i] = i * fb
            return da == db, da != db

        # Same keys and values
        self.assertEqual(foo(10, 10, 3, 3), (True, False))
        # Same keys and diff values
        self.assertEqual(foo(10, 10, 3, 3.1), (False, True))
        # LHS has more keys
        self.assertEqual(foo(11, 10, 3, 3), (False, True))
        # RHS has more keys
        self.assertEqual(foo(10, 11, 3, 3), (False, True))

    def test_dict_equality_more(self):
        """
        Exercise dict.__eq__
        """
        @njit
        def foo(ak, av, bk, bv):
            # The key-value types are different in the two dictionaries
            da = dictobject.new_dict(int32, float64)
            db = dictobject.new_dict(int64, float32)
            for i in range(len(ak)):
                da[ak[i]] = av[i]
            for i in range(len(bk)):
                db[bk[i]] = bv[i]
            return da == db

        # Simple equal case
        ak = [1, 2, 3]
        av = [2, 3, 4]
        bk = [1, 2, 3]
        bv = [2, 3, 4]
        self.assertTrue(foo(ak, av, bk, bv))

        # Equal with replacement
        ak = [1, 2, 3]
        av = [2, 3, 4]
        bk = [1, 2, 2, 3]
        bv = [2, 1, 3, 4]
        self.assertTrue(foo(ak, av, bk, bv))

        # Diff values
        ak = [1, 2, 3]
        av = [2, 3, 4]
        bk = [1, 2, 3]
        bv = [2, 1, 4]
        self.assertFalse(foo(ak, av, bk, bv))

        # Diff keys
        ak = [0, 2, 3]
        av = [2, 3, 4]
        bk = [1, 2, 3]
        bv = [2, 3, 4]
        self.assertFalse(foo(ak, av, bk, bv))

    def test_dict_equality_diff_type(self):
        """
        Exercise dict.__eq__
        """
        @njit
        def foo(na, b):
            da = dictobject.new_dict(int32, float64)
            for i in range(na):
                da[i] = i
            return da == b

        # dict != int
        self.assertFalse(foo(10, 1))
        # dict != tuple[int]
        self.assertFalse(foo(10, (1,)))

    def test_dict_to_from_meminfo(self):
        """
        Exercise dictobject.{_as_meminfo, _from_meminfo}
        """
        @njit
        def make_content(nelem):
            for i in range(nelem):
                yield i, i + (i + 1) / 100

        @njit
        def boxer(nelem):
            d = dictobject.new_dict(int32, float64)
            for k, v in make_content(nelem):
                d[k] = v
            return dictobject._as_meminfo(d)

        dcttype = types.DictType(int32, float64)

        @njit
        def unboxer(mi):
            d = dictobject._from_meminfo(mi, dcttype)
            return list(d.items())

        mi = boxer(10)
        self.assertEqual(mi.refcount, 1)

        got = unboxer(mi)
        expected = list(make_content.py_func(10))
        self.assertEqual(got, expected)

    def test_001_cannot_downcast_key(self):
        @njit
        def foo(n):
            d = dictobject.new_dict(int32, float64)
            for i in range(n):
                d[i] = i + 1
            # bad key type
            z = d.get(1j)
            return z

        with self.assertRaises(TypingError) as raises:
            foo(10)
        self.assertIn(
            'cannot safely cast complex128 to int32',
            str(raises.exception),
        )

    def test_002_cannot_downcast_default(self):
        @njit
        def foo(n):
            d = dictobject.new_dict(int32, float64)
            for i in range(n):
                d[i] = i + 1
            # bad default type
            z = d.get(2 * n, 1j)
            return z

        with self.assertRaises(TypingError) as raises:
            foo(10)
        self.assertIn(
            'cannot safely cast complex128 to float64',
            str(raises.exception),
        )

    def test_003_cannot_downcast_key(self):
        @njit
        def foo(n):
            d = dictobject.new_dict(int32, float64)
            for i in range(n):
                d[i] = i + 1
            # bad cast!?
            z = d.get(2.4)
            return z

        # should raise
        with self.assertRaises(TypingError) as raises:
            foo(10)
        self.assertIn(
            'cannot safely cast float64 to int32',
            str(raises.exception),
        )

    def test_004_cannot_downcast_key(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            # should raise TypingError
            d[1j] = 7.

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            'cannot safely cast complex128 to int32',
            str(raises.exception),
        )

    def test_005_cannot_downcast_value(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            # should raise TypingError
            d[1] = 1j

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            'cannot safely cast complex128 to float64',
            str(raises.exception),
        )

    def test_006_cannot_downcast_key(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            # raise TypingError
            d[11.5]

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            'cannot safely cast float64 to int32',
            str(raises.exception),
        )

    @unittest.skipUnless(sys.maxsize > 2 ** 32, "64 bit test only")
    def test_007_collision_checks(self):
        # this checks collisions in real life for 64bit systems
        @njit
        def foo(v1, v2):
            d = dictobject.new_dict(int64, float64)
            c1 = np.uint64(2 ** 61 - 1)
            c2 = np.uint64(0)
            assert hash(c1) == hash(c2)
            d[c1] = v1
            d[c2] = v2
            return (d[c1], d[c2])

        a, b = 10., 20.
        x, y = foo(a, b)
        self.assertEqual(x, a)
        self.assertEqual(y, b)

    def test_008_lifo_popitem(self):
        # check that (keys, vals) are LIFO .popitem()
        @njit
        def foo(n):
            d = dictobject.new_dict(int32, float64)
            for i in range(n):
                d[i] = i + 1
            keys = []
            vals = []
            for i in range(n):
                tmp = d.popitem()
                keys.append(tmp[0])
                vals.append(tmp[1])
            return keys, vals

        z = 10
        gk, gv = foo(z)

        self.assertEqual(gk, [x for x in reversed(range(z))])
        self.assertEqual(gv, [x + 1 for x in reversed(range(z))])

    def test_010_cannot_downcast_default(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            d[0] = 6.
            d[1] = 7.
            # pop'd default must have same type as value
            d.pop(11, 12j)

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "cannot safely cast complex128 to float64",
            str(raises.exception),
        )

    def test_011_cannot_downcast_key(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            d[0] = 6.
            d[1] = 7.
            # pop'd key must have same type as key
            d.pop(11j)

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "cannot safely cast complex128 to int32",
            str(raises.exception),
        )

    def test_012_cannot_downcast_key(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            d[0] = 6.
            # invalid key type
            return 1j in d

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "cannot safely cast complex128 to int32",
            str(raises.exception),
        )

    def test_013_contains_empty_dict(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            # contains on empty dict
            return 1 in d

        self.assertFalse(foo())

    def test_014_not_contains_empty_dict(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            # not contains empty dict
            return 1 not in d

        self.assertTrue(foo())

    def test_015_dict_clear(self):
        @njit
        def foo(n):
            d = dictobject.new_dict(int32, float64)
            for i in range(n):
                d[i] = i + 1
            x = len(d)
            d.clear()
            y = len(d)
            return x, y

        m = 10
        self.assertEqual(foo(m), (m, 0))

    def test_016_cannot_downcast_key(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            # key is wrong type
            d.setdefault(1j, 12.)

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "cannot safely cast complex128 to int32",
            str(raises.exception),
        )

    def test_017_cannot_downcast_default(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            # default value is wrong type
            d.setdefault(1, 12.j)

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "cannot safely cast complex128 to float64",
            str(raises.exception),
        )

    def test_018_keys_iter_are_views(self):
        # this is broken somewhere in llvmlite, intent of test is to check if
        # keys behaves like a view or not
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            d[11] = 12.
            k1 = d.keys()
            d[22] = 9.
            k2 = d.keys()
            rk1 = [x for x in k1]
            rk2 = [x for x in k2]
            return rk1, rk2

        a, b = foo()
        self.assertEqual(a, b)
        self.assertEqual(a, [11, 22])

    # Not implemented yet
    @unittest.expectedFailure
    def test_019(self):
        # should keys/vals be set-like?
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            d[11] = 12.
            d[22] = 9.
            k2 = d.keys() & {12, }
            return k2

        print(foo())

    def test_020_string_key(self):
        @njit
        def foo():
            d = dictobject.new_dict(types.unicode_type, float64)
            d['a'] = 1.
            d['b'] = 2.
            d['c'] = 3.
            d['d'] = 4.
            out = []
            for x in d.items():
                out.append(x)
            return out, d['a']

        items, da = foo()
        self.assertEqual(items, [('a', 1.), ('b', 2.), ('c', 3.), ('d', 4)])
        self.assertEqual(da, 1.)

    def test_021_long_str_key(self):
        @njit
        def foo():
            d = dictobject.new_dict(types.unicode_type, float64)
            tmp = []
            for i in range(10000):
                tmp.append('a')
            s = ''.join(tmp)
            d[s] = 1.
            out = list(d.items())
            return out
        self.assertEqual(foo(), [('a' * 10000, 1)])

    def test_022_references_juggle(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)
            e = d
            d[1] = 12.
            e[2] = 14.
            e = dictobject.new_dict(int32, float64)
            e[1] = 100.
            e[2] = 1000.
            f = d
            d = e

            k1 = [x for x in d.items()]
            k2 = [x for x in e.items()]
            k3 = [x for x in f.items()]

            return k1, k2, k3

        k1, k2, k3 = foo()
        self.assertEqual(k1, [(1, 100.0), (2, 1000.0)])
        self.assertEqual(k2, [(1, 100.0), (2, 1000.0)])
        self.assertEqual(k3, [(1, 12), (2, 14)])

    def test_023_closure(self):
        @njit
        def foo():
            d = dictobject.new_dict(int32, float64)

            def bar():
                d[1] = 12.
                d[2] = 14.
            bar()
            return [x for x in d.keys()]

        self.assertEqual(foo(), [1, 2])

    def test_024_unicode_getitem_keys(self):
        # See issue #6135
        @njit
        def foo():
            s = 'a\u1234'
            d = {s[0] : 1}
            return d['a']

        self.assertEqual(foo(), foo.py_func())

        @njit
        def foo():
            s = 'abc\u1234'
            d = {s[:1] : 1}
            return d['a']

        self.assertEqual(foo(), foo.py_func())

    def test_issue6570_alignment_padding(self):
        # Create a key type that is 12-bytes long on a 8-byte aligned system
        # so that the a 4-byte padding is needed.
        # If the 4-byte padding is not zero-filled, it will have garbage data
        # that affects key matching in the lookup.
        keyty = types.Tuple([types.uint64, types.float32])

        @njit
        def foo():
            d = dictobject.new_dict(keyty, float64)
            t1 = np.array([3], dtype=np.uint64)
            t2 = np.array([5.67], dtype=np.float32)
            v1 = np.array([10.23], dtype=np.float32)
            d[(t1[0], t2[0])] = v1[0]
            return (t1[0], t2[0]) in d

        self.assertTrue(foo())

    def test_dict_update(self):
        """
        Tests dict.update works with various dictionaries.
        """
        n = 10

        def f1(n):
            """
            Test update with a regular dictionary.
            """
            d1 = {i: i + 1 for i in range(n)}
            d2 = {3 * i: i for i in range(n)}
            d1.update(d2)
            return d1

        py_func = f1
        cfunc = njit()(f1)
        a = py_func(n)
        b = cfunc(n)
        self.assertEqual(a, b)

        def f2(n):
            """
            Test update where one of the dictionaries
            is created as a Python literal.
            """
            d1 = {
                1: 2,
                3: 4,
                5: 6
            }
            d2 = {3 * i: i for i in range(n)}
            d1.update(d2)
            return d1

        py_func = f2
        cfunc = njit()(f2)
        a = py_func(n)
        b = cfunc(n)
        self.assertEqual(a, b)


class TestDictTypeCasting(TestCase):
    def check_good(self, fromty, toty):
        _sentry_safe_cast(fromty, toty)

    def check_bad(self, fromty, toty):
        with self.assertRaises(TypingError) as raises:
            _sentry_safe_cast(fromty, toty)
        self.assertIn(
            'cannot safely cast {fromty} to {toty}'.format(**locals()),
            str(raises.exception),
        )

    def test_cast_int_to(self):
        self.check_good(types.int32, types.float32)
        self.check_good(types.int32, types.float64)
        self.check_good(types.int32, types.complex128)
        self.check_good(types.int64, types.complex128)
        self.check_bad(types.int32, types.complex64)
        self.check_good(types.int8, types.complex64)

    def test_cast_float_to(self):
        self.check_good(types.float32, types.float64)
        self.check_good(types.float32, types.complex64)
        self.check_good(types.float64, types.complex128)

    def test_cast_bool_to(self):
        self.check_good(types.boolean, types.int32)
        self.check_good(types.boolean, types.float64)
        self.check_good(types.boolean, types.complex128)


class TestTypedDict(MemoryLeakMixin, TestCase):
    def test_basic(self):
        d = Dict.empty(int32, float32)
        # len
        self.assertEqual(len(d), 0)
        # setitems
        d[1] = 1
        d[2] = 2.3
        d[3] = 3.4
        self.assertEqual(len(d), 3)
        # keys
        self.assertEqual(list(d.keys()), [1, 2, 3])
        # values
        for x, y in zip(list(d.values()), [1, 2.3, 3.4]):
            self.assertAlmostEqual(x, y, places=4)
        # getitem
        self.assertAlmostEqual(d[1], 1)
        self.assertAlmostEqual(d[2], 2.3, places=4)
        self.assertAlmostEqual(d[3], 3.4, places=4)
        # deltiem
        del d[2]
        self.assertEqual(len(d), 2)
        # get
        self.assertIsNone(d.get(2))
        # setdefault
        d.setdefault(2, 100)
        d.setdefault(3, 200)
        self.assertEqual(d[2], 100)
        self.assertAlmostEqual(d[3], 3.4, places=4)
        # update
        d.update({4: 5, 5: 6})
        self.assertAlmostEqual(d[4], 5)
        self.assertAlmostEqual(d[5], 6)
        # contains
        self.assertTrue(4 in d)
        # items
        pyd = dict(d.items())
        self.assertEqual(len(pyd), len(d))
        # pop
        self.assertAlmostEqual(d.pop(4), 5)
        # popitem
        nelem = len(d)
        k, v = d.popitem()
        self.assertEqual(len(d), nelem - 1)
        self.assertTrue(k not in d)
        # __eq__ & copy
        copied = d.copy()
        self.assertEqual(copied, d)
        self.assertEqual(list(copied.items()), list(d.items()))

    def test_copy_from_dict(self):
        expect = {k: float(v) for k, v in zip(range(10), range(10, 20))}
        nbd = Dict.empty(int32, float64)
        for k, v in expect.items():
            nbd[k] = v
        got = dict(nbd)
        self.assertEqual(got, expect)

    def test_compiled(self):
        @njit
        def producer():
            d = Dict.empty(int32, float64)
            d[1] = 1.23
            return d

        @njit
        def consumer(d):
            return d[1]

        d = producer()
        val = consumer(d)
        self.assertEqual(val, 1.23)

    def test_gh7908(self):
        d = Dict.empty(
            key_type=types.Tuple([types.uint32,
                                  types.uint32]),
            value_type=int64)

        d[(1, 1)] = 12345
        self.assertEqual(d[(1, 1)], d.get((1, 1)))

    def check_stringify(self, strfn, prefix=False):
        nbd = Dict.empty(int32, int32)
        d = {}
        nbd[1] = 2
        d[1] = 2
        checker = self.assertIn if prefix else self.assertEqual
        checker(strfn(d), strfn(nbd))
        nbd[2] = 3
        d[2] = 3
        checker(strfn(d), strfn(nbd))
        for i in range(10, 20):
            nbd[i] = i + 1
            d[i] = i + 1
        checker(strfn(d), strfn(nbd))
        if prefix:
            self.assertTrue(strfn(nbd).startswith('DictType'))

    def test_repr(self):
        self.check_stringify(repr, prefix=True)

    def test_str(self):
        self.check_stringify(str)


class DictIterableCtor:

    def test_iterable_type_constructor(self):
        # https://docs.python.org/3/library/stdtypes.html#dict
        @njit
        def func1(a, b):
            d = Dict(zip(a, b))
            return d

        @njit
        def func2(a_, b):
            a = range(3)
            return Dict(zip(a, b))

        @njit
        def func3(a_, b):
            a = [0, 1, 2]
            return Dict(zip(a, b))

        @njit
        def func4(a, b):
            c = zip(a, b)
            return Dict(zip(a, zip(c, a)))

        @njit
        def func5(a, b):
            return Dict(zip(zip(a, b), b))

        @njit
        def func6(items):
            return Dict(items)

        @njit
        def func7(k, v):
            return Dict({k: v})  # mapping - not supported

        @njit
        def func8(k, v):
            d = Dict()
            d[k] = v
            return d

        def _get_dict(py_dict):
            d = Dict()
            for k, v in py_dict.items():
                d[k] = v
            return d

        vals = (
            (func1, [(0, 1, 2), 'abc'], _get_dict({0: 'a', 1: 'b', 2: 'c'})),
            (func2, [(0, 1, 2), 'abc'], _get_dict({0: 'a', 1: 'b', 2: 'c'})),
            (func3, [(0, 1, 2), 'abc'], _get_dict({0: 'a', 1: 'b', 2: 'c'})),
            (func4, [(0, 1, 2), 'abc'], _get_dict(
                {0: ((0, 'a'), 0), 1: ((1, 'b'), 1), 2: ((2, 'c'), 2)})),
            (func5, [(0, 1, 2), 'abc'], _get_dict(
                {(0, 'a'): 'a', (1, 'b'): 'b', (2, 'c'): 'c'})),
            # (func6, [(),], Dict({})),
            (func6, [((1, 'a'), (3, 'b')),], _get_dict({1: 'a', 3: 'b'})),
            (func1, ['key', _get_dict({1: 'abc'})], _get_dict({'k': 1})),
            (func8, ['key', _get_dict({1: 'abc'})], _get_dict(
                {'key': _get_dict({1: 'abc'})})),
            (func8, ['key', List([1, 2, 3])], _get_dict(
                {'key': List([1, 2, 3])})),
        )

        for func, args, expected in vals:
            if self.jit_enabled:
                got = func(*args)
            else:
                got = func.py_func(*args)
            self.assertPreciseEqual(expected, got)


class TestDictIterableCtorJit(TestCase, DictIterableCtor):

    def setUp(self):
        self.jit_enabled = True

    def test_exception_no_iterable_arg(self):
        @njit
        def ctor():
            return Dict(3)

        msg = ".*No implementation of function.*"
        with self.assertRaisesRegex(TypingError, msg):
            ctor()

    def test_exception_dict_mapping(self):
        @njit
        def ctor():
            return Dict({1: 2, 3: 4})

        msg = ".*No implementation of function.*"
        with self.assertRaisesRegex(TypingError, msg):
            ctor()

    def test_exception_setitem(self):
        @njit
        def ctor():
            return Dict(((1, 'a'), (2, 'b', 3)))

        msg = ".*No implementation of function.*"
        with self.assertRaisesRegex(TypingError, msg):
            ctor()


class TestDictIterableCtorNoJit(TestCase, DictIterableCtor):

    def setUp(self):
        self.jit_enabled = False

    def test_exception_nargs(self):
        msg = 'Dict expect at most 1 argument, got 2'
        with self.assertRaisesRegex(TypingError, msg):
            Dict(1, 2)

    def test_exception_mapping_ctor(self):
        msg = r'.*dict\(mapping\) is not supported.*'  # noqa: W605
        with self.assertRaisesRegex(TypingError, msg):
            Dict({1: 2})

    def test_exception_non_iterable_arg(self):
        msg = '.*object is not iterable.*'
        with self.assertRaisesRegex(TypingError, msg):
            Dict(3)

    def test_exception_setitem(self):
        msg = ".*dictionary update sequence element #1 has length 3.*"
        with self.assertRaisesRegex(ValueError, msg):
            Dict(((1, 'a'), (2, 'b', 3)))


class TestDictRefctTypes(MemoryLeakMixin, TestCase):

    def test_str_key(self):
        @njit
        def foo():
            d = Dict.empty(
                key_type=types.unicode_type,
                value_type=types.int32,
            )
            d["123"] = 123
            d["321"] = 321
            return d

        d = foo()
        self.assertEqual(d['123'], 123)
        self.assertEqual(d['321'], 321)
        expect = {'123': 123, '321': 321}
        self.assertEqual(dict(d), expect)
        # Test insert replacement
        d['123'] = 231
        expect['123'] = 231
        self.assertEqual(d['123'], 231)
        self.assertEqual(dict(d), expect)
        # Test dictionary growth
        nelem = 100
        for i in range(nelem):
            d[str(i)] = i
            expect[str(i)] = i
        for i in range(nelem):
            self.assertEqual(d[str(i)], i)
        self.assertEqual(dict(d), expect)

    def test_str_val(self):
        @njit
        def foo():
            d = Dict.empty(
                key_type=types.int32,
                value_type=types.unicode_type,
            )
            d[123] = "123"
            d[321] = "321"
            return d

        d = foo()
        self.assertEqual(d[123], '123')
        self.assertEqual(d[321], '321')
        expect = {123: '123', 321: '321'}
        self.assertEqual(dict(d), expect)
        # Test insert replacement
        d[123] = "231"
        expect[123] = "231"
        self.assertEqual(dict(d), expect)
        # Test dictionary growth
        nelem = 1
        for i in range(nelem):
            d[i] = str(i)
            expect[i] = str(i)
        for i in range(nelem):
            self.assertEqual(d[i], str(i))
        self.assertEqual(dict(d), expect)

    def test_str_key_array_value(self):
        np.random.seed(123)
        d = Dict.empty(
            key_type=types.unicode_type,
            value_type=types.float64[:],
        )
        expect = []
        expect.append(np.random.random(10))
        d['mass'] = expect[-1]
        expect.append(np.random.random(20))
        d['velocity'] = expect[-1]
        for i in range(100):
            expect.append(np.random.random(i))
            d[str(i)] = expect[-1]
        self.assertEqual(len(d), len(expect))
        self.assertPreciseEqual(d['mass'], expect[0])
        self.assertPreciseEqual(d['velocity'], expect[1])
        # Ordering is kept
        for got, exp in zip(d.values(), expect):
            self.assertPreciseEqual(got, exp)

        # Try deleting
        self.assertTrue('mass' in d)
        self.assertTrue('velocity' in d)
        del d['mass']
        self.assertFalse('mass' in d)
        del d['velocity']
        self.assertFalse('velocity' in d)
        del expect[0:2]

        for i in range(90):
            k, v = d.popitem()
            w = expect.pop()
            self.assertPreciseEqual(v, w)

        # Trigger a resize
        expect.append(np.random.random(10))
        d["last"] = expect[-1]

        # Ordering is kept
        for got, exp in zip(d.values(), expect):
            self.assertPreciseEqual(got, exp)

    def test_dict_of_dict_int_keyval(self):
        def inner_numba_dict():
            d = Dict.empty(
                key_type=types.intp,
                value_type=types.intp,
            )
            return d

        d = Dict.empty(
            key_type=types.intp,
            value_type=types.DictType(types.intp, types.intp),
        )

        def usecase(d, make_inner_dict):
            for i in range(100):
                mid = make_inner_dict()
                for j in range(i + 1):
                    mid[j] = j * 10000
                d[i] = mid
            return d

        got = usecase(d, inner_numba_dict)
        expect = usecase({}, dict)

        self.assertIsInstance(expect, dict)

        self.assertEqual(dict(got), expect)

        # Delete items
        for where in [12, 3, 6, 8, 10]:
            del got[where]
            del expect[where]
            self.assertEqual(dict(got), expect)

    def test_dict_of_dict_npm(self):
        inner_dict_ty = types.DictType(types.intp, types.intp)

        @njit
        def inner_numba_dict():
            d = Dict.empty(
                key_type=types.intp,
                value_type=types.intp,
            )
            return d

        @njit
        def foo(count):
            d = Dict.empty(
                key_type=types.intp,
                value_type=inner_dict_ty,
            )
            for i in range(count):
                d[i] = inner_numba_dict()
                for j in range(i + 1):
                    d[i][j] = j

            return d

        d = foo(100)
        ct = 0
        for k, dd in d.items():
            ct += 1
            self.assertEqual(len(dd), k + 1)
            for kk, vv in dd.items():
                self.assertEqual(kk, vv)

        self.assertEqual(ct, 100)

    def test_delitem(self):
        d = Dict.empty(types.int64, types.unicode_type)
        d[1] = 'apple'

        @njit
        def foo(x, k):
            del x[1]

        foo(d, 1)
        self.assertEqual(len(d), 0)
        self.assertFalse(d)

    def test_getitem_return_type(self):
        # Dict.__getitem__ must return non-optional type.
        d = Dict.empty(types.int64, types.int64[:])
        d[1] = np.arange(10, dtype=np.int64)

        @njit
        def foo(d):
            d[1] += 100
            return d[1]

        foo(d)
        # Return type is an array, not optional
        retty = foo.nopython_signatures[0].return_type
        self.assertIsInstance(retty, types.Array)
        self.assertNotIsInstance(retty, types.Optional)
        # Value is correctly updated
        self.assertPreciseEqual(d[1], np.arange(10, dtype=np.int64) + 100)

    def test_storage_model_mismatch(self):
        # https://github.com/numba/numba/issues/4520
        # check for storage model mismatch in refcount ops generation
        dct = Dict()
        ref = [
            ("a", True, "a"),
            ("b", False, "b"),
            ("c", False, "c"),
        ]
        # populate
        for x in ref:
            dct[x] = x
        # test
        for i, x in enumerate(ref):
            self.assertEqual(dct[x], x)


class TestDictForbiddenTypes(TestCase):
    def assert_disallow(self, expect, callable):
        with self.assertRaises(TypingError) as raises:
            callable()
        msg = str(raises.exception)
        self.assertIn(expect, msg)

    def assert_disallow_key(self, ty):
        msg = '{} as key is forbidden'.format(ty)
        self.assert_disallow(msg, lambda: Dict.empty(ty, types.intp))

        @njit
        def foo():
            Dict.empty(ty, types.intp)
        self.assert_disallow(msg, foo)

    def assert_disallow_value(self, ty):
        msg = '{} as value is forbidden'.format(ty)
        self.assert_disallow(msg, lambda: Dict.empty(types.intp, ty))

        @njit
        def foo():
            Dict.empty(types.intp, ty)
        self.assert_disallow(msg, foo)

    def test_disallow_list(self):
        self.assert_disallow_key(types.List(types.intp))
        self.assert_disallow_value(types.List(types.intp))

    def test_disallow_set(self):
        self.assert_disallow_key(types.Set(types.intp))
        self.assert_disallow_value(types.Set(types.intp))


class TestDictInferred(TestCase):
    def test_simple_literal(self):
        @njit
        def foo():
            d = Dict()
            d[123] = 321
            return d

        k, v = 123, 321
        d = foo()
        self.assertEqual(dict(d), {k: v})
        self.assertEqual(typeof(d).key_type, typeof(k))
        self.assertEqual(typeof(d).value_type, typeof(v))

    def test_simple_args(self):
        @njit
        def foo(k, v):
            d = Dict()
            d[k] = v
            return d

        k, v = 123, 321
        d = foo(k, v)
        self.assertEqual(dict(d), {k: v})
        self.assertEqual(typeof(d).key_type, typeof(k))
        self.assertEqual(typeof(d).value_type, typeof(v))

    def test_simple_upcast(self):
        @njit
        def foo(k, v, w):
            d = Dict()
            d[k] = v
            d[k] = w
            return d

        k, v, w = 123, 32.1, 321
        d = foo(k, v, w)
        self.assertEqual(dict(d), {k: w})
        self.assertEqual(typeof(d).key_type, typeof(k))
        self.assertEqual(typeof(d).value_type, typeof(v))

    def test_conflicting_value_type(self):
        @njit
        def foo(k, v, w):
            d = Dict()
            d[k] = v
            d[k] = w
            return d

        k, v, w = 123, 321, 32.1
        with self.assertRaises(TypingError) as raises:
            foo(k, v, w)
        self.assertIn(
            'cannot safely cast float64 to {}'.format(typeof(v)),
            str(raises.exception),
        )

    def test_conflicting_key_type(self):
        @njit
        def foo(k, h, v):
            d = Dict()
            d[k] = v
            d[h] = v
            return d

        k, h, v = 123, 123.1, 321
        with self.assertRaises(TypingError) as raises:
            foo(k, h, v)
        self.assertIn(
            'cannot safely cast float64 to {}'.format(typeof(v)),
            str(raises.exception),
        )

    def test_conflict_key_type_non_number(self):
        # Allow non-number types to cast unsafely
        @njit
        def foo(k1, v1, k2):
            d = Dict()
            d[k1] = v1
            return d, d[k2]

        # k2 will unsafely downcast typeof(k1)
        k1 = (np.int8(1), np.int8(2))
        k2 = (np.int32(1), np.int32(2))
        v1 = np.intp(123)

        with warnings.catch_warnings(record=True) as w:
            d, dk2 = foo(k1, v1, k2)
        self.assertEqual(len(w), 1)
        # Make sure the warning is about unsafe cast
        msg = 'unsafe cast from UniTuple(int32 x 2) to UniTuple(int8 x 2)'
        self.assertIn(msg, str(w[0]))

        keys = list(d.keys())
        self.assertEqual(keys[0], (1, 2))
        self.assertEqual(dk2, d[(np.int32(1), np.int32(2))])

    def test_ifelse_filled_both_branches(self):
        @njit
        def foo(k, v):
            d = Dict()
            if k:
                d[k] = v
            else:
                d[0xdead] = v + 1

            return d

        k, v = 123, 321
        d = foo(k, v)
        self.assertEqual(dict(d), {k: v})
        k, v = 0, 0
        d = foo(k, v)
        self.assertEqual(dict(d), {0xdead: v + 1})

    def test_ifelse_empty_one_branch(self):
        @njit
        def foo(k, v):
            d = Dict()
            if k:
                d[k] = v

            return d

        k, v = 123, 321
        d = foo(k, v)
        self.assertEqual(dict(d), {k: v})
        k, v = 0, 0
        d = foo(k, v)
        self.assertEqual(dict(d), {})
        self.assertEqual(typeof(d).key_type, typeof(k))
        self.assertEqual(typeof(d).value_type, typeof(v))

    def test_loop(self):
        @njit
        def foo(ks, vs):
            d = Dict()
            for k, v in zip(ks, vs):
                d[k] = v
            return d

        vs = list(range(4))
        ks = list(map(lambda x : x + 100, vs))
        d = foo(ks, vs)
        self.assertEqual(dict(d), dict(zip(ks, vs)))

    def test_unused(self):
        @njit
        def foo():
            d = Dict()
            return d

        with self.assertRaises(TypingError) as raises:
            foo()
        self.assertIn(
            "imprecise type",
            str(raises.exception)
        )

    def test_define_after_use(self):
        @njit
        def foo(define):
            d = Dict()
            ct = len(d)
            for k, v in d.items():
                ct += v

            if define:
                # This will set the type
                d[1] = 2
            return ct, d, len(d)

        ct, d, n = foo(True)
        self.assertEqual(ct, 0)
        self.assertEqual(n, 1)
        self.assertEqual(dict(d), {1: 2})

        ct, d, n = foo(False)
        self.assertEqual(ct, 0)
        self.assertEqual(dict(d), {})
        self.assertEqual(n, 0)

    def test_dict_of_dict(self):
        @njit
        def foo(k1, k2, v):
            d = Dict()
            z1 = Dict()
            z1[k1 + 1] = v + k1
            z2 = Dict()
            z2[k2 + 2] = v + k2
            d[k1] = z1
            d[k2] = z2
            return d

        k1, k2, v = 100, 200, 321
        d = foo(k1, k2, v)
        self.assertEqual(
            dict(d),
            {
                k1: {k1 + 1: k1 + v},
                k2: {k2 + 2: k2 + v},
            },
        )

    def test_comprehension_basic(self):
        @njit
        def foo():
            return {i: 2 * i for i in range(10)}

        self.assertEqual(foo(), foo.py_func())

    def test_comprehension_basic_mixed_type(self):
        @njit
        def foo():
            return {i: float(j) for i, j in zip(range(10), range(10, 0, -1))}

        self.assertEqual(foo(), foo.py_func())

    def test_comprehension_involved(self):
        @njit
        def foo():
            a = {0: 'A', 1: 'B', 2: 'C'}
            return {3 + i: a[i] for i in range(3)}

        self.assertEqual(foo(), foo.py_func())

    def test_comprehension_fail_mixed_type(self):
        @njit
        def foo():
            a = {0: 'A', 1: 'B', 2: 1j}
            return {3 + i: a[i] for i in range(3)}

        with self.assertRaises(TypingError) as e:
            foo()

        excstr = str(e.exception)
        self.assertIn("Cannot cast complex128 to unicode_type", excstr)


class TestNonCompiledInfer(TestCase):
    def test_check_untyped_dict_ops(self):
        # Check operation on untyped dictionary
        d = Dict()
        self.assertFalse(d._typed)
        self.assertEqual(len(d), 0)
        self.assertEqual(str(d), str({}))
        self.assertEqual(list(iter(d)), [])
        # Test __getitem__
        with self.assertRaises(KeyError) as raises:
            d[1]
        self.assertEqual(str(raises.exception), str(KeyError(1)))
        # Test __delitem__
        with self.assertRaises(KeyError) as raises:
            del d[1]
        self.assertEqual(str(raises.exception), str(KeyError(1)))
        # Test .pop
        with self.assertRaises(KeyError):
            d.pop(1)
        self.assertEqual(str(raises.exception), str(KeyError(1)))
        # Test .pop
        self.assertIs(d.pop(1, None), None)
        # Test .get
        self.assertIs(d.get(1), None)
        # Test .popitem
        with self.assertRaises(KeyError) as raises:
            d.popitem()
        self.assertEqual(str(raises.exception),
                         str(KeyError('dictionary is empty')))
        # Test setdefault(k)
        with self.assertRaises(TypeError) as raises:
            d.setdefault(1)
        self.assertEqual(
            str(raises.exception),
            str(TypeError('invalid operation on untyped dictionary')),
        )
        # Test __contains__
        self.assertFalse(1 in d)
        # It's untyped
        self.assertFalse(d._typed)

    def test_getitem(self):
        # Test __getitem__
        d = Dict()
        d[1] = 2
        # It's typed now
        self.assertTrue(d._typed)
        self.assertEqual(d[1], 2)

    def test_setdefault(self):
        # Test setdefault(k, d)
        d = Dict()
        d.setdefault(1, 2)
        # It's typed now
        self.assertTrue(d._typed)
        self.assertEqual(d[1], 2)


@jitclass(spec=[('a', types.intp)])
class Bag(object):
    def __init__(self, a):
        self.a = a

    def __hash__(self):
        return hash(self.a)


class TestDictWithJitclass(TestCase):
    def test_jitclass_as_value(self):
        @njit
        def foo(x):
            d = Dict()
            d[0] = x
            d[1] = Bag(101)
            return d

        d = foo(Bag(a=100))
        self.assertEqual(d[0].a, 100)
        self.assertEqual(d[1].a, 101)


class TestNoJit(TestCase):
    """Exercise dictionary creation with JIT disabled. """

    def test_dict_create_no_jit_using_new_dict(self):
        with override_config('DISABLE_JIT', True):
            with forbid_codegen():
                d = dictobject.new_dict(int32, float32)
                self.assertEqual(type(d), dict)

    def test_dict_create_no_jit_using_Dict(self):
        with override_config('DISABLE_JIT', True):
            with forbid_codegen():
                d = Dict()
                self.assertEqual(type(d), dict)

    def test_dict_create_no_jit_using_empty(self):
        with override_config('DISABLE_JIT', True):
            with forbid_codegen():
                d = Dict.empty(types.int32, types.float32)
                self.assertEqual(type(d), dict)


class TestDictIterator(TestCase):
    def test_dict_iterator(self):
        @njit
        def fun1():
            dd = Dict.empty(key_type=types.intp,
                            value_type=types.intp)
            dd[0] = 10
            dd[1] = 20
            dd[2] = 30

            return list(dd.keys()), list(dd.values())

        @njit
        def fun2():
            dd = Dict.empty(key_type=types.intp,
                            value_type=types.intp)
            dd[4] = 77
            dd[5] = 88
            dd[6] = 99

            return list(dd.keys()), list(dd.values())
        res1 = fun1()
        res2 = fun2()

        self.assertEqual([0,1,2], res1[0])
        self.assertEqual([10,20,30], res1[1])
        self.assertEqual([4,5,6], res2[0])
        self.assertEqual([77,88,99], res2[1])


class TestTypedDictInitialValues(MemoryLeakMixin, TestCase):
    """Tests that typed dictionaries carry their initial value if present"""

    def test_homogeneous_and_literal(self):
        def bar(d):
            ...

        @overload(bar)
        def ol_bar(d):
            if d.initial_value is None:
                return lambda d: literally(d)
            self.assertTrue(isinstance(d, types.DictType))
            self.assertEqual(d.initial_value, {'a': 1, 'b': 2, 'c': 3})
            self.assertEqual(hasattr(d, 'literal_value'), False)
            return lambda d: d

        @njit
        def foo():
            # keys and values all have literal representation
            x = {'a': 1, 'b': 2, 'c': 3}
            bar(x)

        foo()

    def test_heterogeneous_but_castable_to_homogeneous(self):
        def bar(d):
            ...

        @overload(bar)
        def ol_bar(d):
            self.assertTrue(isinstance(d, types.DictType))
            self.assertEqual(d.initial_value, None)
            self.assertEqual(hasattr(d, 'literal_value'), False)
            return lambda d: d

        @njit
        def foo():
            # This dictionary will be typed based on 1j, i.e. complex128
            # as the values are not all literals, there's no "initial_value"
            # available irrespective of whether it's possible to rip this
            # information out of the bytecode.
            x = {'a': 1j, 'b': 2, 'c': 3}
            bar(x)

        foo()

    def test_heterogeneous_but_not_castable_to_homogeneous(self):
        def bar(d):
            ...

        @overload(bar)
        def ol_bar(d):
            a = {'a': 1, 'b': 2j, 'c': 3}

            def specific_ty(z):
                return types.literal(z) if types.maybe_literal(z) else typeof(z)
            expected = {types.literal(x): specific_ty(y) for x, y in a.items()}
            self.assertTrue(isinstance(d, types.LiteralStrKeyDict))
            self.assertEqual(d.literal_value, expected)
            self.assertEqual(hasattr(d, 'initial_value'), False)
            return lambda d: d

        @njit
        def foo():
            # This dictionary will be typed based on 1, i.e. intp, as the values
            # cannot all be cast to this type, but the keys are literal strings
            # this is a LiteralStrKey[Dict], there's no initial_value but there
            # is a literal_value.
            x = {'a': 1, 'b': 2j, 'c': 3}
            bar(x)

        foo()

    def test_mutation_not_carried(self):
        def bar(d):
            ...

        @overload(bar)
        def ol_bar(d):
            if d.initial_value is None:
                return lambda d: literally(d)
            self.assertTrue(isinstance(d, types.DictType))
            self.assertEqual(d.initial_value, {'a': 1, 'b': 2, 'c': 3})
            return lambda d: d

        @njit
        def foo():
            # This dictionary is mutated, check the initial_value carries
            # correctly and is not mutated
            x = {'a': 1, 'b': 2, 'c': 3}
            x['d'] = 4
            bar(x)

        foo()

    def test_mutation_not_carried_single_function(self):
        # this is another pattern for using literally

        @njit
        def nop(*args):
            pass

        for fn, iv in (nop, None), (literally, {'a': 1, 'b': 2, 'c': 3}):
            @njit
            def baz(x):
                pass

            def bar(z):
                pass

            @overload(bar)
            def ol_bar(z):
                def impl(z):
                    fn(z)
                    baz(z)
                return impl

            @njit
            def foo():
                x = {'a': 1, 'b': 2, 'c': 3}
                bar(x)
                x['d'] = 4
                return x

            foo()
            # baz should be specialised based on literally being invoked and
            # the literal/unliteral arriving at the call site
            larg = baz.signatures[0][0]
            self.assertEqual(larg.initial_value, iv)

    def test_unify_across_function_call(self):

        @njit
        def bar(x):
            o = {1: 2}
            if x:
                o = {2: 3}
            return o

        @njit
        def foo(x):
            if x:
                d = {3: 4}
            else:
                d = bar(x)
            return d

        e1 = Dict()
        e1[3] = 4
        e2 = Dict()
        e2[1] = 2
        self.assertEqual(foo(True), e1)
        self.assertEqual(foo(False), e2)


class TestLiteralStrKeyDict(MemoryLeakMixin, TestCase):
    """ Tests for dictionaries with string keys that can map to anything!"""

    def test_basic_const_lowering_boxing(self):
        @njit
        def foo():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            return (ld['a'], ld['b'], ld['c'])

        self.assertEqual(foo(), (1, 2j, 'd'))

    def test_basic_nonconst_in_scope(self):
        @njit
        def foo(x):
            y = x + 5
            e = True if y > 2 else False
            ld = {'a': 1, 'b': 2j, 'c': 'd', 'non_const': e}
            return ld['non_const']

        # Recall that key non_const has a value of a known type, bool, and it's
        # value is stuffed in at run time, this is permitted as the dictionary
        # is immutable in type
        self.assertTrue(foo(34))
        self.assertFalse(foo(-100))

    def test_basic_nonconst_freevar(self):
        e = 5

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            self.assertEqual(x.literal_value,
                             {types.literal('a'): types.literal(1),
                              types.literal('b'): typeof(2j),
                              types.literal('c'): types.literal('d'),
                              types.literal('d'): types.literal(5)})

            def impl(x):
                pass
            return impl

        @njit
        def foo():
            ld = {'a': 1, 'b': 2j, 'c': 'd', 'd': e}
            bar(ld)

        foo()

    def test_literal_value(self):

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            self.assertEqual(x.literal_value,
                             {types.literal('a'): types.literal(1),
                              types.literal('b'): typeof(2j),
                              types.literal('c'): types.literal('d')})

            def impl(x):
                pass
            return impl

        @njit
        def foo():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            bar(ld)

        foo()

    def test_list_and_array_as_value(self):

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            self.assertEqual(x.literal_value,
                             {types.literal('a'): types.literal(1),
                              types.literal('b'):
                              types.List(types.intp, initial_value=[1,2,3]),
                              types.literal('c'): typeof(np.zeros(5))})

            def impl(x):
                pass
            return impl

        @njit
        def foo():
            b = [1, 2, 3]
            ld = {'a': 1, 'b': b, 'c': np.zeros(5)}
            bar(ld)

        foo()

    def test_repeated_key_literal_value(self):

        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            # order is important, 'a' was seen first, but updated later
            self.assertEqual(x.literal_value,
                             {types.literal('a'): types.literal('aaaa'),
                              types.literal('b'): typeof(2j),
                              types.literal('c'): types.literal('d')})

            def impl(x):
                pass
            return impl

        @njit
        def foo():
            ld = {'a': 1, 'a': 10, 'b': 2j, 'c': 'd', 'a': 'aaaa'} # noqa #F601
            bar(ld)

        foo()

    def test_read_only(self):

        def _len():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            return len(ld)

        def static_getitem():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            return ld['b']

        def contains():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            return 'b' in ld, 'f' in ld

        def copy():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            new = ld.copy()
            return ld == new

        rdonlys = (_len, static_getitem, contains, copy)

        for test in rdonlys:
            with self.subTest(test.__name__):
                self.assertPreciseEqual(njit(test)(), test())

    def test_mutation_failure(self):

        def setitem():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            ld['a'] = 12

        def delitem():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            del ld['a']

        def popitem():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            ld.popitem()

        def pop():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            ld.pop()

        def clear():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            ld.clear()

        def setdefault():
            ld = {'a': 1, 'b': 2j, 'c': 'd'}
            ld.setdefault('f', 1)

        illegals = (setitem, delitem, popitem, pop, clear, setdefault)

        for test in illegals:
            with self.subTest(test.__name__):
                with self.assertRaises(TypingError) as raises:
                    njit(test)()
                expect = "Cannot mutate a literal dictionary"
                self.assertIn(expect, str(raises.exception))

    def test_get(self):

        @njit
        def get(x):
            ld = {'a': 2j, 'c': 'd'}
            return ld.get(x)

        @njit
        def getitem(x):
            ld = {'a': 2j, 'c': 'd'}
            return ld[x]

        for test in (get, getitem):
            with self.subTest(test.__name__):
                with self.assertRaises(TypingError) as raises:
                    test('a')
                expect = "Cannot get{item}() on a literal dictionary"
                self.assertIn(expect, str(raises.exception))

    def test_dict_keys(self):

        @njit
        def foo():
            ld = {'a': 2j, 'c': 'd'}
            return [x for x in ld.keys()]

        self.assertEqual(foo(), ['a', 'c'])

    def test_dict_values(self):

        @njit
        def foo():
            ld = {'a': 2j, 'c': 'd'}
            return ld.values()

        self.assertEqual(foo(), (2j, 'd'))

    def test_dict_items(self):
        @njit
        def foo():
            ld = {'a': 2j, 'c': 'd', 'f': np.zeros((5))}
            return ld.items()

        self.assertPreciseEqual(foo(),
                                (('a', 2j), ('c', 'd'), ('f', np.zeros((5)))))

    def test_dict_return(self):

        @njit
        def foo():
            ld = {'a': 2j, 'c': 'd'}
            return ld

        # escaping heterogeneous dictionary is not supported
        with self.assertRaises(TypeError) as raises:
            foo()

        excstr = str(raises.exception)
        self.assertIn("cannot convert native LiteralStrKey", excstr)

    def test_dict_unify(self):
        @njit
        def foo(x):
            if x + 7 > 4:
                a = {'a': 2j, 'c': 'd', 'e': np.zeros(4)}
            else:
                # Note the use of a different literal str for key 'c'
                a = {'a': 5j, 'c': 'CAT', 'e': np.zeros((5,))}
            return a['c']

        self.assertEqual(foo(100), 'd')
        self.assertEqual(foo(-100), 'CAT')
        self.assertEqual(foo(100), foo.py_func(100))
        self.assertEqual(foo(-100), foo.py_func(-100))

    def test_dict_not_unify(self):

        @njit
        def key_mismatch(x):
            if x + 7 > 4:
                a = {'BAD_KEY': 2j, 'c': 'd', 'e': np.zeros(4)}
            else:
                a = {'a': 5j, 'c': 'CAT', 'e': np.zeros((5,))}
            # prevents inline of return on py310
            py310_defeat1 = 1  # noqa
            py310_defeat2 = 2  # noqa
            py310_defeat3 = 3  # noqa
            py310_defeat4 = 4  # noqa
            return a['a']

        with self.assertRaises(TypingError) as raises:
            key_mismatch(100)

        self.assertIn("Cannot unify LiteralStrKey", str(raises.exception))

        @njit
        def value_type_mismatch(x):
            if x + 7 > 4:
                a = {'a': 2j, 'c': 'd', 'e': np.zeros((4, 3))}
            else:
                a = {'a': 5j, 'c': 'CAT', 'e': np.zeros((5,))}
            # prevents inline of return on py310
            py310_defeat1 = 1  # noqa
            py310_defeat2 = 2  # noqa
            py310_defeat3 = 3  # noqa
            py310_defeat4 = 4  # noqa
            return a['a']

        with self.assertRaises(TypingError) as raises:
            value_type_mismatch(100)

        self.assertIn("Cannot unify LiteralStrKey", str(raises.exception))

    def test_dict_value_coercion(self):
        # checks that things coerce or not!

        p = {# safe and no conversion: TypedDict
             (np.int32, np.int32): types.DictType,
             # safe and convertible: TypedDict
             (np.int32, np.int8): types.DictType,
             # safe convertible: TypedDict
             (np.complex128, np.int32): types.DictType,
             # unsafe not convertible: LiteralStrKey
             (np.int32, np.complex128): types.LiteralStrKeyDict,
             # unsafe not convertible: LiteralStrKey
             (np.int32, np.array): types.LiteralStrKeyDict,
             # unsafe not convertible: LiteralStrKey
             (np.array, np.int32): types.LiteralStrKeyDict,
             # unsafe not convertible: LiteralStrKey
             (np.int8, np.int32): types.LiteralStrKeyDict,
             # unsafe not convertible: LiteralStrKey (issue #6420 case)
             (np.int64, np.float64): types.LiteralStrKeyDict,}

        def bar(x):
            pass

        for dts, container in p.items():
            @overload(bar)
            def ol_bar(x):
                self.assertTrue(isinstance(x, container))

                def impl(x):
                    pass
                return impl

            ty1, ty2 = dts

            @njit
            def foo():
                d = {'a': ty1(1), 'b': ty2(2)}
                bar(d)

            foo()

    def test_build_map_op_code(self):
        # tests building dictionaries via `build_map`, which, for statically
        # determinable str key->things cases is just a single key:value
        # any other build_map would either end up as being non-const str keys
        # or keys of some non-string type and therefore not considered.
        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            def impl(x):
                pass
            return impl

        @njit
        def foo():
            a = {'a': {'b1': 10, 'b2': 'string'}}
            bar(a)

        foo()

    def test_dict_as_arg(self):
        @njit
        def bar(fake_kwargs=None):
            if fake_kwargs is not None:
                # Add 10 to array in key 'd'
                fake_kwargs['d'][:] += 10

        @njit
        def foo():
            a = 1
            b = 2j
            c = 'string'
            d = np.zeros(3)
            e = {'a': a, 'b': b, 'c': c, 'd': d}
            bar(fake_kwargs=e)
            return e['d']

        np.testing.assert_allclose(foo(), np.ones(3) * 10)

    def test_dict_with_single_literallist_value(self):
        #see issue #6094
        @njit
        def foo():
            z = {"A": [lambda a: 2 * a, "B"]}
            return z["A"][0](5)

        self.assertPreciseEqual(foo(), foo.py_func())

    def test_tuple_not_in_mro(self):
        # Related to #6094, make sure that LiteralStrKey does not inherit from
        # types.BaseTuple as this breaks isinstance checks.
        def bar(x):
            pass

        @overload(bar)
        def ol_bar(x):
            self.assertFalse(isinstance(x, types.BaseTuple))
            self.assertTrue(isinstance(x, types.LiteralStrKeyDict))
            return lambda x: ...

        @njit
        def foo():
            d = {'a': 1, 'b': 'c'}
            bar(d)

        foo()

    def test_const_key_not_in_dict(self):

        @njit
        def foo():
            a = {'not_a': 2j, 'c': 'd', 'e': np.zeros(4)}
            return a['a']

        with self.assertRaises(TypingError) as raises:
            foo()

        self.assertIn("Key 'a' is not in dict.", str(raises.exception))

    def test_uncommon_identifiers(self):
        # Tests uncommon identifiers like numerical values and operators in
        # the key fields. See #6518 and #7416.

        # Numerical values in keys
        @njit
        def foo():
            d = {'0': np.ones(5), '1': 4}
            return len(d)

        self.assertPreciseEqual(foo(), foo.py_func())

        # operators in keys
        @njit
        def bar():
            d = {'+': np.ones(5), 'x--': 4}
            return len(d)

        self.assertPreciseEqual(bar(), bar.py_func())

    def test_update_error(self):
        # Tests that dict.update produces a reasonable
        # error with a LiteralStrKeyDict input.
        @njit
        def foo():

            d1 = {
                'a': 2,
                'b': 4,
                'c': 'a'
            }
            d1.update({'x': 3})
            return d1

        with self.assertRaises(TypingError) as raises:
            foo()

        self.assertIn(
            "Cannot mutate a literal dictionary",
            str(raises.exception)
        )


if __name__ == '__main__':
    unittest.main()
