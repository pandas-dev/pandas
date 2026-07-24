"""
Testing C implementation of the numba sets
"""

import ctypes

from numba.tests.support import TestCase
from numba import _helperlib


ENTRY_PRESENT = 1
OK = 0
ERR_KEY_NOT_FOUND = -1
ERR_SET_MUTATED = -2
ERR_ITER_EXHAUSTED = -3
ERR_SET_EMPTY = -4
ERR_CMP_FAILED = -5


def to_bytes(key, key_size):
    if isinstance(key, str):
        key_bytes = (key_size - len(key)) * b'\0' + key.encode()
    else:
        key_bytes = key.to_bytes(key_size, 'big')
    return key_bytes


class SetIter(object):
    """An iterator for the `Set`.
    """
    def __init__(self, parent):
        self.parent = parent
        itsize = self.parent.tc.numba_set_iter_sizeof()
        self.it_state_buf = (ctypes.c_char_p * itsize)(0)
        self.it = ctypes.cast(self.it_state_buf, ctypes.c_void_p)
        self.parent.set_iter(self.it)

    def __iter__(self):
        return self

    def __next__(self):
        return self.parent.set_iter_next(self.it)


class Set(object):
    """A wrapper around the C-API to provide a minimal set object for
    testing.
    """
    def __init__(self, tc, key_size, allocated):
        """
        Parameters
        ----------
        tc : TestCase instance
        key_size : int
            byte size for the items
        allocated : int
            number of items to allocate for
        """
        self.tc = tc
        self.key_size = key_size
        self.setp = self.set_new(key_size, allocated)

    # The following methods implement part of the set API

    def __len__(self):
        return self.set_length()

    def __contains__(self, key):
        key_bytes = to_bytes(key, self.key_size)
        return self.set_contains(key_bytes=key_bytes)

    def __iter__(self):
        return SetIter(self)

    def __del__(self):
        self.tc.numba_set_free(self.setp)

    def add(self, key):
        key_bytes = to_bytes(key, self.key_size)
        self.set_add(key_bytes)

    def discard(self, key):
        key_bytes = to_bytes(key, self.key_size)
        self.set_discard(key_bytes)

    # The methods below are higher-level wrappers for the C-API wrappers

    def set_new(self, key_size, allocated):
        setp = ctypes.c_void_p()
        status = self.tc.numba_set_new(
            ctypes.byref(setp), key_size, allocated,
        )
        self.tc.assertEqual(status, OK)
        return setp

    def set_length(self):
        return self.tc.numba_set_length(self.setp)

    def set_add(self, key):
        status = self.tc.numba_set_add(self.setp, key, hash(key))
        self.tc.assertEqual(status, OK)

    def set_contains(self, key_bytes):
        hash_key = hash(key_bytes)
        ix = self.tc.numba_set_contains(
            self.setp, key_bytes, hash_key
        )
        return ix

    def set_discard(self, key_bytes):
        hash_key = hash(key_bytes)
        ix = self.tc.numba_set_discard(
            self.setp, key_bytes, hash_key
        )
        return ix

    def set_iter(self, itptr):
        self.tc.numba_set_iter(itptr, self.setp)

    def set_iter_next(self, itptr):
        bi = ctypes.c_void_p(0)
        status = self.tc.numba_set_iter_next(
            itptr, ctypes.byref(bi),
        )
        if status == ERR_SET_MUTATED:
            raise ValueError('set mutated')
        elif status == ERR_ITER_EXHAUSTED:
            raise StopIteration
        else:
            self.tc.assertGreaterEqual(status, 0)
            item = (ctypes.c_char * self.key_size).from_address(bi.value)
            return item.value


class TestSetImpl(TestCase):
    def setUp(self):
        """Bind to the c_helper library and provide the ctypes wrapper.
        """
        set_ty = ctypes.c_void_p
        hash_ty = ctypes.c_ssize_t
        setiter_ty = ctypes.c_void_p

        def wrap(name, restype, argtypes=()):
            proto = ctypes.CFUNCTYPE(restype, *argtypes)
            return proto(_helperlib.c_helpers[name])

        # numba_set_new(NB_set *setp, Py_ssize_t key_size, Py_ssize_t allocated)
        self.numba_set_new = wrap(
            'set_new',
            ctypes.c_int,
            [ctypes.POINTER(set_ty), ctypes.c_ssize_t, ctypes.c_ssize_t],
        )

        # numba_test_set()
        self.numba_test_set = wrap(
            'test_set',
            ctypes.c_int,
        )
        # numba_set_length(NB_set *setp)
        self.numba_set_length = wrap(
            'set_length',
            ctypes.c_int,
            [set_ty],
        )
        # numba_set_add(NB_set *setp, char *key, Py_ssize_t hash)
        self.numba_set_add = wrap(
            'set_add',
            ctypes.c_int,
            [set_ty, ctypes.c_char_p, hash_ty],
        )
        # numba_set_contains(NB_Set *setp, char *key, Py_ssize_t hash)
        self.numba_set_contains = wrap(
            'set_contains',
            ctypes.c_int,
            [set_ty, ctypes.c_char_p, hash_ty],
        )
        # numba_set_discard(NB_Set *setp, char *key, Py_ssize_t hash)
        self.numba_set_discard = wrap(
            'set_discard',
            ctypes.c_int,
            [set_ty, ctypes.c_char_p, hash_ty],
        )
        # numba_set_iter_sizeof()
        self.numba_set_iter_sizeof = wrap(
            'set_iter_sizeof',
            ctypes.c_size_t,
        )
        # numba_set_iter(NB_SetIter *it, NB_Set *setp)
        self.numba_set_iter = wrap(
            'set_iter',
            None,
            [
                setiter_ty,
                set_ty,
            ],
        )
        # numba_set_iter_next(NB_SetIter *it, const char **key_ptr)
        self.numba_set_iter_next = wrap(
            'set_iter_next',
            ctypes.c_int,
            [
                setiter_ty,                             # it
                ctypes.POINTER(ctypes.c_void_p),    # item_ptr
            ],
        )
        # numba_set_free(NB_Set *setp)
        self.numba_set_free = wrap(
            'set_free',
            None,
            [set_ty],
        )

    def test_simple_c_test(self):
        # Runs the basic test in C.
        ret = self.numba_test_set()
        self.assertEqual(ret, OK)

    def test_init_len(self):
        # Checks length at initialization
        s = Set(self, 8, 8)

        # length == zero
        self.assertEqual(len(s), 0)

    def test_add_check(self):
        # Test adding keys to a set

        for key in [3, 'a', "abc"]:
            s = Set(self, 8, 8)
            # Add key to the set
            s.add(key)
            # Length of set == 1
            self.assertEqual(len(s), 1)
            # Check key exists in set
            has_entry = s.__contains__(key)
            self.assertEqual(has_entry, True)

    def test_remove_check(self):
        # Test removal of keys from a set

        for key in [3, 'a', "abc"]:
            s = Set(self, 8, 8)
            # Add key to the set
            s.add(key)
            # Check key exists in set
            has_entry = s.__contains__(key)
            self.assertEqual(has_entry, True)
            # Discard the key
            s.discard(key)
            # Check key no longer exists in set
            has_entry = s.__contains__(key)
            self.assertEqual(has_entry, False)

    def test_set_iter(self):
        # Test iteration
        s = Set(self, 3, 8)
        nmax = 1000

        # Add elements to the set
        for i in range(100, nmax):
            s.add(f'{i}')

        keys = []
        # Iterate over the set
        for key in s:
            keys.append(key)

        # Check every key was present in the iteration
        self.assertEqual(len(s), len(keys))
        for i in range(100, nmax):
            self.assertIn(to_bytes(f'{i}', 3), keys)
