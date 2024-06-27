"""
Testing C implementation of the numba dictionary
"""

import ctypes
import random

from numba.tests.support import TestCase
from numba import _helperlib, jit, typed, types
from numba.core.config import IS_32BITS
from numba.core.datamodel.models import UniTupleModel
from numba.extending import register_model, typeof_impl, unbox, overload


DKIX_EMPTY = -1


ALIGN = 4 if IS_32BITS else 8


class Dict(object):
    """A wrapper around the C-API to provide a minimal dictionary object for
    testing.
    """
    def __init__(self, tc, keysize, valsize):
        """
        Parameters
        ----------
        tc : TestCase instance
        keysize : int
            byte size for the key
        valsize : int
            byte size for the value
        """
        self.tc = tc
        self.keysize = keysize
        self.valsize = valsize
        self.dp = self.dict_new_minsize(keysize, valsize)

    def __del__(self):
        self.tc.numba_dict_free(self.dp)

    def __len__(self):
        return self.dict_length()

    def __setitem__(self, k, v):
        bk = bytes(k.encode())
        bv = bytes(v.encode())
        self.tc.assertEqual(len(bk), self.keysize)
        self.tc.assertEqual(len(bv), self.valsize)
        self.dict_insert(bk, bv)

    def __getitem__(self, k):
        bk = bytes(k.encode())
        self.tc.assertEqual(len(bk), self.keysize)
        ix, old = self.dict_lookup(bk)
        if ix == DKIX_EMPTY:
            raise KeyError
        else:
            return old.decode()

    def __delitem__(self, k):
        bk = bytes(k.encode())
        self.tc.assertEqual(len(bk), self.keysize)
        if not self.dict_delitem(bk):
            raise KeyError(k)

    def get(self, k):
        try:
            return self[k]
        except KeyError:
            return

    def items(self):
        return DictIter(self)

    def popitem(self):
        k, v = self.dict_popitem()
        return k.decode(), v.decode()

    #
    # The methods below are higher-level wrappers for the C-API
    #

    def dict_new_minsize(self, key_size, val_size):
        dp = ctypes.c_void_p()
        status = self.tc.numba_dict_new_sized(
            ctypes.byref(dp), 0, key_size, val_size,
        )
        self.tc.assertEqual(status, 0)
        return dp

    def dict_length(self):
        return self.tc.numba_dict_length(self.dp)

    def dict_insert(self, key_bytes, val_bytes):
        hashval = hash(key_bytes)
        status = self.tc.numba_dict_insert_ez(
            self.dp, key_bytes, hashval, val_bytes,
        )
        self.tc.assertGreaterEqual(status, 0)

    def dict_lookup(self, key_bytes):
        hashval = hash(key_bytes)
        oldval_bytes = ctypes.create_string_buffer(self.valsize)
        ix = self.tc.numba_dict_lookup(
            self.dp, key_bytes, hashval, oldval_bytes,
        )
        self.tc.assertGreaterEqual(ix, DKIX_EMPTY)
        return ix, oldval_bytes.value

    def dict_delitem(self, key_bytes):
        ix, oldval = self.dict_lookup(key_bytes)
        if ix == DKIX_EMPTY:
            return False
        hashval = hash(key_bytes)
        status = self.tc.numba_dict_delitem(self.dp, hashval, ix)
        self.tc.assertEqual(status, 0)
        return True

    def dict_popitem(self):
        key_bytes = ctypes.create_string_buffer(self.keysize)
        val_bytes = ctypes.create_string_buffer(self.valsize)
        status = self.tc.numba_dict_popitem(self.dp, key_bytes, val_bytes)
        if status != 0:
            if status == -4:
                raise KeyError('popitem(): dictionary is empty')
            else:
                self.tc._fail('Unknown')
        return key_bytes.value, val_bytes.value

    def dict_iter(self, itptr):
        self.tc.numba_dict_iter(itptr, self.dp)

    def dict_iter_next(self, itptr):
        bk = ctypes.c_void_p(0)
        bv = ctypes.c_void_p(0)
        status = self.tc.numba_dict_iter_next(
            itptr, ctypes.byref(bk), ctypes.byref(bv),
        )
        if status == -2:
            raise ValueError('dictionary mutated')
        elif status == -3:
            return
        else:
            self.tc.assertGreaterEqual(status, 0)

            # Check the alignment of the key-value in the entries.
            # We know we are getting the pointers to data in the entries.

            self.tc.assertEqual(bk.value % ALIGN, 0, msg='key not aligned')
            self.tc.assertEqual(bv.value % ALIGN, 0, msg='val not aligned')

            key = (ctypes.c_char * self.keysize).from_address(bk.value)
            val = (ctypes.c_char * self.valsize).from_address(bv.value)
            return key.value, val.value


class DictIter(object):
    """A iterator for the `Dict.items()`.

    Only the `.items()` is needed.  `.keys` and `.values` can be trivially
    implemented on the `.items` iterator.
    """
    def __init__(self, parent):
        self.parent = parent
        itsize = self.parent.tc.numba_dict_iter_sizeof()
        self.it_state_buf = (ctypes.c_char_p * itsize)(0)
        self.it = ctypes.cast(self.it_state_buf, ctypes.c_void_p)
        self.parent.dict_iter(self.it)

    def __iter__(self):
        return self

    def __next__(self):
        out = self.parent.dict_iter_next(self.it)
        if out is None:
            raise StopIteration
        else:
            k, v = out
            return k.decode(), v.decode()

    next = __next__    # needed for py2 only


class Parametrized(tuple):
    """supporting type for TestDictImpl.test_parametrized_types
    needs to be global to be cacheable"""
    def __init__(self, tup):
        assert all(isinstance(v, str) for v in tup)


class ParametrizedType(types.Type):
    """this is essentially UniTuple(unicode_type, n)
    BUT type name is the same for all n"""

    def __init__(self, value):
        super(ParametrizedType, self).__init__('ParametrizedType')
        self.dtype = types.unicode_type
        self.n = len(value)

    @property
    def key(self):
        return self.n

    def __len__(self):
        return self.n


class TestDictImpl(TestCase):
    def setUp(self):
        """Bind to the c_helper library and provide the ctypes wrapper.
        """
        dict_t = ctypes.c_void_p
        iter_t = ctypes.c_void_p
        hash_t = ctypes.c_ssize_t

        def wrap(name, restype, argtypes=()):
            proto = ctypes.CFUNCTYPE(restype, *argtypes)
            return proto(_helperlib.c_helpers[name])

        # numba_test_dict()
        self.numba_test_dict = wrap(
            'test_dict',
            ctypes.c_int,
        )
        # numba_dict_new_sized(
        #    NB_Dict **out,
        #    Py_ssize_t n_keys,
        #    Py_ssize_t key_size,
        #    Py_ssize_t val_size
        # )
        self.numba_dict_new_sized = wrap(
            'dict_new_sized',
            ctypes.c_int,
            [
                ctypes.POINTER(dict_t),  # out
                ctypes.c_ssize_t,        # n_keys
                ctypes.c_ssize_t,        # key_size
                ctypes.c_ssize_t,        # val_size
            ],
        )
        # numba_dict_free(NB_Dict *d)
        self.numba_dict_free = wrap(
            'dict_free',
            None,
            [dict_t],
        )
        # numba_dict_length(NB_Dict *d)
        self.numba_dict_length = wrap(
            'dict_length',
            ctypes.c_ssize_t,
            [dict_t],
        )
        # numba_dict_insert_ez(
        #     NB_Dict    *d,
        #     const char *key_bytes,
        #     Py_hash_t   hash,
        #     const char *val_bytes,
        #     )
        self.numba_dict_insert_ez = wrap(
            'dict_insert_ez',
            ctypes.c_int,
            [
                dict_t,             # d
                ctypes.c_char_p,    # key_bytes
                hash_t,             # hash
                ctypes.c_char_p,    # val_bytes
            ],
        )
        # numba_dict_lookup(
        #       NB_Dict *d,
        #       const char *key_bytes,
        #       Py_hash_t hash,
        #       char *oldval_bytes
        # )
        self.numba_dict_lookup = wrap(
            'dict_lookup',
            ctypes.c_ssize_t,
            [
                dict_t,             # d
                ctypes.c_char_p,    # key_bytes
                hash_t,             # hash
                ctypes.c_char_p,    # oldval_bytes
            ],
        )
        # numba_dict_delitem(
        #     NB_Dict *d,
        #     Py_hash_t hash,
        #     Py_ssize_t ix
        # )
        self.numba_dict_delitem = wrap(
            'dict_delitem',
            ctypes.c_int,
            [
                dict_t,             # d
                hash_t,             # hash
                ctypes.c_ssize_t,   # ix
            ],
        )
        # numba_dict_popitem(
        #   NB_Dict *d,
        #   char *key_bytes,
        #   char *val_bytes
        # )
        self.numba_dict_popitem = wrap(
            'dict_popitem',
            ctypes.c_int,
            [
                dict_t,             # d
                ctypes.c_char_p,    # key_bytes
                ctypes.c_char_p,    # val_bytes
            ],
        )
        # numba_dict_iter_sizeof()
        self.numba_dict_iter_sizeof = wrap(
            'dict_iter_sizeof',
            ctypes.c_size_t,
        )
        # numba_dict_iter(
        #     NB_DictIter *it,
        #     NB_Dict     *d
        # )
        self.numba_dict_iter = wrap(
            'dict_iter',
            None,
            [
                iter_t,
                dict_t,
            ],
        )
        # numba_dict_iter_next(
        #     NB_DictIter *it,
        #     const char **key_ptr,
        #     const char **val_ptr
        # )
        self.numba_dict_iter_next = wrap(
            'dict_iter_next',
            ctypes.c_int,
            [
                iter_t,                             # it
                ctypes.POINTER(ctypes.c_void_p),    # key_ptr
                ctypes.POINTER(ctypes.c_void_p),    # val_ptr
            ],
        )

    def test_simple_c_test(self):
        # Runs the basic test in C.
        ret = self.numba_test_dict()
        self.assertEqual(ret, 0)

    def test_insertion_small(self):
        # Tests insertion and lookup for a small dict.
        d = Dict(self, 4, 8)
        self.assertEqual(len(d), 0)
        self.assertIsNone(d.get('abcd'))

        # First key
        d['abcd'] = 'beefcafe'
        self.assertEqual(len(d), 1)
        self.assertIsNotNone(d.get('abcd'))
        self.assertEqual(d['abcd'], 'beefcafe')

        # Duplicated key replaces
        d['abcd'] = 'cafe0000'
        self.assertEqual(len(d), 1)
        self.assertEqual(d['abcd'], 'cafe0000')

        # Second key
        d['abce'] = 'cafe0001'
        self.assertEqual(len(d), 2)
        self.assertEqual(d['abcd'], 'cafe0000')
        self.assertEqual(d['abce'], 'cafe0001')

        # Third key
        d['abcf'] = 'cafe0002'
        self.assertEqual(len(d), 3)
        self.assertEqual(d['abcd'], 'cafe0000')
        self.assertEqual(d['abce'], 'cafe0001')
        self.assertEqual(d['abcf'], 'cafe0002')

    def check_insertion_many(self, nmax):
        # Helper to test insertion/lookup/resize
        d = Dict(self, 8, 8)

        def make_key(v):
            return "key_{:04}".format(v)

        def make_val(v):
            return "val_{:04}".format(v)

        # Check insert
        for i in range(nmax):
            d[make_key(i)] = make_val(i)
            self.assertEqual(len(d), i + 1)

        # Check lookup
        for i in range(nmax):
            self.assertEqual(d[make_key(i)], make_val(i))

    def test_insertion_many(self):
        # Test insertion for differently sized dict
        # Around minsize
        self.check_insertion_many(nmax=7)
        self.check_insertion_many(nmax=8)
        self.check_insertion_many(nmax=9)
        # Around nmax = 32
        self.check_insertion_many(nmax=31)
        self.check_insertion_many(nmax=32)
        self.check_insertion_many(nmax=33)
        # Around nmax = 1024
        self.check_insertion_many(nmax=1023)
        self.check_insertion_many(nmax=1024)
        self.check_insertion_many(nmax=1025)
        # Around nmax = 4096
        self.check_insertion_many(nmax=4095)
        self.check_insertion_many(nmax=4096)
        self.check_insertion_many(nmax=4097)

    def test_deletion_small(self):
        # Test deletion
        d = Dict(self, 4, 8)
        self.assertEqual(len(d), 0)
        self.assertIsNone(d.get('abcd'))

        d['abcd'] = 'cafe0000'
        d['abce'] = 'cafe0001'
        d['abcf'] = 'cafe0002'

        self.assertEqual(len(d), 3)
        self.assertEqual(d['abcd'], 'cafe0000')
        self.assertEqual(d['abce'], 'cafe0001')
        self.assertEqual(d['abcf'], 'cafe0002')
        self.assertEqual(len(d), 3)

        # Delete first item
        del d['abcd']
        self.assertIsNone(d.get('abcd'))
        self.assertEqual(d['abce'], 'cafe0001')
        self.assertEqual(d['abcf'], 'cafe0002')
        self.assertEqual(len(d), 2)

        # Delete first item again
        with self.assertRaises(KeyError):
            del d['abcd']

        # Delete third
        del d['abcf']
        self.assertIsNone(d.get('abcd'))
        self.assertEqual(d['abce'], 'cafe0001')
        self.assertIsNone(d.get('abcf'))
        self.assertEqual(len(d), 1)

        # Delete second
        del d['abce']
        self.assertIsNone(d.get('abcd'))
        self.assertIsNone(d.get('abce'))
        self.assertIsNone(d.get('abcf'))
        self.assertEqual(len(d), 0)

    def check_delete_randomly(self, nmax, ndrop, nrefill, seed=0):
        # Helper to test deletion
        random.seed(seed)

        d = Dict(self, 8, 8)
        keys = {}

        def make_key(v):
            return "k_{:06x}".format(v)

        def make_val(v):
            return "v_{:06x}".format(v)

        for i in range(nmax):
            d[make_key(i)] = make_val(i)

        # Fill to nmax
        for i in range(nmax):
            k = make_key(i)
            v = make_val(i)
            keys[k] = v
            self.assertEqual(d[k], v)

        self.assertEqual(len(d), nmax)

        # Randomly drop
        droplist = random.sample(list(keys), ndrop)
        remain = keys.copy()
        for i, k in enumerate(droplist, start=1):
            del d[k]
            del remain[k]
            self.assertEqual(len(d), nmax - i)
        self.assertEqual(len(d), nmax - ndrop)

        # Make sure everything dropped is gone
        for k in droplist:
            self.assertIsNone(d.get(k))

        # Make sure everything else is still here
        for k in remain:
            self.assertEqual(d[k], remain[k])

        # Refill
        for i in range(nrefill):
            k = make_key(nmax + i)
            v = make_val(nmax + i)
            remain[k] = v
            d[k] = v

        self.assertEqual(len(remain), len(d))

        # Make sure everything is here
        for k in remain:
            self.assertEqual(d[k], remain[k])

    def test_delete_randomly(self):
        # Test deletion for differently sized dict
        self.check_delete_randomly(nmax=8, ndrop=2, nrefill=2)
        self.check_delete_randomly(nmax=13, ndrop=10, nrefill=31)
        self.check_delete_randomly(nmax=100, ndrop=50, nrefill=200)
        self.check_delete_randomly(nmax=100, ndrop=99, nrefill=100)
        self.check_delete_randomly(nmax=100, ndrop=100, nrefill=100)
        self.check_delete_randomly(nmax=1024, ndrop=999, nrefill=1)
        self.check_delete_randomly(nmax=1024, ndrop=999, nrefill=2048)

    def test_delete_randomly_large(self):
        # Go beyond 2^16 to exercise large indices.
        # Internally, size of index changes as the hashtable size changes.
        # Size of index can be 8, 16, 32 or 64 bytes (on 64-bit).
        # We are not inserting >2^32 elements because of limitation of time.
        self.check_delete_randomly(nmax=2**17, ndrop=2**16, nrefill=2**10)

    def test_popitem(self):
        nmax = 10
        d = Dict(self, 8, 8)

        def make_key(v):
            return "k_{:06x}".format(v)

        def make_val(v):
            return "v_{:06x}".format(v)

        for i in range(nmax):
            d[make_key(i)] = make_val(i)

        self.assertEqual(len(d), nmax)
        k, v = d.popitem()
        self.assertEqual(len(d), nmax - 1)
        self.assertEqual(k, make_key(len(d)))
        self.assertEqual(v, make_val(len(d)))

        while len(d):
            n = len(d)
            k, v = d.popitem()
            self.assertEqual(len(d), n - 1)
            self.assertEqual(k, make_key(len(d)))
            self.assertEqual(v, make_val(len(d)))

        self.assertEqual(len(d), 0)
        with self.assertRaises(KeyError) as raises:
            d.popitem()
        self.assertIn(
            'popitem(): dictionary is empty',
            str(raises.exception),
        )

    def test_iter_items(self):
        # Test .items iteration
        d = Dict(self, 4, 4)
        nmax = 1000

        def make_key(v):
            return "{:04}".format(v)

        def make_val(v):
            return "{:04}".format(v + nmax)

        for i in range(nmax):
            d[make_key(i)] = make_val(i)

        # Check that the everything is ordered
        for i, (k, v) in enumerate(d.items()):
            self.assertEqual(make_key(i), k)
            self.assertEqual(make_val(i), v)

    def check_sizing(self, key_size, val_size, nmax):
        # Helper to verify different key/value sizes.
        d = Dict(self, key_size, val_size)

        def make_key(v):
            return "{:0{}}".format(v, key_size)[:key_size]

        def make_val(v):
            return "{:0{}}".format(nmax - v - 1, val_size)[:val_size]

        for i in range(nmax):
            d[make_key(i)] = make_val(i)

        # Check that the everything is ordered
        for i, (k, v) in enumerate(d.items()):
            self.assertEqual(make_key(i), k)
            self.assertEqual(make_val(i), v)

    def test_sizing(self):
        # Check different sizes of the key & value.
        for i in range(1, 8):
            self.check_sizing(key_size=i, val_size=i, nmax=2**i)

    def test_parameterized_types(self):
        """https://github.com/numba/numba/issues/6401"""

        register_model(ParametrizedType)(UniTupleModel)

        @typeof_impl.register(Parametrized)
        def typeof_unit(val, c):
            return ParametrizedType(val)

        @unbox(ParametrizedType)
        def unbox_parametrized(typ, obj, context):
            return context.unbox(types.UniTuple(typ.dtype, len(typ)), obj)

        def dict_vs_cache_vs_parametrized(v):
            assert 0

        @overload(dict_vs_cache_vs_parametrized)
        def ol_dict_vs_cache_vs_parametrized(v):
            typ = v

            def objmode_vs_cache_vs_parametrized_impl(v):
                # typed.List shows same behaviour after fix for #6397
                d = typed.Dict.empty(types.unicode_type, typ)
                d['data'] = v

            return objmode_vs_cache_vs_parametrized_impl

        @jit(nopython=True, cache=True)
        def set_parametrized_data(x, y):
            # Has had a tendency to segfault when the compiled function
            # was loaded from cache in a different process than the one
            # it was originally compiled in.
            # The new process is simulated below by resetting the dispatchers
            # and the target context
            dict_vs_cache_vs_parametrized(x)
            dict_vs_cache_vs_parametrized(y)

        x, y = Parametrized(('a', 'b')), Parametrized(('a',))
        set_parametrized_data(x, y)

        # reset dispatchers and targetctx to force re-load from cache as if a
        # new process would jit the function
        set_parametrized_data._make_finalizer()()
        set_parametrized_data._reset_overloads()
        set_parametrized_data.targetctx.init()

        for ii in range(50):  # <- sometimes works a few times
            self.assertIsNone(set_parametrized_data(x, y))
