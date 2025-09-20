"""
Testing C implementation of the numba typed-list
"""

import ctypes
import struct

from numba.tests.support import TestCase
from numba import _helperlib


LIST_OK = 0
LIST_ERR_INDEX = -1
LIST_ERR_NO_MEMORY = -2
LIST_ERR_MUTATED = -3
LIST_ERR_ITER_EXHAUSTED = -4
LIST_ERR_IMMUTABLE = -5


class List(object):
    """A wrapper around the C-API to provide a minimal list object for
    testing.
    """
    def __init__(self, tc, item_size, allocated):
        """
        Parameters
        ----------
        tc : TestCase instance
        item_size : int
            byte size for the items
        allocated : int
            number of items to allocate for
        """
        self.tc = tc
        self.item_size = item_size
        self.lp = self.list_new(item_size, allocated)

    # The following methods implement part of the list API

    def __del__(self):
        self.tc.numba_list_free(self.lp)

    def __len__(self):
        return self.list_length()

    def __setitem__(self, i, item):
        return self.list_setitem(i, item)

    def __getitem__(self, i):
        return self.list_getitem(i)

    def __iter__(self):
        return ListIter(self)

    def __delitem__(self, i):
        self.list_delitem(i)

    def handle_index(self, i):
        # handling negative indices is done at the compiler level, so we only
        # support -1 to be last element of the list here
        if i < -1 or len(self) == 0:
            IndexError("list index out of range")
        elif i == -1:
            i = len(self) - 1
        return i

    @property
    def allocated(self):
        return self.list_allocated()

    @property
    def is_mutable(self):
        return self.list_is_mutable()

    def set_mutable(self):
        return self.list_set_is_mutable(1)

    def set_immutable(self):
        return self.list_set_is_mutable(0)

    def append(self, item):
        self.list_append(item)

    def pop(self, i=-1):
        return self.list_pop(i)

    # The methods below are higher-level wrappers for the C-API wrappers

    def list_new(self, item_size, allocated):
        lp = ctypes.c_void_p()
        status = self.tc.numba_list_new(
            ctypes.byref(lp), item_size, allocated,
        )
        self.tc.assertEqual(status, LIST_OK)
        return lp

    def list_length(self):
        return self.tc.numba_list_length(self.lp)

    def list_allocated(self):
        return self.tc.numba_list_allocated(self.lp)

    def list_is_mutable(self):
        return self.tc.numba_list_is_mutable(self.lp)

    def list_set_is_mutable(self, is_mutable):
        return self.tc.numba_list_set_is_mutable(self.lp, is_mutable)

    def list_setitem(self, i, item):
        status = self.tc.numba_list_setitem(self.lp, i, item)
        if status == LIST_ERR_INDEX:
            raise IndexError("list index out of range")
        elif status == LIST_ERR_IMMUTABLE:
            raise ValueError("list is immutable")
        else:
            self.tc.assertEqual(status, LIST_OK)

    def list_getitem(self, i):
        i = self.handle_index(i)
        item_out_buffer = ctypes.create_string_buffer(self.item_size)
        status = self.tc.numba_list_getitem(self.lp, i, item_out_buffer)
        if status == LIST_ERR_INDEX:
            raise IndexError("list index out of range")
        else:
            self.tc.assertEqual(status, LIST_OK)
            return item_out_buffer.raw

    def list_append(self, item):
        status = self.tc.numba_list_append(self.lp, item)
        if status == LIST_ERR_IMMUTABLE:
            raise ValueError("list is immutable")
        self.tc.assertEqual(status, LIST_OK)

    def list_pop(self, i):
        # pop is getitem and delitem
        i = self.handle_index(i)
        item = self.list_getitem(i)
        self.list_delitem(i)
        return item

    def list_delitem(self, i):
        # special case slice
        if isinstance(i, slice):
            status = self.tc.numba_list_delete_slice(self.lp,
                                                     i.start,
                                                     i.stop,
                                                     i.step)
            if status == LIST_ERR_IMMUTABLE:
                raise ValueError("list is immutable")
            self.tc.assertEqual(status, LIST_OK)
        # must be an integer, defer to delitem
        else:
            i = self.handle_index(i)
            status = self.tc.numba_list_delitem(self.lp, i)
            if status == LIST_ERR_INDEX:
                raise IndexError("list index out of range")
            elif status == LIST_ERR_IMMUTABLE:
                raise ValueError("list is immutable")
            self.tc.assertEqual(status, LIST_OK)

    def list_iter(self, itptr):
        self.tc.numba_list_iter(itptr, self.lp)

    def list_iter_next(self, itptr):
        bi = ctypes.c_void_p(0)
        status = self.tc.numba_list_iter_next(
            itptr, ctypes.byref(bi),
        )
        if status == LIST_ERR_MUTATED:
            raise ValueError('list mutated')
        elif status == LIST_ERR_ITER_EXHAUSTED:
            raise StopIteration
        else:
            self.tc.assertGreaterEqual(status, 0)
            item = (ctypes.c_char * self.item_size).from_address(bi.value)
            return item.value


class ListIter(object):
    """An iterator for the `List`.
    """
    def __init__(self, parent):
        self.parent = parent
        itsize = self.parent.tc.numba_list_iter_sizeof()
        self.it_state_buf = (ctypes.c_char_p * itsize)(0)
        self.it = ctypes.cast(self.it_state_buf, ctypes.c_void_p)
        self.parent.list_iter(self.it)

    def __iter__(self):
        return self

    def __next__(self):
        return self.parent.list_iter_next(self.it)

    next = __next__    # needed for py2 only


class TestListImpl(TestCase):
    def setUp(self):
        """Bind to the c_helper library and provide the ctypes wrapper.
        """
        list_t = ctypes.c_void_p
        iter_t = ctypes.c_void_p

        def wrap(name, restype, argtypes=()):
            proto = ctypes.CFUNCTYPE(restype, *argtypes)
            return proto(_helperlib.c_helpers[name])

        # numba_test_list()
        self.numba_test_list = wrap(
            'test_list',
            ctypes.c_int,
        )

        # numba_list_new(NB_List *l, Py_ssize_t item_size, Py_ssize_t allocated)
        self.numba_list_new = wrap(
            'list_new',
            ctypes.c_int,
            [ctypes.POINTER(list_t), ctypes.c_ssize_t, ctypes.c_ssize_t],
        )
        # numba_list_free(NB_List *l)
        self.numba_list_free = wrap(
            'list_free',
            None,
            [list_t],
        )
        # numba_list_length(NB_List *l)
        self.numba_list_length = wrap(
            'list_length',
            ctypes.c_int,
            [list_t],
        )
        # numba_list_allocated(NB_List *l)
        self.numba_list_allocated = wrap(
            'list_allocated',
            ctypes.c_int,
            [list_t],
        )
        # numba_list_is_mutable(NB_List *lp)
        self.numba_list_is_mutable = wrap(
            'list_is_mutable',
            ctypes.c_int,
            [list_t],
        )
        # numba_list_set_is_mutable(NB_List *lp, int is_mutable)
        self.numba_list_set_is_mutable = wrap(
            'list_set_is_mutable',
            None,
            [list_t, ctypes.c_int],
        )
        # numba_list_setitem(NB_List *l, Py_ssize_t i, const char *item)
        self.numba_list_setitem = wrap(
            'list_setitem',
            ctypes.c_int,
            [list_t, ctypes.c_ssize_t, ctypes.c_char_p],
        )
        # numba_list_append(NB_List *l, const char *item)
        self.numba_list_append = wrap(
            'list_append',
            ctypes.c_int,
            [list_t, ctypes.c_char_p],
        )
        # numba_list_getitem(NB_List *l,  Py_ssize_t i, char *out)
        self.numba_list_getitem = wrap(
            'list_getitem',
            ctypes.c_int,
            [list_t, ctypes.c_ssize_t, ctypes.c_char_p],
        )
        # numba_list_delitem(NB_List *l,  Py_ssize_t i)
        self.numba_list_delitem = wrap(
            'list_delitem',
            ctypes.c_int,
            [list_t, ctypes.c_ssize_t],
        )
        # numba_list_delete_slice(NB_List *l,
        #                         Py_ssize_t start,
        #                         Py_ssize_t stop,
        #                         Py_ssize_t step)
        self.numba_list_delete_slice = wrap(
            'list_delete_slice',
            ctypes.c_int,
            [list_t, ctypes.c_ssize_t, ctypes.c_ssize_t, ctypes.c_ssize_t],
        )
        # numba_list_iter_sizeof()
        self.numba_list_iter_sizeof = wrap(
            'list_iter_sizeof',
            ctypes.c_size_t,
        )
        # numba_list_iter(NB_ListIter *it, NB_List *l)
        self.numba_list_iter = wrap(
            'list_iter',
            None,
            [
                iter_t,
                list_t,
            ],
        )
        # numba_list_iter_next(NB_ListIter *it, const char **item_ptr)
        self.numba_list_iter_next = wrap(
            'list_iter_next',
            ctypes.c_int,
            [
                iter_t,                             # it
                ctypes.POINTER(ctypes.c_void_p),    # item_ptr
            ],
        )

    def test_simple_c_test(self):
        # Runs the basic test in C.
        ret = self.numba_test_list()
        self.assertEqual(ret, 0)

    def test_length(self):
        l = List(self, 8, 0)
        self.assertEqual(len(l), 0)

    def test_allocation(self):
        for i in range(16):
            l = List(self, 8, i)
            self.assertEqual(len(l), 0)
            self.assertEqual(l.allocated, i)

    def test_append_get_string(self):
        l = List(self, 8, 1)
        l.append(b"abcdefgh")
        self.assertEqual(len(l), 1)
        r = l[0]
        self.assertEqual(r, b"abcdefgh")

    def test_append_get_int(self):
        l = List(self, 8, 1)
        l.append(struct.pack("q", 1))
        self.assertEqual(len(l), 1)
        r = struct.unpack("q", l[0])[0]
        self.assertEqual(r, 1)

    def test_append_get_string_realloc(self):
        l = List(self, 8, 1)
        l.append(b"abcdefgh")
        self.assertEqual(len(l), 1)
        l.append(b"hijklmno")
        self.assertEqual(len(l), 2)
        r = l[1]
        self.assertEqual(r, b"hijklmno")

    def test_set_item_getitem_index_error(self):
        l = List(self, 8, 0)
        with self.assertRaises(IndexError):
            l[0]
        with self.assertRaises(IndexError):
            l[0] = b"abcdefgh"

    def test_iter(self):
        l = List(self, 1, 0)
        values = [b'a', b'b', b'c', b'd', b'e', b'f', b'g', b'h']
        for i in values:
            l.append(i)
        received = []
        for j in l:
            received.append(j)
        self.assertEqual(values, received)

    def test_pop(self):
        l = List(self, 1, 0)
        values = [b'a', b'b', b'c', b'd', b'e', b'f', b'g', b'h']
        for i in values:
            l.append(i)
        self.assertEqual(len(l), 8)

        received = l.pop()
        self.assertEqual(b'h', received)
        self.assertEqual(len(l), 7)
        received = [j for j in l]
        self.assertEqual(received, values[:-1])

        received = l.pop(0)
        self.assertEqual(b'a', received)
        self.assertEqual(len(l), 6)

        received = l.pop(2)
        self.assertEqual(b'd', received)
        self.assertEqual(len(l), 5)

        expected = [b'b', b'c', b'e', b'f', b'g']
        received = [j for j in l]
        self.assertEqual(received, expected)

    def test_pop_index_error(self):
        l = List(self, 8, 0)
        with self.assertRaises(IndexError):
            l.pop()

    def test_pop_byte(self):
        l = List(self, 4, 0)
        values = [b'aaaa', b'bbbb', b'cccc', b'dddd',
                  b'eeee', b'ffff', b'gggg', b'hhhhh']
        for i in values:
            l.append(i)
        self.assertEqual(len(l), 8)

        received = l.pop()
        self.assertEqual(b'hhhh', received)
        self.assertEqual(len(l), 7)
        received = [j for j in l]
        self.assertEqual(received, values[:-1])

        received = l.pop(0)
        self.assertEqual(b'aaaa', received)
        self.assertEqual(len(l), 6)

        received = l.pop(2)
        self.assertEqual(b'dddd', received)
        self.assertEqual(len(l), 5)

        expected = [b'bbbb', b'cccc', b'eeee', b'ffff', b'gggg']
        received = [j for j in l]
        self.assertEqual(received, expected)

    def test_delitem(self):
        l = List(self, 1, 0)
        values = [b'a', b'b', b'c', b'd', b'e', b'f', b'g', b'h']
        for i in values:
            l.append(i)
        self.assertEqual(len(l), 8)

        # delete first item
        del l[0]
        self.assertEqual(len(l), 7)
        self.assertEqual(list(l), values[1:])
        # delete last item
        del l[-1]
        self.assertEqual(len(l), 6)
        self.assertEqual(list(l), values[1:-1])
        # delete item from middle
        del l[2]
        self.assertEqual(len(l), 5)
        self.assertEqual(list(l), [b'b', b'c', b'e', b'f', b'g'])

    def test_delete_slice(self):
        l = List(self, 1, 0)
        values = [b'a', b'b', b'c', b'd', b'e', b'f', b'g', b'h']
        for i in values:
            l.append(i)
        self.assertEqual(len(l), 8)

        # delete every second item
        # no slice default normalization here, be explicit about start anb stop
        del l[0:8:2]
        self.assertEqual(len(l), 4)
        self.assertEqual(list(l), values[1:8:2])

        # delete first item
        del l[0:1:1]
        self.assertEqual(len(l), 3)
        self.assertEqual(list(l), [b'd', b'f', b'h'])

        # delete last item
        del l[2:3:1]
        self.assertEqual(len(l), 2)
        self.assertEqual(list(l), [b'd', b'f'])

        # delete all left items
        del l[0:2:1]
        self.assertEqual(len(l), 0)
        self.assertEqual(list(l), [])

    def check_sizing(self, item_size, nmax):
        # Helper to verify different item_sizes
        l = List(self, item_size, 0)

        def make_item(v):
            tmp = "{:0{}}".format(nmax - v - 1, item_size).encode("latin-1")
            return tmp[:item_size]

        for i in range(nmax):
            l.append(make_item(i))

        self.assertEqual(len(l), nmax)

        for i in range(nmax):
            self.assertEqual(l[i], make_item(i))

    def test_sizing(self):
        # Check different sizes of the key & value.
        for i in range(1, 16):
            self.check_sizing(item_size=i, nmax=2**i)

    def test_mutability(self):
        # setup and populate a singleton
        l = List(self, 8, 1)
        one = struct.pack("q", 1)
        l.append(one)
        self.assertTrue(l.is_mutable)
        self.assertEqual(len(l), 1)
        r = struct.unpack("q", l[0])[0]
        self.assertEqual(r, 1)

        # set to immutable and test guards
        l.set_immutable()
        self.assertFalse(l.is_mutable)
        # append
        with self.assertRaises(ValueError) as raises:
            l.append(one)
        self.assertIn("list is immutable", str(raises.exception))
        # setitem
        with self.assertRaises(ValueError) as raises:
            l[0] = one
        self.assertIn("list is immutable", str(raises.exception))
        # pop
        with self.assertRaises(ValueError) as raises:
            l.pop()
        self.assertIn("list is immutable", str(raises.exception))
        # delitem with index
        with self.assertRaises(ValueError) as raises:
            del l[0]
        self.assertIn("list is immutable", str(raises.exception))
        # delitem with slice
        with self.assertRaises(ValueError) as raises:
            del l[0:1:1]
        self.assertIn("list is immutable", str(raises.exception))
        l.set_mutable()

        # check that nothing has changed
        self.assertTrue(l.is_mutable)
        self.assertEqual(len(l), 1)
        r = struct.unpack("q", l[0])[0]
        self.assertEqual(r, 1)
