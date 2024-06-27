import array

import numpy as np

from numba import jit
from numba.tests.support import TestCase, compile_function, MemoryLeakMixin
import unittest


@jit(nopython=True)
def len_usecase(buf):
    return len(buf)


@jit(nopython=True)
def getitem_usecase(buf, i):
    return buf[i]


@jit(nopython=True)
def getslice_usecase(buf, i, j):
    s = buf[i:j]
    return s[0] + 2 * s[-1]


@jit(nopython=True)
def setitem_usecase(buf, i, v):
    buf[i] = v


@jit(nopython=True)
def iter_usecase(buf):
    res = 0.0
    for i, x in enumerate(buf):
        res += x
        res *= i + 1
    return res


def attrgetter(attr):
    code = """def func(x):
        return x.%(attr)s
""" % locals()
    pyfunc = compile_function("func", code, globals())
    return jit(nopython=True)(pyfunc)


contiguous_usecase = attrgetter("contiguous")
c_contiguous_usecase = attrgetter("c_contiguous")
f_contiguous_usecase = attrgetter("f_contiguous")
itemsize_usecase = attrgetter("itemsize")
nbytes_usecase = attrgetter("nbytes")
ndim_usecase = attrgetter("ndim")
readonly_usecase = attrgetter("readonly")
shape_usecase = attrgetter("shape")
strides_usecase = attrgetter("strides")


class TestBufferProtocol(MemoryLeakMixin, TestCase):
    """
    Test operations on buffer-providing objects.
    """

    def _arrays(self):
        n = 10
        for letter, offset in [
            ('b', -3),
            ('B', 0),
            ('h', -5000),
            ('H', 40000),
            ('i', -100000),
            ('I', 1000000),
            ('l', -100000),
            ('L', 1000000),
            ('q', -2**60),
            ('Q', 2**63 + 1),
            ('f', 1.5),
            ('d', -1.5),
            ]:
            yield array.array(letter, [i + offset for i in range(n)])

    def _memoryviews(self):
        n = 10
        yield memoryview(bytearray(b"abcdefghi"))
        yield memoryview(b"abcdefghi")
        # Different item types
        for dtype, start, stop in [
            ('int8', -10, 10),
            ('uint8', 0, 10),
            ('int16', -5000, 1000),
            ('uint16', 40000, 50000),
            ('int32', -100000, 100000),
            ('uint32', 0, 1000000),
            ('int64', -2**60, 10),
            ('uint64', 0, 2**64 - 10),
            ('float32', 1.5, 3.5),
            ('float64', 1.5, 3.5),
            ('complex64', -8j, 12 + 5j),
            ('complex128', -8j, 12 + 5j),
            ]:
            yield memoryview(np.linspace(start, stop, n).astype(dtype))
        # Different layouts
        arr = np.arange(12).reshape((3, 4))
        assert arr.flags.c_contiguous and not arr.flags.f_contiguous
        yield memoryview(arr)
        arr = arr.T
        assert arr.flags.f_contiguous and not arr.flags.c_contiguous
        yield memoryview(arr)
        arr = arr[::2]
        assert not arr.flags.f_contiguous and not arr.flags.c_contiguous
        yield memoryview(arr)

    def _readonlies(self):
        yield b"xyz"
        yield memoryview(b"abcdefghi")
        arr = np.arange(5)
        arr.setflags(write=False)
        yield memoryview(arr)

    def _check_unary(self, jitfunc, *args):
        pyfunc = jitfunc.py_func
        self.assertPreciseEqual(jitfunc(*args), pyfunc(*args))

    def check_len(self, obj):
        self._check_unary(len_usecase, obj)

    def check_iter(self, obj):
        self._check_unary(iter_usecase, obj)

    def check_getitem(self, obj):
        # Be careful to index all dimensions, since we don't support
        # partial indexing yet.
        def yield_indices(obj):
            try:
                shape = obj.shape
            except AttributeError:
                shape = len(obj),
            for tup in np.ndindex(shape):
                # Simple 1d buffer-providing objects usually don't support
                # tuple indexing.
                if len(tup) == 1:
                    yield tup[0]
                else:
                    yield tup

        for i in yield_indices(obj):
            try:
                expected = obj[i]
            except (NotImplementedError, TypeError):
                if isinstance(obj, memoryview):
                    # The memoryview object doesn't support all codes yet,
                    # fall back on the underlying object.
                    expected = obj.obj[i]
                else:
                    raise
            self.assertPreciseEqual(getitem_usecase(obj, i), expected)

    def check_setitem(self, obj):
        for i in range(len(obj)):
            orig = list(obj)
            val = obj[i] // 2 + 1
            setitem_usecase(obj, i, val)
            self.assertEqual(obj[i], val)
            for j, val in enumerate(orig):
                if j != i:
                    self.assertEqual(obj[j], val)

    def check_getslice(self, obj):
        self._check_unary(getslice_usecase, obj, 1, len(obj) - 1)

    def test_len(self):
        self.check_len(bytearray(5))
        self.check_len(b"xyz")
        for mem in self._memoryviews():
            self.check_len(mem)
        for arr in self._arrays():
            self.check_len(arr)
        for buf in self._readonlies():
            self.check_getitem(buf)

    def test_getitem(self):
        self.check_getitem(bytearray(b"abc"))
        self.check_getitem(b"xyz")
        for mem in self._memoryviews():
            self.check_getitem(mem)
        for arr in self._arrays():
            self.check_getitem(arr)
        for buf in self._readonlies():
            self.check_getitem(buf)

    def test_getslice(self):
        with self.assertTypingError():
            self.check_getslice(bytearray(b"abcde"))
        self.check_getslice(b"xyzuvw")
        self.check_getslice(memoryview(b"xyzuvw"))
        with self.assertTypingError():
            self.check_getslice(array.array('i', range(10)))
        for buf in self._readonlies():
            self.check_getitem(buf)

    def test_setitem(self):
        self.check_setitem(bytearray(b"abcdefghi"))
        for arr in self._arrays():
            self.check_setitem(arr)
        for mem in self._memoryviews():
            self.check_getitem(mem)
        # Read-only buffers
        for buf in self._readonlies():
            with self.assertTypingError():
                self.check_setitem(buf)

    def test_iter(self):
        self.check_iter(bytearray(b"abc"))
        self.check_iter(b"xyz")
        self.check_iter(memoryview(b"xyz"))
        for arr in self._arrays():
            self.check_iter(arr)
        for buf in self._readonlies():
            self.check_getitem(buf)


class TestMemoryView(MemoryLeakMixin, TestCase):
    """
    Test memoryview-specific attributes and operations.
    """

    def _arrays(self):
        arr = np.arange(12)
        yield arr
        arr = arr.reshape((3, 4))
        yield arr
        yield arr.T
        yield arr[::2]
        arr.setflags(write=False)
        yield arr
        arr = np.zeros(())
        assert arr.ndim == 0
        yield arr

    def test_ndim(self):
        for arr in self._arrays():
            m = memoryview(arr)
            self.assertPreciseEqual(ndim_usecase(m), arr.ndim)

    def test_shape(self):
        for arr in self._arrays():
            m = memoryview(arr)
            self.assertPreciseEqual(shape_usecase(m), arr.shape)

    def test_strides(self):
        for arr in self._arrays():
            m = memoryview(arr)
            self.assertPreciseEqual(strides_usecase(m), arr.strides)

    def test_itemsize(self):
        for arr in self._arrays():
            m = memoryview(arr)
            self.assertPreciseEqual(itemsize_usecase(m), arr.itemsize)

    def test_nbytes(self):
        for arr in self._arrays():
            m = memoryview(arr)
            self.assertPreciseEqual(nbytes_usecase(m), arr.size * arr.itemsize)

    def test_readonly(self):
        for arr in self._arrays():
            m = memoryview(arr)
            self.assertIs(readonly_usecase(m), not arr.flags.writeable)
        m = memoryview(b"xyz")
        self.assertIs(readonly_usecase(m), True)
        m = memoryview(bytearray(b"xyz"))
        self.assertIs(readonly_usecase(m), False)

    def test_contiguous(self):
        m = memoryview(bytearray(b"xyz"))
        self.assertIs(contiguous_usecase(m), True)
        self.assertIs(c_contiguous_usecase(m), True)
        self.assertIs(f_contiguous_usecase(m), True)
        for arr in self._arrays():
            m = memoryview(arr)
            # Note `arr.flags.contiguous` is wrong (it mimics c_contiguous)
            self.assertIs(contiguous_usecase(m),
                          arr.flags.f_contiguous or arr.flags.c_contiguous)
            self.assertIs(c_contiguous_usecase(m), arr.flags.c_contiguous)
            self.assertIs(f_contiguous_usecase(m), arr.flags.f_contiguous)


if __name__ == '__main__':
    unittest.main()
