import numpy as np

import unittest
from numba import jit, from_dtype
from numba.core import types
from numba.typed import Dict
from numba.tests.support import (TestCase, skip_ppc64le_issue4563)


def getitem(x, i):
    return x[i]


def getitem2(x, i, j):
    return x[i][j]


def setitem(x, i, v):
    x[i] = v
    return x


def setitem2(x, i, y, j):
    x[i] = y[j]
    return x


def setitem_literal(x, i):
    x[i] = '123'
    return x


def getitem_key(x, y, j):
    x[y[j]] = 123


def return_len(x, i):
    return len(x[i])


def return_bool(x, i):
    return bool(x[i])


def equal_getitem(x, i, j):
    return x[i] == x[j]


def notequal_getitem(x, i, j):
    return x[i] != x[j]


def lessthan_getitem(x, i, j):
    return x[i] < x[j]


def greaterthan_getitem(x, i, j):
    return x[i] > x[j]


def lessequal_getitem(x, i, j):
    return x[i] <= x[j]


def greaterequal_getitem(x, i, j):
    return x[i] >= x[j]


def contains_getitem2(x, i, y, j):
    return x[i] in y[j]


def equal_getitem_value(x, i, v):
    r1 = x[i] == v
    r2 = v == x[i]
    if r1 == r2:
        return r1
    raise ValueError('x[i] == v and v == x[i] are unequal')


def notequal_getitem_value(x, i, v):
    r1 = x[i] != v
    r2 = v != x[i]
    if r1 == r2:
        return r1
    raise ValueError('x[i] != v and v != x[i] are unequal')


def return_isascii(x, i):
    return x[i].isascii()


def return_isupper(x, i):
    return x[i].isupper()


def return_upper(x, i):
    return x[i].upper()


def return_str(x, i):
    return str(x[i])


def return_bytes(x, i):
    return bytes(x[i])


def return_hash(x, i):
    return hash(x[i])


def return_find(x, i, y, j):
    return x[i].find(y[j])


def return_rfind(x, i, y, j):
    return x[i].rfind(y[j])


def return_startswith(x, i, y, j):
    return x[i].startswith(y[j])


def return_endswith(x, i, y, j):
    return x[i].endswith(y[j])


def return_split1(x, i):
    return x[i].split()


def return_split2(x, i, y, j):
    return x[i].split(y[j])


def return_split3(x, i, y, j, maxsplit):
    return x[i].split(sep=y[j], maxsplit=maxsplit)


# NOT IMPLEMENTED: tests for rsplit


def return_center1(x, i, w):
    return x[i].center(w)


def return_center2(x, i, w, y, j):
    return x[i].center(w, y[j])


def return_ljust1(x, i, w):
    return x[i].ljust(w)


def return_ljust2(x, i, w, y, j):
    return x[i].ljust(w, y[j])


def return_rjust1(x, i, w):
    return x[i].rjust(w)


def return_rjust2(x, i, w, y, j):
    return x[i].rjust(w, y[j])


def return_join(x, i, y, j, z, k):
    return x[i].join([y[j], z[k]])


def return_zfill(x, i, w):
    return x[i].zfill(w)


def return_lstrip1(x, i):
    return x[i].lstrip()


def return_lstrip2(x, i, y, j):
    return x[i].lstrip(y[j])


def return_rstrip1(x, i):
    return x[i].rstrip()


def return_rstrip2(x, i, y, j):
    return x[i].rstrip(y[j])


def return_strip1(x, i):
    return x[i].strip()


def return_strip2(x, i, y, j):
    return x[i].strip(y[j])


def return_add(x, i, y, j):
    return x[i] + y[j]


def return_iadd(x, i, y, j):
    x[i] += y[j]
    return x[i]


def return_mul(x, i, y, j):
    return x[i] * y[j]  # one of the operants must be integer


def return_not(x, i):
    return not x[i]


def join_string_array(str_arr):
    return ",".join(str_arr)


@skip_ppc64le_issue4563
class TestUnicodeArray(TestCase):

    def _test(self, pyfunc, cfunc, *args, **kwargs):
        expected = pyfunc(*args, **kwargs)
        self.assertPreciseEqual(cfunc(*args, **kwargs), expected)

    def test_getitem2(self):
        cgetitem2 = jit(nopython=True)(getitem2)

        arr = np.array(b'12')
        self.assertPreciseEqual(cgetitem2(arr, (), 0), getitem2(arr, (), 0))
        with self.assertRaisesRegex(IndexError, 'index out of range'):
            cgetitem2(arr, (), 2)

        arr = np.array('12')
        self.assertPreciseEqual(cgetitem2(arr, (), 0), getitem2(arr, (), 0))
        with self.assertRaisesRegex(IndexError, 'index out of range'):
            cgetitem2(arr, (), 2)

        arr = np.array([b'12', b'3'])
        self.assertPreciseEqual(cgetitem2(arr, 0, 0), getitem2(arr, 0, 0))
        self.assertPreciseEqual(cgetitem2(arr, 0, 1), getitem2(arr, 0, 1))
        self.assertPreciseEqual(cgetitem2(arr, 1, 0), getitem2(arr, 1, 0))
        with self.assertRaisesRegex(IndexError, 'index out of range'):
            cgetitem2(arr, 1, 1)

        arr = np.array(['12', '3'])
        self.assertPreciseEqual(cgetitem2(arr, 0, 0), getitem2(arr, 0, 0))
        self.assertPreciseEqual(cgetitem2(arr, 0, 1), getitem2(arr, 0, 1))
        self.assertPreciseEqual(cgetitem2(arr, 1, 0), getitem2(arr, 1, 0))
        with self.assertRaisesRegex(IndexError, 'index out of range'):
            cgetitem2(arr, 1, 1)

    def test_getitem(self):
        pyfunc = getitem
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, b'12', 1)
        self._test(pyfunc, cfunc, np.array(b'12'), ())
        self._test(pyfunc, cfunc, np.array([b'12', b'3']), 0)
        self._test(pyfunc, cfunc, np.array([b'12', b'3']), 1)

        self._test(pyfunc, cfunc, '12', 1)
        self._test(pyfunc, cfunc, np.array('12'), ())
        self._test(pyfunc, cfunc, np.array(['12', '3']), 0)
        self._test(pyfunc, cfunc, np.array(['12', '3']), 1)

    def test_getitem_key(self):
        pyfunc = getitem_key
        cfunc = jit(nopython=True)(pyfunc)

        for x, i in [
                (np.array('123'), ()),
                (np.array(['123']), 0),
                (np.array(b'123'), ()),
                (np.array([b'123']), 0)
        ]:
            d1 = {}
            d2 = Dict.empty(from_dtype(x.dtype), types.int64)
            pyfunc(d1, x, i)
            cfunc(d2, x, i)
            self.assertEqual(d1, d2)
            # check for charseq to str conversion:
            str(d2)

    def test_setitem(self):
        pyfunc = setitem
        cfunc = jit(nopython=True)(pyfunc)

        x = np.array(12)
        self._test(pyfunc, cfunc, x, (), 34)

        x1 = np.array(b'123')
        x2 = np.array(b'123')
        y1 = pyfunc(x1, (), b'34')
        y2 = cfunc(x2, (), b'34')
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        x1 = np.array(['123'])
        x2 = np.array(['123'])
        y1 = pyfunc(x1, 0, '34')
        y2 = cfunc(x2, 0, '34')
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

    def test_setitem2(self):
        pyfunc = setitem2
        cfunc = jit(nopython=True)(pyfunc)

        x1 = np.array(['123', 'ABC'])
        x2 = np.array(['123', 'ABC'])
        y1 = pyfunc(x1, 0, x1, 1)
        y2 = cfunc(x2, 0, x2, 1)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        x1 = np.array([b'123', b'ABC'])
        x2 = np.array([b'123', b'ABC'])
        y1 = pyfunc(x1, 0, x1, 1)
        y2 = cfunc(x2, 0, x2, 1)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        x1 = np.array('123')
        x2 = np.array('123')
        z1 = np.array('ABC')
        z2 = np.array('ABC')
        y1 = pyfunc(x1, (), z1, ())
        y2 = cfunc(x2, (), z2, ())
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        x1 = np.array(123)
        x2 = np.array(123)
        z1 = 456,
        z2 = 456,
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        # bytes
        x1 = np.array(b'123')
        x2 = np.array(b'123')
        z1 = b'ABC',
        z2 = b'ABC',
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        # UTF-8
        x1 = np.array('123')
        x2 = np.array('123')
        z1 = 'ABC',
        z2 = 'ABC',
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        # UTF-16
        x1 = np.array('123')
        x2 = np.array('123')
        z1 = 'AB\u01e9',
        z2 = 'AB\u01e9',
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        # UTF-32
        x1 = np.array('123')
        x2 = np.array('123')
        z1 = 'AB\U00108a0e',
        z2 = 'AB\U00108a0e',
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        # UTF-8, assign longer value (truncates as in numpy)
        x1 = np.array('123')
        x2 = np.array('123')
        z1 = 'ABCD',
        z2 = 'ABCD',
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        # UTF-8, assign shorter value
        x1 = np.array('123')
        x2 = np.array('123')
        z1 = 'AB',
        z2 = 'AB',
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        # bytes, assign longer value (truncates as in numpy)
        x1 = np.array(b'123')
        x2 = np.array(b'123')
        z1 = b'ABCD',
        z2 = b'ABCD',
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        # bytes, assign shorter value
        x1 = np.array(b'123')
        x2 = np.array(b'123')
        z1 = b'AB',
        z2 = b'AB',
        y1 = pyfunc(x1, (), z1, 0)
        y2 = cfunc(x2, (), z2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

    def test_setitem_literal(self):
        pyfunc = setitem_literal
        cfunc = jit(nopython=True)(pyfunc)

        x1 = np.array('ABC')
        x2 = np.array('ABC')
        y1 = pyfunc(x1, ())
        y2 = cfunc(x2, ())
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        x1 = np.array(['ABC', '5678'])
        x2 = np.array(['ABC', '5678'])
        y1 = pyfunc(x1, 0)
        y2 = cfunc(x2, 0)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

        x1 = np.array(['ABC', '5678'])
        x2 = np.array(['ABC', '5678'])
        y1 = pyfunc(x1, 1)
        y2 = cfunc(x2, 1)
        self.assertPreciseEqual(x1, x2)
        self.assertPreciseEqual(y1, y2)

    def test_return_len(self):
        pyfunc = return_len
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array(''), ())
        self._test(pyfunc, cfunc, np.array(b''), ())
        self._test(pyfunc, cfunc, np.array(b'12'), ())
        self._test(pyfunc, cfunc, np.array('12'), ())
        self._test(pyfunc, cfunc, np.array([b'12', b'3']), 0)
        self._test(pyfunc, cfunc, np.array(['12', '3']), 0)
        self._test(pyfunc, cfunc, np.array([b'12', b'3']), 1)
        self._test(pyfunc, cfunc, np.array(['12', '3']), 1)

    def test_return_bool(self):
        pyfunc = return_bool
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array(''), ())
        self._test(pyfunc, cfunc, np.array(b''), ())
        self._test(pyfunc, cfunc, np.array(b'12'), ())
        self._test(pyfunc, cfunc, np.array('12'), ())
        self._test(pyfunc, cfunc, np.array([b'12', b'']), 0)
        self._test(pyfunc, cfunc, np.array(['12', '']), 0)
        self._test(pyfunc, cfunc, np.array([b'12', b'']), 1)
        self._test(pyfunc, cfunc, np.array(['12', '']), 1)

    def _test_op_getitem(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array([1, 2]), 0, 1)
        self._test(pyfunc, cfunc, '12', 0, 1)
        self._test(pyfunc, cfunc, b'12', 0, 1)
        self._test(pyfunc, cfunc, np.array(b'12'), (), ())
        self._test(pyfunc, cfunc, np.array('1234'), (), ())

        self._test(pyfunc, cfunc, np.array([b'1', b'2']), 0, 0)
        self._test(pyfunc, cfunc, np.array([b'1', b'2']), 0, 1)
        self._test(pyfunc, cfunc, np.array([b'12', b'3']), 0, 0)
        self._test(pyfunc, cfunc, np.array([b'12', b'3']), 1, 1)
        self._test(pyfunc, cfunc, np.array([b'12', b'3']), 0, 1)
        self._test(pyfunc, cfunc, np.array([b'12', b'3']), 1, 0)

        self._test(pyfunc, cfunc, np.array(['1', '2']), 0, 0)
        self._test(pyfunc, cfunc, np.array(['1', '2']), 0, 1)
        self._test(pyfunc, cfunc, np.array(['12', '3']), 0, 0)
        self._test(pyfunc, cfunc, np.array(['12', '3']), 1, 1)
        self._test(pyfunc, cfunc, np.array(['12', '3']), 0, 1)
        self._test(pyfunc, cfunc, np.array(['12', '3']), 1, 0)

    def test_equal_getitem(self):
        self._test_op_getitem(equal_getitem)

    def test_notequal_getitem(self):
        self._test_op_getitem(notequal_getitem)

    def test_lessthan_getitem(self):
        self._test_op_getitem(lessthan_getitem)

    def test_greaterthan_getitem(self):
        self._test_op_getitem(greaterthan_getitem)

    def test_lessequal_getitem(self):
        self._test_op_getitem(lessequal_getitem)

    def test_greaterequal_getitem(self):
        self._test_op_getitem(greaterequal_getitem)

    def _test_op_getitem_value(self, pyfunc):
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array([1, 2]), 0, 1)
        self._test(pyfunc, cfunc, '12', 0, '1')
        self._test(pyfunc, cfunc, '12', 1, '3')
        self._test(pyfunc, cfunc, np.array('1234'), (), '1234')
        self._test(pyfunc, cfunc, np.array(['1234']), 0, '1234')
        self._test(pyfunc, cfunc, np.array(['1234']), 0, 'abc')
        # fails: No conversion from array(bool, 1d, C) to bool
        #self._test(pyfunc, cfunc, b'12', 0, b'1')
        self._test(pyfunc, cfunc, np.array(b'12'), (), b'12')
        self._test(pyfunc, cfunc, np.array([b'12']), 0, b'12')
        self._test(pyfunc, cfunc, np.array([b'12']), 0, b'a')

    def test_equal_getitem_value(self):
        self._test_op_getitem_value(equal_getitem_value)

    def test_notequal_getitem_value(self):
        self._test_op_getitem_value(notequal_getitem_value)

    def test_contains_getitem2(self):
        pyfunc = contains_getitem2
        cfunc = jit(nopython=True)(pyfunc)

        x = np.array('123')
        y = np.array('12345')
        self._test(pyfunc, cfunc, x, (), y, ())
        self._test(pyfunc, cfunc, y, (), x, ())

        x = np.array(b'123')
        y = np.array(b'12345')
        self._test(pyfunc, cfunc, x, (), y, ())
        self._test(pyfunc, cfunc, y, (), x, ())

        x = ('123',)
        y = np.array('12345')
        self._test(pyfunc, cfunc, x, 0, y, ())
        self._test(pyfunc, cfunc, y, (), x, 0)

        x = (b'123',)
        y = np.array(b'12345')
        self._test(pyfunc, cfunc, x, 0, y, ())
        self._test(pyfunc, cfunc, y, (), x, 0)

    def test_return_isascii(self):
        pyfunc = return_isascii
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array('1234'), ())
        self._test(pyfunc, cfunc, np.array(['1234']), 0)
        self._test(pyfunc, cfunc, np.array('1234\u00e9'), ())
        self._test(pyfunc, cfunc, np.array(['1234\u00e9']), 0)

    def test_return_isupper(self):
        pyfunc = return_isupper
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array('abc'), ())
        self._test(pyfunc, cfunc, np.array(['abc']), 0)

        self._test(pyfunc, cfunc, np.array(b'abc'), ())
        self._test(pyfunc, cfunc, np.array([b'abc']), 0)

    def test_return_str(self):
        pyfunc = return_str
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array('1234'), ())
        self._test(pyfunc, cfunc, np.array(['1234']), 0)

    def test_return_bytes(self):
        pyfunc = return_bytes
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array(b'1234'), ())
        self._test(pyfunc, cfunc, np.array([b'1234']), 0)

    def test_return_upper(self):
        pyfunc = return_upper
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array('abc'), ())
        self._test(pyfunc, cfunc, np.array(['abc']), 0)

        self._test(pyfunc, cfunc, np.array(b'abc'), ())
        self._test(pyfunc, cfunc, np.array([b'abc']), 0)

    def test_hash(self):
        pyfunc = return_hash
        cfunc = jit(nopython=True)(pyfunc)

        hash1 = pyfunc(np.array('123'), ())
        hash2 = hash('123')
        hash3 = hash(np.array('123')[()])
        self.assertTrue(hash1 == hash2 == hash3)

        self._test(pyfunc, cfunc, np.array('1234'), ())
        self._test(pyfunc, cfunc, np.array(['1234']), 0)

        self._test(pyfunc, cfunc, np.array('1234\u00e9'), ())
        self._test(pyfunc, cfunc, np.array(['1234u00e9']), 0)

        self._test(pyfunc, cfunc, np.array('1234\U00108a0e'), ())
        self._test(pyfunc, cfunc, np.array(['1234\U00108a0e']), 0)

        self._test(pyfunc, cfunc, np.array(b'1234'), ())
        self._test(pyfunc, cfunc, np.array([b'1234']), 0)

    def test_return_find(self):
        pyfunc = return_find
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array('1234'), (), np.array('23'), ())
        self._test(pyfunc, cfunc, np.array('1234'), (), ('23',), 0)
        self._test(pyfunc, cfunc, ('1234',), 0, np.array('23'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), (), np.array(b'23'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), (), (b'23',), 0)
        self._test(pyfunc, cfunc, (b'1234',), 0, np.array(b'23'), ())

    def test_return_rfind(self):
        pyfunc = return_rfind
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array('1234'), (), np.array('23'), ())
        self._test(pyfunc, cfunc, np.array('1234'), (), ('23',), 0)
        self._test(pyfunc, cfunc, ('1234',), 0, np.array('23'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), (), np.array(b'23'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), (), (b'23',), 0)
        self._test(pyfunc, cfunc, (b'1234',), 0, np.array(b'23'), ())

    def test_return_startswith(self):
        pyfunc = return_startswith
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1234'), (), np.array('23'), ())
        self._test(pyfunc, cfunc, np.array('1234'), (), ('23',), 0)
        self._test(pyfunc, cfunc, ('1234',), 0, np.array('23'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), (), np.array(b'23'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), (), (b'23',), 0)
        self._test(pyfunc, cfunc, (b'1234',), 0, np.array(b'23'), ())

    def test_return_endswith(self):
        pyfunc = return_endswith
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1234'), (), np.array('23'), ())
        self._test(pyfunc, cfunc, np.array('1234'), (), ('23',), 0)
        self._test(pyfunc, cfunc, ('1234',), 0, np.array('23'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), (), np.array(b'23'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), (), (b'23',), 0)
        self._test(pyfunc, cfunc, (b'1234',), 0, np.array(b'23'), ())

    def test_return_split1(self):
        pyfunc = return_split1
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('12 34'), ())
        self._test(pyfunc, cfunc, np.array(b'1234'), ())

    def test_return_split2(self):
        pyfunc = return_split2
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('12 34'), (), np.array(' '), ())
        self._test(pyfunc, cfunc, np.array('12 34'), (), (' ',), 0)
        self._test(pyfunc, cfunc, ('12 34',), 0, np.array(' '), ())
        self._test(pyfunc, cfunc, np.array(b'12 34'), (), np.array(b' '), ())
        self._test(pyfunc, cfunc, np.array(b'12 34'), (), (b' ',), 0)
        self._test(pyfunc, cfunc, (b'12 34',), 0, np.array(b' '), ())

    def test_return_split3(self):
        pyfunc = return_split3
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (),
                   np.array(' '), (), 2)
        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), (' ',), 0, 2)
        self._test(pyfunc, cfunc, ('1 2 3 4',), 0, np.array(' '), (), 2)
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (),
                   np.array(b' '), (), 2)
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), (b' ',), 0, 2)
        self._test(pyfunc, cfunc, (b'1 2 3 4',), 0, np.array(b' '), (), 2)

    def test_return_ljust1(self):
        pyfunc = return_ljust1
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40)
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40)

    def test_return_ljust2(self):
        pyfunc = return_ljust2
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40,
                   np.array('='), ())
        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40, ('=',), 0)
        self._test(pyfunc, cfunc, ('1 2 3 4',), 0, 40, np.array('='), ())
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40,
                   np.array(b'='), ())
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40, (b'=',), 0)
        self._test(pyfunc, cfunc, (b'1 2 3 4',), 0, 40, np.array(b'='), ())

    def test_return_rjust1(self):
        pyfunc = return_rjust1
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40)
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40)

    def test_return_rjust2(self):
        pyfunc = return_rjust2
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40,
                   np.array('='), ())
        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40, ('=',), 0)
        self._test(pyfunc, cfunc, ('1 2 3 4',), 0, 40, np.array('='), ())
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40,
                   np.array(b'='), ())
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40, (b'=',), 0)
        self._test(pyfunc, cfunc, (b'1 2 3 4',), 0, 40, np.array(b'='), ())

    def test_return_center1(self):
        pyfunc = return_center1
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40)
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40)

    def test_return_center2(self):
        pyfunc = return_center2
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40,
                   np.array('='), ())
        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40, ('=',), 0)
        self._test(pyfunc, cfunc, ('1 2 3 4',), 0, 40, np.array('='), ())
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40,
                   np.array(b'='), ())
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40, (b'=',), 0)
        self._test(pyfunc, cfunc, (b'1 2 3 4',), 0, 40, np.array(b'='), ())

    def test_return_join(self):
        pyfunc = return_join
        cfunc = jit(nopython=True)(pyfunc)
        self._test(pyfunc, cfunc, np.array(','), (), np.array('abc'), (),
                   np.array('123'), ())
        self._test(pyfunc, cfunc, np.array(','), (), np.array('abc'), (),
                   ('123',), 0)
        self._test(pyfunc, cfunc, (',',), 0, np.array('abc'), (),
                   np.array('123'), ())
        self._test(pyfunc, cfunc, (',',), 0, np.array('abc'), (),
                   ('123',), 0)
        self._test(pyfunc, cfunc, np.array(b','), (), np.array(b'abc'), (),
                   np.array(b'123'), ())
        self._test(pyfunc, cfunc, np.array(b','), (), np.array(b'abc'), (),
                   (b'123',), 0)
        self._test(pyfunc, cfunc, (b',',), 0, np.array(b'abc'), (),
                   np.array(b'123'), ())
        self._test(pyfunc, cfunc, (b',',), 0, np.array(b'abc'), (),
                   (b'123',), 0)

    def test_return_zfill(self):
        pyfunc = return_zfill
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('1 2 3 4'), (), 40)
        self._test(pyfunc, cfunc, np.array(b'1 2 3 4'), (), 40)

    def test_return_lstrip1(self):
        pyfunc = return_lstrip1
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('  123  '), ())
        self._test(pyfunc, cfunc, np.array(b'  123  '), ())

    def test_return_lstrip2(self):
        pyfunc = return_lstrip2
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('  123  '), (), np.array(' '), ())
        self._test(pyfunc, cfunc, np.array('  123  '), (), (' ',), 0)
        self._test(pyfunc, cfunc, ('  123  ',), 0, np.array(' '), ())

        self._test(pyfunc, cfunc, np.array(b'  123  '), (), np.array(b' '), ())
        self._test(pyfunc, cfunc, np.array(b'  123  '), (), (b' ',), 0)
        self._test(pyfunc, cfunc, (b'  123  ',), 0, np.array(b' '), ())

    def test_return_rstrip1(self):
        pyfunc = return_rstrip1
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('  123  '), ())
        self._test(pyfunc, cfunc, np.array(b'  123  '), ())

    def test_return_rstrip2(self):
        pyfunc = return_rstrip2
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('  123  '), (), np.array(' '), ())
        self._test(pyfunc, cfunc, np.array('  123  '), (), (' ',), 0)
        self._test(pyfunc, cfunc, ('  123  ',), 0, np.array(' '), ())

        self._test(pyfunc, cfunc, np.array(b'  123  '), (), np.array(b' '), ())
        self._test(pyfunc, cfunc, np.array(b'  123  '), (), (b' ',), 0)
        self._test(pyfunc, cfunc, (b'  123  ',), 0, np.array(b' '), ())

    def test_return_strip1(self):
        pyfunc = return_strip1
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('  123  '), ())
        self._test(pyfunc, cfunc, np.array(b'  123  '), ())

    def test_return_strip2(self):
        pyfunc = return_strip2
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('  123  '), (), np.array(' '), ())
        self._test(pyfunc, cfunc, np.array('  123  '), (), (' ',), 0)
        self._test(pyfunc, cfunc, ('  123  ',), 0, np.array(' '), ())

        self._test(pyfunc, cfunc, np.array(b'  123  '), (), np.array(b' '), ())
        self._test(pyfunc, cfunc, np.array(b'  123  '), (), (b' ',), 0)
        self._test(pyfunc, cfunc, (b'  123  ',), 0, np.array(b' '), ())

    def test_return_add(self):
        pyfunc = return_add
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('ab'), (), np.array('cd'), ())
        self._test(pyfunc, cfunc, np.array('ab'), (), ('cd',), 0)
        self._test(pyfunc, cfunc, ('ab',), 0, np.array('cd'), ())

        self._test(pyfunc, cfunc, np.array(b'ab'), (), np.array(b'cd'), ())
        self._test(pyfunc, cfunc, np.array(b'ab'), (), (b'cd',), 0)
        self._test(pyfunc, cfunc, (b'ab',), 0, np.array(b'cd'), ())

    def test_return_iadd(self):
        pyfunc = return_iadd
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('ab'), (), np.array('cd'), ())
        self._test(pyfunc, cfunc, np.array('ab'), (), ('cd',), 0)
        expected = pyfunc(['ab'], 0, np.array('cd'), ())
        result = pyfunc(['ab'], 0, np.array('cd'), ())
        self.assertPreciseEqual(result, expected)

        self._test(pyfunc, cfunc, np.array(b'ab'), (), np.array(b'cd'), ())
        self._test(pyfunc, cfunc, np.array(b'ab'), (), (b'cd',), 0)
        expected = pyfunc([b'ab'], 0, np.array(b'cd'), ())
        result = pyfunc([b'ab'], 0, np.array(b'cd'), ())
        self.assertPreciseEqual(result, expected)

    def test_return_mul(self):
        pyfunc = return_mul
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('ab'), (), (5,), 0)
        self._test(pyfunc, cfunc, (5,), 0, np.array('ab'), ())
        self._test(pyfunc, cfunc, np.array(b'ab'), (), (5,), 0)
        self._test(pyfunc, cfunc, (5,), 0, np.array(b'ab'), ())

    def test_return_not(self):
        pyfunc = return_not
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array('ab'), ())
        self._test(pyfunc, cfunc, np.array(b'ab'), ())
        self._test(pyfunc, cfunc, (b'ab',), 0)

        self._test(pyfunc, cfunc, np.array(''), ())
        self._test(pyfunc, cfunc, np.array(b''), ())
        self._test(pyfunc, cfunc, (b'',), 0)

    def test_join(self):
        pyfunc = join_string_array
        cfunc = jit(nopython=True)(pyfunc)

        self._test(pyfunc, cfunc, np.array(["hi", "there"]))


if __name__ == '__main__':
    unittest.main()
