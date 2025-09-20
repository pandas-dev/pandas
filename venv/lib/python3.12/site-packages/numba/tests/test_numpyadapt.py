from ctypes import *

import numpy as np

import unittest
from numba import _helperlib


class ArrayStruct3D(Structure):
    # Mimic the structure defined in numba.targets.arrayobj's make_array()
    _fields_ = [
        ("meminfo", c_void_p),
        ("parent", c_void_p),
        ("nitems", c_ssize_t),
        ("itemsize", c_ssize_t),
        ("data", c_void_p),
        ("shape", (c_ssize_t * 3)),
        ("strides", (c_ssize_t * 3)),
    ]


class TestArrayAdaptor(unittest.TestCase):
    def test_array_adaptor(self):
        arystruct = ArrayStruct3D()

        adaptorptr = _helperlib.c_helpers['adapt_ndarray']
        adaptor = PYFUNCTYPE(c_int, py_object, c_void_p)(adaptorptr)

        ary = np.arange(60).reshape(2, 3, 10)
        status = adaptor(ary, byref(arystruct))
        self.assertEqual(status, 0)
        self.assertEqual(arystruct.data, ary.ctypes.data)
        self.assertNotEqual(arystruct.meminfo, 0)
        self.assertEqual(arystruct.parent, id(ary))
        self.assertEqual(arystruct.nitems, 60)
        self.assertEqual(arystruct.itemsize, ary.itemsize)
        for i in range(3):
            self.assertEqual(arystruct.shape[i], ary.ctypes.shape[i])
            self.assertEqual(arystruct.strides[i], ary.ctypes.strides[i])


if __name__ == '__main__':
    unittest.main()
