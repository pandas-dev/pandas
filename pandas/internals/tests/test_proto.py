# -*- coding: utf-8 -*-

import gc
import sys

import numpy as np

import pandas.util.testing as tm
import pandas.native as lib


class TestLibPandas(tm.TestCase):

    def test_convert_numpy_int_dtypes(self):
        cases = [
            ('i1', lib.INT8, 'int8'),
            ('i2', lib.INT16, 'int16'),
            ('i4', lib.INT32, 'int32'),
            ('i8', lib.INT64, 'int64'),
            ('u1', lib.UINT8, 'uint8'),
            ('u2', lib.UINT16, 'uint16'),
            ('u4', lib.UINT32, 'uint32'),
            ('u8', lib.UINT64, 'uint64'),
            ('f4', lib.FLOAT, 'float'),
            ('f8', lib.DOUBLE, 'double'),
        ]

        for np_type, pd_type, name in cases:
            dtype = np.dtype(np_type)

            result = lib.convert_numpy_dtype(dtype)
            expected = lib.primitive_type(pd_type)

            assert result.equals(expected)
            assert name in repr(result)

    def test_libpandas_decrefs(self):
        arr = np.array([1, 2, 3, 4, 5])
        before = sys.getrefcount(arr)

        result = lib.to_array(arr)
        assert sys.getrefcount(arr) == (before + 1)

        result = None  # noqa
        gc.collect()
        assert sys.getrefcount(arr) == before

    def test_convert_integer_arrays(self):
        arr = np.array([1, 2, 3, 4, 5])
        result = lib.to_array(arr)

        ex_type = lib.primitive_type(lib.INT8)

        assert len(result) == 5
        # assert result.dtype.equals(ex_type)
