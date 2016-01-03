# -*- coding: utf-8 -*-

import numpy as np

import pandas.util.testing as tm
import pandas.native as lib


class TestLibPandas(tm.TestCase):

    def test_convert_numpy_int_dtypes(self):
        cases = [
            'i1', lib.INT8
        ]

        for np_type, pd_type in cases:
            dtype = np.dtype(np_type)

            result = lib.convert_numpy_dtype(dtype)
            expected = lib.primitive_type(pd_type)

            assert result.equals(expected)
