import re

import numpy as np
import pytest

from pandas.core.dtypes.common import pandas_dtype

from pandas import Float64Index, Index, Int64Index
import pandas._testing as tm


class TestAstype:
    def test_astype(self):

        float_index = Float64Index([0.0, 2.5, 5.0, 7.5, 10.0])
        result = float_index.astype(object)
        assert result.equals(float_index)
        assert float_index.equals(result)
        assert isinstance(result, Index) and not isinstance(result, Float64Index)

        # mixed int-float
        i = Float64Index([1.5, 2, 3, 4, 5])
        i.name = "foo"
        result = i.astype(object)
        assert result.equals(i)
        assert i.equals(result)
        assert isinstance(result, Index) and not isinstance(result, Float64Index)

        # GH#12881
        # a float astype int
        for dtype in ["int16", "int32", "int64"]:
            i = Float64Index([0, 1, 2])
            result = i.astype(dtype)
            expected = Int64Index([0, 1, 2])
            tm.assert_index_equal(result, expected)

            i = Float64Index([0, 1.1, 2])
            result = i.astype(dtype)
            expected = Int64Index([0, 1, 2])
            tm.assert_index_equal(result, expected)

        for dtype in ["float32", "float64"]:
            i = Float64Index([0, 1, 2])
            result = i.astype(dtype)
            expected = i
            tm.assert_index_equal(result, expected)

            i = Float64Index([0, 1.1, 2])
            result = i.astype(dtype)
            expected = Index(i.values.astype(dtype))
            tm.assert_index_equal(result, expected)

        # invalid
        for dtype in ["M8[ns]", "m8[ns]"]:
            msg = (
                f"Cannot convert Float64Index to dtype {pandas_dtype(dtype)}; "
                f"integer values are required for conversion"
            )
            with pytest.raises(TypeError, match=re.escape(msg)):
                i.astype(dtype)

        # GH#13149
        for dtype in ["int16", "int32", "int64"]:
            i = Float64Index([0, 1.1, np.NAN])
            msg = r"Cannot convert non-finite values \(NA or inf\) to integer"
            with pytest.raises(ValueError, match=msg):
                i.astype(dtype)

    def test_cannot_cast_inf_to_int(self):
        idx = Float64Index([1, 2, np.inf])

        msg = r"Cannot convert non-finite values \(NA or inf\) to integer"
        with pytest.raises(ValueError, match=msg):
            idx.astype(int)

    def test_astype_from_object(self):
        index = Index([1.0, np.nan, 0.2], dtype="object")
        result = index.astype(float)
        expected = Float64Index([1.0, np.nan, 0.2])
        assert result.dtype == expected.dtype
        tm.assert_index_equal(result, expected)
