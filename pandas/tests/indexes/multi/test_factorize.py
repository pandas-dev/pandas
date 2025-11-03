"""
Tests for MultiIndex.factorize method
"""

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestMultiIndexFactorize:
    def test_factorize_extension_dtype_int32(self):
        # GH#62337: factorize should preserve Int32 extension dtype
        df = pd.DataFrame({"col": pd.Series([1, None, 2], dtype="Int32")})
        mi = pd.MultiIndex.from_frame(df)

        codes, uniques = mi.factorize()

        result_dtype = uniques.to_frame().iloc[:, 0].dtype
        expected_dtype = pd.Int32Dtype()
        assert result_dtype == expected_dtype

        # Verify codes are correct
        expected_codes = np.array([0, 1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(codes, expected_codes)

    @pytest.mark.parametrize("dtype", ["Int32", "Int64", "string", "boolean"])
    def test_factorize_extension_dtypes(self, dtype):
        # GH#62337: factorize should preserve various extension dtypes
        if dtype == "boolean":
            values = [True, None, False]
        elif dtype == "string":
            values = ["a", None, "b"]
        else:  # Int32, Int64
            values = [1, None, 2]

        df = pd.DataFrame({"col": pd.Series(values, dtype=dtype)})
        mi = pd.MultiIndex.from_frame(df)

        codes, uniques = mi.factorize()
        result_dtype = uniques.to_frame().iloc[:, 0].dtype

        assert str(result_dtype) == dtype

    def test_factorize_multiple_extension_dtypes(self):
        # GH#62337: factorize with multiple columns having extension dtypes
        df = pd.DataFrame(
            {
                "int_col": pd.Series([1, 2, 1], dtype="Int64"),
                "str_col": pd.Series(["a", "b", "a"], dtype="string"),
            }
        )
        mi = pd.MultiIndex.from_frame(df)

        codes, uniques = mi.factorize()

        result_frame = uniques.to_frame()
        assert result_frame.iloc[:, 0].dtype == pd.Int64Dtype()
        assert result_frame.iloc[:, 1].dtype == pd.StringDtype()

        # Should have 2 unique combinations: (1,'a') and (2,'b')
        assert len(uniques) == 2

    def test_factorize_preserves_names(self):
        # GH#62337: factorize should preserve MultiIndex names
        df = pd.DataFrame(
            {
                "level_1": pd.Series([1, 2], dtype="Int32"),
                "level_2": pd.Series(["a", "b"], dtype="string"),
            }
        )
        mi = pd.MultiIndex.from_frame(df)

        codes, uniques = mi.factorize()

        tm.assert_index_equal(pd.Index(uniques.names), pd.Index(mi.names))

    def test_factorize_extension_dtype_with_sort(self):
        # GH#62337: factorize with sort=True should preserve extension dtypes
        df = pd.DataFrame({"col": pd.Series([2, None, 1], dtype="Int32")})
        mi = pd.MultiIndex.from_frame(df)

        codes, uniques = mi.factorize(sort=True)

        result_dtype = uniques.to_frame().iloc[:, 0].dtype
        assert result_dtype == pd.Int32Dtype()

    def test_factorize_empty_extension_dtype(self):
        # GH#62337: factorize on empty MultiIndex with extension dtype
        df = pd.DataFrame({"col": pd.Series([], dtype="Int32")})
        mi = pd.MultiIndex.from_frame(df)

        codes, uniques = mi.factorize()

        assert len(codes) == 0
        assert len(uniques) == 0
        assert uniques.to_frame().iloc[:, 0].dtype == pd.Int32Dtype()

    def test_factorize_regular_dtypes_unchanged(self):
        # Ensure regular dtypes still work as before
        df = pd.DataFrame({"int_col": [1, 2, 1], "float_col": [1.1, 2.2, 1.1]})
        mi = pd.MultiIndex.from_frame(df)

        codes, uniques = mi.factorize()

        result_frame = uniques.to_frame()
        assert result_frame.iloc[:, 0].dtype == np.dtype("int64")
        assert result_frame.iloc[:, 1].dtype == np.dtype("float64")

        # Should have 2 unique combinations
        assert len(uniques) == 2

    def test_factorize_mixed_extension_regular_dtypes(self):
        # Mix of extension and regular dtypes
        df = pd.DataFrame(
            {
                "ext_col": pd.Series([1, 2, 1], dtype="Int64"),
                "reg_col": [1.1, 2.2, 1.1],  # regular float64
            }
        )
        mi = pd.MultiIndex.from_frame(df)

        codes, uniques = mi.factorize()

        result_frame = uniques.to_frame()
        assert result_frame.iloc[:, 0].dtype == pd.Int64Dtype()
        assert result_frame.iloc[:, 1].dtype == np.dtype("float64")
