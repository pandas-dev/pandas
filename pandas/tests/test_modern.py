import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

# Ensure accessor is registered
import pandas.core.modern


class TestModernAccessor:
    def test_optimize(self):
        df = pd.DataFrame({
            "A": ["a", "b", "c", "a"] * 100,
            "B": [1, 2, 3, 4] * 100,
            "C": [1.1, 2.2, 3.3, 4.4] * 100
        })
        # Force types to expensive ones
        df["A"] = df["A"].astype(object)
        df["B"] = df["B"].astype("int64")
        df["C"] = df["C"].astype("float64")

        opt_df = df.modern.optimize(verbose=False)

        assert opt_df["A"].dtype == "category"
        assert str(opt_df["B"].dtype) in ["int8", "int16"]
        assert str(opt_df["C"].dtype) == "float32"

    def test_clean_smart(self):
        df = pd.DataFrame({
            "id": [" 1 ", "2", "3"],
            "val": ["na", "N/A", "100"]
        })
        clean_df = df.modern.clean_smart()
        
        # Check standard NaNs
        assert pd.isna(clean_df["val"].iloc[0])
        assert pd.isna(clean_df["val"].iloc[1])
        # Check stripping
        assert clean_df["id"].iloc[0] == "1"

    def test_expect_validate(self):
        df = pd.DataFrame({"id": [1, 2, 2], "val": [10, None, 30]})

        # Test failure (Unique)
        with pytest.raises(ValueError, match="is not unique"):
            df.modern.expect.unique("id").validate()

        # Test failure (No Nulls)
        with pytest.raises(ValueError, match="Nulls found"):
            df.modern.expect.no_nulls("val").validate()
        
        # Test Success
        df_good = pd.DataFrame({"id": [1, 2], "val": [10, 20]})
        res = df_good.modern.expect.unique("id").no_nulls().validate()
        assert res is df_good

    def test_ai_mock(self):
        df = pd.DataFrame({"a": [1, 2]})
        # Just check it runs and prints (simulated) without error
        res = df.modern.ai.ask("Plot this")
        assert res == "AI_Response_Object"

