import numpy as np
import pytest
import pandas as pd
from pandas import Series

@pytest.mark.parametrize("backend", ["numpy_nullable", "pyarrow"])
def test_to_numeric_nan_mask_sync(backend):
    # GH 63732
    ser = Series([1.0, np.nan], dtype=np.float64)
    result = pd.to_numeric(ser, dtype_backend=backend)
    
    # Verify the mask is correctly synced
    assert result.isna().sum() == 1
    assert result.dropna().size == 1
    assert result.iloc[1] is pd.NA or pd.isna(result.iloc[1])
