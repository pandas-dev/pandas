import numpy as np

import pandas as pd
from pandas import Series


def test_to_numeric_nan_mask_sync():
    # GH 63732
    ser = Series([1.0, np.nan], dtype=np.float64)
    result = pd.to_numeric(ser, dtype_backend="numpy_nullable")

    # Verify the mask is correctly synced
    assert result.isna().sum() == 1
