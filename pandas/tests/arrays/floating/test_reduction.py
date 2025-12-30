

def test_skipna_with_nan_values():
    # GH#59965 - skipna=True should skip both masked values and NaN in data
    import pandas as pd
    import numpy as np

    # Create FloatingArrays with np.NaN in the data
    s1 = pd.Series({"a": 0.0, "b": 1, "c": 1, "d": 0})
    s2 = pd.Series({"a": 0.0, "b": 2, "c": 2, "d": 2})
    
    # Division creates NaN in the data (0.0 / 0.0 = NaN)
    s4 = s1.convert_dtypes() / s2.convert_dtypes()
    
    # mean with skipna=True should skip the NaN value
    result = s4.mean(skipna=True)
    expected = 0.3333333333333333  # mean of [0.5, 0.5, 0]
    assert abs(result - expected) < 1e-10, f"Expected {expected}, got {result}"
    
    # Also test with explicit NA
    s5 = pd.Series([None, 0.5, 0.5, 0]).convert_dtypes()
    result = s5.mean(skipna=True)
    assert abs(result - expected) < 1e-10, f"Expected {expected}, got {result}"
    
    # Test sum as well
    result_sum = s4.sum(skipna=True)
    expected_sum = 1.0  # sum of [0.5, 0.5, 0]
    assert abs(result_sum - expected_sum) < 1e-10, f"Expected {expected_sum}, got {result_sum}"
