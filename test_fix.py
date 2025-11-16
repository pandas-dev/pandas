#!/usr/bin/env python3
"""
Test to verify the fix for the NDArrayBacked.__setstate__ issue with datetime64[ns] in MultiIndex.
"""
import numpy as np
import pandas as pd
import pickle
import sys
from pandas._libs.arrays import NDArrayBacked


def test_ndarraybacked_setstate_with_tuple():
    """Test that NDArrayBacked.__setstate__ handles 3-element tuple with (data, dtype, (dtype, array)) format."""
    # Create a mock NDArrayBacked instance to test the __setstate__ method directly
    class MockNDArrayBacked(NDArrayBacked):
        def __init__(self):
            # We'll manually set _ndarray and _dtype when needed
            pass
    
    # Test the problematic state format: a 3-element tuple where the third element 
    # is another (dtype, array) tuple rather than an attributes dict
    arr = np.array(['2026-04-05T16:07:45.133961216'], dtype='datetime64[ns]')
    dtype = arr.dtype
    problematic_state = (arr, dtype, (dtype, arr))
    
    print(f"Testing state: {problematic_state}")
    
    try:
        obj = MockNDArrayBacked()
        obj.__setstate__(problematic_state)
        print("SUCCESS: No NotImplementedError raised!")
        return True
    except NotImplementedError as e:
        print(f"FAILED: NotImplementedError still raised: {e}")
        return False
    except Exception as e:
        print(f"FAILED: Other exception raised: {e}")
        return False


def test_original_scenario():
    """Test the original scenario that triggered the issue."""
    try:
        print("Creating DataFrame with MultiIndex containing datetime64[ns]...")
        df = pd.DataFrame({
            'date': [np.datetime64("20110101", 'ns')], 
            'id': [1], 
            'val': [1]
        }).set_index(["date", "id"])
        
        print("DataFrame created successfully:", df.index)
        
        # Try to pickle and unpickle it (this mimics what joblib does internally)
        print("Testing pickle/unpickle...")
        pickled = pickle.dumps(df)
        unpickled_df = pickle.loads(pickled)
        print("Pickle/unpickle successful:", unpickled_df.index)
        return True
    except Exception as e:
        print(f"FAILED: Error in original scenario: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("Testing NDArrayBacked __setstate__ fix...")
    
    test1_result = test_ndarraybacked_setstate_with_tuple()
    print()
    test2_result = test_original_scenario()
    
    if test1_result and test2_result:
        print("\nAll tests passed! The fix should work.")
        sys.exit(0)
    else:
        print("\nSome tests failed. The fix needs more work.")
        sys.exit(1)