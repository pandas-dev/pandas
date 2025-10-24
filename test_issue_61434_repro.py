"""
Pandas Issue #61434 - Reproduction Test

Issue: When attempting to merge a pandas DataFrame with a polars DataFrame,
the error message is unhelpful. 

Current behavior: Generic error about missing attributes or type errors
Expected behavior: Clear message saying "other must be pandas.DataFrame, 
                   received: polars.DataFrame"

Snippet from issue #61434:
https://github.com/pandas-dev/pandas/issues/61434
"""

import pandas as pd

# Try to import polars for testing
try:
    import polars as pl
    POLARS_AVAILABLE = True
except ImportError:
    POLARS_AVAILABLE = False
    print("Warning: polars not installed. Install with: pip install polars")


def test_merge_with_polars():
    """
    Reproduce the issue: Merging pandas DataFrame with polars DataFrame.
    
    Before fix: Generic/confusing error message
    After fix: Clear message about type mismatch
    """
    if not POLARS_AVAILABLE:
        print("Skipping test - polars not available")
        return False
    
    print("=" * 70)
    print("Test: Merging pandas DataFrame with polars DataFrame")
    print("=" * 70)
    
    # Create pandas DataFrame
    pdf = pd.DataFrame({
        'key': ['a', 'b', 'c'],
        'value_x': [1, 2, 3]
    })
    
    # Create polars DataFrame
    plf = pl.DataFrame({
        'key': ['a', 'b', 'c'],
        'value_y': [10, 20, 30]
    })
    
    print(f"\nPandas DataFrame type: {type(pdf)}")
    print(f"Polars DataFrame type: {type(plf)}")
    print("\nAttempting merge...")
    
    try:
        result = pd.merge(pdf, plf, on='key')
        print(f"✗ Unexpected: merge succeeded with result type {type(result)}")
        return False
    except TypeError as e:
        error_msg = str(e)
        print(f"\nError caught: {type(e).__name__}")
        print(f"Error message: {error_msg}")
        
        # Check if error message is helpful
        if "polars" in error_msg.lower() and "pandas" in error_msg.lower():
            print("\n✓ GOOD: Error message mentions both polars and pandas")
            print("✓ GOOD: User knows what went wrong")
            return True
        elif "must be" in error_msg.lower() or "expected" in error_msg.lower():
            print("\n✓ GOOD: Error message explains what's expected")
            return True
        else:
            print(f"\n✗ BAD: Error message is not helpful enough")
            print(f"   Expected something like:")
            print(f"   'other must be pandas.DataFrame, received: polars.DataFrame'")
            print(f"   But got: {error_msg}")
            return False
    except Exception as e:
        print(f"\n✗ Unexpected error type: {type(e).__name__}")
        print(f"   {e}")
        return False


def test_merge_pandas_baseline():
    """
    Baseline test: merge two pandas DataFrames should work.
    """
    print("\n" + "=" * 70)
    print("Test: Merging two pandas DataFrames (baseline)")
    print("=" * 70)
    
    df1 = pd.DataFrame({
        'key': ['a', 'b', 'c'],
        'value_x': [1, 2, 3]
    })
    
    df2 = pd.DataFrame({
        'key': ['a', 'b', 'c'],
        'value_y': [10, 20, 30]
    })
    
    try:
        result = pd.merge(df1, df2, on='key')
        print(f"✓ Merge succeeded")
        print(f"  Result shape: {result.shape}")
        print(f"  Result columns: {list(result.columns)}")
        return True
    except Exception as e:
        print(f"✗ Baseline test failed: {e}")
        return False


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("PANDAS ISSUE #61434 - REPRODUCTION TEST")
    print("=" * 70)
    print()
    
    baseline_ok = test_merge_pandas_baseline()
    polars_test_ok = test_merge_with_polars()
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Baseline (pandas merge): {'✓ PASS' if baseline_ok else '✗ FAIL'}")
    print(f"Polars test (error msg): {'✓ GOOD' if polars_test_ok else '✗ NEEDS FIX'}")
    print()
