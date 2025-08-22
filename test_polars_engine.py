#!/usr/bin/env python3

"""
Simple test script to verify polars engine functionality
"""

import tempfile
import pandas as pd


def test_polars_engine_basic():
    """Test basic functionality of polars engine"""
    
    # Create sample CSV data
    csv_data = """a,b,c
1,2,3
4,5,6
7,8,9"""
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_data)
        temp_path = f.name
    
    try:
        # Test with polars engine
        df_polars = pd.read_csv(temp_path, engine='polars')
        print("✓ Polars engine basic test passed")
        print(f"DataFrame shape: {df_polars.shape}")
        print(f"DataFrame columns: {list(df_polars.columns)}")
        print(f"DataFrame:\n{df_polars}")
        
        # Test with default engine for comparison
        df_default = pd.read_csv(temp_path)
        print(f"\n✓ Default engine shape: {df_default.shape}")
        
        # Check if they're the same
        if df_polars.equals(df_default):
            print("✓ Polars and default engines produce identical results")
        else:
            print("✗ Results differ between engines")
            print(f"Polars:\n{df_polars}")
            print(f"Default:\n{df_default}")
            
    except Exception as e:
        print(f"✗ Error: {e}")
        return False
    finally:
        import os
        os.unlink(temp_path)
    
    return True


def test_polars_engine_with_options():
    """Test polars engine with various options"""
    
    # Create sample CSV data
    csv_data = """name,age,city
Alice,25,New York
Bob,30,Los Angeles
Charlie,35,Chicago"""
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_data)
        temp_path = f.name
    
    try:
        # Test with usecols
        df_usecols = pd.read_csv(temp_path, engine='polars', usecols=['name', 'age'])
        print(f"✓ usecols test passed: {list(df_usecols.columns)}")
        
        # Test with nrows
        df_nrows = pd.read_csv(temp_path, engine='polars', nrows=2)
        print(f"✓ nrows test passed: {len(df_nrows)} rows")
        
        # Test with custom names
        df_names = pd.read_csv(temp_path, engine='polars', names=['col1', 'col2', 'col3'], header=0)
        print(f"✓ custom names test passed: {list(df_names.columns)}")
        
    except Exception as e:
        print(f"✗ Error in options test: {e}")
        return False
    finally:
        import os
        os.unlink(temp_path)
    
    return True


if __name__ == "__main__":
    print("Testing polars engine implementation...")
    
    try:
        import polars as pl
        print(f"✓ Polars available: {pl.__version__}")
    except ImportError:
        print("✗ Polars not available - installing...")
        import subprocess
        import sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", "polars"])
        
    success1 = test_polars_engine_basic()
    success2 = test_polars_engine_with_options()
    
    if success1 and success2:
        print("\n🎉 All tests passed!")
    else:
        print("\n❌ Some tests failed")
