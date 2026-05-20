import pytest
import pandas as pd
import pandas._testing as tm

def assert_cow_inplace_safe(df: pd.DataFrame, operation_func):
    """
    Ensures that an inplace operation does not corrupt aliases 
    if it routes through an external engine (NumExpr/Cython).
    """
    expected = df.copy()
    alias = df[:]
    
    operation_func(df)
    
    df.iloc[0, df.columns.get_loc("a")] = 999
    
    tm.assert_frame_equal(alias, expected)

def test_query_inplace_cow_semantics():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    
    def run_query(data):
        data.query("a > 0", inplace=True)
        
    assert_cow_inplace_safe(df, run_query)

def test_where_inplace_cow_semantics():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    
    def run_where(data):
        data.where(data > 0, inplace=True)
        
    assert_cow_inplace_safe(df, run_where)