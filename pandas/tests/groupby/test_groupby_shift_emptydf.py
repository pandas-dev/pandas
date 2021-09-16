import numpy as np
import pytest

from pandas import (
    DataFrame   
)

import pandas._testing as tm

def test_empty_shift_with_fill():
    #This test is designed to replicate the extra named indexes in 41264 caused by having a fill_value parameter
    
    empty_df = pd.DataFrame(columns=["a", "b", "c"])
    
    shifted = empty_df.groupby(['a']).shift(1)
    shifted_with_fill = empty_df.groupby(['a']).shift(1, fill_value=0)
    
    tm.assert_frame_equal(shifted, shifted_with_fill)
    tm.assert_index_equal(shifted.index, shifted_with_fill.index)
    
def test_multindex_empty_shift_with_fill():
    #This test is designed to replicate the extra named indexes in 41264 caused by having a fill_value parameter for multindexes
    
    empty_df = pd.DataFrame(columns=["a", "b", "c"])
    
    shifted = empty_df.groupby(['a','b']).shift(1)
    shifted_with_fill = empty_df.groupby(['a','b']).shift(1, fill_value=0)
        
    tm.assert_frame_equal(shifted, shifted_with_fill)
    tm.assert_index_equal(shifted.index, shifted_with_fill.index)
    
def test_multindex_shift_with_fill():
    #This test is designed to replicate the extra named indexes in 41264. 
    #Both numbered dataframes and empty dataframes should return FrozenList([None])
    
    num_df = pd.DataFrame({'a':[1,2,3], 'b':[4,5,6], 'c':[7,8,9]})
    empty_df = pd.DataFrame(columns=["a", "b", "c"])
    
    shifted = empty_df.groupby(['a','b']).shift(1)
    shifted_with_fill = num_df.groupby(['a','b']).shift(1, fill_value=0)
        
    tm.assert_class_equal(shifted.index.names, shifted_with_fill.index.names)
    

   