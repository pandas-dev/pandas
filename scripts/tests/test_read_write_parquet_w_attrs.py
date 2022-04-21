import pytest
from typing import Union
import pandas as pd
from os import remove

@pytest.fixture
def dataframe(attr_param: Union[list, None, int, float, str]) -> pd.DataFrame:
    df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
    df.attrs['test_attr'] = attr_param
    return df

@pytest.mark.parametrize('attr_param', [[1,2,3], ['a', 'b', 'c'],None, 'test_attr_str', 1, 1.0])
def test_read_write_pyarrow(dataframe):
    dataframe.to_parquet('test_read_write_parquet_w_attrs.parquet', 
                         engine='pyarrow')
    dataframe_read = pd.read_parquet('test_read_write_parquet_w_attrs.parquet', 
                                     engine='pyarrow')
    remove('test_read_write_parquet_w_attrs.parquet')
    assert dataframe.attrs == dataframe_read.attrs
    
@pytest.mark.parametrize('attr_param', [[1,2,3], ['a', 'b', 'c'],None, 'test_attr_str', 1, 1.0])
def test_read_write_fastparquet(dataframe):
    dataframe.to_parquet('test_read_write_parquet_w_attrs.parquet', 
                         engine='fastparquet')
    dataframe_read = pd.read_parquet('test_read_write_parquet_w_attrs.parquet', 
                                     engine='fastparquet')
    remove('test_read_write_prquet_w_attrs.parquet')
    assert dataframe.attrs == dataframe_read.attrs
    