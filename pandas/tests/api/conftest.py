import pytest
import pandas as pd
from pandas.api.exchange.implementation import _from_dataframe

@pytest.fixture(scope='package')
def df_from_dict():
    def maker(dct, is_categorical=False):
        df = pd.DataFrame(dct)
        return df.astype('category') if is_categorical else df
    return maker

@pytest.fixture(scope='package')
def df_from_xchg():
    def maker(xchg):
        return _from_dataframe(xchg)
    return maker
