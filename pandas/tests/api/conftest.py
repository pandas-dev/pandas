import pytest
import pandas as pd

@pytest.fixture(scope='package')
def create_df_from_dict():
    def maker(dct):
        return pd.DataFrame(dct)
    return maker
