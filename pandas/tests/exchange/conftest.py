import pytest

import pandas as pd
from pandas.core.exchange.from_dataframe import _from_dataframe


@pytest.fixture(scope="package")
def df_from_dict():
    def maker(dct, is_categorical=False):
        df = pd.DataFrame(dct)
        return df.astype("category") if is_categorical else df

    return maker
