import pytest

import pandas as pd


@pytest.fixture(name="df_from_dict")
def fixture_df_from_dict():
    def maker(dct, is_categorical=False):
        df = pd.DataFrame(dct)
        return df.astype("category") if is_categorical else df

    return maker
