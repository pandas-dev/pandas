import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame


@td.skip_if_installed("bodo")
def test_bodo_not_installed_df_apply():
    "Test that importing bodo when not installed results in ImportError."

    df = DataFrame({"A": [1, 2, 3, 4, 5]})

    def f(x):
        return 1

    with pytest.raises(ImportError, match="Missing optional"):
        df.apply(f, engine="bodo")
