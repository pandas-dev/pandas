import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame

pytestmark = pytest.mark.bodo_udf_engine


def test_bodo_not_installed_df_apply():
    "Test that importing bodo when not installed results in ImportError."
    bodo_installed = bool(td.import_optional_dependency("bodo", errors="ignore"))
    if bodo_installed:
        pytest.skip("bodo is installed.")

    df = DataFrame({"A": [1, 2, 3, 4, 5]})

    def f(x):
        return 1

    with pytest.raises(ImportError, match="Missing optional"):
        df.apply(f, engine="bodo", axis=1)
