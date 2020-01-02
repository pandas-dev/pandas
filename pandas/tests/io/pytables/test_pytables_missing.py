import pytest

import pandas.util._test_decorators as td
import pandas.util._testing as tm

import pandas as pd


@td.skip_if_installed("tables")
def test_pytables_raises():
    df = pd.DataFrame({"A": [1, 2]})
    with pytest.raises(ImportError, match="tables"):
        with tm.ensure_clean("foo.h5") as path:
            df.to_hdf(path, "df")
