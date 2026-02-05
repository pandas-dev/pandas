import pytest

import pandas.util._test_decorators as td

import pandas as pd


@td.skip_if_installed("tables")
def test_pytables_raises(temp_h5_path):
    df = pd.DataFrame({"A": [1, 2]})
    with pytest.raises(ImportError, match="tables"):
        df.to_hdf(temp_h5_path, key="df")
