import pytest

import pandas.util._test_decorators as td

import pandas as pd


@td.skip_if_installed("tables")
def test_pytables_raises(tmp_path):
    df = pd.DataFrame({"A": [1, 2]})
    path = tmp_path / "foo.h5"
    with pytest.raises(ImportError, match="tables"):
        df.to_hdf(path, key="df")
