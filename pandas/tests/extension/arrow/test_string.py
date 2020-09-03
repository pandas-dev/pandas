import pytest

import pandas as pd

pytest.importorskip("pyarrow", minversion="0.13.0")


def test_constructor_from_list():
    # GH 27673
    result = pd.Series(["E"], dtype=pd.StringDtype(storage="pyarrow"))
    assert isinstance(result.dtype, pd.StringDtype)
    assert result.dtype.storage == "pyarrow"
