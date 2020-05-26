import pytest

import pandas as pd

pytest.importorskip("pyarrow", minversion="0.13.0")

from .arrays import ArrowStringDtype  # isort:skip


def test_constructor_from_list():
    # GH 27673
    result = pd.Series(["E"], dtype=ArrowStringDtype())
    assert isinstance(result.dtype, ArrowStringDtype)
