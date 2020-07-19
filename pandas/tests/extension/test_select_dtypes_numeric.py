import pytest
import numpy as np

import pandas as pd
from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.dtypes import ExtensionDtype


class DummyDtype(ExtensionDtype):
    type = int
    _numeric = False

    @property
    def name(self):
        return "Dummy"

    @property
    def _is_numeric(self):
        return self._numeric


class DummyArray(ExtensionArray):
    _dtype = DummyDtype()

    def __init__(self, data):
        self.data = data

    def __array__(self, dtype):
        return self.data

    @property
    def dtype(self):
        return self._dtype

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, item):
        pass


@pytest.mark.parametrize("numeric", [True, False])
def test_select_dtypes_numeric(numeric):
    da = DummyArray([1, 2])
    da._dtype._numeric = numeric
    df = pd.DataFrame(da)
    if numeric:
        assert df.select_dtypes(np.number).shape == df.shape
    else:
        assert df.select_dtypes(np.number).shape != df.shape
