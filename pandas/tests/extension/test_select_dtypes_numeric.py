import numpy as np

from pandas.core.dtypes.dtypes import ExtensionDtype

import pandas as pd
from pandas.core.arrays import ExtensionArray


class DummyDtype(ExtensionDtype):
    type = int

    def __init__(self, numeric):
        self._numeric = numeric

    @property
    def name(self):
        return "Dummy"

    @property
    def _is_numeric(self):
        return self._numeric


class DummyArray(ExtensionArray):
    def __init__(self, data, dtype):
        self.data = data
        self._dtype = dtype

    def __array__(self, dtype):
        return self.data

    @property
    def dtype(self):
        return self._dtype

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, item):
        pass


def test_select_dtypes_numeric():
    da = DummyArray([1, 2], dtype=DummyDtype(numeric=True))
    df = pd.DataFrame(da)
    assert df.select_dtypes(np.number).shape == df.shape

    da = DummyArray([1, 2], dtype=DummyDtype(numeric=False))
    df = pd.DataFrame(da)
    assert df.select_dtypes(np.number).shape != df.shape
