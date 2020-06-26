import numpy as np

import pandas as pd

from .base import BaseExtensionTests


class BaseNumpyArrayFunctionTests(BaseExtensionTests):
    def test_tile(self, data):
        expected = pd.array(list(data) * 3, dtype=data.dtype)
        result = np.tile(data, 3)
        self.assert_extension_array_equal(result, expected)
