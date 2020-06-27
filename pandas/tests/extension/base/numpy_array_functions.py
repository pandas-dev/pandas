import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd

from .base import BaseExtensionTests


@td.skip_if_no_nep18
class BaseNumpyArrayFunctionTests(BaseExtensionTests):
    def test_tile(self, data):
        expected = pd.array(list(data) * 3, dtype=data.dtype)

        result = np.tile(data, 3)
        self.assert_extension_array_equal(result, expected)

        result = np.tile(data, (3,))
        self.assert_extension_array_equal(result, expected)

        expected = np.array([data.to_numpy()] * 3)

        result = np.tile(data, (3, 1))
        self.assert_numpy_array_equal(result, expected)

    def test_concatenate(self, data):
        expected = pd.array(list(data) * 3, dtype=data.dtype)

        result = np.concatenate([data] * 3)
        self.assert_extension_array_equal(result, expected)

        result = np.concatenate([data] * 3, axis=0)
        self.assert_extension_array_equal(result, expected)

        msg = "axis 1 is out of bounds for array of dimension 1"
        with pytest.raises(np.AxisError, match=msg):
            np.concatenate([data] * 3, axis=1)

    def test_ndim(self, data):
        assert np.ndim(data) == 1
