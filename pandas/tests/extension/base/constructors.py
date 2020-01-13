import numpy as np
import pytest

import pandas as pd
from pandas.core.internals import ExtensionBlock

from .base import BaseExtensionTests


class BaseConstructorsTests(BaseExtensionTests):
    def test_from_sequence_from_cls(self, data):
        result = type(data)._from_sequence(data, dtype=data.dtype)
        self.assert_extension_array_equal(result, data)

        data = data[:0]
        result = type(data)._from_sequence(data, dtype=data.dtype)
        self.assert_extension_array_equal(result, data)

    def test_array_from_scalars(self, data):
        scalars = [data[0], data[1], data[2]]
        result = data._from_sequence(scalars)
        assert isinstance(result, type(data))

    def test_series_constructor(self, data):
        result = pd.Series(data)
        assert result.dtype == data.dtype
        assert len(result) == len(data)
        assert isinstance(result._data.blocks[0], ExtensionBlock)
        assert result._data.blocks[0].values is data

        # Series[EA] is unboxed / boxed correctly
        result2 = pd.Series(result)
        assert result2.dtype == data.dtype
        assert isinstance(result2._data.blocks[0], ExtensionBlock)

    @pytest.mark.parametrize("from_series", [True, False])
    def test_dataframe_constructor_from_dict(self, data, from_series):
        if from_series:
            data = pd.Series(data)
        result = pd.DataFrame({"A": data})
        assert result.dtypes["A"] == data.dtype
        assert result.shape == (len(data), 1)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    def test_dataframe_from_series(self, data):
        result = pd.DataFrame(pd.Series(data))
        assert result.dtypes[0] == data.dtype
        assert result.shape == (len(data), 1)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    def test_series_given_mismatched_index_raises(self, data):
        msg = "Length of passed values is 3, index implies 5"
        with pytest.raises(ValueError, match=msg):
            pd.Series(data[:3], index=[0, 1, 2, 3, 4])

    def test_from_dtype(self, data):
        # construct from our dtype & string dtype
        dtype = data.dtype

        expected = pd.Series(data)
        result = pd.Series(list(data), dtype=dtype)
        self.assert_series_equal(result, expected)

        result = pd.Series(list(data), dtype=str(dtype))
        self.assert_series_equal(result, expected)

        # gh-30280

        expected = pd.DataFrame(data).astype(dtype)
        result = pd.DataFrame(list(data), dtype=dtype)
        self.assert_frame_equal(result, expected)

        result = pd.DataFrame(list(data), dtype=str(dtype))
        self.assert_frame_equal(result, expected)

    def test_pandas_array(self, data):
        # pd.array(extension_array) should be idempotent...
        result = pd.array(data)
        self.assert_extension_array_equal(result, data)

    def test_pandas_array_dtype(self, data):
        # ... but specifying dtype will override idempotency
        result = pd.array(data, dtype=np.dtype(object))
        expected = pd.arrays.PandasArray(np.asarray(data, dtype=object))
        self.assert_equal(result, expected)
