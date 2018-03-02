import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.core.internals import ExtensionBlock

from .base import BaseExtensionTests


class BaseConstructorsTests(BaseExtensionTests):

    def test_array_from_scalars(self, data):
        scalars = [data[0], data[1], data[2]]
        result = data._constructor_from_sequence(scalars)
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
        assert result.dtypes['A'] == data.dtype
        assert result.shape == (len(data), 1)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    def test_dataframe_from_series(self, data):
        result = pd.DataFrame(pd.Series(data))
        assert result.dtypes[0] == data.dtype
        assert result.shape == (len(data), 1)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    @pytest.mark.xfail(reason="GH-19342")
    def test_series_given_mismatched_index_raises(self, data):
        msg = 'Wrong number of items passed 3, placement implies 4'
        with tm.assert_raises_regex(ValueError, None) as m:
            pd.Series(data[:3], index=[0, 1, 2, 3, 4])

        assert m.match(msg)
