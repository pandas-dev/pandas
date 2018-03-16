import pytest

import pandas as pd
from pandas.core.internals import ExtensionBlock

from .base import BaseExtensionTests


class BaseReshapingTests(BaseExtensionTests):
    """Tests for reshaping and concatenation."""
    @pytest.mark.parametrize('in_frame', [True, False])
    def test_concat(self, data, in_frame):
        wrapped = pd.Series(data)
        if in_frame:
            wrapped = pd.DataFrame(wrapped)
        result = pd.concat([wrapped, wrapped], ignore_index=True)

        assert len(result) == len(data) * 2

        if in_frame:
            dtype = result.dtypes[0]
        else:
            dtype = result.dtype

        assert dtype == data.dtype
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    @pytest.mark.parametrize('in_frame', [True, False])
    def test_concat_all_na_block(self, data_missing, in_frame):
        valid_block = pd.Series(data_missing.take([1, 1]), index=[0, 1])
        na_block = pd.Series(data_missing.take([0, 0]), index=[2, 3])
        if in_frame:
            valid_block = pd.DataFrame({"a": valid_block})
            na_block = pd.DataFrame({"a": na_block})
        result = pd.concat([valid_block, na_block])
        if in_frame:
            expected = pd.DataFrame({"a": data_missing.take([1, 1, 0, 0])})
            self.assert_frame_equal(result, expected)
        else:
            expected = pd.Series(data_missing.take([1, 1, 0, 0]))
            self.assert_series_equal(result, expected)

    def test_align(self, data, na_value):
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.Series(a).align(pd.Series(b, index=[1, 2, 3]))

        # Assumes that the ctor can take a list of scalars of the type
        e1 = pd.Series(type(data)(list(a) + [na_value]))
        e2 = pd.Series(type(data)([na_value] + list(b)))
        self.assert_series_equal(r1, e1)
        self.assert_series_equal(r2, e2)

    def test_align_frame(self, data, na_value):
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.DataFrame({'A': a}).align(
            pd.DataFrame({'A': b}, index=[1, 2, 3])
        )

        # Assumes that the ctor can take a list of scalars of the type
        e1 = pd.DataFrame({'A': type(data)(list(a) + [na_value])})
        e2 = pd.DataFrame({'A': type(data)([na_value] + list(b))})
        self.assert_frame_equal(r1, e1)
        self.assert_frame_equal(r2, e2)

    def test_set_frame_expand_regular_with_extension(self, data):
        df = pd.DataFrame({"A": [1] * len(data)})
        df['B'] = data
        expected = pd.DataFrame({"A": [1] * len(data), "B": data})
        self.assert_frame_equal(df, expected)

    def test_set_frame_expand_extension_with_regular(self, data):
        df = pd.DataFrame({'A': data})
        df['B'] = [1] * len(data)
        expected = pd.DataFrame({"A": data, "B": [1] * len(data)})
        self.assert_frame_equal(df, expected)
