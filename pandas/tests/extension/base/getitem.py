import numpy as np

import pandas as pd

from .base import BaseExtensionTests


class BaseGetitemTests(BaseExtensionTests):
    """Tests for ExtensionArray.__getitem__."""

    def test_iloc_series(self, data):
        ser = pd.Series(data)
        result = ser.iloc[:4]
        expected = pd.Series(data[:4])
        self.assert_series_equal(result, expected)

        result = ser.iloc[[0, 1, 2, 3]]
        self.assert_series_equal(result, expected)

    def test_iloc_frame(self, data):
        df = pd.DataFrame({"A": data, 'B':
                           np.arange(len(data), dtype='int64')})
        expected = pd.DataFrame({"A": data[:4]})

        # slice -> frame
        result = df.iloc[:4, [0]]
        self.assert_frame_equal(result, expected)

        # sequence -> frame
        result = df.iloc[[0, 1, 2, 3], [0]]
        self.assert_frame_equal(result, expected)

        expected = pd.Series(data[:4], name='A')

        # slice -> series
        result = df.iloc[:4, 0]
        self.assert_series_equal(result, expected)

        # sequence -> series
        result = df.iloc[:4, 0]
        self.assert_series_equal(result, expected)

    def test_loc_series(self, data):
        ser = pd.Series(data)
        result = ser.loc[:3]
        expected = pd.Series(data[:4])
        self.assert_series_equal(result, expected)

        result = ser.loc[[0, 1, 2, 3]]
        self.assert_series_equal(result, expected)

    def test_loc_frame(self, data):
        df = pd.DataFrame({"A": data,
                           'B': np.arange(len(data), dtype='int64')})
        expected = pd.DataFrame({"A": data[:4]})

        # slice -> frame
        result = df.loc[:3, ['A']]
        self.assert_frame_equal(result, expected)

        # sequence -> frame
        result = df.loc[[0, 1, 2, 3], ['A']]
        self.assert_frame_equal(result, expected)

        expected = pd.Series(data[:4], name='A')

        # slice -> series
        result = df.loc[:3, 'A']
        self.assert_series_equal(result, expected)

        # sequence -> series
        result = df.loc[:3, 'A']
        self.assert_series_equal(result, expected)

    def test_getitem_scalar(self, data):
        result = data[0]
        assert isinstance(result, data.dtype.type)

        result = pd.Series(data)[0]
        assert isinstance(result, data.dtype.type)

    def test_getitem_scalar_na(self, data_missing, na_cmp, na_value):
        result = data_missing[0]
        assert na_cmp(result, na_value)

    def test_getitem_mask(self, data):
        # Empty mask, raw array
        mask = np.zeros(len(data), dtype=bool)
        result = data[mask]
        assert len(result) == 0
        assert isinstance(result, type(data))

        # Empty mask, in series
        mask = np.zeros(len(data), dtype=bool)
        result = pd.Series(data)[mask]
        assert len(result) == 0
        assert result.dtype == data.dtype

        # non-empty mask, raw array
        mask[0] = True
        result = data[mask]
        assert len(result) == 1
        assert isinstance(result, type(data))

        # non-empty mask, in series
        result = pd.Series(data)[mask]
        assert len(result) == 1
        assert result.dtype == data.dtype

    def test_getitem_slice(self, data):
        # getitem[slice] should return an array
        result = data[slice(0)]  # empty
        assert isinstance(result, type(data))

        result = data[slice(1)]  # scalar
        assert isinstance(result, type(data))

    def test_take_sequence(self, data):
        result = pd.Series(data)[[0, 1, 3]]
        assert result.iloc[0] == data[0]
        assert result.iloc[1] == data[1]
        assert result.iloc[2] == data[3]
