import pytest

import numpy as np
import pandas as pd
import pandas.util.testing as tm
from pandas.core.dtypes.dtypes import PeriodDtype
from pandas.core.period import PeriodArray


class TestArray:

    def test_init(self):
        arr = PeriodArray([2017, 2018], freq='A')
        assert isinstance(arr, PeriodArray)

    def test_concat(self):
        p1 = PeriodArray([2017, 2018], freq='A')
        p2 = PeriodArray([2019, 2020], freq='A')
        result = pd.concat([pd.Series(p1), pd.Series(p2)], ignore_index=True)
        expected = pd.Series(PeriodArray([2017, 2018, 2019, 2020], freq='A'))
        tm.assert_series_equal(result, expected)

    def test_equals(self):
        p1 = PeriodArray([2017, 2018], freq='A')
        p2 = PeriodArray([2017, 2018], freq='A')
        assert p1.equals(p2)

    @pytest.mark.parametrize('other', [
        2017,
        [2017, 2018],
        PeriodArray([2016, 2017], freq='A'),
        PeriodArray([2017, 2018], freq='A-JAN'),
        PeriodArray([2017, 2018, 2019], freq='A'),
    ])
    def test_equals_unequal(self, other):
        p1 = PeriodArray([2017, 2018], freq='A')
        assert not p1.equals(other)

    def test_getitem(self):
        p1 = PeriodArray([2017, 2018, 2019], freq='A')
        result = p1[0]
        expected = pd.Period(2017, freq='A')
        assert result == expected

        result = p1[[0, 1]]
        expected = PeriodArray([2017, 2018], freq='A')
        assert result.equals(expected)

        result = p1[slice(2)]
        assert result.equals(expected)

        result = p1[np.array([True, True, False])]
        assert result.equals(expected)

    def test_isna(self):
        result = PeriodArray(['2018', 'NaT'], freq='D').isna()
        expected = np.array([False, True])
        tm.assert_numpy_array_equal(result, expected)


class TestInContainers:

    def test_series_constructor(self):
        result = pd.Series(PeriodArray([2017, 2018], freq='D'))
        assert len(result) == 2
        assert result.dtype == PeriodDtype('D')

    def test_slice(self):
        ser = pd.Series(PeriodArray(pd.period_range(2000, periods=100)))
        result = ser.iloc[:4]
        expected = pd.Series(PeriodArray(['2000-01-01', '2000-01-02',
                                          '2000-01-03', '2000-01-04'],
                                         freq='D'))
        tm.assert_series_equal(result, expected)

        result = ser.loc[[0, 1, 2, 3]]
        tm.assert_series_equal(result, expected)
