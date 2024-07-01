import numpy as np
import pytest

from pandas import (
    Series,
)

class TestSeriesIsConstant:
    def test_isconstant(self):
        s = Series([2, 2, 2, 2])
        result = s.isconstant()
        assert result

        s = Series([2, 2, 2, 3])
        result = s.isconstant()
        assert not result

        s = Series([])
        result = s.isconstant()
        assert result

        s = Series([5])
        result = s.isconstant()
        assert result

    def test_isconstant_with_nan(self):
        s = Series([np.nan])
        result = s.isconstant()
        assert result

        s = Series([np.nan,np.nan])
        result = s.isconstant()
        assert result

        s = Series([np.nan,1])
        result = s.isconstant()
        assert not result

        s = Series([np.nan, np.nan, 1])
        result = s.isconstant()
        assert not result

    def test_isconstant_with_nan_dropna(self):
        s = Series([np.nan])
        result = s.isconstant(True)
        assert result

        s = Series([np.nan,np.nan])
        result = s.isconstant(True)
        assert result

        s = Series([np.nan,1])
        result = s.isconstant(True)
        assert result

        s = Series([np.nan, np.nan, 1])
        result = s.isconstant(True)
        assert result

    def test_isconstant_mixed_types(self):
        s = Series([2, '2', 2])
        result = s.isconstant()
        assert not result

        s = Series([2, 2.0, 2])
        result = s.isconstant()
        assert result