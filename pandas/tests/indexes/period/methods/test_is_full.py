import pytest

from pandas.errors import Pandas4Warning

from pandas import PeriodIndex
import pandas._testing as tm

msg = "PeriodIndex.is_full is deprecated"


class TestIsFull:
    def test_is_full(self):
        index = PeriodIndex([2005, 2007, 2009], freq="Y")
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            assert not index.is_full

        index = PeriodIndex([2005, 2006, 2007], freq="Y")
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            assert index.is_full

        index = PeriodIndex([2005, 2005, 2007], freq="Y")
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            assert not index.is_full

        index = PeriodIndex([2005, 2005, 2006], freq="Y")
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            assert index.is_full

    def test_is_full_not_monotonic(self):
        index = PeriodIndex([2006, 2005, 2005], freq="Y")
        with pytest.raises(ValueError, match="Index is not monotonic"):
            with tm.assert_produces_warning(Pandas4Warning, match=msg):
                index.is_full

    def test_is_full_empty(self):
        index = PeriodIndex([], freq="Y")
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            assert index.is_full
