import pytest

from pandas.errors import Pandas4Warning

from pandas import PeriodIndex
import pandas._testing as tm

msg = "PeriodIndex.is_full is deprecated"


@pytest.mark.parametrize(
    "years, expected",
    [
        (["2005", "2007", "2009"], False),
        (["2005", "2006", "2007"], True),
        (["2005", "2005", "2007"], False),
        (["2005", "2005", "2006"], True),
        ([], True),
    ],
)
def test_is_full(years, expected):
    index = PeriodIndex(years, freq="Y")
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = index.is_full
    assert result is expected


def test_is_full_not_monotonic():
    index = PeriodIndex(["2006", "2005", "2005"], freq="Y")
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        with pytest.raises(ValueError, match="Index is not monotonic"):
            index.is_full
