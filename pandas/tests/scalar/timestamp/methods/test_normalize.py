import pytest

from pandas._libs.tslibs import Timestamp
from pandas._libs.tslibs.dtypes import NpyDatetimeUnit
from pandas._libs.tslibs.np_datetime import OutOfBoundsDatetime


class TestTimestampNormalize:
    @pytest.mark.parametrize("arg", ["2013-11-30", "2013-11-30 12:00:00"])
    def test_normalize(self, tz_naive_fixture, arg, unit):
        tz = tz_naive_fixture
        ts = Timestamp(arg, tz=tz).as_unit(unit)
        result = ts.normalize()
        expected = Timestamp("2013-11-30", tz=tz)
        assert result == expected
        assert result._creso == getattr(NpyDatetimeUnit, f"NPY_FR_{unit}").value

    def test_normalize_pre_epoch_dates(self):
        # GH: 36294
        result = Timestamp("1969-01-01 09:00:00").normalize()
        expected = Timestamp("1969-01-01 00:00:00")
        assert result == expected

    def test_normalize_edge_cases(self):
        # GH: 60583
        expected_msg = (
            r"Cannot normalize 1677-09-21 00:12:43\.145224193 to midnight "
            "without overflow"
        )
        with pytest.raises(OutOfBoundsDatetime, match=expected_msg):
            Timestamp.min.normalize()

        result = Timestamp.max.normalize()
        excepted = Timestamp("2262-04-11 00:00:00")
        assert result == excepted
