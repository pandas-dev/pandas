import pytest

from pandas._libs.tslibs import Timestamp
from pandas._libs.tslibs.dtypes import NpyDatetimeUnit


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

    def test_normalize_overflow_raises(self):
        # GH#60583
        ts = Timestamp.min

        msg = "Cannot normalize Timestamp without integer overflow"
        with pytest.raises(ValueError, match=msg):
            ts.normalize()
