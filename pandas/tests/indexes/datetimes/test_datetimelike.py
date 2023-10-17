""" generic tests from the Datetimelike class """
from pandas import date_range
import pandas._testing as tm


class TestDatetimeIndex:
    def test_format(self):
        # GH35439
        idx = date_range("20130101", periods=5)
        expected = [f"{x:%Y-%m-%d}" for x in idx]
        msg = r"DatetimeIndex\.format is deprecated"
        with tm.assert_produces_warning(FutureWarning, match=msg):
            assert idx.format() == expected
