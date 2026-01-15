from datetime import timedelta

from pandas import (
    DatetimeIndex,
    Index,
    Timestamp,
    date_range,
    isna,
)


class TestAsOf:
    def test_asof_partial(self):
        index = date_range("2010-01-01", periods=2, freq="ME")
        expected = Timestamp("2010-02-28")
        result = index.asof("2010-02")
        assert result == expected
        assert not isinstance(result, Index)

    def test_asof(self):
        index = date_range("2020-01-01", periods=10)

        dt = index[0]
        assert index.asof(dt) == dt
        assert isna(index.asof(dt - timedelta(1)))

        dt = index[-1]
        assert index.asof(dt + timedelta(1)) == dt

        dt = index[0].to_pydatetime()
        assert isinstance(index.asof(dt), Timestamp)

    def test_asof_datetime_string(self):
        # GH#50946

        dti = date_range("2021-08-05", "2021-08-10", freq="1D")

        key = "2021-08-09"
        res = dti.asof(key)
        exp = dti[4]
        assert res == exp

        # add a non-midnight time caused a bug
        dti2 = DatetimeIndex([*list(dti), "2021-08-11 00:00:01"])
        res = dti2.asof(key)
        assert res == exp
