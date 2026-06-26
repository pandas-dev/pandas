import pytest

from pandas import (
    DatetimeIndex,
    Index,
    MultiIndex,
    Period,
    date_range,
)
import pandas._testing as tm


class TestMap:
    def test_map(self):
        rng = date_range("1/1/2000", periods=10)

        f = lambda x: x.strftime("%Y%m%d")
        result = rng.map(f)
        exp = Index([f(x) for x in rng])
        tm.assert_index_equal(result, exp)

    def test_map_fallthrough(self, capsys):
        # GH#22067, check we don't get warnings about silently ignored errors
        dti = date_range("2017-01-01", "2018-01-01", freq="B")

        dti.map(lambda x: Period(year=x.year, month=x.month, freq="M"))

        captured = capsys.readouterr()
        assert captured.err == ""

    def test_map_bug_1677(self):
        index = DatetimeIndex(["2012-04-25 09:30:00.393000"])
        f = index.asof

        result = index.map(f)
        expected = Index([f(index[0])])
        tm.assert_index_equal(result, expected)

    def test_map_tz_aware_to_tz_naive(self):
        # GH#57192 a function that strips the tz should give a tz-naive result,
        #  not silently re-localize back to the original tz
        dti = date_range("2024-02-01", "2024-02-03", freq="8h", tz="UTC")

        result = dti.map(lambda ts: ts.normalize().tz_convert(None))
        expected = DatetimeIndex(
            [ts.normalize().tz_convert(None) for ts in dti], dtype="M8[us]"
        )
        tm.assert_index_equal(result, expected)
        assert result.tz is None

    @pytest.mark.parametrize("name", [None, "name"])
    def test_index_map(self, name):
        # see GH#20990
        count = 6
        index = date_range("2018-01-01", periods=count, freq="ME", name=name).map(
            lambda x: (x.year, x.month)
        )
        exp_index = MultiIndex.from_product(((2018,), range(1, 7)), names=[name, name])
        tm.assert_index_equal(index, exp_index)
