from pandas import DataFrame, date_range
import pandas._testing as tm


class TestTZLocalize:
    # See also:
    # test_tz_convert_and_localize in test_tz_convert

    def test_frame_tz_localize(self):
        rng = date_range("1/1/2011", periods=100, freq="H")

        df = DataFrame({"a": 1}, index=rng)
        result = df.tz_localize("utc")
        expected = DataFrame({"a": 1}, rng.tz_localize("UTC"))
        assert result.index.tz.zone == "UTC"
        tm.assert_frame_equal(result, expected)

        df = df.T
        result = df.tz_localize("utc", axis=1)
        assert result.columns.tz.zone == "UTC"
        tm.assert_frame_equal(result, expected.T)
