import numpy as np
import pytest

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

    @pytest.mark.parametrize("copy", [True, False])
    def test_tz_localize_copy_inplace_mutate(self, copy, frame_or_series):
        # GH#6326
        obj = frame_or_series(
            np.arange(0, 5), index=date_range("20131027", periods=5, freq="1H", tz=None)
        )
        orig = obj.copy()
        result = obj.tz_localize("UTC", copy=copy)
        expected = frame_or_series(
            np.arange(0, 5),
            index=date_range("20131027", periods=5, freq="1H", tz="UTC"),
        )
        tm.assert_equal(result, expected)
        tm.assert_equal(obj, orig)
        assert result.index is not obj.index
        assert result is not obj
