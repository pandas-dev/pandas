from datetime import datetime

import numpy as np
import pytest
import pytz

from pandas import NA, DatetimeIndex, Index, NaT, Timestamp, date_range
import pandas._testing as tm


class TestInsert:
    @pytest.mark.parametrize("null", [None, np.nan, np.datetime64("NaT"), NaT, NA])
    @pytest.mark.parametrize("tz", [None, "UTC", "US/Eastern"])
    def test_insert_nat(self, tz, null):
        # GH#16537, GH#18295 (test missing)
        idx = DatetimeIndex(["2017-01-01"], tz=tz)
        expected = DatetimeIndex(["NaT", "2017-01-01"], tz=tz)
        res = idx.insert(0, null)
        tm.assert_index_equal(res, expected)

    @pytest.mark.parametrize("tz", [None, "UTC", "US/Eastern"])
    def test_insert_invalid_na(self, tz):
        idx = DatetimeIndex(["2017-01-01"], tz=tz)
        msg = "value should be a 'Timestamp' or 'NaT'. Got 'timedelta64' instead."
        with pytest.raises(TypeError, match=msg):
            idx.insert(0, np.timedelta64("NaT"))

    def test_insert_empty_preserves_freq(self, tz_naive_fixture):
        # GH#33573
        tz = tz_naive_fixture
        dti = DatetimeIndex([], tz=tz, freq="D")
        item = Timestamp("2017-04-05").tz_localize(tz)

        result = dti.insert(0, item)
        assert result.freq == dti.freq

        # But not when we insert an item that doesnt conform to freq
        dti = DatetimeIndex([], tz=tz, freq="W-THU")
        result = dti.insert(0, item)
        assert result.freq is None

    def test_insert(self):
        idx = DatetimeIndex(["2000-01-04", "2000-01-01", "2000-01-02"], name="idx")

        result = idx.insert(2, datetime(2000, 1, 5))
        exp = DatetimeIndex(
            ["2000-01-04", "2000-01-01", "2000-01-05", "2000-01-02"], name="idx"
        )
        tm.assert_index_equal(result, exp)

        # insertion of non-datetime should coerce to object index
        result = idx.insert(1, "inserted")
        expected = Index(
            [
                datetime(2000, 1, 4),
                "inserted",
                datetime(2000, 1, 1),
                datetime(2000, 1, 2),
            ],
            name="idx",
        )
        assert not isinstance(result, DatetimeIndex)
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

        idx = date_range("1/1/2000", periods=3, freq="M", name="idx")

        # preserve freq
        expected_0 = DatetimeIndex(
            ["1999-12-31", "2000-01-31", "2000-02-29", "2000-03-31"],
            name="idx",
            freq="M",
        )
        expected_3 = DatetimeIndex(
            ["2000-01-31", "2000-02-29", "2000-03-31", "2000-04-30"],
            name="idx",
            freq="M",
        )

        # reset freq to None
        expected_1_nofreq = DatetimeIndex(
            ["2000-01-31", "2000-01-31", "2000-02-29", "2000-03-31"],
            name="idx",
            freq=None,
        )
        expected_3_nofreq = DatetimeIndex(
            ["2000-01-31", "2000-02-29", "2000-03-31", "2000-01-02"],
            name="idx",
            freq=None,
        )

        cases = [
            (0, datetime(1999, 12, 31), expected_0),
            (-3, datetime(1999, 12, 31), expected_0),
            (3, datetime(2000, 4, 30), expected_3),
            (1, datetime(2000, 1, 31), expected_1_nofreq),
            (3, datetime(2000, 1, 2), expected_3_nofreq),
        ]

        for n, d, expected in cases:
            result = idx.insert(n, d)
            tm.assert_index_equal(result, expected)
            assert result.name == expected.name
            assert result.freq == expected.freq

        # reset freq to None
        result = idx.insert(3, datetime(2000, 1, 2))
        expected = DatetimeIndex(
            ["2000-01-31", "2000-02-29", "2000-03-31", "2000-01-02"],
            name="idx",
            freq=None,
        )
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name
        assert result.freq is None

        # see gh-7299
        idx = date_range("1/1/2000", periods=3, freq="D", tz="Asia/Tokyo", name="idx")
        with pytest.raises(TypeError, match="Cannot compare tz-naive and tz-aware"):
            idx.insert(3, Timestamp("2000-01-04"))
        with pytest.raises(TypeError, match="Cannot compare tz-naive and tz-aware"):
            idx.insert(3, datetime(2000, 1, 4))
        with pytest.raises(ValueError, match="Timezones don't match"):
            idx.insert(3, Timestamp("2000-01-04", tz="US/Eastern"))
        with pytest.raises(ValueError, match="Timezones don't match"):
            idx.insert(3, datetime(2000, 1, 4, tzinfo=pytz.timezone("US/Eastern")))

        for tz in ["US/Pacific", "Asia/Singapore"]:
            idx = date_range("1/1/2000 09:00", periods=6, freq="H", tz=tz, name="idx")
            # preserve freq
            expected = date_range(
                "1/1/2000 09:00", periods=7, freq="H", tz=tz, name="idx"
            )
            for d in [
                Timestamp("2000-01-01 15:00", tz=tz),
                pytz.timezone(tz).localize(datetime(2000, 1, 1, 15)),
            ]:

                result = idx.insert(6, d)
                tm.assert_index_equal(result, expected)
                assert result.name == expected.name
                assert result.freq == expected.freq
                assert result.tz == expected.tz

            expected = DatetimeIndex(
                [
                    "2000-01-01 09:00",
                    "2000-01-01 10:00",
                    "2000-01-01 11:00",
                    "2000-01-01 12:00",
                    "2000-01-01 13:00",
                    "2000-01-01 14:00",
                    "2000-01-01 10:00",
                ],
                name="idx",
                tz=tz,
                freq=None,
            )
            # reset freq to None
            for d in [
                Timestamp("2000-01-01 10:00", tz=tz),
                pytz.timezone(tz).localize(datetime(2000, 1, 1, 10)),
            ]:
                result = idx.insert(6, d)
                tm.assert_index_equal(result, expected)
                assert result.name == expected.name
                assert result.tz == expected.tz
                assert result.freq is None

    @pytest.mark.parametrize(
        "item", [0, np.int64(0), np.float64(0), np.array(0), np.timedelta64(456)]
    )
    def test_insert_mismatched_types_raises(self, tz_aware_fixture, item):
        # GH#33703 dont cast these to dt64
        tz = tz_aware_fixture
        dti = date_range("2019-11-04", periods=9, freq="-1D", name=9, tz=tz)

        msg = "value should be a 'Timestamp' or 'NaT'. Got '.*' instead"
        with pytest.raises(TypeError, match=msg):
            dti.insert(1, item)

    def test_insert_object_casting(self, tz_aware_fixture):
        # GH#33703
        tz = tz_aware_fixture
        dti = date_range("2019-11-04", periods=3, freq="-1D", name=9, tz=tz)

        # ATM we treat this as a string, but we could plausibly wrap it in Timestamp
        value = "2019-11-05"
        result = dti.insert(0, value)
        expected = Index(["2019-11-05"] + list(dti), dtype=object, name=9)
        tm.assert_index_equal(result, expected)
