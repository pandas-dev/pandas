import numpy as np
import pytest

from pandas import (
    DataFrame,
    DatetimeIndex,
)


class TestResampleDSTAfricaCairo:
    """DST transition tests for Africa/Cairo timezone."""

    def test_resample_across_dst_transition(self):
        df = DataFrame(
            {"value": [1, 2]},
            index=DatetimeIndex(
                [
                    "2024-04-26 01:00:00",
                    "2024-04-27 00:00:00",
                ]
            ).tz_localize("Africa/Cairo", nonexistent="shift_forward"),
        )

        result = df.resample("D").mean()

        assert len(result) == 2
        assert isinstance(result.index, DatetimeIndex)
        assert result.index.tz is not None
        assert not result.isna().any().any()

    def test_resample_before_dst_boundary(self):
        df = DataFrame(
            {"value": [76.0, 42.0]},
            index=DatetimeIndex(
                [
                    "2024-04-24 00:00:00",
                    "2024-04-25 00:00:00",
                ]
            ).tz_localize("Africa/Cairo"),
        )

        result = df.resample("D").mean()

        assert len(result) == 2
        assert isinstance(result.index, DatetimeIndex)
        assert "Africa/Cairo" in str(result.index.tz)
        assert result.iloc[0, 0] == 76.0
        assert result.iloc[1, 0] == 42.0

    @pytest.mark.parametrize("freq", ["2h", "6h", "12h"])
    def test_resample_various_freq(self, freq):
        df = DataFrame(
            {"value": [1, 2, 3, 4, 5]},
            index=DatetimeIndex(
                [
                    "2024-04-25 22:00:00",
                    "2024-04-25 23:00:00",
                    "2024-04-26 01:00:00",
                    "2024-04-26 02:00:00",
                    "2024-04-26 03:00:00",
                ]
            ).tz_localize("Africa/Cairo", nonexistent="shift_forward"),
        )

        result = df.resample(freq).mean()

        assert isinstance(result, DataFrame)
        assert len(result) > 0
        assert not result.isna().all().any()

    def test_resample_closed_label_combinations(self):
        df = DataFrame(
            {"value": [1, 2]},
            index=DatetimeIndex(
                [
                    "2024-04-26 01:00:00",
                    "2024-04-27 00:00:00",
                ]
            ).tz_localize("Africa/Cairo", nonexistent="shift_forward"),
        )

        for closed in ["left", "right"]:
            for label in ["left", "right"]:
                result = df.resample("D", closed=closed, label=label).mean()
                assert len(result) >= 1
                assert not result.isna().all().any()

    def test_resample_nonexistent_times(self):
        timestamps = [
            "2024-04-25 23:00:00",
            "2024-04-26 00:30:00",
            "2024-04-26 01:00:00",
        ]

        df = DataFrame(
            {"value": [1, 2, 3]},
            index=DatetimeIndex(timestamps).tz_localize(
                "Africa/Cairo", nonexistent="shift_forward"
            ),
        )

        result = df.resample("h").mean()

        assert len(result) > 0
        assert isinstance(result, DataFrame)

    def test_resample_empty_dataframe(self):
        df = DataFrame({"value": []}, index=DatetimeIndex([], tz="Africa/Cairo"))

        result = df.resample("D").mean()

        assert len(result) == 0
        assert isinstance(result.index, DatetimeIndex)

    def test_resample_single_point(self):
        df = DataFrame(
            {"value": [42.0]},
            index=DatetimeIndex(["2024-04-26 12:00:00"]).tz_localize(
                "Africa/Cairo", nonexistent="shift_forward"
            ),
        )

        result = df.resample("D").mean()

        assert len(result) == 1
        assert result.iloc[0, 0] == 42.0


class TestResampleDSTMultipleTimezones:
    """DST handling across multiple timezones."""

    def test_resample_multiple_timezones(self):
        timezones = [
            ("Africa/Cairo", "2024-04-26 01:00:00", "2024-04-27 00:00:00"),
            ("Europe/London", "2024-03-31 01:00:00", "2024-04-01 00:00:00"),
            ("America/New_York", "2024-03-10 01:00:00", "2024-03-11 00:00:00"),
        ]

        for tz, start, end in timezones:
            df = DataFrame(
                {"value": [1, 2]},
                index=DatetimeIndex([start, end]).tz_localize(
                    tz, nonexistent="shift_forward", ambiguous=True
                ),
            )

            result = df.resample("D").mean()

            assert len(result) >= 1
            assert isinstance(result.index, DatetimeIndex)
            assert result.index.tz is not None


class TestResampleDSTEdgeCases:
    """Edge cases around DST transitions."""

    def test_resample_multiple_dst_days(self):
        df = DataFrame(
            {"value": [1, 2, 3, 4]},
            index=DatetimeIndex(
                [
                    "2024-04-25 23:00:00",
                    "2024-04-26 01:00:00",
                    "2024-04-27 00:00:00",
                    "2024-04-28 00:00:00",
                ]
            ).tz_localize("Africa/Cairo", nonexistent="shift_forward"),
        )

        result = df.resample("D").mean()

        assert len(result) >= 3

    def test_resample_microsecond_precision(self):
        df = DataFrame(
            {"value": [1.1, 2.2]},
            index=DatetimeIndex(
                [
                    "2024-04-26 01:00:00.123456",
                    "2024-04-27 00:00:00.654321",
                ]
            ).tz_localize("Africa/Cairo", nonexistent="shift_forward"),
        )

        result = df.resample("D").mean()

        assert len(result) == 2

    def test_resample_with_na_values(self):
        df = DataFrame(
            {"value": [1.0, np.nan, 3.0]},
            index=DatetimeIndex(
                [
                    "2024-04-25 23:00:00",
                    "2024-04-26 01:00:00",
                    "2024-04-26 02:00:00",
                ]
            ).tz_localize("Africa/Cairo", nonexistent="shift_forward"),
        )

        result = df.resample("h").mean()

        assert len(result) > 0
        assert isinstance(result, DataFrame)


class TestResampleDSTOriginalIssues:
    """Tests reproducing the originally reported issues."""

    def test_original_issue_1(self):
        df = DataFrame(
            {"value": [1, 2]},
            index=DatetimeIndex(
                [
                    "2024-04-26 01:00:00",
                    "2024-04-27 00:00:00",
                ]
            ).tz_localize("Africa/Cairo", nonexistent="shift_forward"),
        )

        result = df.resample("D").mean()

        assert len(result) > 0
        assert not result.isna().any().any()

    def test_original_issue_2(self):
        df = DataFrame(
            {"value": [76.0, 42.0]},
            index=DatetimeIndex(
                [
                    "2024-04-24 00:00:00",
                    "2024-04-25 00:00:00",
                ]
            ).tz_localize("Africa/Cairo"),
        )

        result = df.resample("D").mean()

        assert len(result) > 0
        assert not result.isna().any().any()
