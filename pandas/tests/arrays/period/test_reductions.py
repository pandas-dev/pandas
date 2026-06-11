import pandas as pd


class TestReductions:
    def test_min_max(self):
        arr = pd.PeriodIndex(
            [
                "2000-01-03",
                "2000-01-03",
                "NaT",
                "2000-01-02",
                "2000-01-05",
                "2000-01-04",
            ],
            freq="D",
        ).array

        result = arr.min()
        expected = pd.Period("2000-01-02", freq="D")
        assert result == expected

        result = arr.max()
        expected = pd.Period("2000-01-05", freq="D")
        assert result == expected

        result = arr.min(skipna=False)
        assert result is pd.NaT

        result = arr.max(skipna=False)
        assert result is pd.NaT

    def test_min_max_empty(self, skipna):
        arr = pd.PeriodIndex([], freq="D").array
        result = arr.min(skipna=skipna)
        assert result is pd.NaT

        result = arr.max(skipna=skipna)
        assert result is pd.NaT
