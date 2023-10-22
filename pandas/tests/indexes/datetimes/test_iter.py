import dateutil.tz
import numpy as np
import pytest

from pandas import (
    DatetimeIndex,
    date_range,
    to_datetime,
)


class TestDatetimeIndexIteration:
    @pytest.mark.parametrize(
        "tz", [None, "UTC", "US/Central", dateutil.tz.tzoffset(None, -28800)]
    )
    def test_iteration_preserves_nanoseconds(self, tz):
        # GH#19603
        index = DatetimeIndex(
            ["2018-02-08 15:00:00.168456358", "2018-02-08 15:00:00.168456359"], tz=tz
        )
        for i, ts in enumerate(index):
            assert ts == index[i]  # pylint: disable=unnecessary-list-index-lookup

    def test_iter_readonly(self):
        # GH#28055 ints_to_pydatetime with readonly array
        arr = np.array([np.datetime64("2012-02-15T12:00:00.000000000")])
        arr.setflags(write=False)
        dti = to_datetime(arr)
        list(dti)

    def test_iteration_preserves_tz(self):
        # see GH#8890
        index = date_range("2012-01-01", periods=3, freq="h", tz="US/Eastern")

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]  # pylint: disable=unnecessary-list-index-lookup
            assert result == expected

    def test_iteration_preserves_tz2(self):
        index = date_range(
            "2012-01-01", periods=3, freq="h", tz=dateutil.tz.tzoffset(None, -28800)
        )

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]  # pylint: disable=unnecessary-list-index-lookup
            assert result._repr_base == expected._repr_base
            assert result == expected

    def test_iteration_preserves_tz3(self):
        # GH#9100
        index = DatetimeIndex(
            ["2014-12-01 03:32:39.987000-08:00", "2014-12-01 04:12:34.987000-08:00"]
        )
        for i, ts in enumerate(index):
            result = ts
            expected = index[i]  # pylint: disable=unnecessary-list-index-lookup
            assert result._repr_base == expected._repr_base
            assert result == expected

    @pytest.mark.parametrize("periods", [0, 9999, 10000, 10001])
    def test_iteration_over_chunksize(self, periods):
        # GH#21012

        index = date_range("2000-01-01 00:00:00", periods=periods, freq="min")
        num = 0
        for stamp in index:
            assert index[num] == stamp
            num += 1
        assert num == len(index)
