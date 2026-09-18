# NB: This is for the Timestamp.timestamp *method* specifically, not
# the Timestamp class in general.
from datetime import (
    UTC,
    datetime,
)

import pytest

from pandas._libs.tslibs import Timestamp
from pandas.compat import WASM
from pandas.errors import Pandas4Warning
import pandas.util._test_decorators as td

import pandas._testing as tm

naive_msg = "Timestamp.timestamp treating a tz-naive Timestamp as UTC"


class TestTimestampMethod:
    @td.skip_if_windows
    @pytest.mark.skipif(WASM, reason="tzset is not available on WASM")
    def test_timestamp(self, fixed_now_ts):
        # GH#17329
        # tz-naive --> treat it as if it were UTC for purposes of timestamp()
        ts = fixed_now_ts
        uts = ts.replace(tzinfo=UTC)
        with tm.assert_produces_warning(Pandas4Warning, match=naive_msg):
            result = ts.timestamp()
        assert result == uts.timestamp()

        tsc = Timestamp("2014-10-11 11:00:01.12345678", tz="US/Central")
        utsc = tsc.tz_convert("UTC")

        # utsc is a different representation of the same time
        assert tsc.timestamp() == utsc.timestamp()

        # datetime.timestamp() converts in the local timezone
        with tm.set_timezone("UTC"):
            # should agree with datetime.timestamp method
            dt = ts.to_pydatetime()
            with tm.assert_produces_warning(Pandas4Warning, match=naive_msg):
                result = ts.timestamp()
            assert dt.timestamp() == result

    def test_timestamp_naive_deprecated(self):
        # GH#50298
        ts = Timestamp("2020-03-14T15:32:52.192548")
        with tm.assert_produces_warning(Pandas4Warning, match=naive_msg):
            result = ts.timestamp()

        # the suggested alternative preserves the value
        assert result == ts.tz_localize("UTC").timestamp()

    @pytest.mark.parametrize("tz", ["UTC", "US/Central", "+01:00"])
    def test_timestamp_aware_not_deprecated(self, tz):
        # GH#50298
        ts = Timestamp("2020-03-14T15:32:52.192548", tz=tz)
        with tm.assert_produces_warning(None):
            ts.timestamp()

    @td.skip_if_windows
    @pytest.mark.skipif(WASM, reason="tzset is not available on WASM")
    def test_timestamp_naive_diverges_from_stdlib(self):
        # GH#50298 stdlib reads a naive value as local time, so
        #  Timestamp.fromtimestamp(x).timestamp() != x
        with tm.set_timezone("US/Central"):
            epoch = 1584199972.192548
            ts = Timestamp.fromtimestamp(epoch)
            assert ts.to_pydatetime() == datetime.fromtimestamp(epoch)

            with tm.assert_produces_warning(Pandas4Warning, match=naive_msg):
                result = ts.timestamp()
            assert result != epoch
            assert result == ts.tz_localize("UTC").timestamp()
            assert ts.to_pydatetime().timestamp() == epoch
