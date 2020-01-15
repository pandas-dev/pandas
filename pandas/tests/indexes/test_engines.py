import re

import pytest

import pandas as pd


class TestDatetimeEngine:
    @pytest.mark.parametrize(
        "scalar",
        [
            pd.Timedelta(pd.Timestamp("2016-01-01").asm8.view("m8[ns]")),
            pd.Timestamp("2016-01-01").value,
            pd.Timestamp("2016-01-01").to_pydatetime(),
            pd.Timestamp("2016-01-01").to_datetime64(),
        ],
    )
    def test_not_contains_requires_timestamp(self, scalar):
        dti1 = pd.date_range("2016-01-01", periods=3)
        dti2 = dti1.insert(1, pd.NaT)  # non-monotonic
        dti3 = dti1.insert(3, dti1[0])  # non-unique
        dti4 = pd.date_range("2016-01-01", freq="ns", periods=2_000_000)
        dti5 = dti4.insert(0, dti4[0])  # over size threshold, not unique

        msg = "|".join([re.escape(str(scalar)), re.escape(repr(scalar))])
        for dti in [dti1, dti2, dti3, dti4, dti5]:
            with pytest.raises(TypeError, match=msg):
                scalar in dti._engine

            with pytest.raises(KeyError, match=msg):
                dti._engine.get_loc(scalar)


class TestTimedeltaEngine:
    @pytest.mark.parametrize(
        "scalar",
        [
            pd.Timestamp(pd.Timedelta(days=42).asm8.view("datetime64[ns]")),
            pd.Timedelta(days=42).value,
            pd.Timedelta(days=42).to_pytimedelta(),
            pd.Timedelta(days=42).to_timedelta64(),
        ],
    )
    def test_not_contains_requires_timestamp(self, scalar):
        tdi1 = pd.timedelta_range("42 days", freq="9h", periods=1234)
        tdi2 = tdi1.insert(1, pd.NaT)  # non-monotonic
        tdi3 = tdi1.insert(3, tdi1[0])  # non-unique
        tdi4 = pd.timedelta_range("42 days", freq="ns", periods=2_000_000)
        tdi5 = tdi4.insert(0, tdi4[0])  # over size threshold, not unique

        msg = "|".join([re.escape(str(scalar)), re.escape(repr(scalar))])
        for tdi in [tdi1, tdi2, tdi3, tdi4, tdi5]:
            with pytest.raises(TypeError, match=msg):
                scalar in tdi._engine

            with pytest.raises(KeyError, match=msg):
                tdi._engine.get_loc(scalar)
