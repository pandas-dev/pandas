import pandas as pd


def test_is_monotonic_with_nat():
    # GH#31437
    # PeriodIndex.is_monotonic_increasing should behave analogously to DatetimeIndex,
    #  in particular never be monotonic when we have NaT
    dti = pd.date_range("2016-01-01", periods=3)
    pi = dti.to_period("D")
    tdi = pd.Index(dti.view("timedelta64[ns]"))

    for obj in [pi, pi._engine, dti, dti._engine, tdi, tdi._engine]:
        if isinstance(obj, pd.Index):
            # i.e. not Engines
            assert obj.is_monotonic_increasing
        assert obj.is_monotonic_increasing
        assert not obj.is_monotonic_decreasing
        assert obj.is_unique

    dti1 = dti.insert(0, pd.NaT)
    pi1 = dti1.to_period("D")
    tdi1 = pd.Index(dti1.view("timedelta64[ns]"))

    for obj in [pi1, pi1._engine, dti1, dti1._engine, tdi1, tdi1._engine]:
        if isinstance(obj, pd.Index):
            # i.e. not Engines
            assert not obj.is_monotonic_increasing
        assert not obj.is_monotonic_increasing
        assert not obj.is_monotonic_decreasing
        assert obj.is_unique

    dti2 = dti.insert(3, pd.NaT)
    pi2 = dti2.to_period("h")
    tdi2 = pd.Index(dti2.view("timedelta64[ns]"))

    for obj in [pi2, pi2._engine, dti2, dti2._engine, tdi2, tdi2._engine]:
        if isinstance(obj, pd.Index):
            # i.e. not Engines
            assert not obj.is_monotonic_increasing
        assert not obj.is_monotonic_increasing
        assert not obj.is_monotonic_decreasing
        assert obj.is_unique
