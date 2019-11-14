import warnings

import pandas as pd

try:
    import pandas.tseries.holiday
except ImportError:
    pass

hcal = pd.tseries.holiday.USFederalHolidayCalendar()
# These offsets currently raise a NotImplimentedError with .apply_index()
non_apply = [
    pd.offsets.Day(),
    pd.offsets.BYearEnd(),
    pd.offsets.BYearBegin(),
    pd.offsets.BQuarterEnd(),
    pd.offsets.BQuarterBegin(),
    pd.offsets.BMonthEnd(),
    pd.offsets.BMonthBegin(),
    pd.offsets.CustomBusinessDay(),
    pd.offsets.CustomBusinessDay(calendar=hcal),
    pd.offsets.CustomBusinessMonthBegin(calendar=hcal),
    pd.offsets.CustomBusinessMonthEnd(calendar=hcal),
    pd.offsets.CustomBusinessMonthEnd(calendar=hcal),
]
other_offsets = [
    pd.offsets.YearEnd(),
    pd.offsets.YearBegin(),
    pd.offsets.QuarterEnd(),
    pd.offsets.QuarterBegin(),
    pd.offsets.MonthEnd(),
    pd.offsets.MonthBegin(),
    pd.offsets.DateOffset(months=2, days=2),
    pd.offsets.BusinessDay(),
    pd.offsets.SemiMonthEnd(),
    pd.offsets.SemiMonthBegin(),
]
offsets = non_apply + other_offsets


class ApplyIndex:

    params = other_offsets
    param_names = ["offset"]

    def setup(self, offset):
        N = 10000
        self.rng = pd.date_range(start="1/1/2000", periods=N, freq="T")

    def time_apply_index(self, offset):
        offset.apply_index(self.rng)


class OffsetSeriesArithmetic:

    params = offsets
    param_names = ["offset"]

    def setup(self, offset):
        N = 1000
        rng = pd.date_range(start="1/1/2000", periods=N, freq="T")
        self.data = pd.Series(rng)

    def time_add_offset(self, offset):
        with warnings.catch_warnings(record=True):
            self.data + offset


class OffsetDatetimeIndexArithmetic:

    params = offsets
    param_names = ["offset"]

    def setup(self, offset):
        N = 1000
        self.data = pd.date_range(start="1/1/2000", periods=N, freq="T")

    def time_add_offset(self, offset):
        with warnings.catch_warnings(record=True):
            self.data + offset
