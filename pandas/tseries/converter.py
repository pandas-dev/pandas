# flake8: noqa
import warnings

from pandas.plotting._converter import (
    DatetimeConverter, MilliSecondLocator, PandasAutoDateFormatter,
    PandasAutoDateLocator, PeriodConverter, TimeConverter, TimeFormatter,
    TimeSeries_DateFormatter, TimeSeries_DateLocator, get_datevalue,
    get_finder, time2num)


def register():
    from pandas.plotting._converter import register as register_
    msg = ("'pandas.tseries.converter.register' has been moved and renamed to "
           "'pandas.plotting.register_matplotlib_converters'. ")
    warnings.warn(msg, FutureWarning, stacklevel=2)
    register_()
