# flake8: noqa
import warnings

# TODO `_matplotlib` module should be private, so the plotting backend
# can be change. Decide whether all these should be public and exponsed
# in `pandas.plotting`, or remove from here (I guess they are here for
# legacy reasons
from pandas.plotting._matplotlib.converter import (
    DatetimeConverter,
    MilliSecondLocator,
    PandasAutoDateFormatter,
    PandasAutoDateLocator,
    PeriodConverter,
    TimeConverter,
    TimeFormatter,
    TimeSeries_DateFormatter,
    TimeSeries_DateLocator,
    get_datevalue,
    get_finder,
    time2num,
)


def register():
    from pandas.plotting import register_matplotlib_converters

    msg = (
        "'pandas.tseries.converter.register' has been moved and renamed to "
        "'pandas.plotting.register_matplotlib_converters'. "
    )
    warnings.warn(msg, FutureWarning, stacklevel=2)
    register_matplotlib_converters()
