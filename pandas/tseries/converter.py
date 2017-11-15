# flake8: noqa

from pandas.plotting._converter import (time2num,
                                        TimeConverter, TimeFormatter,
                                        PeriodConverter, get_datevalue,
                                        DatetimeConverter,
                                        PandasAutoDateFormatter,
                                        PandasAutoDateLocator,
                                        MilliSecondLocator, get_finder,
                                        TimeSeries_DateLocator,
                                        TimeSeries_DateFormatter)


def register():
    import warnings

    msg = ("'pandas.tseries.converter' has been moved to 'pandas.plotting'. "
           "Update your import to 'from pandas.plotting import converter"')
    warnings.warn(msg, FutureWarning, stacklevel=2)
    from pandas.plotting import _converter
    return _converter.register()
