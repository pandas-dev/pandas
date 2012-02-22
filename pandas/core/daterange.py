# pylint: disable=E1101,E1103

from pandas.core.index import DatetimeIndex
import pandas.core.datetools as datetools

__all__ = ['DateRange']

#-----------------------------------------------------------------------------
# DateRange class

class DateRange(DatetimeIndex):
    def __new__(cls, start=None, end=None, periods=None,
                offset=datetools.bday, time_rule=None,
                tzinfo=None, name=None, **kwds):

        import warnings
        warnings.warn("DateRange is deprecated, use DatetimeIndex instead",
                       FutureWarning)

        # use old mapping
        if time_rule is not None:
            offset = datetools._offsetMap[time_rule]

        return super(DateRange, cls).__new__(cls, start=start, end=end,
                periods=periods, offset=offset, tzinfo=tzinfo, name=name,
                _deprecated=True, **kwds)
