# pylint: disable=E1101,E1103

from pandas.core.index import DatetimeIndex, Index
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
        elif 'timeRule' in kwds and kwds['timeRule'] is not None:
            offset = datetools._offsetMap[kwds['timeRule']]

        return super(DateRange, cls).__new__(cls, start=start, end=end,
                periods=periods, offset=offset, tzinfo=tzinfo, name=name,
                _deprecated=True, **kwds)


def date_range(start=None, end=None, periods=None, freq='D', tz=None):
    """
    Return a fixed frequency datetime index, with day (calendar) as the default
    frequency


    Parameters
    ----------
    start :
    end :

    Returns
    -------

    """
    return DatetimeIndex(start=start, end=end, periods=periods,
                         freq=freq, tz=tz)


def bdate_range(start=None, end=None, periods=None, freq='B', tz=None):
    """
    Return a fixed frequency datetime index, with business day as the default
    frequency

    Parameters
    ----------


    Returns
    -------
    date_range : DatetimeIndex

    """

    return DatetimeIndex(start=start, end=end, periods=periods,
                         freq=freq, tz=tz)
