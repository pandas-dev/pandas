# pylint: disable=E1101,E1103

from pandas.core.index import Index
from pandas.tseries.index import DatetimeIndex
import pandas.core.datetools as datetools


#-----------------------------------------------------------------------------
# DateRange class

class DateRange(Index):

    """Deprecated
    """

    offset = tzinfo = None

    def __new__(cls, start=None, end=None, periods=None,
                offset=datetools.bday, time_rule=None,
                tzinfo=None, name=None, **kwds):

        import warnings
        warnings.warn("DateRange is deprecated, use DatetimeIndex instead",
                      FutureWarning)

        if time_rule is None:
            time_rule = kwds.get('timeRule')
        if time_rule is not None:
            offset = datetools.get_offset(time_rule)

        return DatetimeIndex(start=start, end=end,
                             periods=periods, freq=offset,
                             tzinfo=tzinfo, name=name, **kwds)

    def __setstate__(self, aug_state):
        """Necessary for making this object picklable"""
        index_state = aug_state[:1]
        offset = aug_state[1]

        # for backwards compatibility
        if len(aug_state) > 2:
            tzinfo = aug_state[2]
        else:  # pragma: no cover
            tzinfo = None

        self.offset = offset
        self.tzinfo = tzinfo
        Index.__setstate__(self, *index_state)
