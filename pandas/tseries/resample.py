import numpy as np

from pandas.core.groupby import BinGrouper
from pandas.tseries.frequencies import to_offset
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.offsets import DateOffset, Tick
from pandas.tseries.period import PeriodIndex
from pandas.util.decorators import cache_readonly
import pandas.core.common as com

from pandas._tseries import Timestamp
import pandas._tseries as lib

class TimeGrouper(BinGrouper):
    """
    Custom groupby class for time-interval grouping

    Parameters
    ----------
    rule : pandas offset string or object for identifying bin edges
    closed : closed end of interval; left (default) or right
    label : interval boundary to use for labeling; left (default) or right
    begin : optional, timestamp-like
    end : optional, timestamp-like
    nperiods : optional, integer

    Notes
    -----
    Use begin, end, nperiods to generate intervals that cannot be derived
    directly from the associated object
    """

    axis = None
    bins = None
    binlabels = None
    begin = None
    end = None
    nperiods = None
    binner = None

    def __init__(self, offset='Min', closed='left', label='left',
                 begin=None, end=None, nperiods=None, axis=None,
                 kind=None):
        if isinstance(offset, basestring):
            offset = to_offset(offset)

        if not isinstance(offset, DateOffset):
            raise ValueError("Rule not a recognized offset")

        self.offset = offset
        self.closed = closed
        self.label = label
        self.begin = begin
        self.end = end
        self.nperiods = None

        if axis is not None:
            self.set_axis(axis)

    def set_axis(self, axis):
        """
        Injects the axisect we'll act on, which we use to initialize grouper
        """
        if id(self.axis) == id(axis):
            return

        self.axis = axis

        if len(self.axis) < 1:
            self.bins = []
            self.binlabels = []
            return

        if isinstance(self.axis, DatetimeIndex):
            self.binner = _generate_time_binner(self.axis, self.offset,
                                                self.begin, self.end,
                                                self.nperiods)

            int_axis = self.axis.asi8
            int_binner = com._ensure_int64(self.binner)

            # general version, knowing nothing about relative frequencies
            bins, labels = lib.generate_bins_dt64(int_axis, int_binner,
                                                  self.closed, self.label)

            self.bins = bins
            self.binlabels = labels.view('M8[us]')
        elif isinstance(self.axis, PeriodIndex):
            pass
        else:
            raise ValueError('Invalid index: %s' % type(self.axis))


    @property
    def names(self):
        return [self.axis.name]

    @property
    def levels(self):
        return [self.binlabels]

    @cache_readonly
    def ngroups(self):
        return len(self.binlabels)

    @cache_readonly
    def result_index(self):
        return self.binlabels


def _generate_time_binner(dtindex, offset, begin=None, end=None, nperiods=None):
    if isinstance(offset, basestring):
        offset = to_offset(offset)

    if begin is None:
        first = Timestamp(dtindex[0] - offset)
    else:
        first = Timestamp(offset.rollback(begin))

    if end is None:
        last = Timestamp(dtindex[-1] + offset)
    else:
        last = Timestamp(offset.rollforward(end))

    if isinstance(offset, Tick):
        return np.arange(first.value, last.value+1, offset.us_stride(),
                         dtype=np.int64)

    result = DatetimeIndex(freq=offset, start=first, end=last,
                           periods=nperiods)
    return result.asi8


def _generate_period_binner(dtindex, offset, begin=None, end=None,
                            nperiods=None):
    if isinstance(offset, basestring):
        offset = to_offset(offset)

    if begin is None:
        first = Timestamp(dtindex[0] - offset)
    else:
        first = Timestamp(offset.rollback(begin))

    if end is None:
        last = Timestamp(dtindex[-1] + offset)
    else:
        last = Timestamp(offset.rollforward(end))

    if isinstance(offset, Tick):
        return np.arange(first.value, last.value+1, offset.us_stride(),
                         dtype=np.int64)

    return DatetimeIndex(freq=offset, start=first, end=last, periods=nperiods)

