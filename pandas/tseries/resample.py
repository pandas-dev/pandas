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

    def __init__(self, offset='Min', closed='left', label='left',
                 begin=None, end=None, nperiods=None, axis=None,
                 kind=None):
        self.offset = offset
        self.closed = closed
        self.label = label
        self.begin = begin
        self.end = end
        self.nperiods = None
        self.kind = kind

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
            # TODO: Should we be a bit more careful here?
            self.bins = []
            self.binlabels = []
            return

        if isinstance(self.axis, DatetimeIndex):
            self.bins, self.binlabels = self._group_timestamps()
        elif isinstance(self.axis, PeriodIndex):
            self.bins, self.binlabels = self._group_periods()
        else:
            raise ValueError('Invalid index: %s' % type(self.axis))

    def _group_timestamps(self):
        if self.kind is None or self.kind == 'timestamp':
            binner = _generate_time_binner(self.axis, self.offset,
                                           self.begin, self.end,
                                           self.nperiods)

            int_axis = self.axis.asi8
            int_binner = com._ensure_int64(binner)

            # general version, knowing nothing about relative frequencies
            bins, labels = lib.generate_bins_dt64(int_axis, int_binner,
                                                  self.closed, self.label)
            return bins, labels.view('M8[us]')
        elif self.kind == 'period':
            index = PeriodIndex(start=self.axis[0], end=self.axis[-1],
                                freq=self.offset)

            end_stamps = (index + 1).asfreq('D', 's').to_timestamp()
            bins = self.axis.searchsorted(end_stamps, side='left')

            return bins, index

    def _group_periods(self):
        raise NotImplementedError

    def _generate_time_binner(self):
        offset = self.offset
        if isinstance(offset, basestring):
            offset = to_offset(offset)

        if not isinstance(offset, DateOffset):
            raise ValueError("Rule not a recognized offset")

        if self.begin is None:
            first = Timestamp(self.axis[0] - self.offset)
        else:
            first = Timestamp(self.offset.rollback(self.begin))

        if self.end is None:
            last = Timestamp(self.axis[-1] + self.offset)
        else:
            last = Timestamp(self.offset.rollforward(self.end))

        if isinstance(offset, Tick):
            return np.arange(first.value, last.value + 1,
                             self.offset.us_stride(), dtype=np.int64)

        result = DatetimeIndex(freq=offset, start=first, end=last,
                               periods=self.nperiods)
        return result.asi8


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
    # if isinstance(offset, basestring):
    #     offset = to_offset(offset)

    first = dtindex[0]
    last = dtindex[-1]
    # if isinstance(offset, Tick):
    #     return np.arange(first.value, last.value+1, offset.us_stride(),
    #                      dtype=np.int64)

    return PeriodIndex(freq=offset, start=first, end=last, periods=nperiods)

