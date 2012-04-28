from datetime import timedelta

import numpy as np

from pandas.core.groupby import BinGrouper, CustomGrouper
from pandas.tseries.frequencies import to_offset, is_subperiod, is_superperiod
from pandas.tseries.index import DatetimeIndex, date_range
from pandas.tseries.offsets import DateOffset
from pandas.tseries.period import PeriodIndex, period_range
from pandas.util.decorators import cache_readonly
import pandas.core.common as com

from pandas._tseries import Timestamp
import pandas._tseries as lib


class TimeGrouper(CustomGrouper):
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
    convention : {'start', 'end', 'e', 's'}
        If axis is PeriodIndex

    Notes
    -----
    Use begin, end, nperiods to generate intervals that cannot be derived
    directly from the associated object
    """
    def __init__(self, freq='Min', closed='right', label='right', how='mean',
                 begin=None, end=None, nperiods=None, axis=0,
                 fill_method=None, limit=None, loffset=None, kind=None,
                 convention=None):
        self.freq = freq
        self.closed = closed
        self.label = label
        self.begin = begin
        self.end = end
        self.nperiods = nperiods
        self.kind = kind
        self.convention = convention or 'E'
        self.axis = axis
        self.loffset = loffset
        self.how = how
        self.fill_method = fill_method
        self.limit = limit

        if axis != 0:
            raise NotImplementedError

    def resample(self, obj):
        axis = obj._get_axis(self.axis)
        if isinstance(axis, DatetimeIndex):
            return self._resample_timestamps(obj)
        elif isinstance(axis, PeriodIndex):
            return self._resample_periods(obj)
        else:
            raise TypeError('Only valid with DatetimeIndex or PeriodIndex')

    def get_grouper(self, obj):
        # Only return grouper
        return self._get_time_grouper(obj)[1]

    def _get_time_grouper(self, obj):
        axis = obj._get_axis(self.axis)

        if self.kind is None or self.kind == 'timestamp':
            binner, bins, binlabels = self._get_time_bins(axis)
            grouper = BinGrouper(bins, binlabels)
        else:
            index = binner = PeriodIndex(start=axis[0], end=axis[-1],
                                         freq=self.freq)

            end_stamps = (index + 1).asfreq('D', 's').to_timestamp()
            bins = axis.searchsorted(end_stamps, side='left')

            grouper = BinGrouper(bins, index)

        return binner, grouper

    def _get_time_bins(self, axis):
        return _make_time_bins(axis, self.freq, begin=self.begin,
                               end=self.end, nperiods=self.nperiods,
                               closed=self.closed, label=self.label)

    def _resample_timestamps(self, obj):
        axis = obj._get_axis(self.axis)

        binner, grouper = self._get_time_grouper(obj)

        # downsamples
        if len(grouper.binlabels) < len(axis):
            grouped  = obj.groupby(grouper, axis=axis)
            result = grouped.agg(self.how)
        else:
            assert(axis == 0)
            # upsampling

            # this is sort of a hack
            result = obj.reindex(binner[1:], method=self.fill_method)

        if isinstance(self.loffset, (DateOffset, timedelta)):
            if (isinstance(result.index, DatetimeIndex)
                and len(result.index) > 0):

                result.index = result.index + self.loffset

        return result

    def _resample_periods(self, obj):
        axis = obj._get_axis(self.axis)

        # Start vs. end of period
        memb = axis.asfreq(self.freq, how=self.convention)

        if is_subperiod(self.axis.freq, self.freq):
            # Downsampling
            if len(memb) > 1:
                rng = np.arange(memb.values[0], memb.values[-1])
                bins = memb.searchsorted(rng, side='right')
            else:
                bins = np.array([], dtype=np.int32)

            index = period_range(memb[0], memb[-1], freq=self.freq)
            grouper = BinGrouper(bins, index)

            grouped = obj.groupby(grouper, axis=axis)
            return grouped.agg(self.how)
        elif is_superperiod(self.axis.freq, self.freq):
            # Generate full range
            new_index = period_range(memb[0], memb[-1], freq=self.freq)

            # Get the fill indexer
            indexer = memb.get_indexer(new_index, method=self.fill_method,
                                       limit=self.limit)

            return _take_new_index(obj, indexer, new_index, axis=self.axis)
        else:
            raise ValueError('Frequency %s cannot be resampled to %s'
                             % (self.axis.freq, self.freq))


def _take_new_index(obj, indexer, new_index, axis=0):
    from pandas.core.api import Series, DataFrame
    from pandas.core.internals import BlockManager

    if isinstance(obj, Series):
        new_values = com.take_1d(obj.values, indexer)
        return Series(new_values, index=new_index, name=obj.name)
    elif isinstance(obj, DataFrame):
        if axis == 1:
            raise NotImplementedError
        data = obj._data

        new_blocks = [b.take(indexer, axis=1) for b in data.blocks]
        new_axes = list(data.axes)
        new_axes[1] = new_index
        new_data = BlockManager(new_blocks, new_axes)
        return DataFrame(new_data)
    else:
        raise NotImplementedError


def _make_period_bins(axis, freq):
    index = PeriodIndex(start=axis[0], end=axis[-1], freq=freq)
    end_stamps = (index + 1).asfreq('D', 's').to_timestamp()
    bins = axis.searchsorted(end_stamps, side='left')

    return index, bins, index


def _make_time_bins(axis, freq, begin=None, end=None, nperiods=None,
                    closed='right', label='right'):
    assert(isinstance(axis, DatetimeIndex))

    if len(axis) == 0:
        # TODO: Should we be a bit more careful here?
        return [], [], []

    if isinstance(freq, basestring):
        freq = to_offset(freq)

    if not isinstance(freq, DateOffset):
        raise ValueError("Rule not a recognized offset")

    first, last = _get_range_edges(axis, begin, end, freq, closed=closed)
    binner = DatetimeIndex(freq=freq, start=first, end=last,
                           periods=nperiods)

    # a little hack
    trimmed = False
    if len(binner) > 2 and binner[-2] == axis[-1]:
        binner = binner[:-1]
        trimmed = True

    # general version, knowing nothing about relative frequencies
    bins = lib.generate_bins_dt64(axis.asi8, binner.asi8, closed)

    if label == 'right':
        labels = binner[1:]
    elif not trimmed:
        labels = binner[:-1]
    else:
        labels = binner

    return binner, bins, labels

def _get_range_edges(axis, begin, end, offset, closed='left'):
    if begin is None:
        if closed == 'left':
            first = Timestamp(offset.rollback(axis[0]))
        else:
            first = Timestamp(axis[0] - offset)
    else:
        first = Timestamp(offset.rollback(begin))

    if end is None:
        last = Timestamp(axis[-1] + offset)
        # last = Timestamp(offset.rollforward(axis[-1]))
    else:
        last = Timestamp(offset.rollforward(end))

    return first, last

def asfreq(obj, freq, method=None, how=None):
    """
    Utility frequency conversion method for Series/DataFrame
    """
    if isinstance(obj.index, PeriodIndex):
        if method is not None:
            raise NotImplementedError

        if how is None:
            how = 'E'

        new_index = obj.index.asfreq(freq, how=how)
        new_obj = obj.copy()
        new_obj.index = new_index
        return new_obj
    else:
        if len(obj.index) == 0:
            return obj.copy()
        dti = date_range(obj.index[0], obj.index[-1], freq=freq)
        return obj.reindex(dti, method=method)

def values_at_time(obj, time, tz=None, asof=False):
    """
    Select values at particular time of day (e.g. 9:30AM)

    Parameters
    ----------
    time : datetime.time or string
    tz : string or pytz.timezone
        Time zone for time. Corresponding timestamps would be converted to
        time zone of the TimeSeries

    Returns
    -------
    values_at_time : TimeSeries
    """
    from dateutil.parser import parse

    if asof:
        raise NotImplementedError
    if tz:
        raise NotImplementedError

    if not isinstance(obj.index, DatetimeIndex):
        raise NotImplementedError

    if isinstance(time, basestring):
        time = parse(time).time()

    # TODO: time object with tzinfo?

    mus = _time_to_microsecond(time)
    indexer = lib.values_at_time(obj.index.asi8, mus)
    indexer = com._ensure_platform_int(indexer)
    return obj.take(indexer)

def _time_to_microsecond(time):
    seconds = time.hour * 60 * 60 + 60 * time.minute + time.second
    return 1000000 * seconds + time.microsecond
