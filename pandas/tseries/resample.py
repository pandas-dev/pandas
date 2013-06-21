from datetime import timedelta

import numpy as np

from pandas.core.groupby import BinGrouper, CustomGrouper
from pandas.tseries.frequencies import to_offset, is_subperiod, is_superperiod
from pandas.tseries.index import DatetimeIndex, date_range
from pandas.tseries.offsets import DateOffset, Tick, _delta_to_nanoseconds
from pandas.tseries.period import PeriodIndex, period_range
import pandas.tseries.tools as tools
import pandas.core.common as com

from pandas.lib import Timestamp
import pandas.lib as lib


_DEFAULT_METHOD = 'mean'


class TimeGrouper(CustomGrouper):
    """
    Custom groupby class for time-interval grouping

    Parameters
    ----------
    freq : pandas date offset or offset alias for identifying bin edges
    closed : closed end of interval; left or right
    label : interval boundary to use for labeling; left or right
    nperiods : optional, integer
    convention : {'start', 'end', 'e', 's'}
        If axis is PeriodIndex

    Notes
    -----
    Use begin, end, nperiods to generate intervals that cannot be derived
    directly from the associated object
    """
    def __init__(self, freq='Min', closed=None, label=None, how='mean',
                 nperiods=None, axis=0,
                 fill_method=None, limit=None, loffset=None, kind=None,
                 convention=None, base=0):
        self.freq = to_offset(freq)

        end_types = set(['M', 'A', 'Q', 'BM', 'BA', 'BQ', 'W'])
        rule = self.freq.rule_code
        if (rule in end_types or
                ('-' in rule and rule[:rule.find('-')] in end_types)):
            if closed is None:
                closed = 'right'
            if label is None:
                label = 'right'
        else:
            if closed is None:
                closed = 'left'
            if label is None:
                label = 'left'

        self.closed = closed
        self.label = label
        self.nperiods = nperiods
        self.kind = kind

        self.convention = convention or 'E'
        self.convention = self.convention.lower()

        self.axis = axis
        self.loffset = loffset
        self.how = how
        self.fill_method = fill_method
        self.limit = limit
        self.base = base

    def resample(self, obj):
        axis = obj._get_axis(self.axis)

        if not axis.is_monotonic:
            try:
                obj = obj.sort_index(axis=self.axis)
            except TypeError:
                obj = obj.sort_index()

        if isinstance(axis, DatetimeIndex):
            rs = self._resample_timestamps(obj)
        elif isinstance(axis, PeriodIndex):
            offset = to_offset(self.freq)
            if offset.n > 1:
                if self.kind == 'period':  # pragma: no cover
                    print ('Warning: multiple of frequency -> timestamps')
                # Cannot have multiple of periods, convert to timestamp
                self.kind = 'timestamp'

            if self.kind is None or self.kind == 'period':
                rs = self._resample_periods(obj)
            else:
                obj = obj.to_timestamp(how=self.convention)
                rs = self._resample_timestamps(obj)
        elif len(axis) == 0:
            return obj
        else:  # pragma: no cover
            raise TypeError('Only valid with DatetimeIndex or PeriodIndex')

        rs_axis = rs._get_axis(self.axis)
        rs_axis.name = axis.name
        return rs

    def get_grouper(self, obj):
        # Only return grouper
        return self._get_time_grouper(obj)[1]

    def _get_time_grouper(self, obj):
        axis = obj._get_axis(self.axis)

        if self.kind is None or self.kind == 'timestamp':
            binner, bins, binlabels = self._get_time_bins(axis)
        else:
            binner, bins, binlabels = self._get_time_period_bins(axis)

        grouper = BinGrouper(bins, binlabels)
        return binner, grouper

    def _get_time_bins(self, axis):
        if not (isinstance(axis, DatetimeIndex)):
            raise AssertionError()

        if len(axis) == 0:
            binner = labels = DatetimeIndex(data=[], freq=self.freq)
            return binner, [], labels

        first, last = _get_range_edges(axis, self.freq, closed=self.closed,
                                       base=self.base)
        tz = axis.tz
        binner = labels = DatetimeIndex(freq=self.freq,
                                        start=first.replace(tzinfo=None),
                                        end=last.replace(tzinfo=None), tz=tz)

        # a little hack
        trimmed = False
        if (len(binner) > 2 and binner[-2] == axis[-1] and
                self.closed == 'right'):

            binner = binner[:-1]
            trimmed = True

        ax_values = axis.asi8
        binner, bin_edges = self._adjust_bin_edges(binner, ax_values)

        # general version, knowing nothing about relative frequencies
        bins = lib.generate_bins_dt64(ax_values, bin_edges, self.closed)

        if self.closed == 'right':
            labels = binner
            if self.label == 'right':
                labels = labels[1:]
            elif not trimmed:
                labels = labels[:-1]
        else:
            if self.label == 'right':
                labels = labels[1:]
            elif not trimmed:
                labels = labels[:-1]

        return binner, bins, labels

    def _adjust_bin_edges(self, binner, ax_values):
        # Some hacks for > daily data, see #1471, #1458, #1483

        bin_edges = binner.asi8

        if self.freq != 'D' and is_superperiod(self.freq, 'D'):
            day_nanos = _delta_to_nanoseconds(timedelta(1))
            if self.closed == 'right':
                bin_edges = bin_edges + day_nanos - 1

            # intraday values on last day
            if bin_edges[-2] > ax_values[-1]:
                bin_edges = bin_edges[:-1]
                binner = binner[:-1]

        return binner, bin_edges

    def _get_time_period_bins(self, axis):
        if not(isinstance(axis, DatetimeIndex)):
            raise AssertionError()

        if len(axis) == 0:
            binner = labels = PeriodIndex(data=[], freq=self.freq)
            return binner, [], labels

        labels = binner = PeriodIndex(start=axis[0], end=axis[-1],
                                      freq=self.freq)

        end_stamps = (labels + 1).asfreq('D', 's').to_timestamp()
        bins = axis.searchsorted(end_stamps, side='left')

        return binner, bins, labels

    @property
    def _agg_method(self):
        return self.how if self.how else _DEFAULT_METHOD

    def _resample_timestamps(self, obj):
        axlabels = obj._get_axis(self.axis)

        binner, grouper = self._get_time_grouper(obj)

        # Determine if we're downsampling
        if axlabels.freq is not None or axlabels.inferred_freq is not None:
            if len(grouper.binlabels) < len(axlabels) or self.how is not None:
                grouped = obj.groupby(grouper, axis=self.axis)
                result = grouped.aggregate(self._agg_method)
            else:
                # upsampling shortcut
                if not (self.axis == 0):
                    raise AssertionError()

                if self.closed == 'right':
                    res_index = binner[1:]
                else:
                    res_index = binner[:-1]

                result = obj.reindex(res_index, method=self.fill_method,
                                     limit=self.limit)
        else:
            # Irregular data, have to use groupby
            grouped = obj.groupby(grouper, axis=self.axis)
            result = grouped.aggregate(self._agg_method)

            if self.fill_method is not None:
                result = result.fillna(method=self.fill_method,
                                       limit=self.limit)

        loffset = self.loffset
        if isinstance(loffset, basestring):
            loffset = to_offset(self.loffset)

        if isinstance(loffset, (DateOffset, timedelta)):
            if (isinstance(result.index, DatetimeIndex)
                    and len(result.index) > 0):

                result.index = result.index + loffset

        return result

    def _resample_periods(self, obj):
        axlabels = obj._get_axis(self.axis)

        if len(axlabels) == 0:
            new_index = PeriodIndex(data=[], freq=self.freq)
            return obj.reindex(new_index)
        else:
            start = axlabels[0].asfreq(self.freq, how=self.convention)
            end = axlabels[-1].asfreq(self.freq, how='end')

            new_index = period_range(start, end, freq=self.freq)

        # Start vs. end of period
        memb = axlabels.asfreq(self.freq, how=self.convention)

        if is_subperiod(axlabels.freq, self.freq) or self.how is not None:
            # Downsampling
            rng = np.arange(memb.values[0], memb.values[-1] + 1)
            bins = memb.searchsorted(rng, side='right')
            grouper = BinGrouper(bins, new_index)

            grouped = obj.groupby(grouper, axis=self.axis)
            return grouped.aggregate(self._agg_method)
        elif is_superperiod(axlabels.freq, self.freq):
            # Get the fill indexer
            indexer = memb.get_indexer(new_index, method=self.fill_method,
                                       limit=self.limit)

            return _take_new_index(obj, indexer, new_index, axis=self.axis)
        else:
            raise ValueError('Frequency %s cannot be resampled to %s'
                             % (axlabels.freq, self.freq))


def _take_new_index(obj, indexer, new_index, axis=0):
    from pandas.core.api import Series, DataFrame
    from pandas.core.internals import BlockManager

    if isinstance(obj, Series):
        new_values = com.take_1d(obj.values, indexer)
        return Series(new_values, index=new_index, name=obj.name)
    elif isinstance(obj, DataFrame):
        if axis == 1:
            raise NotImplementedError
        return DataFrame(obj._data.take(indexer,new_index=new_index,axis=1))
    else:
        raise NotImplementedError


def _get_range_edges(axis, offset, closed='left', base=0):
    if isinstance(offset, basestring):
        offset = to_offset(offset)

    if isinstance(offset, Tick):
        day_nanos = _delta_to_nanoseconds(timedelta(1))
        # #1165
        if (day_nanos % offset.nanos) == 0:
            return _adjust_dates_anchored(axis[0], axis[-1], offset,
                                          closed=closed, base=base)

    first, last = axis[0], axis[-1]
    if not isinstance(offset, Tick):  # and first.time() != last.time():
        # hack!
        first = tools.normalize_date(first)
        last = tools.normalize_date(last)

    if closed == 'left':
        first = Timestamp(offset.rollback(first))
    else:
        first = Timestamp(first - offset)

    last = Timestamp(last + offset)

    return first, last


def _adjust_dates_anchored(first, last, offset, closed='right', base=0):
    from pandas.tseries.tools import normalize_date

    start_day_nanos = Timestamp(normalize_date(first)).value
    last_day_nanos = Timestamp(normalize_date(last)).value

    base_nanos = (base % offset.n) * offset.nanos // offset.n
    start_day_nanos += base_nanos
    last_day_nanos += base_nanos

    foffset = (first.value - start_day_nanos) % offset.nanos
    loffset = (last.value - last_day_nanos) % offset.nanos

    if closed == 'right':
        if foffset > 0:
            # roll back
            fresult = first.value - foffset
        else:
            fresult = first.value - offset.nanos

        if loffset > 0:
            # roll forward
            lresult = last.value + (offset.nanos - loffset)
        else:
            # already the end of the road
            lresult = last.value
    else:  # closed == 'left'
        if foffset > 0:
            fresult = first.value - foffset
        else:
            # start of the road
            fresult = first.value

        if loffset > 0:
            # roll forward
            lresult = last.value + (offset.nanos - loffset)
        else:
            lresult = last.value + offset.nanos

    return (Timestamp(fresult, tz=first.tz),
            Timestamp(lresult, tz=last.tz))


def asfreq(obj, freq, method=None, how=None, normalize=False):
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
        rs = obj.reindex(dti, method=method)
        if normalize:
            rs.index = rs.index.normalize()
        return rs
