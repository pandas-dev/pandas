from datetime import timedelta
import numpy as np
from pandas.core.groupby import BinGrouper, Grouper
from pandas.tseries.frequencies import to_offset, is_subperiod, is_superperiod
from pandas.tseries.index import DatetimeIndex, date_range
from pandas.tseries.tdi import TimedeltaIndex
from pandas.tseries.offsets import DateOffset, Tick, Day, _delta_to_nanoseconds
from pandas.tseries.period import PeriodIndex, period_range
import pandas.core.common as com
import pandas.compat as compat

from pandas.lib import Timestamp
import pandas.lib as lib
import pandas.tslib as tslib


_DEFAULT_METHOD = 'mean'


class TimeGrouper(Grouper):
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
                 convention=None, base=0, **kwargs):
        freq = to_offset(freq)

        end_types = set(['M', 'A', 'Q', 'BM', 'BA', 'BQ', 'W'])
        rule = freq.rule_code
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

        self.loffset = loffset
        self.how = how
        self.fill_method = fill_method
        self.limit = limit
        self.base = base

        # always sort time groupers
        kwargs['sort'] = True

        super(TimeGrouper, self).__init__(freq=freq, axis=axis, **kwargs)

    def resample(self, obj):
        self._set_grouper(obj, sort=True)
        ax = self.grouper

        if isinstance(ax, DatetimeIndex):
            rs = self._resample_timestamps()
        elif isinstance(ax, PeriodIndex):
            offset = to_offset(self.freq)
            if offset.n > 1:
                if self.kind == 'period':  # pragma: no cover
                    print('Warning: multiple of frequency -> timestamps')
                # Cannot have multiple of periods, convert to timestamp
                self.kind = 'timestamp'

            if self.kind is None or self.kind == 'period':
                rs = self._resample_periods()
            else:
                obj = self.obj.to_timestamp(how=self.convention)
                self._set_grouper(obj)
                rs = self._resample_timestamps()
        elif isinstance(ax, TimedeltaIndex):
            rs = self._resample_timestamps(kind='timedelta')
        elif len(ax) == 0:
            return self.obj
        else:  # pragma: no cover
            raise TypeError('Only valid with DatetimeIndex, TimedeltaIndex or PeriodIndex')

        rs_axis = rs._get_axis(self.axis)
        rs_axis.name = ax.name
        return rs

    def _get_grouper(self, obj):
        self._set_grouper(obj)
        return self._get_binner_for_resample()

    def _get_binner_for_resample(self, kind=None):
        # create the BinGrouper
        # assume that self.set_grouper(obj) has already been called

        ax = self.ax
        if kind is None:
            kind = self.kind
        if kind is None or kind == 'timestamp':
            self.binner, bins, binlabels = self._get_time_bins(ax)
        elif kind == 'timedelta':
            self.binner, bins, binlabels = self._get_time_delta_bins(ax)
        else:
            self.binner, bins, binlabels = self._get_time_period_bins(ax)

        self.grouper = BinGrouper(bins, binlabels)
        return self.binner, self.grouper, self.obj

    def _get_binner_for_grouping(self, obj):
        # return an ordering of the transformed group labels,
        # suitable for multi-grouping, e.g the labels for
        # the resampled intervals
        ax = self._set_grouper(obj)
        self._get_binner_for_resample()

        # create the grouper
        binner = self.binner
        l = []
        for key, group in self.grouper.get_iterator(ax):
            l.extend([key]*len(group))
        grouper = binner.__class__(l,freq=binner.freq,name=binner.name)

        # since we may have had to sort
        # may need to reorder groups here
        if self.indexer is not None:
            indexer = self.indexer.argsort(kind='quicksort')
            grouper = grouper.take(indexer)
        return grouper

    def _get_time_bins(self, ax):
        if not isinstance(ax, DatetimeIndex):
            raise TypeError('axis must be a DatetimeIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        if len(ax) == 0:
            binner = labels = DatetimeIndex(data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        first, last = ax.min(), ax.max()
        first, last = _get_range_edges(first, last, self.freq, closed=self.closed,
                                       base=self.base)
        tz = ax.tz
        binner = labels = DatetimeIndex(freq=self.freq,
                                        start=first.replace(tzinfo=None),
                                        end=last.replace(tzinfo=None),
                                        tz=tz,
                                        name=ax.name)

        # a little hack
        trimmed = False
        if (len(binner) > 2 and binner[-2] == last and
                self.closed == 'right'):

            binner = binner[:-1]
            trimmed = True

        ax_values = ax.asi8
        binner, bin_edges = self._adjust_bin_edges(binner, ax_values)

        # general version, knowing nothing about relative frequencies
        bins = lib.generate_bins_dt64(ax_values, bin_edges, self.closed, hasnans=ax.hasnans)

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

        if ax.hasnans:
            binner = binner.insert(0, tslib.NaT)
            labels = labels.insert(0, tslib.NaT)

        # if we end up with more labels than bins
        # adjust the labels
        # GH4076
        if len(bins) < len(labels):
            labels = labels[:len(bins)]

        return binner, bins, labels

    def _adjust_bin_edges(self, binner, ax_values):
        # Some hacks for > daily data, see #1471, #1458, #1483

        bin_edges = binner.asi8

        if self.freq != 'D' and is_superperiod(self.freq, 'D'):
            day_nanos = _delta_to_nanoseconds(timedelta(1))
            if self.closed == 'right':
                bin_edges = bin_edges + day_nanos - 1

            # intraday values on last day
            if bin_edges[-2] > ax_values.max():
                bin_edges = bin_edges[:-1]
                binner = binner[:-1]

        return binner, bin_edges

    def _get_time_delta_bins(self, ax):
        if not isinstance(ax, TimedeltaIndex):
            raise TypeError('axis must be a TimedeltaIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        if not len(ax):
            binner = labels = TimedeltaIndex(data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        labels = binner = TimedeltaIndex(start=ax[0],
                                         end=ax[-1],
                                         freq=self.freq,
                                         name=ax.name)

        end_stamps = labels + 1
        bins = ax.searchsorted(end_stamps, side='left')

        # Addresses GH #10530
        if self.base > 0:
            labels += type(self.freq)(self.base)

        return binner, bins, labels

    def _get_time_period_bins(self, ax):
        if not isinstance(ax, DatetimeIndex):
            raise TypeError('axis must be a DatetimeIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        if not len(ax):
            binner = labels = PeriodIndex(data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        labels = binner = PeriodIndex(start=ax[0],
                                      end=ax[-1],
                                      freq=self.freq,
                                      name=ax.name)

        end_stamps = (labels + 1).asfreq(self.freq, 's').to_timestamp()
        if ax.tzinfo:
            end_stamps = end_stamps.tz_localize(ax.tzinfo)
        bins = ax.searchsorted(end_stamps, side='left')

        return binner, bins, labels

    @property
    def _agg_method(self):
        return self.how if self.how else _DEFAULT_METHOD

    def _resample_timestamps(self, kind=None):
        # assumes set_grouper(obj) already called
        axlabels = self.ax

        self._get_binner_for_resample(kind=kind)
        grouper = self.grouper
        binner = self.binner
        obj = self.obj

        # Determine if we're downsampling
        if axlabels.freq is not None or axlabels.inferred_freq is not None:

            if len(grouper.binlabels) < len(axlabels) or self.how is not None:
                # downsample
                grouped = obj.groupby(grouper, axis=self.axis)
                result = grouped.aggregate(self._agg_method)
                # GH2073
                if self.fill_method is not None:
                    result = result.fillna(method=self.fill_method,
                                           limit=self.limit)

            else:
                # upsampling shortcut
                if self.axis:
                    raise AssertionError('axis must be 0')

                if self.closed == 'right':
                    res_index = binner[1:]
                else:
                    res_index = binner[:-1]

                # if we have the same frequency as our axis, then we are equal sampling
                # even if how is None
                if self.fill_method is None and self.limit is None and to_offset(
                    axlabels.inferred_freq) == self.freq:
                    result = obj.copy()
                    result.index = res_index
                else:
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
        if isinstance(loffset, compat.string_types):
            loffset = to_offset(self.loffset)

        if isinstance(loffset, (DateOffset, timedelta)):
            if (isinstance(result.index, DatetimeIndex)
                    and len(result.index) > 0):

                result.index = result.index + loffset

        return result

    def _resample_periods(self):
        # assumes set_grouper(obj) already called
        axlabels = self.ax
        obj = self.obj

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

    if isinstance(obj, Series):
        new_values = com.take_1d(obj.values, indexer)
        return Series(new_values, index=new_index, name=obj.name)
    elif isinstance(obj, DataFrame):
        if axis == 1:
            raise NotImplementedError("axis 1 is not supported")
        return DataFrame(obj._data.reindex_indexer(
            new_axis=new_index, indexer=indexer, axis=1))
    else:
        raise ValueError("'obj' should be either a Series or a DataFrame")


def _get_range_edges(first, last, offset, closed='left', base=0):
    if isinstance(offset, compat.string_types):
        offset = to_offset(offset)

    if isinstance(offset, Tick):
        is_day = isinstance(offset, Day)
        day_nanos = _delta_to_nanoseconds(timedelta(1))

        # #1165
        if (is_day and day_nanos % offset.nanos == 0) or not is_day:
            return _adjust_dates_anchored(first, last, offset,
                                          closed=closed, base=base)

    if not isinstance(offset, Tick):  # and first.time() != last.time():
        # hack!
        first = first.normalize()
        last = last.normalize()

    if closed == 'left':
        first = Timestamp(offset.rollback(first))
    else:
        first = Timestamp(first - offset)

    last = Timestamp(last + offset)

    return first, last


def _adjust_dates_anchored(first, last, offset, closed='right', base=0):
#     from pandas.tseries.tools import normalize_date

    # First and last offsets should be calculated from the start day to fix an
    # error cause by resampling across multiple days when a one day period is
    # not a multiple of the frequency.
    #
    # See https://github.com/pydata/pandas/issues/8683

    first_tzinfo = first.tzinfo
    first = first.tz_localize(None)
    last = last.tz_localize(None)
    start_day_nanos = first.normalize().value

    base_nanos = (base % offset.n) * offset.nanos // offset.n
    start_day_nanos += base_nanos

    foffset = (first.value - start_day_nanos) % offset.nanos
    loffset = (last.value - start_day_nanos) % offset.nanos

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

#     return (Timestamp(fresult, tz=first.tz),
#             Timestamp(lresult, tz=last.tz))

    return (Timestamp(fresult).tz_localize(first_tzinfo),
            Timestamp(lresult).tz_localize(first_tzinfo))


def asfreq(obj, freq, method=None, how=None, normalize=False):
    """
    Utility frequency conversion method for Series/DataFrame
    """
    if isinstance(obj.index, PeriodIndex):
        if method is not None:
            raise NotImplementedError("'method' argument is not supported")

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
        dti.name = obj.index.name
        rs = obj.reindex(dti, method=method)
        if normalize:
            rs.index = rs.index.normalize()
        return rs
