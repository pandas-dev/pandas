import copy
from datetime import timedelta
from textwrap import dedent
import warnings

import numpy as np

from pandas._libs import lib
from pandas._libs.tslibs import NaT, Timestamp
from pandas._libs.tslibs.frequencies import is_subperiod, is_superperiod
from pandas._libs.tslibs.period import IncompatibleFrequency
from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError
from pandas.util._decorators import Appender, Substitution

from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries

import pandas as pd
import pandas.core.algorithms as algos
from pandas.core.generic import _shared_docs
from pandas.core.groupby.base import GroupByMixin
from pandas.core.groupby.generic import SeriesGroupBy
from pandas.core.groupby.groupby import (
    GroupBy, _GroupBy, _pipe_template, groupby)
from pandas.core.groupby.grouper import Grouper
from pandas.core.groupby.ops import BinGrouper
from pandas.core.indexes.datetimes import DatetimeIndex, date_range
from pandas.core.indexes.period import PeriodIndex
from pandas.core.indexes.timedeltas import TimedeltaIndex, timedelta_range

from pandas.tseries.frequencies import to_offset
from pandas.tseries.offsets import DateOffset, Day, Nano, Tick

_shared_docs_kwargs = dict()


class Resampler(_GroupBy):
    """
    Class for resampling datetimelike data, a groupby-like operation.
    See aggregate, transform, and apply functions on this object.

    It's easiest to use obj.resample(...) to use Resampler.

    Parameters
    ----------
    obj : pandas object
    groupby : a TimeGrouper object
    axis : int, default 0
    kind : str or None
        'period', 'timestamp' to override default index treatement

    Returns
    -------
    a Resampler of the appropriate type

    Notes
    -----
    After resampling, see aggregate, apply, and transform functions.
    """

    # to the groupby descriptor
    _attributes = ['freq', 'axis', 'closed', 'label', 'convention',
                   'loffset', 'base', 'kind']

    def __init__(self, obj, groupby=None, axis=0, kind=None, **kwargs):
        self.groupby = groupby
        self.keys = None
        self.sort = True
        self.axis = axis
        self.kind = kind
        self.squeeze = False
        self.group_keys = True
        self.as_index = True
        self.exclusions = set()
        self.binner = None
        self.grouper = None

        if self.groupby is not None:
            self.groupby._set_grouper(self._convert_obj(obj), sort=True)

    def __unicode__(self):
        """
        Provide a nice str repr of our rolling object.
        """
        attrs = ("{k}={v}".format(k=k, v=getattr(self.groupby, k))
                 for k in self._attributes if
                 getattr(self.groupby, k, None) is not None)
        return "{klass} [{attrs}]".format(klass=self.__class__.__name__,
                                          attrs=', '.join(attrs))

    def __getattr__(self, attr):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self._attributes:
            return getattr(self.groupby, attr)
        if attr in self.obj:
            return self[attr]

        return object.__getattribute__(self, attr)

    def __iter__(self):
        """
        Resampler iterator.

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group.

        See Also
        --------
        GroupBy.__iter__
        """
        self._set_binner()
        return super(Resampler, self).__iter__()

    @property
    def obj(self):
        return self.groupby.obj

    @property
    def ax(self):
        return self.groupby.ax

    @property
    def _typ(self):
        """
        Masquerade for compat as a Series or a DataFrame.
        """
        if isinstance(self._selected_obj, pd.Series):
            return 'series'
        return 'dataframe'

    @property
    def _from_selection(self):
        """
        Is the resampling from a DataFrame column or MultiIndex level.
        """
        # upsampling and PeriodIndex resampling do not work
        # with selection, this state used to catch and raise an error
        return (self.groupby is not None and
                (self.groupby.key is not None or
                 self.groupby.level is not None))

    def _convert_obj(self, obj):
        """
        Provide any conversions for the object in order to correctly handle.

        Parameters
        ----------
        obj : the object to be resampled

        Returns
        -------
        obj : converted object
        """
        obj = obj._consolidate()
        return obj

    def _get_binner_for_time(self):
        raise AbstractMethodError(self)

    def _set_binner(self):
        """
        Setup our binners.

        Cache these as we are an immutable object
        """
        if self.binner is None:
            self.binner, self.grouper = self._get_binner()

    def _get_binner(self):
        """
        Create the BinGrouper, assume that self.set_grouper(obj)
        has already been called.
        """

        binner, bins, binlabels = self._get_binner_for_time()
        bin_grouper = BinGrouper(bins, binlabels, indexer=self.groupby.indexer)
        return binner, bin_grouper

    def _assure_grouper(self):
        """
        Make sure that we are creating our binner & grouper.
        """
        self._set_binner()

    @Substitution(klass='Resampler',
                  versionadded='.. versionadded:: 0.23.0',
                  examples="""
    >>> df = pd.DataFrame({'A': [1, 2, 3, 4]},
    ...                   index=pd.date_range('2012-08-02', periods=4))
    >>> df
                A
    2012-08-02  1
    2012-08-03  2
    2012-08-04  3
    2012-08-05  4

    To get the difference between each 2-day period's maximum and minimum
    value in one pass, you can do

    >>> df.resample('2D').pipe(lambda x: x.max() - x.min())
                A
    2012-08-02  1
    2012-08-04  1
    """)
    @Appender(_pipe_template)
    def pipe(self, func, *args, **kwargs):
        return super(Resampler, self).pipe(func, *args, **kwargs)

    _agg_see_also_doc = dedent("""
    See Also
    --------
    DataFrame.groupby.aggregate
    DataFrame.resample.transform
    DataFrame.aggregate
    """)

    _agg_examples_doc = dedent("""
    Examples
    --------
    >>> s = pd.Series([1,2,3,4,5],
                      index=pd.date_range('20130101', periods=5,freq='s'))
    2013-01-01 00:00:00    1
    2013-01-01 00:00:01    2
    2013-01-01 00:00:02    3
    2013-01-01 00:00:03    4
    2013-01-01 00:00:04    5
    Freq: S, dtype: int64

    >>> r = s.resample('2s')
    DatetimeIndexResampler [freq=<2 * Seconds>, axis=0, closed=left,
                            label=left, convention=start, base=0]

    >>> r.agg(np.sum)
    2013-01-01 00:00:00    3
    2013-01-01 00:00:02    7
    2013-01-01 00:00:04    5
    Freq: 2S, dtype: int64

    >>> r.agg(['sum','mean','max'])
                         sum  mean  max
    2013-01-01 00:00:00    3   1.5    2
    2013-01-01 00:00:02    7   3.5    4
    2013-01-01 00:00:04    5   5.0    5

    >>> r.agg({'result' : lambda x: x.mean() / x.std(),
               'total' : np.sum})
                         total    result
    2013-01-01 00:00:00      3  2.121320
    2013-01-01 00:00:02      7  4.949747
    2013-01-01 00:00:04      5       NaN
    """)

    @Substitution(see_also=_agg_see_also_doc,
                  examples=_agg_examples_doc,
                  versionadded='',
                  klass='DataFrame',
                  axis='')
    @Appender(_shared_docs['aggregate'])
    def aggregate(self, func, *args, **kwargs):

        self._set_binner()
        result, how = self._aggregate(func, *args, **kwargs)
        if result is None:
            how = func
            grouper = None
            result = self._groupby_and_aggregate(how,
                                                 grouper,
                                                 *args,
                                                 **kwargs)

        result = self._apply_loffset(result)
        return result

    agg = aggregate
    apply = aggregate

    def transform(self, arg, *args, **kwargs):
        """
        Call function producing a like-indexed Series on each group and return
        a Series with the transformed values.

        Parameters
        ----------
        arg : function
            To apply to each group. Should return a Series with the same index.

        Returns
        -------
        transformed : Series

        Examples
        --------
        >>> resampled.transform(lambda x: (x - x.mean()) / x.std())
        """
        return self._selected_obj.groupby(self.groupby).transform(
            arg, *args, **kwargs)

    def _downsample(self, f):
        raise AbstractMethodError(self)

    def _upsample(self, f, limit=None, fill_value=None):
        raise AbstractMethodError(self)

    def _gotitem(self, key, ndim, subset=None):
        """
        Sub-classes to define. Return a sliced object.

        Parameters
        ----------
        key : string / list of selections
        ndim : 1,2
            requested ndim of result
        subset : object, default None
            subset to act on
        """
        self._set_binner()
        grouper = self.grouper
        if subset is None:
            subset = self.obj
        grouped = groupby(subset, by=None, grouper=grouper, axis=self.axis)

        # try the key selection
        try:
            return grouped[key]
        except KeyError:
            return grouped

    def _groupby_and_aggregate(self, how, grouper=None, *args, **kwargs):
        """
        Re-evaluate the obj with a groupby aggregation.
        """

        if grouper is None:
            self._set_binner()
            grouper = self.grouper

        obj = self._selected_obj

        grouped = groupby(obj, by=None, grouper=grouper, axis=self.axis)

        try:
            if isinstance(obj, ABCDataFrame) and callable(how):
                # Check if the function is reducing or not.
                result = grouped._aggregate_item_by_item(how, *args, **kwargs)
            else:
                result = grouped.aggregate(how, *args, **kwargs)
        except Exception:

            # we have a non-reducing function
            # try to evaluate
            result = grouped.apply(how, *args, **kwargs)

        result = self._apply_loffset(result)
        return self._wrap_result(result)

    def _apply_loffset(self, result):
        """
        If loffset is set, offset the result index.

        This is NOT an idempotent routine, it will be applied
        exactly once to the result.

        Parameters
        ----------
        result : Series or DataFrame
            the result of resample
        """

        needs_offset = (
            isinstance(self.loffset, (DateOffset, timedelta,
                                      np.timedelta64)) and
            isinstance(result.index, DatetimeIndex) and
            len(result.index) > 0
        )

        if needs_offset:
            result.index = result.index + self.loffset

        self.loffset = None
        return result

    def _get_resampler_for_grouping(self, groupby, **kwargs):
        """
        Return the correct class for resampling with groupby.
        """
        return self._resampler_for_grouping(self, groupby=groupby, **kwargs)

    def _wrap_result(self, result):
        """
        Potentially wrap any results.
        """
        if isinstance(result, ABCSeries) and self._selection is not None:
            result.name = self._selection

        if isinstance(result, ABCSeries) and result.empty:
            obj = self.obj
            if isinstance(obj.index, PeriodIndex):
                result.index = obj.index.asfreq(self.freq)
            else:
                result.index = obj.index._shallow_copy(freq=self.freq)
            result.name = getattr(obj, 'name', None)

        return result

    def pad(self, limit=None):
        """
        Forward fill the values.

        Parameters
        ----------
        limit : integer, optional
            limit of how many values to fill

        Returns
        -------
        An upsampled Series.

        See Also
        --------
        Series.fillna
        DataFrame.fillna
        """
        return self._upsample('pad', limit=limit)
    ffill = pad

    def nearest(self, limit=None):
        """
        Resample by using the nearest value.

        When resampling data, missing values may appear (e.g., when the
        resampling frequency is higher than the original frequency).
        The `nearest` method will replace ``NaN`` values that appeared in
        the resampled data with the value from the nearest member of the
        sequence, based on the index value.
        Missing values that existed in the original data will not be modified.
        If `limit` is given, fill only this many values in each direction for
        each of the original values.

        Parameters
        ----------
        limit : int, optional
            Limit of how many values to fill.

            .. versionadded:: 0.21.0

        Returns
        -------
        Series or DataFrame
            An upsampled Series or DataFrame with ``NaN`` values filled with
            their nearest value.

        See Also
        --------
        backfill : Backward fill the new missing values in the resampled data.
        pad : Forward fill ``NaN`` values.

        Examples
        --------
        >>> s = pd.Series([1, 2],
        ...               index=pd.date_range('20180101',
        ...                                   periods=2,
        ...                                   freq='1h'))
        >>> s
        2018-01-01 00:00:00    1
        2018-01-01 01:00:00    2
        Freq: H, dtype: int64

        >>> s.resample('15min').nearest()
        2018-01-01 00:00:00    1
        2018-01-01 00:15:00    1
        2018-01-01 00:30:00    2
        2018-01-01 00:45:00    2
        2018-01-01 01:00:00    2
        Freq: 15T, dtype: int64

        Limit the number of upsampled values imputed by the nearest:

        >>> s.resample('15min').nearest(limit=1)
        2018-01-01 00:00:00    1.0
        2018-01-01 00:15:00    1.0
        2018-01-01 00:30:00    NaN
        2018-01-01 00:45:00    2.0
        2018-01-01 01:00:00    2.0
        Freq: 15T, dtype: float64
        """
        return self._upsample('nearest', limit=limit)

    def backfill(self, limit=None):
        """
        Backward fill the new missing values in the resampled data.

        In statistics, imputation is the process of replacing missing data with
        substituted values [1]_. When resampling data, missing values may
        appear (e.g., when the resampling frequency is higher than the original
        frequency). The backward fill will replace NaN values that appeared in
        the resampled data with the next value in the original sequence.
        Missing values that existed in the original data will not be modified.

        Parameters
        ----------
        limit : integer, optional
            Limit of how many values to fill.

        Returns
        -------
        Series, DataFrame
            An upsampled Series or DataFrame with backward filled NaN values.

        See Also
        --------
        bfill : Alias of backfill.
        fillna : Fill NaN values using the specified method, which can be
            'backfill'.
        nearest : Fill NaN values with nearest neighbor starting from center.
        pad : Forward fill NaN values.
        Series.fillna : Fill NaN values in the Series using the
            specified method, which can be 'backfill'.
        DataFrame.fillna : Fill NaN values in the DataFrame using the
            specified method, which can be 'backfill'.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Imputation_(statistics)

        Examples
        --------

        Resampling a Series:

        >>> s = pd.Series([1, 2, 3],
        ...               index=pd.date_range('20180101', periods=3, freq='h'))
        >>> s
        2018-01-01 00:00:00    1
        2018-01-01 01:00:00    2
        2018-01-01 02:00:00    3
        Freq: H, dtype: int64

        >>> s.resample('30min').backfill()
        2018-01-01 00:00:00    1
        2018-01-01 00:30:00    2
        2018-01-01 01:00:00    2
        2018-01-01 01:30:00    3
        2018-01-01 02:00:00    3
        Freq: 30T, dtype: int64

        >>> s.resample('15min').backfill(limit=2)
        2018-01-01 00:00:00    1.0
        2018-01-01 00:15:00    NaN
        2018-01-01 00:30:00    2.0
        2018-01-01 00:45:00    2.0
        2018-01-01 01:00:00    2.0
        2018-01-01 01:15:00    NaN
        2018-01-01 01:30:00    3.0
        2018-01-01 01:45:00    3.0
        2018-01-01 02:00:00    3.0
        Freq: 15T, dtype: float64

        Resampling a DataFrame that has missing values:

        >>> df = pd.DataFrame({'a': [2, np.nan, 6], 'b': [1, 3, 5]},
        ...                   index=pd.date_range('20180101', periods=3,
        ...                                       freq='h'))
        >>> df
                               a  b
        2018-01-01 00:00:00  2.0  1
        2018-01-01 01:00:00  NaN  3
        2018-01-01 02:00:00  6.0  5

        >>> df.resample('30min').backfill()
                               a  b
        2018-01-01 00:00:00  2.0  1
        2018-01-01 00:30:00  NaN  3
        2018-01-01 01:00:00  NaN  3
        2018-01-01 01:30:00  6.0  5
        2018-01-01 02:00:00  6.0  5

        >>> df.resample('15min').backfill(limit=2)
                               a    b
        2018-01-01 00:00:00  2.0  1.0
        2018-01-01 00:15:00  NaN  NaN
        2018-01-01 00:30:00  NaN  3.0
        2018-01-01 00:45:00  NaN  3.0
        2018-01-01 01:00:00  NaN  3.0
        2018-01-01 01:15:00  NaN  NaN
        2018-01-01 01:30:00  6.0  5.0
        2018-01-01 01:45:00  6.0  5.0
        2018-01-01 02:00:00  6.0  5.0
        """
        return self._upsample('backfill', limit=limit)
    bfill = backfill

    def fillna(self, method, limit=None):
        """
        Fill missing values introduced by upsampling.

        In statistics, imputation is the process of replacing missing data with
        substituted values [1]_. When resampling data, missing values may
        appear (e.g., when the resampling frequency is higher than the original
        frequency).

        Missing values that existed in the original data will
        not be modified.

        Parameters
        ----------
        method : {'pad', 'backfill', 'ffill', 'bfill', 'nearest'}
            Method to use for filling holes in resampled data

            * 'pad' or 'ffill': use previous valid observation to fill gap
              (forward fill).
            * 'backfill' or 'bfill': use next valid observation to fill gap.
            * 'nearest': use nearest valid observation to fill gap.

        limit : integer, optional
            Limit of how many consecutive missing values to fill.

        Returns
        -------
        Series or DataFrame
            An upsampled Series or DataFrame with missing values filled.

        See Also
        --------
        backfill : Backward fill NaN values in the resampled data.
        pad : Forward fill NaN values in the resampled data.
        nearest : Fill NaN values in the resampled data
            with nearest neighbor starting from center.
        interpolate : Fill NaN values using interpolation.
        Series.fillna : Fill NaN values in the Series using the
            specified method, which can be 'bfill' and 'ffill'.
        DataFrame.fillna : Fill NaN values in the DataFrame using the
            specified method, which can be 'bfill' and 'ffill'.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Imputation_(statistics)

        Examples
        --------
        Resampling a Series:

        >>> s = pd.Series([1, 2, 3],
        ...               index=pd.date_range('20180101', periods=3, freq='h'))
        >>> s
        2018-01-01 00:00:00    1
        2018-01-01 01:00:00    2
        2018-01-01 02:00:00    3
        Freq: H, dtype: int64

        Without filling the missing values you get:

        >>> s.resample("30min").asfreq()
        2018-01-01 00:00:00    1.0
        2018-01-01 00:30:00    NaN
        2018-01-01 01:00:00    2.0
        2018-01-01 01:30:00    NaN
        2018-01-01 02:00:00    3.0
        Freq: 30T, dtype: float64

        >>> s.resample('30min').fillna("backfill")
        2018-01-01 00:00:00    1
        2018-01-01 00:30:00    2
        2018-01-01 01:00:00    2
        2018-01-01 01:30:00    3
        2018-01-01 02:00:00    3
        Freq: 30T, dtype: int64

        >>> s.resample('15min').fillna("backfill", limit=2)
        2018-01-01 00:00:00    1.0
        2018-01-01 00:15:00    NaN
        2018-01-01 00:30:00    2.0
        2018-01-01 00:45:00    2.0
        2018-01-01 01:00:00    2.0
        2018-01-01 01:15:00    NaN
        2018-01-01 01:30:00    3.0
        2018-01-01 01:45:00    3.0
        2018-01-01 02:00:00    3.0
        Freq: 15T, dtype: float64

        >>> s.resample('30min').fillna("pad")
        2018-01-01 00:00:00    1
        2018-01-01 00:30:00    1
        2018-01-01 01:00:00    2
        2018-01-01 01:30:00    2
        2018-01-01 02:00:00    3
        Freq: 30T, dtype: int64

        >>> s.resample('30min').fillna("nearest")
        2018-01-01 00:00:00    1
        2018-01-01 00:30:00    2
        2018-01-01 01:00:00    2
        2018-01-01 01:30:00    3
        2018-01-01 02:00:00    3
        Freq: 30T, dtype: int64

        Missing values present before the upsampling are not affected.

        >>> sm = pd.Series([1, None, 3],
        ...               index=pd.date_range('20180101', periods=3, freq='h'))
        >>> sm
        2018-01-01 00:00:00    1.0
        2018-01-01 01:00:00    NaN
        2018-01-01 02:00:00    3.0
        Freq: H, dtype: float64

        >>> sm.resample('30min').fillna('backfill')
        2018-01-01 00:00:00    1.0
        2018-01-01 00:30:00    NaN
        2018-01-01 01:00:00    NaN
        2018-01-01 01:30:00    3.0
        2018-01-01 02:00:00    3.0
        Freq: 30T, dtype: float64

        >>> sm.resample('30min').fillna('pad')
        2018-01-01 00:00:00    1.0
        2018-01-01 00:30:00    1.0
        2018-01-01 01:00:00    NaN
        2018-01-01 01:30:00    NaN
        2018-01-01 02:00:00    3.0
        Freq: 30T, dtype: float64

        >>> sm.resample('30min').fillna('nearest')
        2018-01-01 00:00:00    1.0
        2018-01-01 00:30:00    NaN
        2018-01-01 01:00:00    NaN
        2018-01-01 01:30:00    3.0
        2018-01-01 02:00:00    3.0
        Freq: 30T, dtype: float64

        DataFrame resampling is done column-wise. All the same options are
        available.

        >>> df = pd.DataFrame({'a': [2, np.nan, 6], 'b': [1, 3, 5]},
        ...                   index=pd.date_range('20180101', periods=3,
        ...                                       freq='h'))
        >>> df
                               a  b
        2018-01-01 00:00:00  2.0  1
        2018-01-01 01:00:00  NaN  3
        2018-01-01 02:00:00  6.0  5

        >>> df.resample('30min').fillna("bfill")
                               a  b
        2018-01-01 00:00:00  2.0  1
        2018-01-01 00:30:00  NaN  3
        2018-01-01 01:00:00  NaN  3
        2018-01-01 01:30:00  6.0  5
        2018-01-01 02:00:00  6.0  5
        """
        return self._upsample(method, limit=limit)

    @Appender(_shared_docs['interpolate'] % _shared_docs_kwargs)
    def interpolate(self, method='linear', axis=0, limit=None, inplace=False,
                    limit_direction='forward', limit_area=None,
                    downcast=None, **kwargs):
        """
        Interpolate values according to different methods.

        .. versionadded:: 0.18.1
        """
        result = self._upsample(None)
        return result.interpolate(method=method, axis=axis, limit=limit,
                                  inplace=inplace,
                                  limit_direction=limit_direction,
                                  limit_area=limit_area,
                                  downcast=downcast, **kwargs)

    def asfreq(self, fill_value=None):
        """
        Return the values at the new freq, essentially a reindex.

        Parameters
        ----------
        fill_value : scalar, optional
            Value to use for missing values, applied during upsampling (note
            this does not fill NaNs that already were present).

            .. versionadded:: 0.20.0

        See Also
        --------
        Series.asfreq
        DataFrame.asfreq
        """
        return self._upsample('asfreq', fill_value=fill_value)

    def std(self, ddof=1, *args, **kwargs):
        """
        Compute standard deviation of groups, excluding missing values.

        Parameters
        ----------
        ddof : integer, default 1
            Degrees of freedom.
        """
        nv.validate_resampler_func('std', args, kwargs)
        return self._downsample('std', ddof=ddof)

    def var(self, ddof=1, *args, **kwargs):
        """
        Compute variance of groups, excluding missing values.

        Parameters
        ----------
        ddof : integer, default 1
            degrees of freedom
        """
        nv.validate_resampler_func('var', args, kwargs)
        return self._downsample('var', ddof=ddof)

    @Appender(GroupBy.size.__doc__)
    def size(self):
        # It's a special case as higher level does return
        # a copy of 0-len objects. GH14962
        result = self._downsample('size')
        if not len(self.ax) and isinstance(self._selected_obj, ABCDataFrame):
            result = pd.Series([], index=result.index, dtype='int64')
        return result

    def quantile(self, q=0.5, **kwargs):
        """
        Return value at the given quantile.

        .. versionadded:: 0.24.0

        Parameters
        ----------
        q : float or array-like, default 0.5 (50% quantile)

        See Also
        --------
        Series.quantile
        DataFrame.quantile
        DataFrameGroupBy.quantile
        """
        return self._downsample('quantile', q=q, **kwargs)


# downsample methods
for method in ['sum', 'prod']:

    def f(self, _method=method, min_count=0, *args, **kwargs):
        nv.validate_resampler_func(_method, args, kwargs)
        return self._downsample(_method, min_count=min_count)
    f.__doc__ = getattr(GroupBy, method).__doc__
    setattr(Resampler, method, f)


# downsample methods
for method in ['min', 'max', 'first', 'last', 'mean', 'sem',
               'median', 'ohlc']:

    def f(self, _method=method, *args, **kwargs):
        nv.validate_resampler_func(_method, args, kwargs)
        return self._downsample(_method)
    f.__doc__ = getattr(GroupBy, method).__doc__
    setattr(Resampler, method, f)

# groupby & aggregate methods
for method in ['count']:
    def f(self, _method=method):
        return self._downsample(_method)
    f.__doc__ = getattr(GroupBy, method).__doc__
    setattr(Resampler, method, f)

# series only methods
for method in ['nunique']:
    def f(self, _method=method):
        return self._downsample(_method)
    f.__doc__ = getattr(SeriesGroupBy, method).__doc__
    setattr(Resampler, method, f)


def _maybe_process_deprecations(r, how=None, fill_method=None, limit=None):
    """
    Potentially we might have a deprecation warning, show it
    but call the appropriate methods anyhow.
    """

    if how is not None:

        # .resample(..., how='sum')
        if isinstance(how, str):
            method = "{0}()".format(how)

            # .resample(..., how=lambda x: ....)
        else:
            method = ".apply(<func>)"

        # if we have both a how and fill_method, then show
        # the following warning
        if fill_method is None:
            warnings.warn("how in .resample() is deprecated\n"
                          "the new syntax is "
                          ".resample(...).{method}".format(
                              method=method),
                          FutureWarning, stacklevel=3)
        r = r.aggregate(how)

    if fill_method is not None:

        # show the prior function call
        method = '.' + method if how is not None else ''

        args = "limit={0}".format(limit) if limit is not None else ""
        warnings.warn("fill_method is deprecated to .resample()\n"
                      "the new syntax is .resample(...){method}"
                      ".{fill_method}({args})".format(
                          method=method,
                          fill_method=fill_method,
                          args=args),
                      FutureWarning, stacklevel=3)

        if how is not None:
            r = getattr(r, fill_method)(limit=limit)
        else:
            r = r.aggregate(fill_method, limit=limit)

    return r


class _GroupByMixin(GroupByMixin):
    """
    Provide the groupby facilities.
    """
    def __init__(self, obj, *args, **kwargs):

        parent = kwargs.pop('parent', None)
        groupby = kwargs.pop('groupby', None)
        if parent is None:
            parent = obj

        # initialize our GroupByMixin object with
        # the resampler attributes
        for attr in self._attributes:
            setattr(self, attr, kwargs.get(attr, getattr(parent, attr)))

        super(_GroupByMixin, self).__init__(None)
        self._groupby = groupby
        self._groupby.mutated = True
        self._groupby.grouper.mutated = True
        self.groupby = copy.copy(parent.groupby)

    def _apply(self, f, grouper=None, *args, **kwargs):
        """
        Dispatch to _upsample; we are stripping all of the _upsample kwargs and
        performing the original function call on the grouped object.
        """

        def func(x):
            x = self._shallow_copy(x, groupby=self.groupby)

            if isinstance(f, str):
                return getattr(x, f)(**kwargs)

            return x.apply(f, *args, **kwargs)

        result = self._groupby.apply(func)
        return self._wrap_result(result)

    _upsample = _apply
    _downsample = _apply
    _groupby_and_aggregate = _apply


class DatetimeIndexResampler(Resampler):

    @property
    def _resampler_for_grouping(self):
        return DatetimeIndexResamplerGroupby

    def _get_binner_for_time(self):

        # this is how we are actually creating the bins
        if self.kind == 'period':
            return self.groupby._get_time_period_bins(self.ax)
        return self.groupby._get_time_bins(self.ax)

    def _downsample(self, how, **kwargs):
        """
        Downsample the cython defined function.

        Parameters
        ----------
        how : string / cython mapped function
        **kwargs : kw args passed to how function
        """
        self._set_binner()
        how = self._is_cython_func(how) or how
        ax = self.ax
        obj = self._selected_obj

        if not len(ax):
            # reset to the new freq
            obj = obj.copy()
            obj.index.freq = self.freq
            return obj

        # do we have a regular frequency
        if ax.freq is not None or ax.inferred_freq is not None:

            if len(self.grouper.binlabels) > len(ax) and how is None:

                # let's do an asfreq
                return self.asfreq()

        # we are downsampling
        # we want to call the actual grouper method here
        result = obj.groupby(
            self.grouper, axis=self.axis).aggregate(how, **kwargs)

        result = self._apply_loffset(result)
        return self._wrap_result(result)

    def _adjust_binner_for_upsample(self, binner):
        """
        Adjust our binner when upsampling.

        The range of a new index should not be outside specified range
        """
        if self.closed == 'right':
            binner = binner[1:]
        else:
            binner = binner[:-1]
        return binner

    def _upsample(self, method, limit=None, fill_value=None):
        """
        Parameters
        ----------
        method : string {'backfill', 'bfill', 'pad',
            'ffill', 'asfreq'} method for upsampling
        limit : int, default None
            Maximum size gap to fill when reindexing
        fill_value : scalar, default None
            Value to use for missing values

        See Also
        --------
        .fillna

        """
        self._set_binner()
        if self.axis:
            raise AssertionError('axis must be 0')
        if self._from_selection:
            raise ValueError("Upsampling from level= or on= selection"
                             " is not supported, use .set_index(...)"
                             " to explicitly set index to"
                             " datetime-like")

        ax = self.ax
        obj = self._selected_obj
        binner = self.binner
        res_index = self._adjust_binner_for_upsample(binner)

        # if we have the same frequency as our axis, then we are equal sampling
        if limit is None and to_offset(ax.inferred_freq) == self.freq:
            result = obj.copy()
            result.index = res_index
        else:
            result = obj.reindex(res_index, method=method,
                                 limit=limit, fill_value=fill_value)

        result = self._apply_loffset(result)
        return self._wrap_result(result)

    def _wrap_result(self, result):
        result = super(DatetimeIndexResampler, self)._wrap_result(result)

        # we may have a different kind that we were asked originally
        # convert if needed
        if self.kind == 'period' and not isinstance(result.index, PeriodIndex):
            result.index = result.index.to_period(self.freq)
        return result


class DatetimeIndexResamplerGroupby(_GroupByMixin, DatetimeIndexResampler):
    """
    Provides a resample of a groupby implementation

    .. versionadded:: 0.18.1
    """
    @property
    def _constructor(self):
        return DatetimeIndexResampler


class PeriodIndexResampler(DatetimeIndexResampler):

    @property
    def _resampler_for_grouping(self):
        return PeriodIndexResamplerGroupby

    def _get_binner_for_time(self):
        if self.kind == 'timestamp':
            return super(PeriodIndexResampler, self)._get_binner_for_time()
        return self.groupby._get_period_bins(self.ax)

    def _convert_obj(self, obj):
        obj = super(PeriodIndexResampler, self)._convert_obj(obj)

        if self._from_selection:
            # see GH 14008, GH 12871
            msg = ("Resampling from level= or on= selection"
                   " with a PeriodIndex is not currently supported,"
                   " use .set_index(...) to explicitly set index")
            raise NotImplementedError(msg)

        if self.loffset is not None:
            # Cannot apply loffset/timedelta to PeriodIndex -> convert to
            # timestamps
            self.kind = 'timestamp'

        # convert to timestamp
        if self.kind == 'timestamp':
            obj = obj.to_timestamp(how=self.convention)

        return obj

    def _downsample(self, how, **kwargs):
        """
        Downsample the cython defined function.

        Parameters
        ----------
        how : string / cython mapped function
        **kwargs : kw args passed to how function
        """

        # we may need to actually resample as if we are timestamps
        if self.kind == 'timestamp':
            return super(PeriodIndexResampler, self)._downsample(how, **kwargs)

        how = self._is_cython_func(how) or how
        ax = self.ax

        if is_subperiod(ax.freq, self.freq):
            # Downsampling
            return self._groupby_and_aggregate(how, grouper=self.grouper,
                                               **kwargs)
        elif is_superperiod(ax.freq, self.freq):
            if how == 'ohlc':
                # GH #13083
                # upsampling to subperiods is handled as an asfreq, which works
                # for pure aggregating/reducing methods
                # OHLC reduces along the time dimension, but creates multiple
                # values for each period -> handle by _groupby_and_aggregate()
                return self._groupby_and_aggregate(how, grouper=self.grouper)
            return self.asfreq()
        elif ax.freq == self.freq:
            return self.asfreq()

        raise IncompatibleFrequency(
            'Frequency {} cannot be resampled to {}, as they are not '
            'sub or super periods'.format(ax.freq, self.freq))

    def _upsample(self, method, limit=None, fill_value=None):
        """
        Parameters
        ----------
        method : string {'backfill', 'bfill', 'pad', 'ffill'}
            method for upsampling
        limit : int, default None
            Maximum size gap to fill when reindexing
        fill_value : scalar, default None
            Value to use for missing values

        See Also
        --------
        .fillna

        """

        # we may need to actually resample as if we are timestamps
        if self.kind == 'timestamp':
            return super(PeriodIndexResampler, self)._upsample(
                method, limit=limit, fill_value=fill_value)

        self._set_binner()
        ax = self.ax
        obj = self.obj
        new_index = self.binner

        # Start vs. end of period
        memb = ax.asfreq(self.freq, how=self.convention)

        # Get the fill indexer
        indexer = memb.get_indexer(new_index, method=method, limit=limit)
        return self._wrap_result(_take_new_index(
            obj, indexer, new_index, axis=self.axis))


class PeriodIndexResamplerGroupby(_GroupByMixin, PeriodIndexResampler):
    """
    Provides a resample of a groupby implementation.

    .. versionadded:: 0.18.1
    """
    @property
    def _constructor(self):
        return PeriodIndexResampler


class TimedeltaIndexResampler(DatetimeIndexResampler):

    @property
    def _resampler_for_grouping(self):
        return TimedeltaIndexResamplerGroupby

    def _get_binner_for_time(self):
        return self.groupby._get_time_delta_bins(self.ax)

    def _adjust_binner_for_upsample(self, binner):
        """
        Adjust our binner when upsampling.

        The range of a new index is allowed to be greater than original range
        so we don't need to change the length of a binner, GH 13022
        """
        return binner


class TimedeltaIndexResamplerGroupby(_GroupByMixin, TimedeltaIndexResampler):
    """
    Provides a resample of a groupby implementation.

    .. versionadded:: 0.18.1
    """
    @property
    def _constructor(self):
        return TimedeltaIndexResampler


def resample(obj, kind=None, **kwds):
    """
    Create a TimeGrouper and return our resampler.
    """
    tg = TimeGrouper(**kwds)
    return tg._get_resampler(obj, kind=kind)


resample.__doc__ = Resampler.__doc__


def get_resampler_for_grouping(groupby, rule, how=None, fill_method=None,
                               limit=None, kind=None, **kwargs):
    """
    Return our appropriate resampler when grouping as well.
    """

    # .resample uses 'on' similar to how .groupby uses 'key'
    kwargs['key'] = kwargs.pop('on', None)

    tg = TimeGrouper(freq=rule, **kwargs)
    resampler = tg._get_resampler(groupby.obj, kind=kind)
    r = resampler._get_resampler_for_grouping(groupby=groupby)
    return _maybe_process_deprecations(r,
                                       how=how,
                                       fill_method=fill_method,
                                       limit=limit)


class TimeGrouper(Grouper):
    """
    Custom groupby class for time-interval grouping.

    Parameters
    ----------
    freq : pandas date offset or offset alias for identifying bin edges
    closed : closed end of interval; 'left' or 'right'
    label : interval boundary to use for labeling; 'left' or 'right'
    convention : {'start', 'end', 'e', 's'}
        If axis is PeriodIndex
    """
    _attributes = Grouper._attributes + ('closed', 'label', 'how',
                                         'loffset', 'kind', 'convention',
                                         'base')

    def __init__(self, freq='Min', closed=None, label=None, how='mean',
                 axis=0, fill_method=None, limit=None, loffset=None,
                 kind=None, convention=None, base=0, **kwargs):
        # Check for correctness of the keyword arguments which would
        # otherwise silently use the default if misspelled
        if label not in {None, 'left', 'right'}:
            raise ValueError('Unsupported value {} for `label`'.format(label))
        if closed not in {None, 'left', 'right'}:
            raise ValueError('Unsupported value {} for `closed`'.format(
                closed))
        if convention not in {None, 'start', 'end', 'e', 's'}:
            raise ValueError('Unsupported value {} for `convention`'
                             .format(convention))

        freq = to_offset(freq)

        end_types = {'M', 'A', 'Q', 'BM', 'BA', 'BQ', 'W'}
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
        self.kind = kind

        self.convention = convention or 'E'
        self.convention = self.convention.lower()

        if isinstance(loffset, str):
            loffset = to_offset(loffset)
        self.loffset = loffset

        self.how = how
        self.fill_method = fill_method
        self.limit = limit
        self.base = base

        # always sort time groupers
        kwargs['sort'] = True

        super(TimeGrouper, self).__init__(freq=freq, axis=axis, **kwargs)

    def _get_resampler(self, obj, kind=None):
        """
        Return my resampler or raise if we have an invalid axis.

        Parameters
        ----------
        obj : input object
        kind : string, optional
            'period','timestamp','timedelta' are valid

        Returns
        -------
        a Resampler

        Raises
        ------
        TypeError if incompatible axis

        """
        self._set_grouper(obj)

        ax = self.ax
        if isinstance(ax, DatetimeIndex):
            return DatetimeIndexResampler(obj,
                                          groupby=self,
                                          kind=kind,
                                          axis=self.axis)
        elif isinstance(ax, PeriodIndex) or kind == 'period':
            return PeriodIndexResampler(obj,
                                        groupby=self,
                                        kind=kind,
                                        axis=self.axis)
        elif isinstance(ax, TimedeltaIndex):
            return TimedeltaIndexResampler(obj,
                                           groupby=self,
                                           axis=self.axis)

        raise TypeError("Only valid with DatetimeIndex, "
                        "TimedeltaIndex or PeriodIndex, "
                        "but got an instance of %r" % type(ax).__name__)

    def _get_grouper(self, obj, validate=True):
        # create the resampler and return our binner
        r = self._get_resampler(obj)
        r._set_binner()
        return r.binner, r.grouper, r.obj

    def _get_time_bins(self, ax):
        if not isinstance(ax, DatetimeIndex):
            raise TypeError('axis must be a DatetimeIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        if len(ax) == 0:
            binner = labels = DatetimeIndex(
                data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        first, last = _get_timestamp_range_edges(ax.min(), ax.max(),
                                                 self.freq,
                                                 closed=self.closed,
                                                 base=self.base)
        # GH #12037
        # use first/last directly instead of call replace() on them
        # because replace() will swallow the nanosecond part
        # thus last bin maybe slightly before the end if the end contains
        # nanosecond part and lead to `Values falls after last bin` error
        binner = labels = date_range(freq=self.freq,
                                     start=first,
                                     end=last,
                                     tz=ax.tz,
                                     name=ax.name,
                                     ambiguous='infer',
                                     nonexistent='shift_forward')

        ax_values = ax.asi8
        binner, bin_edges = self._adjust_bin_edges(binner, ax_values)

        # general version, knowing nothing about relative frequencies
        bins = lib.generate_bins_dt64(
            ax_values, bin_edges, self.closed, hasnans=ax.hasnans)

        if self.closed == 'right':
            labels = binner
            if self.label == 'right':
                labels = labels[1:]
        elif self.label == 'right':
            labels = labels[1:]

        if ax.hasnans:
            binner = binner.insert(0, NaT)
            labels = labels.insert(0, NaT)

        # if we end up with more labels than bins
        # adjust the labels
        # GH4076
        if len(bins) < len(labels):
            labels = labels[:len(bins)]

        return binner, bins, labels

    def _adjust_bin_edges(self, binner, ax_values):
        # Some hacks for > daily data, see #1471, #1458, #1483

        if self.freq != 'D' and is_superperiod(self.freq, 'D'):
            if self.closed == 'right':
                # GH 21459, GH 9119: Adjust the bins relative to the wall time
                bin_edges = binner.tz_localize(None)
                bin_edges = bin_edges + timedelta(1) - Nano(1)
                bin_edges = bin_edges.tz_localize(binner.tz).asi8
            else:
                bin_edges = binner.asi8

            # intraday values on last day
            if bin_edges[-2] > ax_values.max():
                bin_edges = bin_edges[:-1]
                binner = binner[:-1]
        else:
            bin_edges = binner.asi8
        return binner, bin_edges

    def _get_time_delta_bins(self, ax):
        if not isinstance(ax, TimedeltaIndex):
            raise TypeError('axis must be a TimedeltaIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        if not len(ax):
            binner = labels = TimedeltaIndex(
                data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        start, end = ax.min(), ax.max()
        labels = binner = timedelta_range(start=start,
                                          end=end,
                                          freq=self.freq,
                                          name=ax.name)

        end_stamps = labels + self.freq
        bins = ax.searchsorted(end_stamps, side='left')

        # Addresses GH #10530
        if self.base > 0:
            labels += type(self.freq)(self.base)

        return binner, bins, labels

    def _get_time_period_bins(self, ax):
        if not isinstance(ax, DatetimeIndex):
            raise TypeError('axis must be a DatetimeIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        freq = self.freq

        if not len(ax):
            binner = labels = PeriodIndex(data=[], freq=freq, name=ax.name)
            return binner, [], labels

        labels = binner = pd.period_range(start=ax[0],
                                          end=ax[-1],
                                          freq=freq,
                                          name=ax.name)

        end_stamps = (labels + freq).asfreq(freq, 's').to_timestamp()
        if ax.tzinfo:
            end_stamps = end_stamps.tz_localize(ax.tzinfo)
        bins = ax.searchsorted(end_stamps, side='left')

        return binner, bins, labels

    def _get_period_bins(self, ax):
        if not isinstance(ax, PeriodIndex):
            raise TypeError('axis must be a PeriodIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        memb = ax.asfreq(self.freq, how=self.convention)

        # NaT handling as in pandas._lib.lib.generate_bins_dt64()
        nat_count = 0
        if memb.hasnans:
            nat_count = np.sum(memb._isnan)
            memb = memb[~memb._isnan]

        # if index contains no valid (non-NaT) values, return empty index
        if not len(memb):
            binner = labels = PeriodIndex(
                data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        freq_mult = self.freq.n

        start = ax.min().asfreq(self.freq, how=self.convention)
        end = ax.max().asfreq(self.freq, how='end')
        bin_shift = 0

        # GH 23882
        if self.base:
            # get base adjusted bin edge labels
            p_start, end = _get_period_range_edges(start,
                                                   end,
                                                   self.freq,
                                                   closed=self.closed,
                                                   base=self.base)

            # Get offset for bin edge (not label edge) adjustment
            start_offset = (pd.Period(start, self.freq)
                            - pd.Period(p_start, self.freq))
            bin_shift = start_offset.n % freq_mult
            start = p_start

        labels = binner = pd.period_range(start=start, end=end,
                                          freq=self.freq, name=ax.name)

        i8 = memb.asi8

        # when upsampling to subperiods, we need to generate enough bins
        expected_bins_count = len(binner) * freq_mult
        i8_extend = expected_bins_count - (i8[-1] - i8[0])
        rng = np.arange(i8[0], i8[-1] + i8_extend, freq_mult)
        rng += freq_mult
        # adjust bin edge indexes to account for base
        rng -= bin_shift
        bins = memb.searchsorted(rng, side='left')

        if nat_count > 0:
            # NaT handling as in pandas._lib.lib.generate_bins_dt64()
            # shift bins by the number of NaT
            bins += nat_count
            bins = np.insert(bins, 0, nat_count)
            binner = binner.insert(0, NaT)
            labels = labels.insert(0, NaT)

        return binner, bins, labels


def _take_new_index(obj, indexer, new_index, axis=0):
    from pandas.core.api import Series, DataFrame

    if isinstance(obj, Series):
        new_values = algos.take_1d(obj.values, indexer)
        return Series(new_values, index=new_index, name=obj.name)
    elif isinstance(obj, DataFrame):
        if axis == 1:
            raise NotImplementedError("axis 1 is not supported")
        return DataFrame(obj._data.reindex_indexer(
            new_axis=new_index, indexer=indexer, axis=1))
    else:
        raise ValueError("'obj' should be either a Series or a DataFrame")


def _get_timestamp_range_edges(first, last, offset, closed='left', base=0):
    """
    Adjust the `first` Timestamp to the preceeding Timestamp that resides on
    the provided offset. Adjust the `last` Timestamp to the following
    Timestamp that resides on the provided offset. Input Timestamps that
    already reside on the offset will be adjusted depending on the type of
    offset and the `closed` parameter.

    Parameters
    ----------
    first : pd.Timestamp
        The beginning Timestamp of the range to be adjusted.
    last : pd.Timestamp
        The ending Timestamp of the range to be adjusted.
    offset : pd.DateOffset
        The dateoffset to which the Timestamps will be adjusted.
    closed : {'right', 'left'}, default None
        Which side of bin interval is closed.
    base : int, default 0
        The "origin" of the adjusted Timestamps.

    Returns
    -------
    A tuple of length 2, containing the adjusted pd.Timestamp objects.
    """
    if isinstance(offset, Tick):
        if isinstance(offset, Day):
            # _adjust_dates_anchored assumes 'D' means 24H, but first/last
            # might contain a DST transition (23H, 24H, or 25H).
            # So "pretend" the dates are naive when adjusting the endpoints
            tz = first.tz
            first = first.tz_localize(None)
            last = last.tz_localize(None)

        first, last = _adjust_dates_anchored(first, last, offset,
                                             closed=closed, base=base)
        if isinstance(offset, Day):
            first = first.tz_localize(tz)
            last = last.tz_localize(tz)
        return first, last

    else:
        first = first.normalize()
        last = last.normalize()

    if closed == 'left':
        first = Timestamp(offset.rollback(first))
    else:
        first = Timestamp(first - offset)

    last = Timestamp(last + offset)

    return first, last


def _get_period_range_edges(first, last, offset, closed='left', base=0):
    """
    Adjust the provided `first` and `last` Periods to the respective Period of
    the given offset that encompasses them.

    Parameters
    ----------
    first : pd.Period
        The beginning Period of the range to be adjusted.
    last : pd.Period
        The ending Period of the range to be adjusted.
    offset : pd.DateOffset
        The dateoffset to which the Periods will be adjusted.
    closed : {'right', 'left'}, default None
        Which side of bin interval is closed.
    base : int, default 0
        The "origin" of the adjusted Periods.

    Returns
    -------
    A tuple of length 2, containing the adjusted pd.Period objects.
    """
    if not all(isinstance(obj, pd.Period) for obj in [first, last]):
        raise TypeError("'first' and 'last' must be instances of type Period")

    # GH 23882
    first = first.to_timestamp()
    last = last.to_timestamp()
    adjust_first = not offset.onOffset(first)
    adjust_last = offset.onOffset(last)

    first, last = _get_timestamp_range_edges(first, last, offset,
                                             closed=closed, base=base)

    first = (first + adjust_first * offset).to_period(offset)
    last = (last - adjust_last * offset).to_period(offset)
    return first, last


def _adjust_dates_anchored(first, last, offset, closed='right', base=0):
    # First and last offsets should be calculated from the start day to fix an
    # error cause by resampling across multiple days when a one day period is
    # not a multiple of the frequency.
    #
    # See https://github.com/pandas-dev/pandas/issues/8683

    # GH 10117 & GH 19375. If first and last contain timezone information,
    # Perform the calculation in UTC in order to avoid localizing on an
    # Ambiguous or Nonexistent time.
    first_tzinfo = first.tzinfo
    last_tzinfo = last.tzinfo
    start_day_nanos = first.normalize().value
    if first_tzinfo is not None:
        first = first.tz_convert('UTC')
    if last_tzinfo is not None:
        last = last.tz_convert('UTC')

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
    fresult = Timestamp(fresult)
    lresult = Timestamp(lresult)
    if first_tzinfo is not None:
        fresult = fresult.tz_localize('UTC').tz_convert(first_tzinfo)
    if last_tzinfo is not None:
        lresult = lresult.tz_localize('UTC').tz_convert(last_tzinfo)
    return fresult, lresult


def asfreq(obj, freq, method=None, how=None, normalize=False, fill_value=None):
    """
    Utility frequency conversion method for Series/DataFrame.
    """
    if isinstance(obj.index, PeriodIndex):
        if method is not None:
            raise NotImplementedError("'method' argument is not supported")

        if how is None:
            how = 'E'

        new_obj = obj.copy()
        new_obj.index = obj.index.asfreq(freq, how=how)

    elif len(obj.index) == 0:
        new_obj = obj.copy()
        new_obj.index = obj.index._shallow_copy(freq=to_offset(freq))

    else:
        dti = date_range(obj.index[0], obj.index[-1], freq=freq)
        dti.name = obj.index.name
        new_obj = obj.reindex(dti, method=method, fill_value=fill_value)
        if normalize:
            new_obj.index = new_obj.index.normalize()

    return new_obj
