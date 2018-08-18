# -*- coding: utf-8 -*-
"""
Base and utility classes for tseries type pandas objects.
"""
import warnings

from pandas import compat
from pandas.compat.numpy import function as nv
from pandas.core.tools.timedeltas import to_timedelta

import numpy as np

from pandas._libs import lib, iNaT, NaT
from pandas._libs.tslibs.timestamps import round_ns

from pandas.core.dtypes.common import (
    ensure_int64,
    is_dtype_equal,
    is_float,
    is_integer,
    is_list_like,
    is_scalar,
    is_bool_dtype,
    is_categorical_dtype,
    is_datetime_or_timedelta_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_string_dtype)
from pandas.core.dtypes.generic import (
    ABCIndex, ABCSeries, ABCPeriodIndex, ABCIndexClass)
from pandas.core.dtypes.missing import isna
from pandas.core import common as com, algorithms, ops

import pandas.io.formats.printing as printing

from pandas.core.arrays.datetimelike import DatetimeLikeArrayMixin
from pandas.core.indexes.base import Index, _index_shared_docs
from pandas.util._decorators import Appender, cache_readonly
import pandas.core.dtypes.concat as _concat

import pandas.core.indexes.base as ibase
_index_doc_kwargs = dict(ibase._index_doc_kwargs)


class DatelikeOps(object):
    """ common ops for DatetimeIndex/PeriodIndex, but not TimedeltaIndex """

    def strftime(self, date_format):
        return Index(self.format(date_format=date_format),
                     dtype=compat.text_type)
    strftime.__doc__ = """
    Convert to Index using specified date_format.

    Return an Index of formatted strings specified by date_format, which
    supports the same string format as the python standard library. Details
    of the string format can be found in `python string format doc <{0}>`__

    Parameters
    ----------
    date_format : str
        Date format string (e.g. "%Y-%m-%d").

    Returns
    -------
    Index
        Index of formatted strings

    See Also
    --------
    pandas.to_datetime : Convert the given argument to datetime
    DatetimeIndex.normalize : Return DatetimeIndex with times to midnight.
    DatetimeIndex.round : Round the DatetimeIndex to the specified freq.
    DatetimeIndex.floor : Floor the DatetimeIndex to the specified freq.

    Examples
    --------
    >>> rng = pd.date_range(pd.Timestamp("2018-03-10 09:00"),
    ...                     periods=3, freq='s')
    >>> rng.strftime('%B %d, %Y, %r')
    Index(['March 10, 2018, 09:00:00 AM', 'March 10, 2018, 09:00:01 AM',
           'March 10, 2018, 09:00:02 AM'],
          dtype='object')
    """.format("https://docs.python.org/3/library/datetime.html"
               "#strftime-and-strptime-behavior")


class TimelikeOps(object):
    """ common ops for TimedeltaIndex/DatetimeIndex, but not PeriodIndex """

    _round_doc = (
        """
        {op} the data to the specified `freq`.

        Parameters
        ----------
        freq : str or Offset
            The frequency level to {op} the index to. Must be a fixed
            frequency like 'S' (second) not 'ME' (month end). See
            :ref:`frequency aliases <timeseries.offset_aliases>` for
            a list of possible `freq` values.

        Returns
        -------
        DatetimeIndex, TimedeltaIndex, or Series
            Index of the same type for a DatetimeIndex or TimedeltaIndex,
            or a Series with the same index for a Series.

        Raises
        ------
        ValueError if the `freq` cannot be converted.

        Examples
        --------
        **DatetimeIndex**

        >>> rng = pd.date_range('1/1/2018 11:59:00', periods=3, freq='min')
        >>> rng
        DatetimeIndex(['2018-01-01 11:59:00', '2018-01-01 12:00:00',
                       '2018-01-01 12:01:00'],
                      dtype='datetime64[ns]', freq='T')
        """)

    _round_example = (
        """>>> rng.round('H')
        DatetimeIndex(['2018-01-01 12:00:00', '2018-01-01 12:00:00',
                       '2018-01-01 12:00:00'],
                      dtype='datetime64[ns]', freq=None)

        **Series**

        >>> pd.Series(rng).dt.round("H")
        0   2018-01-01 12:00:00
        1   2018-01-01 12:00:00
        2   2018-01-01 12:00:00
        dtype: datetime64[ns]
        """)

    _floor_example = (
        """>>> rng.floor('H')
        DatetimeIndex(['2018-01-01 11:00:00', '2018-01-01 12:00:00',
                       '2018-01-01 12:00:00'],
                      dtype='datetime64[ns]', freq=None)

        **Series**

        >>> pd.Series(rng).dt.floor("H")
        0   2018-01-01 11:00:00
        1   2018-01-01 12:00:00
        2   2018-01-01 12:00:00
        dtype: datetime64[ns]
        """
    )

    _ceil_example = (
        """>>> rng.ceil('H')
        DatetimeIndex(['2018-01-01 12:00:00', '2018-01-01 12:00:00',
                       '2018-01-01 13:00:00'],
                      dtype='datetime64[ns]', freq=None)

        **Series**

        >>> pd.Series(rng).dt.ceil("H")
        0   2018-01-01 12:00:00
        1   2018-01-01 12:00:00
        2   2018-01-01 13:00:00
        dtype: datetime64[ns]
        """
    )

    def _round(self, freq, rounder):
        # round the local times
        values = _ensure_datetimelike_to_i8(self)
        result = round_ns(values, rounder, freq)
        result = self._maybe_mask_results(result, fill_value=NaT)

        attribs = self._get_attributes_dict()
        if 'freq' in attribs:
            attribs['freq'] = None
        if 'tz' in attribs:
            attribs['tz'] = None
        return self._ensure_localized(
            self._shallow_copy(result, **attribs))

    @Appender((_round_doc + _round_example).format(op="round"))
    def round(self, freq, *args, **kwargs):
        return self._round(freq, np.round)

    @Appender((_round_doc + _floor_example).format(op="floor"))
    def floor(self, freq):
        return self._round(freq, np.floor)

    @Appender((_round_doc + _ceil_example).format(op="ceil"))
    def ceil(self, freq):
        return self._round(freq, np.ceil)


class DatetimeIndexOpsMixin(DatetimeLikeArrayMixin):
    """ common ops mixin to support a unified interface datetimelike Index """

    # DatetimeLikeArrayMixin assumes subclasses are mutable, so these are
    # properties there.  They can be made into cache_readonly for Index
    # subclasses bc they are immutable
    inferred_freq = cache_readonly(DatetimeLikeArrayMixin.inferred_freq.fget)
    _isnan = cache_readonly(DatetimeLikeArrayMixin._isnan.fget)
    hasnans = cache_readonly(DatetimeLikeArrayMixin.hasnans.fget)
    _resolution = cache_readonly(DatetimeLikeArrayMixin._resolution.fget)
    resolution = cache_readonly(DatetimeLikeArrayMixin.resolution.fget)

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if not isinstance(other, ABCIndexClass):
            return False
        elif not isinstance(other, type(self)):
            try:
                other = type(self)(other)
            except Exception:
                return False

        if not is_dtype_equal(self.dtype, other.dtype):
            # have different timezone
            return False

        # ToDo: Remove this when PeriodDtype is added
        elif isinstance(self, ABCPeriodIndex):
            if not isinstance(other, ABCPeriodIndex):
                return False
            if self.freq != other.freq:
                return False

        return np.array_equal(self.asi8, other.asi8)

    @staticmethod
    def _join_i8_wrapper(joinf, dtype, with_indexers=True):
        """ create the join wrapper methods """

        @staticmethod
        def wrapper(left, right):
            if isinstance(left, (np.ndarray, ABCIndex, ABCSeries)):
                left = left.view('i8')
            if isinstance(right, (np.ndarray, ABCIndex, ABCSeries)):
                right = right.view('i8')
            results = joinf(left, right)
            if with_indexers:
                join_index, left_indexer, right_indexer = results
                join_index = join_index.view(dtype)
                return join_index, left_indexer, right_indexer
            return results

        return wrapper

    @Appender(DatetimeLikeArrayMixin._evaluate_compare.__doc__)
    def _evaluate_compare(self, other, op):
        result = DatetimeLikeArrayMixin._evaluate_compare(self, other, op)
        if is_bool_dtype(result):
            return result
        try:
            return Index(result)
        except TypeError:
            return result

    def _ensure_localized(self, result):
        """
        ensure that we are re-localized

        This is for compat as we can then call this on all datetimelike
        indexes generally (ignored for Period/Timedelta)

        Parameters
        ----------
        result : DatetimeIndex / i8 ndarray

        Returns
        -------
        localized DTI
        """

        # reconvert to local tz
        if getattr(self, 'tz', None) is not None:
            if not isinstance(result, ABCIndexClass):
                result = self._simple_new(result)
            result = result.tz_localize(self.tz)
        return result

    def _box_values_as_index(self):
        """
        return object Index which contains boxed values
        """
        from pandas.core.index import Index
        return Index(self._box_values(self.asi8), name=self.name, dtype=object)

    def _format_with_header(self, header, **kwargs):
        return header + list(self._format_native_types(**kwargs))

    @Appender(_index_shared_docs['__contains__'] % _index_doc_kwargs)
    def __contains__(self, key):
        try:
            res = self.get_loc(key)
            return (is_scalar(res) or isinstance(res, slice) or
                    (is_list_like(res) and len(res)))
        except (KeyError, TypeError, ValueError):
            return False

    contains = __contains__

    # Try to run function on index first, and then on elements of index
    # Especially important for group-by functionality
    def map(self, f):
        try:
            result = f(self)

            # Try to use this result if we can
            if isinstance(result, np.ndarray):
                result = Index(result)

            if not isinstance(result, Index):
                raise TypeError('The map function must return an Index object')
            return result
        except Exception:
            return self.astype(object).map(f)

    def sort_values(self, return_indexer=False, ascending=True):
        """
        Return sorted copy of Index
        """
        if return_indexer:
            _as = self.argsort()
            if not ascending:
                _as = _as[::-1]
            sorted_index = self.take(_as)
            return sorted_index, _as
        else:
            sorted_values = np.sort(self._ndarray_values)
            attribs = self._get_attributes_dict()
            freq = attribs['freq']

            if freq is not None and not isinstance(self, ABCPeriodIndex):
                if freq.n > 0 and not ascending:
                    freq = freq * -1
                elif freq.n < 0 and ascending:
                    freq = freq * -1
            attribs['freq'] = freq

            if not ascending:
                sorted_values = sorted_values[::-1]

            return self._simple_new(sorted_values, **attribs)

    @Appender(_index_shared_docs['take'] % _index_doc_kwargs)
    def take(self, indices, axis=0, allow_fill=True,
             fill_value=None, **kwargs):
        nv.validate_take(tuple(), kwargs)
        indices = ensure_int64(indices)

        maybe_slice = lib.maybe_indices_to_slice(indices, len(self))
        if isinstance(maybe_slice, slice):
            return self[maybe_slice]

        taken = self._assert_take_fillable(self.asi8, indices,
                                           allow_fill=allow_fill,
                                           fill_value=fill_value,
                                           na_value=iNaT)

        # keep freq in PeriodIndex, reset otherwise
        freq = self.freq if isinstance(self, ABCPeriodIndex) else None
        return self._shallow_copy(taken, freq=freq)

    _can_hold_na = True

    _na_value = NaT
    """The expected NA value to use with this index."""

    @property
    def asobject(self):
        """Return object Index which contains boxed values.

        .. deprecated:: 0.23.0
            Use ``astype(object)`` instead.

        *this is an internal non-public method*
        """
        warnings.warn("'asobject' is deprecated. Use 'astype(object)'"
                      " instead", FutureWarning, stacklevel=2)
        return self.astype(object)

    def _convert_tolerance(self, tolerance, target):
        tolerance = np.asarray(to_timedelta(tolerance, box=False))
        if target.size != tolerance.size and tolerance.size > 1:
            raise ValueError('list-like tolerance size must match '
                             'target index size')
        return tolerance

    def tolist(self):
        """
        return a list of the underlying data
        """
        return list(self.astype(object))

    def min(self, axis=None, *args, **kwargs):
        """
        Return the minimum value of the Index or minimum along
        an axis.

        See also
        --------
        numpy.ndarray.min
        """
        nv.validate_min(args, kwargs)

        try:
            i8 = self.asi8

            # quick check
            if len(i8) and self.is_monotonic:
                if i8[0] != iNaT:
                    return self._box_func(i8[0])

            if self.hasnans:
                min_stamp = self[~self._isnan].asi8.min()
            else:
                min_stamp = i8.min()
            return self._box_func(min_stamp)
        except ValueError:
            return self._na_value

    def argmin(self, axis=None, *args, **kwargs):
        """
        Returns the indices of the minimum values along an axis.
        See `numpy.ndarray.argmin` for more information on the
        `axis` parameter.

        See also
        --------
        numpy.ndarray.argmin
        """
        nv.validate_argmin(args, kwargs)

        i8 = self.asi8
        if self.hasnans:
            mask = self._isnan
            if mask.all():
                return -1
            i8 = i8.copy()
            i8[mask] = np.iinfo('int64').max
        return i8.argmin()

    def max(self, axis=None, *args, **kwargs):
        """
        Return the maximum value of the Index or maximum along
        an axis.

        See also
        --------
        numpy.ndarray.max
        """
        nv.validate_max(args, kwargs)

        try:
            i8 = self.asi8

            # quick check
            if len(i8) and self.is_monotonic:
                if i8[-1] != iNaT:
                    return self._box_func(i8[-1])

            if self.hasnans:
                max_stamp = self[~self._isnan].asi8.max()
            else:
                max_stamp = i8.max()
            return self._box_func(max_stamp)
        except ValueError:
            return self._na_value

    def argmax(self, axis=None, *args, **kwargs):
        """
        Returns the indices of the maximum values along an axis.
        See `numpy.ndarray.argmax` for more information on the
        `axis` parameter.

        See also
        --------
        numpy.ndarray.argmax
        """
        nv.validate_argmax(args, kwargs)

        i8 = self.asi8
        if self.hasnans:
            mask = self._isnan
            if mask.all():
                return -1
            i8 = i8.copy()
            i8[mask] = 0
        return i8.argmax()

    @property
    def _formatter_func(self):
        raise com.AbstractMethodError(self)

    def _format_attrs(self):
        """
        Return a list of tuples of the (attr,formatted_value)
        """
        attrs = super(DatetimeIndexOpsMixin, self)._format_attrs()
        for attrib in self._attributes:
            if attrib == 'freq':
                freq = self.freqstr
                if freq is not None:
                    freq = "'%s'" % freq
                attrs.append(('freq', freq))
        return attrs

    def _convert_scalar_indexer(self, key, kind=None):
        """
        we don't allow integer or float indexing on datetime-like when using
        loc

        Parameters
        ----------
        key : label of the slice bound
        kind : {'ix', 'loc', 'getitem', 'iloc'} or None
        """

        assert kind in ['ix', 'loc', 'getitem', 'iloc', None]

        # we don't allow integer/float indexing for loc
        # we don't allow float indexing for ix/getitem
        if is_scalar(key):
            is_int = is_integer(key)
            is_flt = is_float(key)
            if kind in ['loc'] and (is_int or is_flt):
                self._invalid_indexer('index', key)
            elif kind in ['ix', 'getitem'] and is_flt:
                self._invalid_indexer('index', key)

        return (super(DatetimeIndexOpsMixin, self)
                ._convert_scalar_indexer(key, kind=kind))

    @classmethod
    def _add_datetimelike_methods(cls):
        """
        add in the datetimelike methods (as we may have to override the
        superclass)
        """

        def __add__(self, other):
            # dispatch to ExtensionArray implementation
            result = super(cls, self).__add__(other)
            return wrap_arithmetic_op(self, other, result)

        cls.__add__ = __add__

        def __radd__(self, other):
            # alias for __add__
            return self.__add__(other)
        cls.__radd__ = __radd__

        def __sub__(self, other):
            # dispatch to ExtensionArray implementation
            result = super(cls, self).__sub__(other)
            return wrap_arithmetic_op(self, other, result)

        cls.__sub__ = __sub__

        def __rsub__(self, other):
            result = super(cls, self).__rsub__(other)
            return wrap_arithmetic_op(self, other, result)

        cls.__rsub__ = __rsub__

    def isin(self, values):
        """
        Compute boolean array of whether each index value is found in the
        passed set of values

        Parameters
        ----------
        values : set or sequence of values

        Returns
        -------
        is_contained : ndarray (boolean dtype)
        """
        if not isinstance(values, type(self)):
            try:
                values = type(self)(values)
            except ValueError:
                return self.astype(object).isin(values)

        return algorithms.isin(self.asi8, values.asi8)

    def repeat(self, repeats, *args, **kwargs):
        """
        Analogous to ndarray.repeat
        """
        nv.validate_repeat(args, kwargs)
        if isinstance(self, ABCPeriodIndex):
            freq = self.freq
        else:
            freq = None
        return self._shallow_copy(self.asi8.repeat(repeats),
                                  freq=freq)

    @Appender(_index_shared_docs['where'] % _index_doc_kwargs)
    def where(self, cond, other=None):
        other = _ensure_datetimelike_to_i8(other)
        values = _ensure_datetimelike_to_i8(self)
        result = np.where(cond, values, other).astype('i8')

        result = self._ensure_localized(result)
        return self._shallow_copy(result,
                                  **self._get_attributes_dict())

    def _summary(self, name=None):
        """
        Return a summarized representation

        Parameters
        ----------
        name : str
            name to use in the summary representation

        Returns
        -------
        String with a summarized representation of the index
        """
        formatter = self._formatter_func
        if len(self) > 0:
            index_summary = ', %s to %s' % (formatter(self[0]),
                                            formatter(self[-1]))
        else:
            index_summary = ''

        if name is None:
            name = type(self).__name__
        result = '%s: %s entries%s' % (printing.pprint_thing(name),
                                       len(self), index_summary)
        if self.freq:
            result += '\nFreq: %s' % self.freqstr

        # display as values, not quoted
        result = result.replace("'", "")
        return result

    def _concat_same_dtype(self, to_concat, name):
        """
        Concatenate to_concat which has the same class
        """
        attribs = self._get_attributes_dict()
        attribs['name'] = name

        if not isinstance(self, ABCPeriodIndex):
            # reset freq
            attribs['freq'] = None

        if getattr(self, 'tz', None) is not None:
            return _concat._concat_datetimetz(to_concat, name)
        else:
            new_data = np.concatenate([c.asi8 for c in to_concat])
        return self._simple_new(new_data, **attribs)

    def astype(self, dtype, copy=True):
        if is_object_dtype(dtype):
            return self._box_values_as_index()
        elif is_string_dtype(dtype) and not is_categorical_dtype(dtype):
            return Index(self.format(), name=self.name, dtype=object)
        elif is_integer_dtype(dtype):
            return Index(self.values.astype('i8', copy=copy), name=self.name,
                         dtype='i8')
        elif (is_datetime_or_timedelta_dtype(dtype) and
              not is_dtype_equal(self.dtype, dtype)) or is_float_dtype(dtype):
            # disallow conversion between datetime/timedelta,
            # and conversions for any datetimelike to float
            msg = 'Cannot cast {name} to dtype {dtype}'
            raise TypeError(msg.format(name=type(self).__name__, dtype=dtype))
        return super(DatetimeIndexOpsMixin, self).astype(dtype, copy=copy)


def _ensure_datetimelike_to_i8(other):
    """ helper for coercing an input scalar or array to i8 """
    if is_scalar(other) and isna(other):
        other = iNaT
    elif isinstance(other, ABCIndexClass):
        # convert tz if needed
        if getattr(other, 'tz', None) is not None:
            other = other.tz_localize(None).asi8
        else:
            other = other.asi8
    else:
        try:
            other = np.array(other, copy=False).view('i8')
        except TypeError:
            # period array cannot be coerces to int
            other = Index(other).asi8
    return other


def wrap_arithmetic_op(self, other, result):
    if result is NotImplemented:
        return NotImplemented

    if not isinstance(result, Index):
        # Index.__new__ will choose appropriate subclass for dtype
        result = Index(result)

    res_name = ops.get_op_result_name(self, other)
    result.name = res_name
    return result
