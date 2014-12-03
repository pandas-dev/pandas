"""
Base and utility classes for tseries type pandas objects.
"""


from datetime import datetime, time, timedelta

from pandas import compat
import numpy as np
from pandas.core import common as com
import pandas.tslib as tslib
import pandas.lib as lib
from pandas.core.index import Index
from pandas.util.decorators import Appender, cache_readonly
from pandas.tseries.frequencies import (
    infer_freq, to_offset, get_period_alias,
    Resolution)
import pandas.algos as _algos

class DatetimeIndexOpsMixin(object):
    """ common ops mixin to support a unified inteface datetimelike Index """

    def __iter__(self):
        return (self._box_func(v) for v in self.asi8)

    @staticmethod
    def _join_i8_wrapper(joinf, dtype, with_indexers=True):
        """ create the join wrapper methods """

        @staticmethod
        def wrapper(left, right):
            if isinstance(left, (np.ndarray, com.ABCIndex, com.ABCSeries)):
                left = left.view('i8')
            if isinstance(right, (np.ndarray, com.ABCIndex, com.ABCSeries)):
                right = right.view('i8')
            results = joinf(left, right)
            if with_indexers:
                join_index, left_indexer, right_indexer = results
                join_index = join_index.view(dtype)
                return join_index, left_indexer, right_indexer
            return results

        return wrapper

    @property
    def _box_func(self):
        """
        box function to get object from internal representation
        """
        raise NotImplementedError

    def _box_values(self, values):
        """
        apply box func to passed values
        """
        return lib.map_infer(values, self._box_func)

    def groupby(self, f):
        objs = self.asobject.values
        return _algos.groupby_object(objs, f)

    def _format_with_header(self, header, **kwargs):
        return header + self._format_native_types(**kwargs)

    def __contains__(self, key):
        try:
            res = self.get_loc(key)
            return np.isscalar(res) or type(res) == slice
        except (KeyError, TypeError):
            return False

    @cache_readonly
    def inferred_freq(self):
        try:
            return infer_freq(self)
        except ValueError:
            return None

    # Try to run function on index first, and then on elements of index
    # Especially important for group-by functionality
    def map(self, f):
        try:
            result = f(self)
            if not isinstance(result, (np.ndarray, Index)):
                raise TypeError
            return result
        except Exception:
            return _algos.arrmap_object(self.asobject.values, f)

    def order(self, return_indexer=False, ascending=True):
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
            sorted_values = np.sort(self.values)
            if not ascending:
                sorted_values = sorted_values[::-1]
            attribs = self._get_attributes_dict()
            attribs['freq'] = None
            return self._simple_new(sorted_values, **attribs)

    def take(self, indices, axis=0):
        """
        Analogous to ndarray.take
        """
        maybe_slice = lib.maybe_indices_to_slice(com._ensure_int64(indices))
        if isinstance(maybe_slice, slice):
            return self[maybe_slice]
        return super(DatetimeIndexOpsMixin, self).take(indices, axis)

    def get_duplicates(self):
        values = Index.get_duplicates(self)
        return self._simple_new(values)

    @cache_readonly
    def hasnans(self):
        """ return if I have any nans; enables various perf speedups """
        return (self.asi8 == tslib.iNaT).any()

    @property
    def asobject(self):
        from pandas.core.index import Index
        return Index(self._box_values(self.asi8), name=self.name, dtype=object)

    def _maybe_mask_results(self, result, fill_value=None, convert=None):
        """
        Parameters
        ----------
        result : a ndarray
        convert : string/dtype or None

        Returns
        -------
        result : ndarray with values replace by the fill_value

        mask the result if needed, convert to the provided dtype if its not None

        This is an internal routine
        """

        if self.hasnans:
            mask = self.asi8 == tslib.iNaT
            if convert:
                result = result.astype(convert)
            if fill_value is None:
                fill_value = np.nan
            result[mask] = fill_value
        return result

    def tolist(self):
        """
        return a list of the underlying data
        """
        return list(self.asobject)

    def min(self, axis=None):
        """
        return the minimum value of the Index

        See also
        --------
        numpy.ndarray.min
        """
        try:
            i8 = self.asi8

            # quick check
            if len(i8) and self.is_monotonic:
                if i8[0] != tslib.iNaT:
                    return self._box_func(i8[0])

            if self.hasnans:
                mask = i8 == tslib.iNaT
                min_stamp = self[~mask].asi8.min()
            else:
                min_stamp = i8.min()
            return self._box_func(min_stamp)
        except ValueError:
            return self._na_value

    def argmin(self, axis=None):
        """
        return a ndarray of the minimum argument indexer

        See also
        --------
        numpy.ndarray.argmin
        """

        i8 = self.asi8
        if self.hasnans:
            mask = i8 == tslib.iNaT
            if mask.all():
                return -1
            i8 = i8.copy()
            i8[mask] = np.iinfo('int64').max
        return i8.argmin()

    def max(self, axis=None):
        """
        return the maximum value of the Index

        See also
        --------
        numpy.ndarray.max
        """
        try:
            i8 = self.asi8

            # quick check
            if len(i8) and self.is_monotonic:
                if i8[-1] != tslib.iNaT:
                    return self._box_func(i8[-1])

            if self.hasnans:
                mask = i8 == tslib.iNaT
                max_stamp = self[~mask].asi8.max()
            else:
                max_stamp = i8.max()
            return self._box_func(max_stamp)
        except ValueError:
            return self._na_value

    def argmax(self, axis=None):
        """
        return a ndarray of the maximum argument indexer

        See also
        --------
        numpy.ndarray.argmax
        """

        i8 = self.asi8
        if self.hasnans:
            mask = i8 == tslib.iNaT
            if mask.all():
                return -1
            i8 = i8.copy()
            i8[mask] = 0
        return i8.argmax()

    @property
    def _formatter_func(self):
        """
        Format function to convert value to representation
        """
        return str

    def _format_footer(self):
        raise NotImplementedError

    def __unicode__(self):
        formatter = self._formatter_func
        summary = str(self.__class__) + '\n'

        n = len(self)
        if n == 0:
            pass
        elif n == 1:
            first = formatter(self[0])
            summary += '[%s]\n' % first
        elif n == 2:
            first = formatter(self[0])
            last = formatter(self[-1])
            summary += '[%s, %s]\n' % (first, last)
        else:
            first = formatter(self[0])
            last = formatter(self[-1])
            summary += '[%s, ..., %s]\n' % (first, last)

        summary += self._format_footer()
        return summary

    @cache_readonly
    def _resolution(self):
        from pandas.tseries.frequencies import Resolution
        return Resolution.get_reso_from_freq(self.freqstr)

    @cache_readonly
    def resolution(self):
        """
        Returns day, hour, minute, second, millisecond or microsecond
        """
        from pandas.tseries.frequencies import get_reso_string
        return get_reso_string(self._resolution)

    def _add_datelike(self, other):
        raise NotImplementedError

    def _sub_datelike(self, other):
        raise NotImplementedError

    @classmethod
    def _add_datetimelike_methods(cls):
        """ add in the datetimelike methods (as we may have to override the superclass) """

        def __add__(self, other):
            from pandas.core.index import Index
            from pandas.tseries.tdi import TimedeltaIndex
            from pandas.tseries.offsets import DateOffset
            if isinstance(other, TimedeltaIndex):
                return self._add_delta(other)
            elif isinstance(self, TimedeltaIndex) and isinstance(other, Index):
                if hasattr(other,'_add_delta'):
                    return other._add_delta(self)
                raise TypeError("cannot add TimedeltaIndex and {typ}".format(typ=type(other)))
            elif isinstance(other, Index):
                return self.union(other)
            elif isinstance(other, (DateOffset, timedelta, np.timedelta64, tslib.Timedelta)):
                return self._add_delta(other)
            elif com.is_integer(other):
                return self.shift(other)
            elif isinstance(other, (tslib.Timestamp, datetime)):
                return self._add_datelike(other)
            else:  # pragma: no cover
                return NotImplemented
        cls.__add__ = __add__
        cls.__radd__ = __add__

        def __sub__(self, other):
            from pandas.core.index import Index
            from pandas.tseries.tdi import TimedeltaIndex
            from pandas.tseries.offsets import DateOffset
            if isinstance(other, TimedeltaIndex):
                return self._add_delta(-other)
            elif isinstance(self, TimedeltaIndex) and isinstance(other, Index):
                if not isinstance(other, TimedeltaIndex):
                    raise TypeError("cannot subtract TimedeltaIndex and {typ}".format(typ=type(other)))
                return self._add_delta(-other)
            elif isinstance(other, Index):
                return self.difference(other)
            elif isinstance(other, (DateOffset, timedelta, np.timedelta64, tslib.Timedelta)):
                return self._add_delta(-other)
            elif com.is_integer(other):
                return self.shift(-other)
            elif isinstance(other, (tslib.Timestamp, datetime)):
                return self._sub_datelike(other)
            else:  # pragma: no cover
                return NotImplemented
        cls.__sub__ = __sub__

        def __rsub__(self, other):
            return -(self - other)
        cls.__rsub__ = __rsub__

        cls.__iadd__ = __add__
        cls.__isub__ = __sub__

    def _add_delta(self, other):
        return NotImplemented

    def _add_delta_td(self, other):
        # add a delta of a timedeltalike
        # return the i8 result view

        inc = tslib._delta_to_nanoseconds(other)
        mask = self.asi8 == tslib.iNaT
        new_values = (self.asi8 + inc).view(self.dtype)
        new_values[mask] = tslib.iNaT
        return new_values.view(self.dtype)

    def _add_delta_tdi(self, other):
        # add a delta of a TimedeltaIndex
        # return the i8 result view

        # delta operation
        if not len(self) == len(other):
            raise ValueError("cannot add indices of unequal length")

        self_i8 = self.asi8
        other_i8 = other.asi8
        mask = (self_i8 == tslib.iNaT) | (other_i8 == tslib.iNaT)
        new_values = self_i8 + other_i8
        new_values[mask] = tslib.iNaT
        return new_values.view(self.dtype)

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
                return self.asobject.isin(values)

        value_set = set(values.asi8)
        return lib.ismember(self.asi8, value_set)

    def shift(self, n, freq=None):
        """
        Specialized shift which produces a DatetimeIndex

        Parameters
        ----------
        n : int
            Periods to shift by
        freq : DateOffset or timedelta-like, optional

        Returns
        -------
        shifted : DatetimeIndex
        """
        if freq is not None and freq != self.freq:
            if isinstance(freq, compat.string_types):
                freq = to_offset(freq)
            result = Index.shift(self, n, freq)

            if hasattr(self,'tz'):
                result.tz = self.tz

            return result

        if n == 0:
            # immutable so OK
            return self

        if self.freq is None:
            raise ValueError("Cannot shift with no freq")

        start = self[0] + n * self.freq
        end = self[-1] + n * self.freq
        attribs = self._get_attributes_dict()
        attribs['start'] = start
        attribs['end'] = end
        return type(self)(**attribs)

    def unique(self):
        """
        Index.unique with handling for DatetimeIndex/PeriodIndex metadata

        Returns
        -------
        result : DatetimeIndex or PeriodIndex
        """
        from pandas.core.index import Int64Index
        result = Int64Index.unique(self)
        return self._simple_new(result, name=self.name, freq=self.freq,
                                tz=getattr(self, 'tz', None))

    def repeat(self, repeats, axis=None):
        """
        Analogous to ndarray.repeat
        """
        return self._simple_new(self.values.repeat(repeats),
                                name=self.name)
