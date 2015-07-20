# pylint: disable=E1101,E1103,W0232
from datetime import datetime, timedelta
import numpy as np
import pandas.tseries.frequencies as frequencies
from pandas.tseries.frequencies import get_freq_code as _gfc
from pandas.tseries.index import DatetimeIndex, Int64Index, Index
from pandas.tseries.base import DatelikeOps, DatetimeIndexOpsMixin
from pandas.tseries.tools import parse_time_string
import pandas.tseries.offsets as offsets

from pandas._period import Period
import pandas._period as period
from pandas._period import (
    get_period_field_arr,
    _validate_end_alias,
    _quarter_to_myear,
)

import pandas.core.common as com
from pandas.core.common import (isnull, _INT64_DTYPE, _maybe_box,
                                _values_from_object, ABCSeries,
                                is_integer, is_float, is_object_dtype,
                                is_float_dtype)
from pandas import compat
from pandas.util.decorators import cache_readonly

from pandas.lib import Timestamp, Timedelta
import pandas.lib as lib
import pandas.tslib as tslib
import pandas.algos as _algos
from pandas.compat import zip, u


def _field_accessor(name, alias, docstring=None):
    def f(self):
        base, mult = _gfc(self.freq)
        return get_period_field_arr(alias, self.values, base)
    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def _get_ordinals(data, freq):
    f = lambda x: Period(x, freq=freq).ordinal
    if isinstance(data[0], Period):
        return period.extract_ordinals(data, freq)
    else:
        return lib.map_infer(data, f)


def dt64arr_to_periodarr(data, freq, tz):
    if data.dtype != np.dtype('M8[ns]'):
        raise ValueError('Wrong dtype: %s' % data.dtype)

    base, mult = _gfc(freq)
    return period.dt64arr_to_periodarr(data.view('i8'), base, tz)

# --- Period index sketch

_DIFFERENT_FREQ_ERROR = "Input has different freq={1} from PeriodIndex(freq={0})"

def _period_index_cmp(opname, nat_result=False):
    """
    Wrap comparison operations to convert datetime-like to datetime64
    """
    def wrapper(self, other):
        if isinstance(other, Period):
            func = getattr(self.values, opname)
            other_base, _ = _gfc(other.freq)
            if other.freq != self.freq:
                msg = _DIFFERENT_FREQ_ERROR.format(self.freqstr, other.freqstr)
                raise ValueError(msg)

            result = func(other.ordinal)
        elif isinstance(other, PeriodIndex):
            if other.freq != self.freq:
                msg = _DIFFERENT_FREQ_ERROR.format(self.freqstr, other.freqstr)
                raise ValueError(msg)

            result = getattr(self.values, opname)(other.values)

            mask = (com.mask_missing(self.values, tslib.iNaT) |
                    com.mask_missing(other.values, tslib.iNaT))
            if mask.any():
                result[mask] = nat_result

            return result
        else:
            other = Period(other, freq=self.freq)
            func = getattr(self.values, opname)
            result = func(other.ordinal)

        if other.ordinal == tslib.iNaT:
            result.fill(nat_result)
        mask = self.values == tslib.iNaT
        if mask.any():
            result[mask] = nat_result

        return result
    return wrapper


class PeriodIndex(DatelikeOps, DatetimeIndexOpsMixin, Int64Index):
    """
    Immutable ndarray holding ordinal values indicating regular periods in
    time such as particular years, quarters, months, etc. A value of 1 is the
    period containing the Gregorian proleptic datetime Jan 1, 0001 00:00:00.
    This ordinal representation is from the scikits.timeseries project.

    For instance,
        # construct period for day 1/1/1 and get the first second
        i = Period(year=1,month=1,day=1,freq='D').asfreq('S', 'S')
        i.ordinal
        ===> 1

    Index keys are boxed to Period objects which carries the metadata (eg,
    frequency information).

    Parameters
    ----------
    data : array-like (1-dimensional), optional
        Optional period-like data to construct index with
    dtype : NumPy dtype (default: i8)
    copy : bool
        Make a copy of input ndarray
    freq : string or period object, optional
        One of pandas period strings or corresponding objects
    start : starting value, period-like, optional
        If data is None, used as the start point in generating regular
        period data.
    periods : int, optional, > 0
        Number of periods to generate, if generating index. Takes precedence
        over end argument
    end : end value, period-like, optional
        If periods is none, generated index will extend to first conforming
        period on or just past end argument
    year : int, array, or Series, default None
    month : int, array, or Series, default None
    quarter : int, array, or Series, default None
    day : int, array, or Series, default None
    hour : int, array, or Series, default None
    minute : int, array, or Series, default None
    second : int, array, or Series, default None
    tz : object, default None
        Timezone for converting datetime64 data to Periods

    Examples
    --------
    >>> idx = PeriodIndex(year=year_arr, quarter=q_arr)

    >>> idx2 = PeriodIndex(start='2000', end='2010', freq='A')
    """
    _box_scalars = True
    _typ = 'periodindex'
    _attributes = ['name','freq']
    _datetimelike_ops = ['year','month','day','hour','minute','second',
                         'weekofyear','week','dayofweek','weekday','dayofyear','quarter', 'qyear', 'freq', 'days_in_month', 'daysinmonth']
    _is_numeric_dtype = False
    freq = None

    __eq__ = _period_index_cmp('__eq__')
    __ne__ = _period_index_cmp('__ne__', nat_result=True)
    __lt__ = _period_index_cmp('__lt__')
    __gt__ = _period_index_cmp('__gt__')
    __le__ = _period_index_cmp('__le__')
    __ge__ = _period_index_cmp('__ge__')

    def __new__(cls, data=None, ordinal=None, freq=None, start=None, end=None,
                periods=None, copy=False, name=None, tz=None, **kwargs):

        if periods is not None:
            if is_float(periods):
                periods = int(periods)
            elif not is_integer(periods):
                raise ValueError('Periods must be a number, got %s' %
                                 str(periods))

        if data is None:
            if ordinal is not None:
                data = np.asarray(ordinal, dtype=np.int64)
            else:
                data, freq = cls._generate_range(start, end, periods,
                                                 freq, kwargs)
        else:
            ordinal, freq = cls._from_arraylike(data, freq, tz)
            data = np.array(ordinal, dtype=np.int64, copy=False)

        return cls._simple_new(data, name=name, freq=freq)

    @classmethod
    def _generate_range(cls, start, end, periods, freq, fields):
        field_count = len(fields)
        if com._count_not_none(start, end) > 0:
            if field_count > 0:
                raise ValueError('Can either instantiate from fields '
                                 'or endpoints, but not both')
            subarr, freq = _get_ordinal_range(start, end, periods, freq)
        elif field_count > 0:
            subarr, freq = _range_from_fields(freq=freq, **fields)
        else:
            raise ValueError('Not enough parameters to construct '
                             'Period range')

        return subarr, freq

    @classmethod
    def _from_arraylike(cls, data, freq, tz):

        if not isinstance(data, (np.ndarray, PeriodIndex, DatetimeIndex, Int64Index)):
            if np.isscalar(data) or isinstance(data, Period):
                raise ValueError('PeriodIndex() must be called with a '
                                 'collection of some kind, %s was passed'
                                 % repr(data))

            # other iterable of some kind
            if not isinstance(data, (list, tuple)):
                data = list(data)

            try:
                data = com._ensure_int64(data)
                if freq is None:
                    raise ValueError('freq not specified')
                data = np.array([Period(x, freq=freq).ordinal for x in data],
                                dtype=np.int64)
            except (TypeError, ValueError):
                data = com._ensure_object(data)

                if freq is None and len(data) > 0:
                    freq = getattr(data[0], 'freq', None)

                if freq is None:
                    raise ValueError('freq not specified and cannot be '
                                     'inferred from first element')

                data = _get_ordinals(data, freq)
        else:
            if isinstance(data, PeriodIndex):
                if freq is None or freq == data.freq:
                    freq = data.freq
                    data = data.values
                else:
                    base1, _ = _gfc(data.freq)
                    base2, _ = _gfc(freq)
                    data = period.period_asfreq_arr(data.values,
                                                    base1, base2, 1)
            else:
                if freq is None and len(data) > 0:
                    freq = getattr(data[0], 'freq', None)

                if freq is None:
                    raise ValueError('freq not specified and cannot be '
                                     'inferred from first element')

                if data.dtype != np.int64:
                    if np.issubdtype(data.dtype, np.datetime64):
                        data = dt64arr_to_periodarr(data, freq, tz)
                    else:
                        try:
                            data = com._ensure_int64(data)
                        except (TypeError, ValueError):
                            data = com._ensure_object(data)
                            data = _get_ordinals(data, freq)

        return data, freq

    @classmethod
    def _simple_new(cls, values, name=None, freq=None, **kwargs):
        if not getattr(values,'dtype',None):
            values = np.array(values,copy=False)
        if is_object_dtype(values):
            return PeriodIndex(values, name=name, freq=freq, **kwargs)

        result = object.__new__(cls)
        result._data = values
        result.name = name
        if freq is None:
            raise ValueError('freq is not specified')
        result.freq = Period._maybe_convert_freq(freq)
        result._reset_identity()
        return result

    def _shallow_copy(self, values=None, infer=False, **kwargs):
        """ we always want to return a PeriodIndex """
        return super(PeriodIndex, self)._shallow_copy(values=values, infer=False, **kwargs)

    def _coerce_scalar_to_index(self, item):
        """
        we need to coerce a scalar to a compat for our index type

        Parameters
        ----------
        item : scalar item to coerce
        """
        return PeriodIndex([item], **self._get_attributes_dict())

    @property
    def _na_value(self):
        return self._box_func(tslib.iNaT)

    def __contains__(self, key):
        if not isinstance(key, Period) or key.freq != self.freq:
            if isinstance(key, compat.string_types):
                try:
                    self.get_loc(key)
                    return True
                except Exception:
                    return False
            return False
        return key.ordinal in self._engine

    def __array_wrap__(self, result, context=None):
        """
        Gets called after a ufunc. Needs additional handling as
        PeriodIndex stores internal data as int dtype

        Replace this to __numpy_ufunc__ in future version
        """
        if isinstance(context, tuple) and len(context) > 0:
            func = context[0]
            if (func is np.add):
                return self._add_delta(context[1][1])
            elif (func is np.subtract):
                return self._add_delta(-context[1][1])
            elif isinstance(func, np.ufunc):
                if 'M->M' not in func.types:
                    msg = "ufunc '{0}' not supported for the PeriodIndex"
                    # This should be TypeError, but TypeError cannot be raised
                    # from here because numpy catches.
                    raise ValueError(msg.format(func.__name__))

        if com.is_bool_dtype(result):
            return result
        return PeriodIndex(result, freq=self.freq, name=self.name)

    @property
    def _box_func(self):
        return lambda x: Period._from_ordinal(ordinal=x, freq=self.freq)

    def _to_embed(self, keep_tz=False):
        """ return an array repr of this object, potentially casting to object """
        return self.asobject.values

    @property
    def _formatter_func(self):
        return lambda x: "'%s'" % x

    def asof_locs(self, where, mask):
        """
        where : array of timestamps
        mask : array of booleans where data is not NA

        """
        where_idx = where
        if isinstance(where_idx, DatetimeIndex):
            where_idx = PeriodIndex(where_idx.values, freq=self.freq)

        locs = self.values[mask].searchsorted(where_idx.values, side='right')

        locs = np.where(locs > 0, locs - 1, 0)
        result = np.arange(len(self))[mask].take(locs)

        first = mask.argmax()
        result[(locs == 0) & (where_idx.values < self.values[first])] = -1

        return result

    def _array_values(self):
        return self.asobject

    def astype(self, dtype):
        dtype = np.dtype(dtype)
        if dtype == np.object_:
            return Index(np.array(list(self), dtype), dtype)
        elif dtype == _INT64_DTYPE:
            return Index(self.values, dtype)
        raise ValueError('Cannot cast PeriodIndex to dtype %s' % dtype)

    def searchsorted(self, key, side='left'):
        if isinstance(key, Period):
            if key.freq != self.freq:
                msg = _DIFFERENT_FREQ_ERROR.format(self.freqstr, key.freqstr)
                raise ValueError(msg)
            key = key.ordinal
        elif isinstance(key, compat.string_types):
            key = Period(key, freq=self.freq).ordinal

        return self.values.searchsorted(key, side=side)

    @property
    def is_all_dates(self):
        return True

    @property
    def is_full(self):
        """
        Returns True if there are any missing periods from start to end
        """
        if len(self) == 0:
            return True
        if not self.is_monotonic:
            raise ValueError('Index is not monotonic')
        values = self.values
        return ((values[1:] - values[:-1]) < 2).all()

    def asfreq(self, freq=None, how='E'):
        """
        Convert the PeriodIndex to the specified frequency `freq`.

        Parameters
        ----------

        freq : str
            a frequency
        how : str {'E', 'S'}
            'E', 'END', or 'FINISH' for end,
            'S', 'START', or 'BEGIN' for start.
            Whether the elements should be aligned to the end
            or start within pa period. January 31st ('END') vs.
            Janury 1st ('START') for example.

        Returns
        -------

        new : PeriodIndex with the new frequency

        Examples
        --------
        >>> pidx = pd.period_range('2010-01-01', '2015-01-01', freq='A')
        >>> pidx
        <class 'pandas.tseries.period.PeriodIndex'>
        [2010, ..., 2015]
        Length: 6, Freq: A-DEC

        >>> pidx.asfreq('M')
        <class 'pandas.tseries.period.PeriodIndex'>
        [2010-12, ..., 2015-12]
        Length: 6, Freq: M

        >>> pidx.asfreq('M', how='S')
        <class 'pandas.tseries.period.PeriodIndex'>
        [2010-01, ..., 2015-01]
        Length: 6, Freq: M
        """
        how = _validate_end_alias(how)

        freq = frequencies.get_standard_freq(freq)

        base1, mult1 = _gfc(self.freq)
        base2, mult2 = _gfc(freq)

        asi8 = self.asi8
        # mult1 can't be negative or 0
        end = how == 'E'
        if end:
            ordinal = asi8 + mult1 - 1
        else:
            ordinal = asi8

        new_data = period.period_asfreq_arr(ordinal, base1, base2, end)

        if self.hasnans:
            mask = asi8 == tslib.iNaT
            new_data[mask] = tslib.iNaT

        return self._simple_new(new_data, self.name, freq=freq)

    def to_datetime(self, dayfirst=False):
        return self.to_timestamp()

    year = _field_accessor('year', 0, "The year of the period")
    month = _field_accessor('month', 3, "The month as January=1, December=12")
    day = _field_accessor('day', 4, "The days of the period")
    hour = _field_accessor('hour', 5, "The hour of the period")
    minute = _field_accessor('minute', 6, "The minute of the period")
    second = _field_accessor('second', 7, "The second of the period")
    weekofyear = _field_accessor('week', 8, "The week ordinal of the year")
    week = weekofyear
    dayofweek = _field_accessor('dayofweek', 10, "The day of the week with Monday=0, Sunday=6")
    weekday = dayofweek
    dayofyear = day_of_year = _field_accessor('dayofyear', 9, "The ordinal day of the year")
    quarter = _field_accessor('quarter', 2, "The quarter of the date")
    qyear = _field_accessor('qyear', 1)
    days_in_month = _field_accessor('days_in_month', 11, "The number of days in the month")
    daysinmonth = days_in_month

    def _get_object_array(self):
        freq = self.freq
        return np.array([ Period._from_ordinal(ordinal=x, freq=freq) for x in self.values], copy=False)

    def _mpl_repr(self):
        # how to represent ourselves to matplotlib
        return self._get_object_array()

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if (not hasattr(other, 'inferred_type') or
                other.inferred_type != 'int64'):
            try:
                other = PeriodIndex(other)
            except:
                return False

        return np.array_equal(self.asi8, other.asi8)

    def to_timestamp(self, freq=None, how='start'):
        """
        Cast to DatetimeIndex

        Parameters
        ----------
        freq : string or DateOffset, default 'D' for week or longer, 'S'
               otherwise
            Target frequency
        how : {'s', 'e', 'start', 'end'}

        Returns
        -------
        DatetimeIndex
        """
        how = _validate_end_alias(how)

        if freq is None:
            base, mult = _gfc(self.freq)
            freq = frequencies.get_to_timestamp_base(base)

        base, mult = _gfc(freq)
        new_data = self.asfreq(freq, how)

        new_data = period.periodarr_to_dt64arr(new_data.values, base)
        return DatetimeIndex(new_data, freq='infer', name=self.name)

    def _maybe_convert_timedelta(self, other):
        if isinstance(other, (timedelta, np.timedelta64, offsets.Tick, Timedelta)):
            offset = frequencies.to_offset(self.freq.rule_code)
            if isinstance(offset, offsets.Tick):
                nanos = tslib._delta_to_nanoseconds(other)
                offset_nanos = tslib._delta_to_nanoseconds(offset)
                if nanos % offset_nanos == 0:
                    return nanos // offset_nanos
        elif isinstance(other, offsets.DateOffset):
            freqstr = frequencies.get_standard_freq(other)
            base = frequencies.get_base_alias(freqstr)
            if base == self.freq.rule_code:
                return other.n
        elif isinstance(other, np.ndarray):
            if com.is_integer_dtype(other):
                return other
            elif com.is_timedelta64_dtype(other):
                offset = frequencies.to_offset(self.freq)
                if isinstance(offset, offsets.Tick):
                    nanos = tslib._delta_to_nanoseconds(other)
                    offset_nanos = tslib._delta_to_nanoseconds(offset)
                    if (nanos % offset_nanos).all() == 0:
                        return nanos // offset_nanos
        msg = "Input has different freq from PeriodIndex(freq={0})"
        raise ValueError(msg.format(self.freqstr))

    def _add_delta(self, other):
        ordinal_delta = self._maybe_convert_timedelta(other)
        return self.shift(ordinal_delta)

    def shift(self, n):
        """
        Specialized shift which produces an PeriodIndex

        Parameters
        ----------
        n : int
            Periods to shift by

        Returns
        -------
        shifted : PeriodIndex
        """
        mask = self.values == tslib.iNaT
        values = self.values + n * self.freq.n
        values[mask] = tslib.iNaT
        return PeriodIndex(data=values, name=self.name, freq=self.freq)

    @cache_readonly
    def dtype_str(self):
        """ return the dtype str of the underlying data """
        return self.inferred_type

    @property
    def inferred_type(self):
        # b/c data is represented as ints make sure we can't have ambiguous
        # indexing
        return 'period'

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """
        s = _values_from_object(series)
        try:
            return _maybe_box(self, super(PeriodIndex, self).get_value(s, key), series, key)
        except (KeyError, IndexError):
            try:
                asdt, parsed, reso = parse_time_string(key, self.freq)
                grp = frequencies.Resolution.get_freq_group(reso)
                freqn = frequencies.get_freq_group(self.freq)

                vals = self.values

                # if our data is higher resolution than requested key, slice
                if grp < freqn:
                    iv = Period(asdt, freq=(grp, 1))
                    ord1 = iv.asfreq(self.freq, how='S').ordinal
                    ord2 = iv.asfreq(self.freq, how='E').ordinal

                    if ord2 < vals[0] or ord1 > vals[-1]:
                        raise KeyError(key)

                    pos = np.searchsorted(self.values, [ord1, ord2])
                    key = slice(pos[0], pos[1] + 1)
                    return series[key]
                elif grp == freqn:
                    key = Period(asdt, freq=self.freq).ordinal
                    return _maybe_box(self, self._engine.get_value(s, key), series, key)
                else:
                    raise KeyError(key)
            except TypeError:
                pass

            key = Period(key, self.freq).ordinal
            return _maybe_box(self, self._engine.get_value(s, key), series, key)

    def get_indexer(self, target, method=None, limit=None, tolerance=None):
        if hasattr(target, 'freq') and target.freq != self.freq:
            raise ValueError('target and index have different freq: '
                             '(%s, %s)' % (target.freq, self.freq))
        return Index.get_indexer(self, target, method, limit, tolerance)

    def get_loc(self, key, method=None, tolerance=None):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int
        """
        try:
            return self._engine.get_loc(key)
        except KeyError:
            if is_integer(key):
                raise

            try:
                asdt, parsed, reso = parse_time_string(key, self.freq)
                key = asdt
            except TypeError:
                pass

            key = Period(key, freq=self.freq)
            try:
                return Index.get_loc(self, key.ordinal, method, tolerance)
            except KeyError:
                raise KeyError(key)

    def _maybe_cast_slice_bound(self, label, side, kind):
        """
        If label is a string or a datetime, cast it to Period.ordinal according to
        resolution.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : string / None

        Returns
        -------
        bound : Period or object

        Notes
        -----
        Value of `side` parameter should be validated in caller.

        """
        if isinstance(label, datetime):
            return Period(label, freq=self.freq)
        elif isinstance(label, compat.string_types):
            try:
                _, parsed, reso = parse_time_string(label, self.freq)
                bounds = self._parsed_string_to_bounds(reso, parsed)
                return bounds[0 if side == 'left' else 1]
            except Exception:
                raise KeyError(label)
        elif is_integer(label) or is_float(label):
            self._invalid_indexer('slice',label)

        return label

    def _parsed_string_to_bounds(self, reso, parsed):
        if reso == 'year':
            t1 = Period(year=parsed.year, freq='A')
        elif reso == 'month':
            t1 = Period(year=parsed.year, month=parsed.month, freq='M')
        elif reso == 'quarter':
            q = (parsed.month - 1) // 3 + 1
            t1 = Period(year=parsed.year, quarter=q, freq='Q-DEC')
        elif reso == 'day':
            t1 = Period(year=parsed.year, month=parsed.month, day=parsed.day,
                        freq='D')
        elif reso == 'hour':
            t1 = Period(year=parsed.year, month=parsed.month, day=parsed.day,
                        hour=parsed.hour, freq='H')
        elif reso == 'minute':
            t1 = Period(year=parsed.year, month=parsed.month, day=parsed.day,
                        hour=parsed.hour, minute=parsed.minute, freq='T')
        elif reso == 'second':
            t1 = Period(year=parsed.year, month=parsed.month, day=parsed.day,
                        hour=parsed.hour, minute=parsed.minute, second=parsed.second,
                        freq='S')
        else:
            raise KeyError(key)
        return (t1.asfreq(self.freq, how='start'),
                t1.asfreq(self.freq, how='end'))

    def _get_string_slice(self, key):
        if not self.is_monotonic:
            raise ValueError('Partial indexing only valid for '
                             'ordered time series')

        key, parsed, reso = parse_time_string(key, self.freq)
        grp = frequencies.Resolution.get_freq_group(reso)
        freqn = frequencies.get_freq_group(self.freq)
        if reso in ['day', 'hour', 'minute', 'second'] and not grp < freqn:
            raise KeyError(key)

        t1, t2 = self._parsed_string_to_bounds(reso, parsed)
        return slice(self.searchsorted(t1.ordinal, side='left'),
                     self.searchsorted(t2.ordinal, side='right'))

    def _convert_tolerance(self, tolerance):
        tolerance = DatetimeIndexOpsMixin._convert_tolerance(self, tolerance)
        return self._maybe_convert_timedelta(tolerance)

    def join(self, other, how='left', level=None, return_indexers=False):
        """
        See Index.join
        """
        self._assert_can_do_setop(other)

        result = Int64Index.join(self, other, how=how, level=level,
                                 return_indexers=return_indexers)

        if return_indexers:
            result, lidx, ridx = result
            return self._apply_meta(result), lidx, ridx
        return self._apply_meta(result)

    def _assert_can_do_setop(self, other):
        super(PeriodIndex, self)._assert_can_do_setop(other)

        if not isinstance(other, PeriodIndex):
            raise ValueError('can only call with other PeriodIndex-ed objects')

        if self.freq != other.freq:
            msg = _DIFFERENT_FREQ_ERROR.format(self.freqstr, other.freqstr)
            raise ValueError(msg)

    def _wrap_union_result(self, other, result):
        name = self.name if self.name == other.name else None
        result = self._apply_meta(result)
        result.name = name
        return result

    def _apply_meta(self, rawarr):
        if not isinstance(rawarr, PeriodIndex):
            rawarr = PeriodIndex(rawarr, freq=self.freq)
        return rawarr

    def __getitem__(self, key):
        getitem = self._data.__getitem__
        if np.isscalar(key):
            val = getitem(key)
            return Period(ordinal=val, freq=self.freq)
        else:
            if com.is_bool_indexer(key):
                key = np.asarray(key)

            result = getitem(key)
            if result.ndim > 1:
                # MPL kludge
                # values = np.asarray(list(values), dtype=object)
                # return values.reshape(result.shape)

                return PeriodIndex(result, name=self.name, freq=self.freq)

            return PeriodIndex(result, name=self.name, freq=self.freq)

    def _format_native_types(self, na_rep=u('NaT'), date_format=None, **kwargs):

        values = np.array(list(self), dtype=object)
        mask = isnull(self.values)
        values[mask] = na_rep
        imask = ~mask

        if date_format:
            formatter = lambda dt: dt.strftime(date_format)
        else:
            formatter = lambda dt: u('%s') % dt
        values[imask] = np.array([formatter(dt) for dt in values[imask]])
        return values

    def take(self, indices, axis=0):
        """
        Analogous to ndarray.take
        """
        indices = com._ensure_platform_int(indices)
        taken = self.asi8.take(indices, axis=axis)
        return self._simple_new(taken, self.name, freq=self.freq)

    def append(self, other):
        """
        Append a collection of Index options together

        Parameters
        ----------
        other : Index or list/tuple of indices

        Returns
        -------
        appended : Index
        """
        name = self.name
        to_concat = [self]

        if isinstance(other, (list, tuple)):
            to_concat = to_concat + list(other)
        else:
            to_concat.append(other)

        for obj in to_concat:
            if isinstance(obj, Index) and obj.name != name:
                name = None
                break

        to_concat = self._ensure_compat_concat(to_concat)

        if isinstance(to_concat[0], PeriodIndex):
            if len(set([x.freq for x in to_concat])) > 1:
                # box
                to_concat = [x.asobject.values for x in to_concat]
            else:
                cat_values = np.concatenate([x.values for x in to_concat])
                return PeriodIndex(cat_values, freq=self.freq, name=name)

        to_concat = [x.values if isinstance(x, Index) else x
                     for x in to_concat]
        return Index(com._concat_compat(to_concat), name=name)

    def repeat(self, n):
        """
        Return a new Index of the values repeated n times.

        See also
        --------
        numpy.ndarray.repeat
        """
        # overwrites method from DatetimeIndexOpsMixin
        return self._shallow_copy(self.values.repeat(n))

    def __setstate__(self, state):
        """Necessary for making this object picklable"""

        if isinstance(state, dict):
            super(PeriodIndex, self).__setstate__(state)

        elif isinstance(state, tuple):

            # < 0.15 compat
            if len(state) == 2:
                nd_state, own_state = state
                data = np.empty(nd_state[1], dtype=nd_state[2])
                np.ndarray.__setstate__(data, nd_state)

                # backcompat
                self.freq = Period._maybe_convert_freq(own_state[1])

            else:  # pragma: no cover
                data = np.empty(state)
                np.ndarray.__setstate__(self, state)

            self._data = data

        else:
            raise Exception("invalid pickle state")

    _unpickle_compat = __setstate__

    def tz_convert(self, tz):
        """
        Convert tz-aware DatetimeIndex from one time zone to another (using pytz/dateutil)

        Parameters
        ----------
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time. Corresponding timestamps would be converted to
            time zone of the TimeSeries.
            None will remove timezone holding UTC time.

        Returns
        -------
        normalized : DatetimeIndex

        Note
        ----
        Not currently implemented for PeriodIndex
        """
        raise NotImplementedError("Not yet implemented for PeriodIndex")

    def tz_localize(self, tz, infer_dst=False):
        """
        Localize tz-naive DatetimeIndex to given time zone (using pytz/dateutil),
        or remove timezone from tz-aware DatetimeIndex

        Parameters
        ----------
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time. Corresponding timestamps would be converted to
            time zone of the TimeSeries.
            None will remove timezone holding local time.
        infer_dst : boolean, default False
            Attempt to infer fall dst-transition hours based on order

        Returns
        -------
        localized : DatetimeIndex

        Note
        ----
        Not currently implemented for PeriodIndex
        """
        raise NotImplementedError("Not yet implemented for PeriodIndex")


PeriodIndex._add_numeric_methods_disabled()
PeriodIndex._add_logical_methods_disabled()
PeriodIndex._add_datetimelike_methods()


def _get_ordinal_range(start, end, periods, freq, mult=1):
    if com._count_not_none(start, end, periods) < 2:
        raise ValueError('Must specify 2 of start, end, periods')

    if freq is not None:
        _, mult = _gfc(freq)

    if start is not None:
        start = Period(start, freq)
    if end is not None:
        end = Period(end, freq)

    is_start_per = isinstance(start, Period)
    is_end_per = isinstance(end, Period)

    if is_start_per and is_end_per and start.freq != end.freq:
        raise ValueError('Start and end must have same freq')
    if ((is_start_per and start.ordinal == tslib.iNaT) or
        (is_end_per and end.ordinal == tslib.iNaT)):
        raise ValueError('Start and end must not be NaT')

    if freq is None:
        if is_start_per:
            freq = start.freq
        elif is_end_per:
            freq = end.freq
        else:  # pragma: no cover
            raise ValueError('Could not infer freq from start/end')

    if periods is not None:
        periods = periods * mult
        if start is None:
            data = np.arange(end.ordinal - periods + mult,
                             end.ordinal + 1, mult,
                             dtype=np.int64)
        else:
            data = np.arange(start.ordinal, start.ordinal + periods, mult,
                             dtype=np.int64)
    else:
        data = np.arange(start.ordinal, end.ordinal + 1, mult, dtype=np.int64)

    return data, freq


def _range_from_fields(year=None, month=None, quarter=None, day=None,
                       hour=None, minute=None, second=None, freq=None):
    if hour is None:
        hour = 0
    if minute is None:
        minute = 0
    if second is None:
        second = 0
    if day is None:
        day = 1

    ordinals = []

    if quarter is not None:
        if freq is None:
            freq = 'Q'
            base = frequencies.FreqGroup.FR_QTR
        else:
            base, mult = _gfc(freq)
            if base != frequencies.FreqGroup.FR_QTR:
                raise AssertionError("base must equal FR_QTR")

        year, quarter = _make_field_arrays(year, quarter)
        for y, q in zip(year, quarter):
            y, m = _quarter_to_myear(y, q, freq)
            val = period.period_ordinal(y, m, 1, 1, 1, 1, 0, 0, base)
            ordinals.append(val)
    else:
        base, mult = _gfc(freq)
        arrays = _make_field_arrays(year, month, day, hour, minute, second)
        for y, mth, d, h, mn, s in zip(*arrays):
            ordinals.append(period.period_ordinal(y, mth, d, h, mn, s, 0, 0, base))

    return np.array(ordinals, dtype=np.int64), freq


def _make_field_arrays(*fields):
    length = None
    for x in fields:
        if isinstance(x, (list, np.ndarray, ABCSeries)):
            if length is not None and len(x) != length:
                raise ValueError('Mismatched Period array lengths')
            elif length is None:
                length = len(x)

    arrays = [np.asarray(x) if isinstance(x, (np.ndarray, list, ABCSeries))
              else np.repeat(x, length) for x in fields]

    return arrays


def pnow(freq=None):
    return Period(datetime.now(), freq=freq)


def period_range(start=None, end=None, periods=None, freq='D', name=None):
    """
    Return a fixed frequency datetime index, with day (calendar) as the default
    frequency


    Parameters
    ----------
    start : starting value, period-like, optional
    end : ending value, period-like, optional
    periods : int, default None
        Number of periods in the index
    freq : str/DateOffset, default 'D'
        Frequency alias
    name : str, default None
        Name for the resulting PeriodIndex

    Returns
    -------
    prng : PeriodIndex
    """
    return PeriodIndex(start=start, end=end, periods=periods,
                       freq=freq, name=name)
