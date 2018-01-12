"""Extension array for Period data
"""
import numpy as np

from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.missing import isna
from pandas.core.dtypes.dtypes import PeriodDtype
from pandas.core import common as com
from pandas.core.extensions import ExtensionArray
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin
from pandas._libs import tslib, iNaT
from pandas._libs.lib import infer_dtype
from pandas._libs.tslibs import period
from pandas._libs.tslibs.period import (
    IncompatibleFrequency,
    Period,
    _quarter_to_myear,
    _validate_end_alias,
)
from pandas.core.dtypes.common import (
    is_datetime64_dtype,
    is_float,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_object_dtype,
    is_period_dtype,
    is_scalar,
    pandas_dtype,
    _ensure_object,
)
import pandas.tseries.frequencies as frequencies
from pandas.tseries.frequencies import get_freq_code as _gfc


def dt64arr_to_periodarr(data, freq, tz):
    # TODO: the reverse is in period. move there?
    if data.dtype != np.dtype('M8[ns]'):
        raise ValueError('Wrong dtype: %s' % data.dtype)

    freq = Period._maybe_convert_freq(freq)
    base, mult = _gfc(freq)
    return period.dt64arr_to_periodarr(data.view('i8'), base, tz)


def to_period(data):
    data = np.asanyarray(data)
    if data.dtype != int:
        raise ValueError(data.dtype)

    return data


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
            ordinals.append(period.period_ordinal(
                y, mth, d, h, mn, s, 0, 0, base))

    return np.array(ordinals, dtype=np.int64), freq


def _get_ordinal_range(start, end, periods, freq, mult=1):
    if com._count_not_none(start, end, periods) != 2:
        raise ValueError('Of the three parameters: start, end, and periods, '
                         'exactly two must be specified')

    if freq is not None:
        _, mult = _gfc(freq)

    if start is not None:
        start = Period(start, freq)
    if end is not None:
        end = Period(end, freq)

    is_start_per = isinstance(start, Period)
    is_end_per = isinstance(end, Period)

    if is_start_per and is_end_per and start.freq != end.freq:
        raise ValueError('start and end must have same freq')
    if (start is tslib.NaT or end is tslib.NaT):
        raise ValueError('start and end must not be NaT')

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


# XXX: We inherit from DatetimeIndexOpsMixin to get comparison, arithmetics
# This should be split into an DatetimeArrayOpsMixin, and then any Index
# version that just does index-stuff


class PeriodArray(DatetimeIndexOpsMixin, ExtensionArray):
    dtype = PeriodDtype()
    ndim = 1
    can_hold_na = True
    _dtype = None

    def __new__(cls, data=None, ordinal=None, freq=None, start=None, end=None,
                periods=None, copy=False, name=None, tz=None, dtype=None,
                **kwargs):
        from pandas.core.indexes.datetimes import DatetimeIndex
        from pandas.core.indexes.numeric import Int64Index
        from pandas.core.indexes.period import PeriodIndex

        if periods is not None:
            if is_float(periods):
                periods = int(periods)
            elif not is_integer(periods):
                msg = 'periods must be a number, got {periods}'
                raise TypeError(msg.format(periods=periods))

        if dtype is not None:
            dtype = pandas_dtype(dtype)
            if not is_period_dtype(dtype):
                raise ValueError('dtype must be PeriodDtype')
            if freq is None:
                freq = dtype.freq
            elif freq != dtype.freq:
                msg = 'specified freq and dtype are different'
                raise IncompatibleFrequency(msg)

        # coerce freq to freq object, otherwise it can be coerced elementwise
        # which is slow
        if freq:
            freq = Period._maybe_convert_freq(freq)

        # TODO: handle copy
        if data is None:
            if ordinal is not None:
                data = np.asarray(ordinal, dtype=np.int64)
            else:
                data, freq = cls._generate_range(start, end, periods,
                                                 freq, kwargs)
            return cls._from_ordinals(data, freq=freq)

        if isinstance(data, PeriodIndex):
            data = data._data

        if isinstance(data, cls):
            if freq is None or freq == data.freq:  # no freq change
                freq = data.freq
                data = data._data
            else:
                base1, _ = _gfc(data.freq)
                base2, _ = _gfc(freq)
                data = period.period_asfreq_arr(data._data,
                                                base1, base2, 1)
            return cls._simple_new(data, freq=freq)

        # not array / index
        if not isinstance(data, (np.ndarray, PeriodIndex,
                                 DatetimeIndex, Int64Index)):
            if is_scalar(data) or isinstance(data, Period):
                cls._scalar_data_error(data)

            # other iterable of some kind
            if not isinstance(data, (list, tuple)):
                data = list(data)

            data = np.asarray(data)

        # datetime other than period
        if is_datetime64_dtype(data.dtype):
            data = dt64arr_to_periodarr(data, freq, tz)
            return cls._from_ordinals(data, freq=freq)

        # check not floats
        if infer_dtype(data) == 'floating' and len(data) > 0:
            raise TypeError("PeriodIndex does not allow "
                            "floating point in construction")

        # anything else, likely an array of strings or periods
        data = _ensure_object(data)
        freq = freq or period.extract_freq(data)
        data = period.extract_ordinals(data, freq)
        return cls._from_ordinals(data, freq=freq)

    @classmethod
    def _generate_range(cls, start, end, periods, freq, fields):
        if freq is not None:
            freq = Period._maybe_convert_freq(freq)

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
    def _simple_new(cls, values, freq=None):
        """
        Values can be any type that can be coerced to Periods.
        Ordinals in an ndarray are fastpath-ed to `_from_ordinals`
        """
        if not is_integer_dtype(values):
            values = np.array(values, copy=False)
            if len(values) > 0 and is_float_dtype(values):
                raise TypeError("PeriodArray can't take floats")
            return cls(values, freq=freq)

        return cls._from_ordinals(values, freq)

    @classmethod
    def _from_ordinals(cls, values, freq=None):
        """
        Values should be int ordinals
        `__new__` & `_simple_new` cooerce to ordinals and call this method
        """

        values = np.array(values, dtype='int64', copy=False)

        result = object.__new__(cls)
        result._data = values
        if freq is None:
            raise ValueError('freq is not specified and cannot be inferred')
        result.freq = Period._maybe_convert_freq(freq)
        return result

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        values = self._format_values()
        return "PeriodArray({}, freq={}, dtype={})".format(
            values, self.freq, self.dtype
        )

    def __getitem__(self, item):
        if is_scalar(item):
            return self._box_func(self._data[item])
        else:
            values = self._data[item]
            return self._simple_new(values, self.freq)

    @property
    def dtype(self):
        if self._dtype is None:
            self._dtype = PeriodDtype(self.freq)
        return self._dtype

    @property
    def shape(self):
        return (len(self),)

    @property
    def values(self):
        return self.astype(object)

    @property
    def asi8(self):
        return self._data.view('i8')

    @property
    def _box_func(self):
        return lambda x: Period._from_ordinal(ordinal=x, freq=self.freq)

    def _format_values(self):
        return np.array(['%s' % x for x in self.values], dtype='object')

    def formatting_values(self):
        return self._format_values()

    def astype(self, dtype, copy=True, how='start'):
        dtype = pandas_dtype(dtype)
        if is_object_dtype(dtype):
            # TODO: perf
            return np.array([Period._from_ordinal(p, self.freq)
                             for p in self], dtype='object')
        else:
            raise ValueError('invalid dtype')

    def copy(self):
        return self._from_ordinals(self._data.copy(), freq=self.freq)

    def isna(self):
        return self.asi8 == iNaT

    def nbytes(self):
        return self._data.nbytes

    def take(self, indexer, allow_fill=True, fill_value=None):
        # XXX: is take supposed to be a view?
        return self._from_ordinals(self._data.take(indexer), self.freq)

    take_nd = take

    @classmethod
    def concat_same_type(cls, to_concat):
        dtype = to_concat[0].dtype
        if not all(other.dtype == dtype for other in to_concat):
            raise TypeError("All frequencies must match")
        values = np.concatenate([other._data for other in to_concat])
        return cls._from_ordinals(values, freq=to_concat[0].freq)

    def get_values(self):
        return self._data

    @classmethod
    def _scalar_data_error(cls, data):
        # TODO: array-mixin
        raise TypeError('{0}(...) must be called with a collection of some '
                        'kind, {1} was passed'.format(cls.__name__,
                                                      repr(data)))

    def _get_attributes_dict(self):
        # TODO: from indexes.base, needed for ops, can remove
        return {}

    def view(self, cls=None):
        return self._data.view(cls)

    def equals(self, other):
        if not isinstance(other, type(self)):
            return False
        return (self.freq == other.freq and
                len(self) == len(other) and
                np.all(self._data == other._data))

    def slice(self, slicer):
        return self._from_ordinals(self._data[slicer], freq=self.freq)

    def asfreq(self, freq=None, how='E'):
        """
        Convert the PeriodArray to the specified frequency `freq`.

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

        new : PeriodArray with the new frequency

        Examples
        --------
        >>> pidx = pd.period_range('2010-01-01', '2015-01-01', freq='A')
        >>> pidx
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010, ..., 2015]
        Length: 6, Freq: A-DEC

        >>> pidx.asfreq('M')
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010-12, ..., 2015-12]
        Length: 6, Freq: M

        >>> pidx.asfreq('M', how='S')
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010-01, ..., 2015-01]
        Length: 6, Freq: M
        """
        how = _validate_end_alias(how)

        freq = Period._maybe_convert_freq(freq)

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

        # XXX: PeriodIndex could cache this. We can't, so this will be slower.
        mask = self.isna()
        if isna(self).any():
            new_data[mask] = tslib.iNaT

        return self._from_ordinals(new_data, freq=freq)

    # Pickling
    def __getnewargs__(self):
        # values, oridinal, freq
        return (None, self._data, self.freq)

    def __getstate__(self):
        return {'ordinal': self._data, 'freq': self.freq}

    def __setstate__(self, state):
        self.__dict__.update(state)


PeriodArray._add_datetimelike_methods()
