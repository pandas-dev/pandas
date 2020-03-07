from datetime import date, datetime, time, timedelta, tzinfo
import operator
from typing import Optional
import warnings

import numpy as np

from pandas._libs import NaT, Period, Timestamp, index as libindex, lib, tslib as libts
from pandas._libs.tslibs import fields, parsing, timezones
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import _NS_DTYPE, is_float, is_integer, is_scalar
from pandas.core.dtypes.dtypes import DatetimeTZDtype
from pandas.core.dtypes.missing import is_valid_nat_for_dtype

from pandas.core.arrays.datetimes import (
    DatetimeArray,
    tz_to_dtype,
    validate_tz_from_dtype,
)
import pandas.core.common as com
from pandas.core.indexes.base import Index, InvalidIndexError, maybe_extract_name
from pandas.core.indexes.datetimelike import DatetimeTimedeltaMixin
from pandas.core.indexes.extension import inherit_names
import pandas.core.tools.datetimes as tools

from pandas.tseries.frequencies import Resolution, to_offset
from pandas.tseries.offsets import prefix_mapping


def _new_DatetimeIndex(cls, d):
    """
    This is called upon unpickling, rather than the default which doesn't
    have arguments and breaks __new__
    """
    if "data" in d and not isinstance(d["data"], DatetimeIndex):
        # Avoid need to verify integrity by calling simple_new directly
        data = d.pop("data")
        result = cls._simple_new(data, **d)
    else:
        with warnings.catch_warnings():
            # TODO: If we knew what was going in to **d, we might be able to
            #  go through _simple_new instead
            warnings.simplefilter("ignore")
            result = cls.__new__(cls, **d)

    return result


@inherit_names(
    ["to_period", "to_perioddelta", "to_julian_date", "strftime"]
    + DatetimeArray._field_ops
    + DatetimeArray._datetimelike_methods,
    DatetimeArray,
    wrap=True,
)
@inherit_names(["_timezone", "is_normalized", "_resolution"], DatetimeArray, cache=True)
@inherit_names(
    [
        "_bool_ops",
        "_object_ops",
        "_field_ops",
        "_datetimelike_ops",
        "_datetimelike_methods",
        "tz",
        "tzinfo",
        "dtype",
        "to_pydatetime",
        "_local_timestamps",
        "_has_same_tz",
        "_format_native_types",
        "date",
        "time",
        "timetz",
    ]
    + DatetimeArray._bool_ops,
    DatetimeArray,
)
class DatetimeIndex(DatetimeTimedeltaMixin):
    """
    Immutable ndarray of datetime64 data, represented internally as int64, and
    which can be boxed to Timestamp objects that are subclasses of datetime and
    carry metadata such as frequency information.

    Parameters
    ----------
    data : array-like (1-dimensional), optional
        Optional datetime-like data to construct index with.
    copy : bool
        Make a copy of input ndarray.
    freq : str or pandas offset object, optional
        One of pandas date offset strings or corresponding objects. The string
        'infer' can be passed in order to set the frequency of the index as the
        inferred frequency upon creation.
    tz : pytz.timezone or dateutil.tz.tzfile
    ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
        When clocks moved backward due to DST, ambiguous times may arise.
        For example in Central European Time (UTC+01), when going from 03:00
        DST to 02:00 non-DST, 02:30:00 local time occurs both at 00:30:00 UTC
        and at 01:30:00 UTC. In such a situation, the `ambiguous` parameter
        dictates how ambiguous times should be handled.

        - 'infer' will attempt to infer fall dst-transition hours based on
          order
        - bool-ndarray where True signifies a DST time, False signifies a
          non-DST time (note that this flag is only applicable for ambiguous
          times)
        - 'NaT' will return NaT where there are ambiguous times
        - 'raise' will raise an AmbiguousTimeError if there are ambiguous times.
    name : object
        Name to be stored in the index.
    dayfirst : bool, default False
        If True, parse dates in `data` with the day first order.
    yearfirst : bool, default False
        If True parse dates in `data` with the year first order.

    Attributes
    ----------
    year
    month
    day
    hour
    minute
    second
    microsecond
    nanosecond
    date
    time
    timetz
    dayofyear
    weekofyear
    week
    dayofweek
    weekday
    quarter
    tz
    freq
    freqstr
    is_month_start
    is_month_end
    is_quarter_start
    is_quarter_end
    is_year_start
    is_year_end
    is_leap_year
    inferred_freq

    Methods
    -------
    normalize
    strftime
    snap
    tz_convert
    tz_localize
    round
    floor
    ceil
    to_period
    to_perioddelta
    to_pydatetime
    to_series
    to_frame
    month_name
    day_name
    mean

    See Also
    --------
    Index : The base pandas Index type.
    TimedeltaIndex : Index of timedelta64 data.
    PeriodIndex : Index of Period data.
    to_datetime : Convert argument to datetime.
    date_range : Create a fixed-frequency DatetimeIndex.

    Notes
    -----
    To learn more about the frequency strings, please see `this link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__.
    """

    _typ = "datetimeindex"

    _engine_type = libindex.DatetimeEngine
    _supports_partial_string_indexing = True

    _comparables = ["name", "freqstr", "tz"]
    _attributes = ["name", "tz", "freq"]

    _is_numeric_dtype = False
    _infer_as_myclass = True

    _data: DatetimeArray
    tz: Optional[tzinfo]

    # --------------------------------------------------------------------
    # Constructors

    def __new__(
        cls,
        data=None,
        freq=None,
        tz=None,
        normalize=False,
        closed=None,
        ambiguous="raise",
        dayfirst=False,
        yearfirst=False,
        dtype=None,
        copy=False,
        name=None,
    ):

        if is_scalar(data):
            raise TypeError(
                f"{cls.__name__}() must be called with a "
                f"collection of some kind, {repr(data)} was passed"
            )

        # - Cases checked above all return/raise before reaching here - #

        name = maybe_extract_name(name, data, cls)

        dtarr = DatetimeArray._from_sequence(
            data,
            dtype=dtype,
            copy=copy,
            tz=tz,
            freq=freq,
            dayfirst=dayfirst,
            yearfirst=yearfirst,
            ambiguous=ambiguous,
        )

        subarr = cls._simple_new(dtarr, name=name)
        return subarr

    @classmethod
    def _simple_new(cls, values, name=None, freq=None, tz=None, dtype=None):
        """
        We require the we have a dtype compat for the values
        if we are passed a non-dtype compat, then coerce using the constructor
        """
        if isinstance(values, DatetimeArray):
            if tz:
                tz = validate_tz_from_dtype(dtype, tz)
                dtype = DatetimeTZDtype(tz=tz)
            elif dtype is None:
                dtype = _NS_DTYPE

            values = DatetimeArray(values, freq=freq, dtype=dtype)
            tz = values.tz
            freq = values.freq
            values = values._data

        dtype = tz_to_dtype(tz)
        dtarr = DatetimeArray._simple_new(values, freq=freq, dtype=dtype)
        assert isinstance(dtarr, DatetimeArray)

        result = object.__new__(cls)
        result._data = dtarr
        result.name = name
        result._no_setting_name = False
        # For groupby perf. See note in indexes/base about _index_data
        result._index_data = dtarr._data
        result._reset_identity()
        return result

    # --------------------------------------------------------------------

    @cache_readonly
    def _is_dates_only(self) -> bool:
        """
        Return a boolean if we are only dates (and don't have a timezone)

        Returns
        -------
        bool
        """
        from pandas.io.formats.format import _is_dates_only

        return _is_dates_only(self.values) and self.tz is None

    def __reduce__(self):

        # we use a special reduce here because we need
        # to simply set the .tz (and not reinterpret it)

        d = dict(data=self._data)
        d.update(self._get_attributes_dict())
        return _new_DatetimeIndex, (type(self), d), None

    def _convert_for_op(self, value):
        """
        Convert value to be insertable to ndarray.
        """
        if self._has_same_tz(value):
            return Timestamp(value).asm8
        raise ValueError("Passed item and index have different timezone")

    # --------------------------------------------------------------------
    # Rendering Methods

    def _mpl_repr(self):
        # how to represent ourselves to matplotlib
        return libts.ints_to_pydatetime(self.asi8, self.tz)

    @property
    def _formatter_func(self):
        from pandas.io.formats.format import _get_format_datetime64

        formatter = _get_format_datetime64(is_dates_only=self._is_dates_only)
        return lambda x: f"'{formatter(x, tz=self.tz)}'"

    # --------------------------------------------------------------------
    # Set Operation Methods

    def union_many(self, others):
        """
        A bit of a hack to accelerate unioning a collection of indexes.
        """
        this = self

        for other in others:
            if not isinstance(this, DatetimeIndex):
                this = Index.union(this, other)
                continue

            if not isinstance(other, DatetimeIndex):
                try:
                    other = DatetimeIndex(other)
                except TypeError:
                    pass

            this, other = this._maybe_utc_convert(other)

            if this._can_fast_union(other):
                this = this._fast_union(other)
            else:
                this = Index.union(this, other)
        return this

    # --------------------------------------------------------------------

    def _get_time_micros(self):
        values = self.asi8
        if self.tz is not None and not timezones.is_utc(self.tz):
            values = self._data._local_timestamps()
        return fields.get_time_micros(values)

    def to_series(self, keep_tz=lib.no_default, index=None, name=None):
        """
        Create a Series with both index and values equal to the index keys
        useful with map for returning an indexer based on an index.

        Parameters
        ----------
        keep_tz : optional, defaults True
            Return the data keeping the timezone.

            If keep_tz is True:

              If the timezone is not set, the resulting
              Series will have a datetime64[ns] dtype.

              Otherwise the Series will have an datetime64[ns, tz] dtype; the
              tz will be preserved.

            If keep_tz is False:

              Series will have a datetime64[ns] dtype. TZ aware
              objects will have the tz removed.

            .. versionchanged:: 1.0.0
                The default value is now True.  In a future version,
                this keyword will be removed entirely.  Stop passing the
                argument to obtain the future behavior and silence the warning.

        index : Index, optional
            Index of resulting Series. If None, defaults to original index.
        name : str, optional
            Name of resulting Series. If None, defaults to name of original
            index.

        Returns
        -------
        Series
        """
        from pandas import Series

        if index is None:
            index = self._shallow_copy()
        if name is None:
            name = self.name

        if keep_tz is not lib.no_default:
            if keep_tz:
                warnings.warn(
                    "The 'keep_tz' keyword in DatetimeIndex.to_series "
                    "is deprecated and will be removed in a future version.  "
                    "You can stop passing 'keep_tz' to silence this warning.",
                    FutureWarning,
                    stacklevel=2,
                )
            else:
                warnings.warn(
                    "Specifying 'keep_tz=False' is deprecated and this "
                    "option will be removed in a future release. If "
                    "you want to remove the timezone information, you "
                    "can do 'idx.tz_convert(None)' before calling "
                    "'to_series'.",
                    FutureWarning,
                    stacklevel=2,
                )
        else:
            keep_tz = True

        if keep_tz and self.tz is not None:
            # preserve the tz & copy
            values = self.copy(deep=True)
        else:
            values = self.values.copy()

        return Series(values, index=index, name=name)

    def snap(self, freq="S"):
        """
        Snap time stamps to nearest occurring frequency.

        Returns
        -------
        DatetimeIndex
        """
        # Superdumb, punting on any optimizing
        freq = to_offset(freq)

        snapped = np.empty(len(self), dtype=_NS_DTYPE)

        for i, v in enumerate(self):
            s = v
            if not freq.is_on_offset(s):
                t0 = freq.rollback(s)
                t1 = freq.rollforward(s)
                if abs(s - t0) < abs(t1 - s):
                    s = t0
                else:
                    s = t1
            snapped[i] = s

        dta = DatetimeArray(snapped, dtype=self.dtype)
        return DatetimeIndex._simple_new(dta, name=self.name)

    def _parsed_string_to_bounds(self, reso: str, parsed: datetime):
        """
        Calculate datetime bounds for parsed time string and its resolution.

        Parameters
        ----------
        reso : str
            Resolution provided by parsed string.
        parsed : datetime
            Datetime from parsed string.

        Returns
        -------
        lower, upper: pd.Timestamp
        """
        valid_resos = {
            "year",
            "month",
            "quarter",
            "day",
            "hour",
            "minute",
            "second",
            "minute",
            "second",
            "microsecond",
        }
        if reso not in valid_resos:
            raise KeyError

        grp = Resolution.get_freq_group(reso)
        per = Period(parsed, freq=(grp, 1))
        start, end = per.start_time, per.end_time

        # GH 24076
        # If an incoming date string contained a UTC offset, need to localize
        # the parsed date to this offset first before aligning with the index's
        # timezone
        if parsed.tzinfo is not None:
            if self.tz is None:
                raise ValueError(
                    "The index must be timezone aware when indexing "
                    "with a date string with a UTC offset"
                )
            start = start.tz_localize(parsed.tzinfo).tz_convert(self.tz)
            end = end.tz_localize(parsed.tzinfo).tz_convert(self.tz)
        elif self.tz is not None:
            start = start.tz_localize(self.tz)
            end = end.tz_localize(self.tz)
        return start, end

    def _validate_partial_date_slice(self, reso: str):
        if (
            self.is_monotonic
            and reso in ["day", "hour", "minute", "second"]
            and self._resolution >= Resolution.get_reso(reso)
        ):
            # These resolution/monotonicity validations came from GH3931,
            # GH3452 and GH2369.

            # See also GH14826
            raise KeyError

        if reso == "microsecond":
            # _partial_date_slice doesn't allow microsecond resolution, but
            # _parsed_string_to_bounds allows it.
            raise KeyError

    def _maybe_promote(self, other):
        if other.inferred_type == "date":
            other = DatetimeIndex(other)
        return self, other

    def get_loc(self, key, method=None, tolerance=None):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int
        """
        if not is_scalar(key):
            raise InvalidIndexError(key)

        orig_key = key
        if is_valid_nat_for_dtype(key, self.dtype):
            key = NaT

        if isinstance(key, self._data._recognized_scalars):
            # needed to localize naive datetimes
            key = self._maybe_cast_for_get_loc(key)

        elif isinstance(key, str):
            try:
                return self._get_string_slice(key)
            except (TypeError, KeyError, ValueError, OverflowError):
                pass

            try:
                key = self._maybe_cast_for_get_loc(key)
            except ValueError as err:
                raise KeyError(key) from err

        elif isinstance(key, timedelta):
            # GH#20464
            raise TypeError(
                f"Cannot index {type(self).__name__} with {type(key).__name__}"
            )

        elif isinstance(key, time):
            if method is not None:
                raise NotImplementedError(
                    "cannot yet lookup inexact labels when key is a time object"
                )
            return self.indexer_at_time(key)

        else:
            # unrecognized type
            raise KeyError(key)

        try:
            return Index.get_loc(self, key, method, tolerance)
        except KeyError as err:
            raise KeyError(orig_key) from err

    def _maybe_cast_for_get_loc(self, key) -> Timestamp:
        # needed to localize naive datetimes
        key = Timestamp(key)
        if key.tzinfo is None:
            key = key.tz_localize(self.tz)
        else:
            key = key.tz_convert(self.tz)
        return key

    def _maybe_cast_slice_bound(self, label, side: str, kind):
        """
        If label is a string, cast it to datetime according to resolution.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : {'loc', 'getitem'} or None

        Returns
        -------
        label : object

        Notes
        -----
        Value of `side` parameter should be validated in caller.
        """
        assert kind in ["loc", "getitem", None]

        if is_float(label) or isinstance(label, time) or is_integer(label):
            self._invalid_indexer("slice", label)

        if isinstance(label, str):
            freq = getattr(self, "freqstr", getattr(self, "inferred_freq", None))
            parsed, reso = parsing.parse_time_string(label, freq)
            lower, upper = self._parsed_string_to_bounds(reso, parsed)
            # lower, upper form the half-open interval:
            #   [parsed, parsed + 1 freq)
            # because label may be passed to searchsorted
            # the bounds need swapped if index is reverse sorted and has a
            # length > 1 (is_monotonic_decreasing gives True for empty
            # and length 1 index)
            if self._is_strictly_monotonic_decreasing and len(self) > 1:
                return upper if side == "left" else lower
            return lower if side == "left" else upper
        else:
            return label

    def _get_string_slice(self, key: str, use_lhs: bool = True, use_rhs: bool = True):
        freq = getattr(self, "freqstr", getattr(self, "inferred_freq", None))
        parsed, reso = parsing.parse_time_string(key, freq)
        loc = self._partial_date_slice(reso, parsed, use_lhs=use_lhs, use_rhs=use_rhs)
        return loc

    def slice_indexer(self, start=None, end=None, step=None, kind=None):
        """
        Return indexer for specified label slice.
        Index.slice_indexer, customized to handle time slicing.

        In addition to functionality provided by Index.slice_indexer, does the
        following:

        - if both `start` and `end` are instances of `datetime.time`, it
          invokes `indexer_between_time`
        - if `start` and `end` are both either string or None perform
          value-based selection in non-monotonic cases.

        """
        # For historical reasons DatetimeIndex supports slices between two
        # instances of datetime.time as if it were applying a slice mask to
        # an array of (self.hour, self.minute, self.seconds, self.microsecond).
        if isinstance(start, time) and isinstance(end, time):
            if step is not None and step != 1:
                raise ValueError("Must have step size of 1 with time slices")
            return self.indexer_between_time(start, end)

        if isinstance(start, time) or isinstance(end, time):
            raise KeyError("Cannot mix time and non-time slice keys")

        # Pandas supports slicing with dates, treated as datetimes at midnight.
        # https://github.com/pandas-dev/pandas/issues/31501
        if isinstance(start, date) and not isinstance(start, datetime):
            start = datetime.combine(start, time(0, 0))
        if isinstance(end, date) and not isinstance(end, datetime):
            end = datetime.combine(end, time(0, 0))

        try:
            return Index.slice_indexer(self, start, end, step, kind=kind)
        except KeyError:
            # For historical reasons DatetimeIndex by default supports
            # value-based partial (aka string) slices on non-monotonic arrays,
            # let's try that.
            if (start is None or isinstance(start, str)) and (
                end is None or isinstance(end, str)
            ):
                mask = True
                if start is not None:
                    start_casted = self._maybe_cast_slice_bound(start, "left", kind)
                    mask = start_casted <= self

                if end is not None:
                    end_casted = self._maybe_cast_slice_bound(end, "right", kind)
                    mask = (self <= end_casted) & mask

                indexer = mask.nonzero()[0][::step]
                if len(indexer) == len(self):
                    return slice(None)
                else:
                    return indexer
            else:
                raise

    # --------------------------------------------------------------------

    def is_type_compatible(self, typ) -> bool:
        return typ == self.inferred_type or typ == "datetime"

    @property
    def inferred_type(self) -> str:
        # b/c datetime is represented as microseconds since the epoch, make
        # sure we can't have ambiguous indexing
        return "datetime64"

    def indexer_at_time(self, time, asof=False):
        """
        Return index locations of index values at particular time of day
        (e.g. 9:30AM).

        Parameters
        ----------
        time : datetime.time or str
            datetime.time or string in appropriate format ("%H:%M", "%H%M",
            "%I:%M%p", "%I%M%p", "%H:%M:%S", "%H%M%S", "%I:%M:%S%p",
            "%I%M%S%p").

        Returns
        -------
        values_at_time : array of integers

        See Also
        --------
        indexer_between_time, DataFrame.at_time
        """
        if asof:
            raise NotImplementedError("'asof' argument is not supported")

        if isinstance(time, str):
            from dateutil.parser import parse

            time = parse(time).time()

        if time.tzinfo:
            if self.tz is None:
                raise ValueError("Index must be timezone aware.")
            time_micros = self.tz_convert(time.tzinfo)._get_time_micros()
        else:
            time_micros = self._get_time_micros()
        micros = _time_to_micros(time)
        return (micros == time_micros).nonzero()[0]

    def indexer_between_time(
        self, start_time, end_time, include_start=True, include_end=True
    ):
        """
        Return index locations of values between particular times of day
        (e.g., 9:00-9:30AM).

        Parameters
        ----------
        start_time, end_time : datetime.time, str
            datetime.time or string in appropriate format ("%H:%M", "%H%M",
            "%I:%M%p", "%I%M%p", "%H:%M:%S", "%H%M%S", "%I:%M:%S%p",
            "%I%M%S%p").
        include_start : bool, default True
        include_end : bool, default True

        Returns
        -------
        values_between_time : array of integers

        See Also
        --------
        indexer_at_time, DataFrame.between_time
        """
        start_time = tools.to_time(start_time)
        end_time = tools.to_time(end_time)
        time_micros = self._get_time_micros()
        start_micros = _time_to_micros(start_time)
        end_micros = _time_to_micros(end_time)

        if include_start and include_end:
            lop = rop = operator.le
        elif include_start:
            lop = operator.le
            rop = operator.lt
        elif include_end:
            lop = operator.lt
            rop = operator.le
        else:
            lop = rop = operator.lt

        if start_time <= end_time:
            join_op = operator.and_
        else:
            join_op = operator.or_

        mask = join_op(lop(start_micros, time_micros), rop(time_micros, end_micros))

        return mask.nonzero()[0]


DatetimeIndex._add_numeric_methods_disabled()
DatetimeIndex._add_logical_methods_disabled()


def date_range(
    start=None,
    end=None,
    periods=None,
    freq=None,
    tz=None,
    normalize=False,
    name=None,
    closed=None,
    **kwargs,
) -> DatetimeIndex:
    """
    Return a fixed frequency DatetimeIndex.

    Parameters
    ----------
    start : str or datetime-like, optional
        Left bound for generating dates.
    end : str or datetime-like, optional
        Right bound for generating dates.
    periods : int, optional
        Number of periods to generate.
    freq : str or DateOffset, default 'D'
        Frequency strings can have multiples, e.g. '5H'. See
        :ref:`here <timeseries.offset_aliases>` for a list of
        frequency aliases.
    tz : str or tzinfo, optional
        Time zone name for returning localized DatetimeIndex, for example
        'Asia/Hong_Kong'. By default, the resulting DatetimeIndex is
        timezone-naive.
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range.
    name : str, default None
        Name of the resulting DatetimeIndex.
    closed : {None, 'left', 'right'}, optional
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None, the default).
    **kwargs
        For compatibility. Has no effect on the result.

    Returns
    -------
    rng : DatetimeIndex

    See Also
    --------
    DatetimeIndex : An immutable container for datetimes.
    timedelta_range : Return a fixed frequency TimedeltaIndex.
    period_range : Return a fixed frequency PeriodIndex.
    interval_range : Return a fixed frequency IntervalIndex.

    Notes
    -----
    Of the four parameters ``start``, ``end``, ``periods``, and ``freq``,
    exactly three must be specified. If ``freq`` is omitted, the resulting
    ``DatetimeIndex`` will have ``periods`` linearly spaced elements between
    ``start`` and ``end`` (closed on both sides).

    To learn more about the frequency strings, please see `this link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__.

    Examples
    --------
    **Specifying the values**

    The next four examples generate the same `DatetimeIndex`, but vary
    the combination of `start`, `end` and `periods`.

    Specify `start` and `end`, with the default daily frequency.

    >>> pd.date_range(start='1/1/2018', end='1/08/2018')
    DatetimeIndex(['2018-01-01', '2018-01-02', '2018-01-03', '2018-01-04',
                   '2018-01-05', '2018-01-06', '2018-01-07', '2018-01-08'],
                  dtype='datetime64[ns]', freq='D')

    Specify `start` and `periods`, the number of periods (days).

    >>> pd.date_range(start='1/1/2018', periods=8)
    DatetimeIndex(['2018-01-01', '2018-01-02', '2018-01-03', '2018-01-04',
                   '2018-01-05', '2018-01-06', '2018-01-07', '2018-01-08'],
                  dtype='datetime64[ns]', freq='D')

    Specify `end` and `periods`, the number of periods (days).

    >>> pd.date_range(end='1/1/2018', periods=8)
    DatetimeIndex(['2017-12-25', '2017-12-26', '2017-12-27', '2017-12-28',
                   '2017-12-29', '2017-12-30', '2017-12-31', '2018-01-01'],
                  dtype='datetime64[ns]', freq='D')

    Specify `start`, `end`, and `periods`; the frequency is generated
    automatically (linearly spaced).

    >>> pd.date_range(start='2018-04-24', end='2018-04-27', periods=3)
    DatetimeIndex(['2018-04-24 00:00:00', '2018-04-25 12:00:00',
                   '2018-04-27 00:00:00'],
                  dtype='datetime64[ns]', freq=None)

    **Other Parameters**

    Changed the `freq` (frequency) to ``'M'`` (month end frequency).

    >>> pd.date_range(start='1/1/2018', periods=5, freq='M')
    DatetimeIndex(['2018-01-31', '2018-02-28', '2018-03-31', '2018-04-30',
                   '2018-05-31'],
                  dtype='datetime64[ns]', freq='M')

    Multiples are allowed

    >>> pd.date_range(start='1/1/2018', periods=5, freq='3M')
    DatetimeIndex(['2018-01-31', '2018-04-30', '2018-07-31', '2018-10-31',
                   '2019-01-31'],
                  dtype='datetime64[ns]', freq='3M')

    `freq` can also be specified as an Offset object.

    >>> pd.date_range(start='1/1/2018', periods=5, freq=pd.offsets.MonthEnd(3))
    DatetimeIndex(['2018-01-31', '2018-04-30', '2018-07-31', '2018-10-31',
                   '2019-01-31'],
                  dtype='datetime64[ns]', freq='3M')

    Specify `tz` to set the timezone.

    >>> pd.date_range(start='1/1/2018', periods=5, tz='Asia/Tokyo')
    DatetimeIndex(['2018-01-01 00:00:00+09:00', '2018-01-02 00:00:00+09:00',
                   '2018-01-03 00:00:00+09:00', '2018-01-04 00:00:00+09:00',
                   '2018-01-05 00:00:00+09:00'],
                  dtype='datetime64[ns, Asia/Tokyo]', freq='D')

    `closed` controls whether to include `start` and `end` that are on the
    boundary. The default includes boundary points on either end.

    >>> pd.date_range(start='2017-01-01', end='2017-01-04', closed=None)
    DatetimeIndex(['2017-01-01', '2017-01-02', '2017-01-03', '2017-01-04'],
                  dtype='datetime64[ns]', freq='D')

    Use ``closed='left'`` to exclude `end` if it falls on the boundary.

    >>> pd.date_range(start='2017-01-01', end='2017-01-04', closed='left')
    DatetimeIndex(['2017-01-01', '2017-01-02', '2017-01-03'],
                  dtype='datetime64[ns]', freq='D')

    Use ``closed='right'`` to exclude `start` if it falls on the boundary.

    >>> pd.date_range(start='2017-01-01', end='2017-01-04', closed='right')
    DatetimeIndex(['2017-01-02', '2017-01-03', '2017-01-04'],
                  dtype='datetime64[ns]', freq='D')
    """
    if freq is None and com.any_none(periods, start, end):
        freq = "D"

    dtarr = DatetimeArray._generate_range(
        start=start,
        end=end,
        periods=periods,
        freq=freq,
        tz=tz,
        normalize=normalize,
        closed=closed,
        **kwargs,
    )
    return DatetimeIndex._simple_new(dtarr, name=name)


def bdate_range(
    start=None,
    end=None,
    periods=None,
    freq="B",
    tz=None,
    normalize=True,
    name=None,
    weekmask=None,
    holidays=None,
    closed=None,
    **kwargs,
) -> DatetimeIndex:
    """
    Return a fixed frequency DatetimeIndex, with business day as the default
    frequency.

    Parameters
    ----------
    start : str or datetime-like, default None
        Left bound for generating dates.
    end : str or datetime-like, default None
        Right bound for generating dates.
    periods : int, default None
        Number of periods to generate.
    freq : str or DateOffset, default 'B' (business daily)
        Frequency strings can have multiples, e.g. '5H'.
    tz : str or None
        Time zone name for returning localized DatetimeIndex, for example
        Asia/Beijing.
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range.
    name : str, default None
        Name of the resulting DatetimeIndex.
    weekmask : str or None, default None
        Weekmask of valid business days, passed to ``numpy.busdaycalendar``,
        only used when custom frequency strings are passed.  The default
        value None is equivalent to 'Mon Tue Wed Thu Fri'.

        .. versionadded:: 0.21.0

    holidays : list-like or None, default None
        Dates to exclude from the set of valid business days, passed to
        ``numpy.busdaycalendar``, only used when custom frequency strings
        are passed.

        .. versionadded:: 0.21.0

    closed : str, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None).
    **kwargs
        For compatibility. Has no effect on the result.

    Returns
    -------
    DatetimeIndex

    Notes
    -----
    Of the four parameters: ``start``, ``end``, ``periods``, and ``freq``,
    exactly three must be specified.  Specifying ``freq`` is a requirement
    for ``bdate_range``.  Use ``date_range`` if specifying ``freq`` is not
    desired.

    To learn more about the frequency strings, please see `this link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__.

    Examples
    --------
    Note how the two weekend days are skipped in the result.

    >>> pd.bdate_range(start='1/1/2018', end='1/08/2018')
    DatetimeIndex(['2018-01-01', '2018-01-02', '2018-01-03', '2018-01-04',
               '2018-01-05', '2018-01-08'],
              dtype='datetime64[ns]', freq='B')
    """
    if freq is None:
        msg = "freq must be specified for bdate_range; use date_range instead"
        raise TypeError(msg)

    if isinstance(freq, str) and freq.startswith("C"):
        try:
            weekmask = weekmask or "Mon Tue Wed Thu Fri"
            freq = prefix_mapping[freq](holidays=holidays, weekmask=weekmask)
        except (KeyError, TypeError) as err:
            msg = f"invalid custom frequency string: {freq}"
            raise ValueError(msg) from err
    elif holidays or weekmask:
        msg = (
            "a custom frequency string is required when holidays or "
            f"weekmask are passed, got frequency {freq}"
        )
        raise ValueError(msg)

    return date_range(
        start=start,
        end=end,
        periods=periods,
        freq=freq,
        tz=tz,
        normalize=normalize,
        name=name,
        closed=closed,
        **kwargs,
    )


def _time_to_micros(time):
    seconds = time.hour * 60 * 60 + 60 * time.minute + time.second
    return 1000000 * seconds + time.microsecond
