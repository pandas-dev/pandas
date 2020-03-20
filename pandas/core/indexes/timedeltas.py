""" implement the TimedeltaIndex """

from pandas._libs import NaT, Timedelta, index as libindex
from pandas._typing import DtypeObj, Label
from pandas.util._decorators import Appender

from pandas.core.dtypes.common import (
    _TD_DTYPE,
    is_float,
    is_integer,
    is_scalar,
    is_timedelta64_dtype,
    is_timedelta64_ns_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.missing import is_valid_nat_for_dtype

from pandas.core.arrays import datetimelike as dtl
from pandas.core.arrays.timedeltas import TimedeltaArray
import pandas.core.common as com
from pandas.core.indexes.base import Index, InvalidIndexError, maybe_extract_name
from pandas.core.indexes.datetimelike import (
    DatetimeIndexOpsMixin,
    DatetimeTimedeltaMixin,
)
from pandas.core.indexes.extension import inherit_names

from pandas.tseries.frequencies import to_offset


@inherit_names(
    [
        "_box_values",
        "__neg__",
        "__pos__",
        "__abs__",
        "total_seconds",
        "round",
        "floor",
        "ceil",
    ]
    + TimedeltaArray._field_ops,
    TimedeltaArray,
    wrap=True,
)
@inherit_names(
    [
        "_bool_ops",
        "_object_ops",
        "_field_ops",
        "_datetimelike_ops",
        "_datetimelike_methods",
        "_other_ops",
        "components",
        "to_pytimedelta",
        "sum",
        "std",
        "median",
        "_format_native_types",
        "freq",
    ],
    TimedeltaArray,
)
class TimedeltaIndex(DatetimeTimedeltaMixin, dtl.TimelikeOps):
    """
    Immutable ndarray of timedelta64 data, represented internally as int64, and
    which can be boxed to timedelta objects.

    Parameters
    ----------
    data  : array-like (1-dimensional), optional
        Optional timedelta-like data to construct index with.
    unit : unit of the arg (D,h,m,s,ms,us,ns) denote the unit, optional
        Which is an integer/float number.
    freq : str or pandas offset object, optional
        One of pandas date offset strings or corresponding objects. The string
        'infer' can be passed in order to set the frequency of the index as the
        inferred frequency upon creation.
    copy  : bool
        Make a copy of input ndarray.
    name : object
        Name to be stored in the index.

    Attributes
    ----------
    days
    seconds
    microseconds
    nanoseconds
    components
    inferred_freq

    Methods
    -------
    to_pytimedelta
    to_series
    round
    floor
    ceil
    to_frame
    mean

    See Also
    --------
    Index : The base pandas Index type.
    Timedelta : Represents a duration between two dates or times.
    DatetimeIndex : Index of datetime64 data.
    PeriodIndex : Index of Period data.
    timedelta_range : Create a fixed-frequency TimedeltaIndex.

    Notes
    -----
    To learn more about the frequency strings, please see `this link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__.
    """

    _typ = "timedeltaindex"

    _engine_type = libindex.TimedeltaEngine

    _comparables = ["name", "freq"]
    _attributes = ["name", "freq"]
    _is_numeric_dtype = True
    _infer_as_myclass = True

    _data: TimedeltaArray

    # -------------------------------------------------------------------
    # Constructors

    def __new__(
        cls,
        data=None,
        unit=None,
        freq=None,
        closed=None,
        dtype=_TD_DTYPE,
        copy=False,
        name=None,
    ):
        name = maybe_extract_name(name, data, cls)

        if is_scalar(data):
            raise TypeError(
                f"{cls.__name__}() must be called with a "
                f"collection of some kind, {repr(data)} was passed"
            )

        if unit in {"Y", "y", "M"}:
            raise ValueError(
                "Units 'M' and 'Y' are no longer supported, as they do not "
                "represent unambiguous timedelta values durations."
            )

        if isinstance(data, TimedeltaArray) and freq is None:
            if copy:
                data = data.copy()
            return cls._simple_new(data, name=name)

        if isinstance(data, TimedeltaIndex) and freq is None and name is None:
            if copy:
                return data.copy()
            else:
                return data._shallow_copy()

        # - Cases checked above all return/raise before reaching here - #

        tdarr = TimedeltaArray._from_sequence(
            data, freq=freq, unit=unit, dtype=dtype, copy=copy
        )
        return cls._simple_new(tdarr, name=name)

    @classmethod
    def _simple_new(cls, values: TimedeltaArray, name: Label = None):
        assert isinstance(values, TimedeltaArray)

        result = object.__new__(cls)
        result._data = values
        result._name = name
        result._cache = {}
        # For groupby perf. See note in indexes/base about _index_data
        result._index_data = values._data

        result._reset_identity()
        return result

    # -------------------------------------------------------------------
    # Rendering Methods

    @property
    def _formatter_func(self):
        from pandas.io.formats.format import _get_format_timedelta64

        return _get_format_timedelta64(self, box=True)

    # -------------------------------------------------------------------

    @Appender(Index.astype.__doc__)
    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)
        if is_timedelta64_dtype(dtype) and not is_timedelta64_ns_dtype(dtype):
            # Have to repeat the check for 'timedelta64' (not ns) dtype
            #  so that we can return a numeric index, since pandas will return
            #  a TimedeltaIndex when dtype='timedelta'
            result = self._data.astype(dtype, copy=copy)
            if self.hasnans:
                return Index(result, name=self.name)
            return Index(result.astype("i8"), name=self.name)
        return DatetimeIndexOpsMixin.astype(self, dtype, copy=copy)

    def _maybe_promote(self, other):
        if other.inferred_type == "timedelta":
            other = TimedeltaIndex(other)
        return self, other

    def _is_comparable_dtype(self, dtype: DtypeObj) -> bool:
        """
        Can we compare values of the given dtype to our own?
        """
        return is_timedelta64_dtype(dtype)

    def get_loc(self, key, method=None, tolerance=None):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int, slice, or ndarray[int]
        """
        if not is_scalar(key):
            raise InvalidIndexError(key)

        if is_valid_nat_for_dtype(key, self.dtype):
            key = NaT

        elif isinstance(key, str):
            try:
                key = Timedelta(key)
            except ValueError as err:
                raise KeyError(key) from err

        elif isinstance(key, self._data._recognized_scalars) or key is NaT:
            key = Timedelta(key)

        else:
            raise KeyError(key)

        return Index.get_loc(self, key, method, tolerance)

    def _maybe_cast_slice_bound(self, label, side: str, kind):
        """
        If label is a string, cast it to timedelta according to resolution.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : {'loc', 'getitem'} or None

        Returns
        -------
        label : object
        """
        assert kind in ["loc", "getitem", None]

        if isinstance(label, str):
            parsed = Timedelta(label)
            lbound = parsed.round(parsed.resolution_string)
            if side == "left":
                return lbound
            else:
                return lbound + to_offset(parsed.resolution_string) - Timedelta(1, "ns")
        elif is_integer(label) or is_float(label):
            self._invalid_indexer("slice", label)

        return label

    def is_type_compatible(self, typ) -> bool:
        return typ == self.inferred_type or typ == "timedelta"

    @property
    def inferred_type(self) -> str:
        return "timedelta64"


TimedeltaIndex._add_logical_methods_disabled()


def timedelta_range(
    start=None, end=None, periods=None, freq=None, name=None, closed=None
) -> TimedeltaIndex:
    """
    Return a fixed frequency TimedeltaIndex, with day as the default
    frequency.

    Parameters
    ----------
    start : str or timedelta-like, default None
        Left bound for generating timedeltas.
    end : str or timedelta-like, default None
        Right bound for generating timedeltas.
    periods : int, default None
        Number of periods to generate.
    freq : str or DateOffset, default 'D'
        Frequency strings can have multiples, e.g. '5H'.
    name : str, default None
        Name of the resulting TimedeltaIndex.
    closed : str, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None).

    Returns
    -------
    rng : TimedeltaIndex

    Notes
    -----
    Of the four parameters ``start``, ``end``, ``periods``, and ``freq``,
    exactly three must be specified. If ``freq`` is omitted, the resulting
    ``TimedeltaIndex`` will have ``periods`` linearly spaced elements between
    ``start`` and ``end`` (closed on both sides).

    To learn more about the frequency strings, please see `this link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__.

    Examples
    --------
    >>> pd.timedelta_range(start='1 day', periods=4)
    TimedeltaIndex(['1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[ns]', freq='D')

    The ``closed`` parameter specifies which endpoint is included.  The default
    behavior is to include both endpoints.

    >>> pd.timedelta_range(start='1 day', periods=4, closed='right')
    TimedeltaIndex(['2 days', '3 days', '4 days'],
                   dtype='timedelta64[ns]', freq='D')

    The ``freq`` parameter specifies the frequency of the TimedeltaIndex.
    Only fixed frequencies can be passed, non-fixed frequencies such as
    'M' (month end) will raise.

    >>> pd.timedelta_range(start='1 day', end='2 days', freq='6H')
    TimedeltaIndex(['1 days 00:00:00', '1 days 06:00:00', '1 days 12:00:00',
                    '1 days 18:00:00', '2 days 00:00:00'],
                   dtype='timedelta64[ns]', freq='6H')

    Specify ``start``, ``end``, and ``periods``; the frequency is generated
    automatically (linearly spaced).

    >>> pd.timedelta_range(start='1 day', end='5 days', periods=4)
    TimedeltaIndex(['1 days 00:00:00', '2 days 08:00:00', '3 days 16:00:00',
                '5 days 00:00:00'],
               dtype='timedelta64[ns]', freq=None)
    """
    if freq is None and com.any_none(periods, start, end):
        freq = "D"

    freq, freq_infer = dtl.maybe_infer_freq(freq)
    tdarr = TimedeltaArray._generate_range(start, end, periods, freq, closed=closed)
    return TimedeltaIndex._simple_new(tdarr, name=name)
