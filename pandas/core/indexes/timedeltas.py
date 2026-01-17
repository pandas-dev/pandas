"""implement the TimedeltaIndex"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    cast,
)

from pandas._libs import (
    index as libindex,
    lib,
)
from pandas._libs.tslibs import (
    Resolution,
    Timedelta,
    to_offset,
)
from pandas._libs.tslibs.dtypes import abbrev_to_npy_unit
from pandas.util._decorators import set_module

from pandas.core.dtypes.common import (
    is_scalar,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import ArrowDtype
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.arrays.timedeltas import TimedeltaArray
import pandas.core.common as com
from pandas.core.indexes.base import (
    Index,
    maybe_extract_name,
)
from pandas.core.indexes.datetimelike import DatetimeTimedeltaMixin
from pandas.core.indexes.extension import inherit_names

if TYPE_CHECKING:
    from pandas._libs import NaTType
    from pandas._libs.tslibs import (
        Day,
        Tick,
    )
    from pandas._typing import (
        DtypeObj,
        TimeUnit,
    )


@inherit_names(
    [
        "__neg__",
        "__pos__",
        "__abs__",
        "total_seconds",
        "round",
        "floor",
        "ceil",
        *TimedeltaArray._field_ops,
    ],
    TimedeltaArray,
    wrap=True,
)
@inherit_names(
    [
        "components",
        "to_pytimedelta",
        "sum",
        "std",
        "median",
    ],
    TimedeltaArray,
)
@set_module("pandas")
class TimedeltaIndex(DatetimeTimedeltaMixin):
    """
    Immutable Index of timedelta64 data.

    Represented internally as int64, and scalars returned Timedelta objects.

    Parameters
    ----------
    data : array-like (1-dimensional), optional
        Optional timedelta-like data to construct index with.
    freq : str or pandas offset object, optional
        One of pandas date offset strings or corresponding objects. The string
        ``'infer'`` can be passed in order to set the frequency of the index as
        the inferred frequency upon creation.
    dtype : numpy.dtype or str, default None
        Valid ``numpy`` dtypes are ``timedelta64[ns]``, ``timedelta64[us]``,
        ``timedelta64[ms]``, and ``timedelta64[s]``.
    copy : bool, default None
        Whether to copy input data, only relevant for array, Series, and Index
        inputs (for other input, e.g. a list, a new array is created anyway).
        Defaults to True for array input and False for Index/Series.
        Set to False to avoid copying array input at your own risk (if you
        know the input data won't be modified elsewhere).
        Set to True to force copying Series/Index input up front.
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
    To learn more about the frequency strings, please see
    :ref:`this link<timeseries.offset_aliases>`.

    Examples
    --------
    >>> pd.TimedeltaIndex(["0 days", "1 days", "2 days", "3 days", "4 days"])
    TimedeltaIndex(['0 days', '1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[us]', freq=None)

    We can also let pandas infer the frequency when possible.

    >>> pd.TimedeltaIndex(np.arange(5) * 24 * 3600 * 1e9, freq="infer")
    TimedeltaIndex(['0 days', '1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[ns]', freq='D')
    """

    _typ = "timedeltaindex"

    _data_cls = TimedeltaArray

    @property
    def _engine_type(self) -> type[libindex.TimedeltaEngine]:
        return libindex.TimedeltaEngine

    _data: TimedeltaArray

    # Use base class method instead of DatetimeTimedeltaMixin._get_string_slice
    _get_string_slice = Index._get_string_slice

    # error: Signature of "_resolution_obj" incompatible with supertype
    # "DatetimeIndexOpsMixin"
    @property
    def _resolution_obj(self) -> Resolution | None:  # type: ignore[override]
        return self._data._resolution_obj

    # -------------------------------------------------------------------
    # Constructors

    def __new__(
        cls,
        data=None,
        freq=lib.no_default,
        dtype=None,
        copy: bool | None = None,
        name=None,
    ):
        name = maybe_extract_name(name, data, cls)

        # GH#63388
        data, copy = cls._maybe_copy_array_input(data, copy, dtype)

        if is_scalar(data):
            cls._raise_scalar_data_error(data)

        if dtype is not None:
            dtype = pandas_dtype(dtype)

        if (
            isinstance(data, TimedeltaArray)
            and freq is lib.no_default
            and (dtype is None or dtype == data.dtype)
        ):
            if copy:
                data = data.copy()
            return cls._simple_new(data, name=name)

        if (
            isinstance(data, TimedeltaIndex)
            and freq is lib.no_default
            and name is None
            and (dtype is None or dtype == data.dtype)
        ):
            if copy:
                return data.copy()
            else:
                return data._view()

        # - Cases checked above all return/raise before reaching here - #

        tdarr = TimedeltaArray._from_sequence_not_strict(
            data, freq=freq, unit=None, dtype=dtype, copy=copy
        )
        refs = None
        if not copy and isinstance(data, (ABCSeries, Index)):
            refs = data._references

        return cls._simple_new(tdarr, name=name, refs=refs)

    # -------------------------------------------------------------------

    def _is_comparable_dtype(self, dtype: DtypeObj) -> bool:
        """
        Can we compare values of the given dtype to our own?
        """
        if isinstance(dtype, ArrowDtype):
            return dtype.kind == "m"
        return lib.is_np_dtype(dtype, "m")  # aka self._data._is_recognized_dtype

    # -------------------------------------------------------------------
    # Indexing Methods

    def get_loc(self, key):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int, slice, or ndarray[int]
        """
        self._check_indexing_error(key)

        try:
            key = self._data._validate_scalar(key, unbox=False)
        except TypeError as err:
            raise KeyError(key) from err

        return Index.get_loc(self, key)

    # error: Return type "tuple[Timedelta | NaTType, Resolution]" of
    # "_parse_with_reso" incompatible with return type
    # "tuple[datetime, Resolution]" in supertype
    # "pandas.core.indexes.datetimelike.DatetimeIndexOpsMixin"
    def _parse_with_reso(self, label: str) -> tuple[Timedelta | NaTType, Resolution]:  # type: ignore[override]
        parsed = Timedelta(label)
        if isinstance(parsed, Timedelta):
            reso = Resolution.get_reso_from_freqstr(parsed.unit)
        else:
            # i.e. pd.NaT
            reso = Resolution.get_reso_from_freqstr("s")
        return parsed, reso

    def _parsed_string_to_bounds(self, reso: Resolution, parsed: Timedelta):
        # reso is unused, included to match signature of DTI/PI
        lbound = parsed.round(parsed.resolution_string)
        rbound = (
            lbound
            + to_offset(parsed.resolution_string)
            - Timedelta(1, unit=self.unit).as_unit(self.unit)
        )
        return lbound, rbound

    # -------------------------------------------------------------------

    @property
    def inferred_type(self) -> str:
        return "timedelta64"


@set_module("pandas")
def timedelta_range(
    start=None,
    end=None,
    periods: int | None = None,
    freq=None,
    name=None,
    closed=None,
    *,
    unit: TimeUnit | None = None,
) -> TimedeltaIndex:
    """
    Return a fixed frequency TimedeltaIndex with day as the default.

    Parameters
    ----------
    start : str or timedelta-like, default None
        Left bound for generating timedeltas.
    end : str or timedelta-like, default None
        Right bound for generating timedeltas.
    periods : int, default None
        Number of periods to generate.
    freq : str, Timedelta, datetime.timedelta, or DateOffset, default 'D'
        Frequency strings can have multiples, e.g. '5h'.
    name : Hashable, default None
        Name of the resulting TimedeltaIndex.
    closed : str, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None).
    unit : {'s', 'ms', 'us', 'ns', None}, default None
        Specify the desired resolution of the result.
        If not specified, this is inferred from the 'start', 'end', and 'freq'
        using the same inference as :class:`Timedelta` taking the highest
        resolution of the three that are provided.

        .. versionadded:: 2.0.0

    Returns
    -------
    TimedeltaIndex
        Fixed frequency, with day as the default.

    See Also
    --------
    date_range : Return a fixed frequency DatetimeIndex.
    period_range : Return a fixed frequency PeriodIndex.

    Notes
    -----
    Of the four parameters ``start``, ``end``, ``periods``, and ``freq``,
    a maximum of three can be specified at once. Of the three parameters
    ``start``, ``end``, and ``periods``, at least two must be specified.
    If ``freq`` is omitted, the resulting ``DatetimeIndex`` will have
    ``periods`` linearly spaced elements between ``start`` and ``end``
    (closed on both sides).

    To learn more about the frequency strings, please see
    :ref:`this link<timeseries.offset_aliases>`.

    Examples
    --------
    >>> pd.timedelta_range(start="1 day", periods=4)
    TimedeltaIndex(['1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[us]', freq='D')

    The ``closed`` parameter specifies which endpoint is included.  The default
    behavior is to include both endpoints.

    >>> pd.timedelta_range(start="1 day", periods=4, closed="right")
    TimedeltaIndex(['2 days', '3 days', '4 days'],
                   dtype='timedelta64[us]', freq='D')

    The ``freq`` parameter specifies the frequency of the TimedeltaIndex.
    Only fixed frequencies can be passed, non-fixed frequencies such as
    'M' (month end) will raise.

    >>> pd.timedelta_range(start="1 day", end="2 days", freq="6h")
    TimedeltaIndex(['1 days 00:00:00', '1 days 06:00:00', '1 days 12:00:00',
                    '1 days 18:00:00', '2 days 00:00:00'],
                   dtype='timedelta64[us]', freq='6h')

    Specify ``start``, ``end``, and ``periods``; the frequency is generated
    automatically (linearly spaced).

    >>> pd.timedelta_range(start="1 day", end="5 days", periods=4)
    TimedeltaIndex(['1 days 00:00:00', '2 days 08:00:00', '3 days 16:00:00',
                    '5 days 00:00:00'],
                   dtype='timedelta64[us]', freq=None)

    **Specify a unit**

    >>> pd.timedelta_range("1 Day", periods=3, freq="100000D", unit="s")
    TimedeltaIndex(['1 days', '100001 days', '200001 days'],
                   dtype='timedelta64[s]', freq='100000D')
    """
    if freq is None and com.any_none(periods, start, end):
        freq = "D"
    freq = to_offset(freq)

    if com.count_not_none(start, end, periods, freq) != 3:
        # This check needs to come before the `unit = start.unit` line below
        raise ValueError(
            "Of the four parameters: start, end, periods, "
            "and freq, exactly three must be specified"
        )

    if unit is None:
        # Infer the unit based on the inputs

        if start is not None and end is not None:
            start = Timedelta(start)
            end = Timedelta(end)
            start = cast(Timedelta, start)
            end = cast(Timedelta, end)
            if abbrev_to_npy_unit(start.unit) > abbrev_to_npy_unit(end.unit):
                unit = cast("TimeUnit", start.unit)
            else:
                unit = cast("TimeUnit", end.unit)
        elif start is not None:
            start = Timedelta(start)
            start = cast(Timedelta, start)
            unit = cast("TimeUnit", start.unit)
        else:
            end = Timedelta(end)
            end = cast(Timedelta, end)
            unit = cast("TimeUnit", end.unit)

        # Last we need to watch out for cases where the 'freq' implies a higher
        #  unit than either start or end
        if freq is not None:
            freq = cast("Tick | Day", freq)
            creso = abbrev_to_npy_unit(unit)
            if freq._creso > creso:  # pyright: ignore[reportAttributeAccessIssue]
                unit = cast("TimeUnit", freq.base.freqstr)

    tdarr = TimedeltaArray._generate_range(
        start, end, periods, freq, closed=closed, unit=unit
    )
    return TimedeltaIndex._simple_new(tdarr, name=name)
