import cython

import time
from typing import Any
import warnings
from cpython.datetime cimport (PyDateTime_IMPORT,
                               PyDateTime_Check,
                               PyDate_Check,
                               PyDelta_Check,
                               datetime, timedelta, date,
                               time as dt_time)
PyDateTime_IMPORT

from dateutil.relativedelta import relativedelta

import numpy as np
cimport numpy as cnp
from numpy cimport int64_t
cnp.import_array()

# TODO: formalize having _libs.properties "above" tslibs in the dependency structure
from pandas._libs.properties import cache_readonly

from pandas._libs.tslibs cimport util
from pandas._libs.tslibs.util cimport is_integer_object, is_datetime64_object

from pandas._libs.tslibs.base cimport ABCTimestamp

from pandas._libs.tslibs.ccalendar import (
    MONTHS, DAYS, MONTH_ALIASES, MONTH_TO_CAL_NUM, weekday_to_int, int_to_weekday,
)
from pandas._libs.tslibs.ccalendar cimport get_days_in_month, dayofweek
from pandas._libs.tslibs.conversion cimport (
    convert_datetime_to_tsobject,
    localize_pydatetime,
)
from pandas._libs.tslibs.nattype cimport NPY_NAT, c_NaT as NaT
from pandas._libs.tslibs.np_datetime cimport (
    npy_datetimestruct, dtstruct_to_dt64, dt64_to_dtstruct)
from pandas._libs.tslibs.timezones cimport utc_pytz as UTC
from pandas._libs.tslibs.tzconversion cimport tz_convert_single

from .timedeltas cimport delta_to_nanoseconds

# ---------------------------------------------------------------------
# Constants

_offset_to_period_map = {
    'WEEKDAY': 'D',
    'EOM': 'M',
    'BM': 'M',
    'BQS': 'Q',
    'QS': 'Q',
    'BQ': 'Q',
    'BA': 'A',
    'AS': 'A',
    'BAS': 'A',
    'MS': 'M',
    'D': 'D',
    'C': 'C',
    'B': 'B',
    'T': 'T',
    'S': 'S',
    'L': 'L',
    'U': 'U',
    'N': 'N',
    'H': 'H',
    'Q': 'Q',
    'A': 'A',
    'W': 'W',
    'M': 'M',
    'Y': 'A',
    'BY': 'A',
    'YS': 'A',
    'BYS': 'A'}

need_suffix = ['QS', 'BQ', 'BQS', 'YS', 'AS', 'BY', 'BA', 'BYS', 'BAS']

for __prefix in need_suffix:
    for _m in MONTHS:
        key = f'{__prefix}-{_m}'
        _offset_to_period_map[key] = _offset_to_period_map[__prefix]

for __prefix in ['A', 'Q']:
    for _m in MONTHS:
        _alias = f'{__prefix}-{_m}'
        _offset_to_period_map[_alias] = _alias

for _d in DAYS:
    _offset_to_period_map[f'W-{_d}'] = f'W-{_d}'


# ---------------------------------------------------------------------
# Misc Helpers

cdef bint is_offset_object(object obj):
    return isinstance(obj, BaseOffset)


cdef bint is_tick_object(object obj):
    return isinstance(obj, Tick)


cdef to_offset(object obj):
    """
    Wrap pandas.tseries.frequencies.to_offset to keep centralize runtime
    imports
    """
    if isinstance(obj, BaseOffset):
        return obj
    from pandas.tseries.frequencies import to_offset
    return to_offset(obj)


def as_datetime(obj: datetime) -> datetime:
    if isinstance(obj, ABCTimestamp):
        return obj.to_pydatetime()
    return obj


cpdef bint is_normalized(datetime dt):
    if dt.hour != 0 or dt.minute != 0 or dt.second != 0 or dt.microsecond != 0:
        # Regardless of whether dt is datetime vs Timestamp
        return False
    if isinstance(dt, ABCTimestamp):
        return dt.nanosecond == 0
    return True


def apply_index_wraps(func):
    # Note: normally we would use `@functools.wraps(func)`, but this does
    # not play nicely with cython class methods
    def wrapper(self, other):

        is_index = not util.is_array(other._data)

        # operate on DatetimeArray
        arr = other._data if is_index else other

        result = func(self, arr)

        if is_index:
            # Wrap DatetimeArray result back to DatetimeIndex
            result = type(other)._simple_new(result, name=other.name)

        if self.normalize:
            result = result.to_period('D').to_timestamp()
        return result

    # do @functools.wraps(func) manually since it doesn't work on cdef funcs
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    try:
        wrapper.__module__ = func.__module__
    except AttributeError:
        # AttributeError: 'method_descriptor' object has no
        # attribute '__module__'
        pass
    return wrapper


def apply_wraps(func):
    # Note: normally we would use `@functools.wraps(func)`, but this does
    # not play nicely with cython class methods

    def wrapper(self, other):
        from pandas import Timestamp

        if other is NaT:
            return NaT
        elif isinstance(other, BaseOffset) or PyDelta_Check(other):
            # timedelta path
            return func(self, other)
        elif is_datetime64_object(other) or PyDate_Check(other):
            # PyDate_Check includes date, datetime
            other = Timestamp(other)
        else:
            # This will end up returning NotImplemented back in __add__
            raise ApplyTypeError

        tz = other.tzinfo
        nano = other.nanosecond

        if self._adjust_dst:
            other = other.tz_localize(None)

        result = func(self, other)

        result = Timestamp(result)
        if self._adjust_dst:
            result = result.tz_localize(tz)

        if self.normalize:
            result = result.normalize()

        # nanosecond may be deleted depending on offset process
        if not self.normalize and nano != 0:
            if result.nanosecond != nano:
                if result.tz is not None:
                    # convert to UTC
                    value = result.tz_localize(None).value
                else:
                    value = result.value
                result = Timestamp(value + nano)

        if tz is not None and result.tzinfo is None:
            result = result.tz_localize(tz)

        return result

    # do @functools.wraps(func) manually since it doesn't work on cdef funcs
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    try:
        wrapper.__module__ = func.__module__
    except AttributeError:
        # AttributeError: 'method_descriptor' object has no
        # attribute '__module__'
        pass
    return wrapper


cdef _wrap_timedelta_result(result):
    """
    Tick operations dispatch to their Timedelta counterparts.  Wrap the result
    of these operations in a Tick if possible.

    Parameters
    ----------
    result : object

    Returns
    -------
    object
    """
    if PyDelta_Check(result):
        # convert Timedelta back to a Tick
        return delta_to_tick(result)

    return result

# ---------------------------------------------------------------------
# Business Helpers

cpdef int get_lastbday(int year, int month) nogil:
    """
    Find the last day of the month that is a business day.

    Parameters
    ----------
    year : int
    month : int

    Returns
    -------
    last_bday : int
    """
    cdef:
        int wkday, days_in_month

    wkday = dayofweek(year, month, 1)
    days_in_month = get_days_in_month(year, month)
    return days_in_month - max(((wkday + days_in_month - 1) % 7) - 4, 0)


cpdef int get_firstbday(int year, int month) nogil:
    """
    Find the first day of the month that is a business day.

    Parameters
    ----------
    year : int
    month : int

    Returns
    -------
    first_bday : int
    """
    cdef:
        int first, wkday

    wkday = dayofweek(year, month, 1)
    first = 1
    if wkday == 5:  # on Saturday
        first = 3
    elif wkday == 6:  # on Sunday
        first = 2
    return first


def _get_calendar(weekmask, holidays, calendar):
    """Generate busdaycalendar"""
    if isinstance(calendar, np.busdaycalendar):
        if not holidays:
            holidays = tuple(calendar.holidays)
        elif not isinstance(holidays, tuple):
            holidays = tuple(holidays)
        else:
            # trust that calendar.holidays and holidays are
            # consistent
            pass
        return calendar, holidays

    if holidays is None:
        holidays = []
    try:
        holidays = holidays + calendar.holidays().tolist()
    except AttributeError:
        pass
    holidays = [to_dt64D(dt) for dt in holidays]
    holidays = tuple(sorted(holidays))

    kwargs = {'weekmask': weekmask}
    if holidays:
        kwargs['holidays'] = holidays

    busdaycalendar = np.busdaycalendar(**kwargs)
    return busdaycalendar, holidays


def to_dt64D(dt):
    # Currently
    # > np.datetime64(dt.datetime(2013,5,1),dtype='datetime64[D]')
    # numpy.datetime64('2013-05-01T02:00:00.000000+0200')
    # Thus astype is needed to cast datetime to datetime64[D]
    if getattr(dt, 'tzinfo', None) is not None:
        # Get the nanosecond timestamp,
        #  equiv `Timestamp(dt).value` or `dt.timestamp() * 10**9`
        nanos = getattr(dt, "nanosecond", 0)
        i8 = convert_datetime_to_tsobject(dt, tz=None, nanos=nanos).value
        dt = tz_convert_single(i8, UTC, dt.tzinfo)
        dt = np.int64(dt).astype('datetime64[ns]')
    else:
        dt = np.datetime64(dt)
    if dt.dtype.name != "datetime64[D]":
        dt = dt.astype("datetime64[D]")
    return dt


# ---------------------------------------------------------------------
# Validation


def _validate_business_time(t_input):
    if isinstance(t_input, str):
        try:
            t = time.strptime(t_input, '%H:%M')
            return dt_time(hour=t.tm_hour, minute=t.tm_min)
        except ValueError:
            raise ValueError("time data must match '%H:%M' format")
    elif isinstance(t_input, dt_time):
        if t_input.second != 0 or t_input.microsecond != 0:
            raise ValueError(
                "time data must be specified only with hour and minute")
        return t_input
    else:
        raise ValueError("time data must be string or datetime.time")


# ---------------------------------------------------------------------
# Constructor Helpers

relativedelta_kwds = {'years', 'months', 'weeks', 'days', 'year', 'month',
                      'day', 'weekday', 'hour', 'minute', 'second',
                      'microsecond', 'nanosecond', 'nanoseconds', 'hours',
                      'minutes', 'seconds', 'microseconds'}


def _determine_offset(kwds):
    # timedelta is used for sub-daily plural offsets and all singular
    # offsets relativedelta is used for plural offsets of daily length or
    # more nanosecond(s) are handled by apply_wraps
    kwds_no_nanos = dict(
        (k, v) for k, v in kwds.items()
        if k not in ('nanosecond', 'nanoseconds')
    )
    # TODO: Are nanosecond and nanoseconds allowed somewhere?

    _kwds_use_relativedelta = ('years', 'months', 'weeks', 'days',
                               'year', 'month', 'week', 'day', 'weekday',
                               'hour', 'minute', 'second', 'microsecond')

    use_relativedelta = False
    if len(kwds_no_nanos) > 0:
        if any(k in _kwds_use_relativedelta for k in kwds_no_nanos):
            offset = relativedelta(**kwds_no_nanos)
            use_relativedelta = True
        else:
            # sub-daily offset - use timedelta (tz-aware)
            offset = timedelta(**kwds_no_nanos)
    else:
        offset = timedelta(1)
    return offset, use_relativedelta


# ---------------------------------------------------------------------
# Mixins & Singletons


class ApplyTypeError(TypeError):
    # sentinel class for catching the apply error to return NotImplemented
    pass


# ---------------------------------------------------------------------
# Base Classes

cdef class BaseOffset:
    """
    Base class for DateOffset methods that are not overridden by subclasses
    and will (after pickle errors are resolved) go into a cdef class.
    """
    _typ = "dateoffset"
    _day_opt = None
    _attributes = frozenset(['n', 'normalize'])
    _use_relativedelta = False
    _adjust_dst = True
    _deprecations = frozenset(["isAnchored", "onOffset"])

    cdef readonly:
        int64_t n
        bint normalize
        dict _cache

    def __init__(self, n=1, normalize=False):
        n = self._validate_n(n)
        self.n = n
        self.normalize = normalize
        self._cache = {}

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            try:
                # GH#23524 if to_offset fails, we are dealing with an
                #  incomparable type so == is False and != is True
                other = to_offset(other)
            except ValueError:
                # e.g. "infer"
                return False
        try:
            return self._params == other._params
        except AttributeError:
            # other is not a DateOffset object
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self._params)

    @cache_readonly
    def _params(self):
        """
        Returns a tuple containing all of the attributes needed to evaluate
        equality between two DateOffset objects.
        """
        # NB: non-cython subclasses override property with cache_readonly
        d = getattr(self, "__dict__", {})
        all_paras = d.copy()
        all_paras["n"] = self.n
        all_paras["normalize"] = self.normalize
        for attr in self._attributes:
            if hasattr(self, attr) and attr not in d:
                # cython attributes are not in __dict__
                all_paras[attr] = getattr(self, attr)

        if 'holidays' in all_paras and not all_paras['holidays']:
            all_paras.pop('holidays')
        exclude = ['kwds', 'name', 'calendar']
        attrs = [(k, v) for k, v in all_paras.items()
                 if (k not in exclude) and (k[0] != '_')]
        attrs = sorted(set(attrs))
        params = tuple([str(type(self))] + attrs)
        return params

    @property
    def kwds(self):
        # for backwards-compatibility
        kwds = {name: getattr(self, name, None) for name in self._attributes
                if name not in ['n', 'normalize']}
        return {name: kwds[name] for name in kwds if kwds[name] is not None}

    @property
    def base(self):
        """
        Returns a copy of the calling offset object with n=1 and all other
        attributes equal.
        """
        return type(self)(n=1, normalize=self.normalize, **self.kwds)

    def __add__(self, other):
        if not isinstance(self, BaseOffset):
            # cython semantics; this is __radd__
            return other.__add__(self)
        try:
            return self.apply(other)
        except ApplyTypeError:
            return NotImplemented

    def __sub__(self, other):
        if PyDateTime_Check(other):
            raise TypeError('Cannot subtract datetime from offset.')
        elif type(other) == type(self):
            return type(self)(self.n - other.n, normalize=self.normalize,
                              **self.kwds)
        elif not isinstance(self, BaseOffset):
            # cython semantics, this is __rsub__
            return (-other).__add__(self)
        else:  # pragma: no cover
            return NotImplemented

    def __call__(self, other):
        warnings.warn(
            "DateOffset.__call__ is deprecated and will be removed in a future "
            "version.  Use `offset + other` instead.",
            FutureWarning,
            stacklevel=1,
        )
        return self.apply(other)

    def __mul__(self, other):
        if util.is_array(other):
            return np.array([self * x for x in other])
        elif is_integer_object(other):
            return type(self)(n=other * self.n, normalize=self.normalize,
                              **self.kwds)
        elif not isinstance(self, BaseOffset):
            # cython semantics, this is __rmul__
            return other.__mul__(self)
        return NotImplemented

    def __neg__(self):
        # Note: we are deferring directly to __mul__ instead of __rmul__, as
        # that allows us to use methods that can go in a `cdef class`
        return self * -1

    def copy(self):
        # Note: we are deferring directly to __mul__ instead of __rmul__, as
        # that allows us to use methods that can go in a `cdef class`
        return self * 1

    # ------------------------------------------------------------------
    # Name and Rendering Methods

    def __repr__(self) -> str:
        className = getattr(self, '_outputName', type(self).__name__)

        if abs(self.n) != 1:
            plural = 's'
        else:
            plural = ''

        n_str = ""
        if self.n != 1:
            n_str = f"{self.n} * "

        out = f'<{n_str}{className}{plural}{self._repr_attrs()}>'
        return out

    def _repr_attrs(self) -> str:
        exclude = {"n", "inc", "normalize"}
        attrs = []
        for attr in sorted(self._attributes):
            if attr.startswith("_") or attr == "kwds" or not hasattr(self, attr):
                # DateOffset may not have some of these attributes
                continue
            elif attr not in exclude:
                value = getattr(self, attr)
                attrs.append(f"{attr}={value}")

        out = ""
        if attrs:
            out += ": " + ", ".join(attrs)
        return out

    @property
    def name(self) -> str:
        return self.rule_code

    @property
    def _prefix(self) -> str:
        raise NotImplementedError("Prefix not defined")

    @property
    def rule_code(self) -> str:
        return self._prefix

    @cache_readonly
    def freqstr(self) -> str:
        try:
            code = self.rule_code
        except NotImplementedError:
            return str(repr(self))

        if self.n != 1:
            fstr = f"{self.n}{code}"
        else:
            fstr = code

        try:
            if self._offset:
                fstr += self._offset_str()
        except AttributeError:
            # TODO: standardize `_offset` vs `offset` naming convention
            pass

        return fstr

    def _offset_str(self) -> str:
        return ""

    # ------------------------------------------------------------------

    @apply_index_wraps
    def apply_index(self, index):
        """
        Vectorized apply of DateOffset to DatetimeIndex,
        raises NotImplementedError for offsets without a
        vectorized implementation.

        Parameters
        ----------
        index : DatetimeIndex

        Returns
        -------
        DatetimeIndex
        """
        raise NotImplementedError(
            f"DateOffset subclass {type(self).__name__} "
            "does not have a vectorized implementation"
        )

    def rollback(self, dt):
        """
        Roll provided date backward to next offset only if not on offset.

        Returns
        -------
        TimeStamp
            Rolled timestamp if not on offset, otherwise unchanged timestamp.
        """
        from pandas import Timestamp
        dt = Timestamp(dt)
        if not self.is_on_offset(dt):
            dt = dt - type(self)(1, normalize=self.normalize, **self.kwds)
        return dt

    def rollforward(self, dt):
        """
        Roll provided date forward to next offset only if not on offset.

        Returns
        -------
        TimeStamp
            Rolled timestamp if not on offset, otherwise unchanged timestamp.
        """
        from pandas import Timestamp
        dt = Timestamp(dt)
        if not self.is_on_offset(dt):
            dt = dt + type(self)(1, normalize=self.normalize, **self.kwds)
        return dt

    def _get_offset_day(self, datetime other):
        # subclass must implement `_day_opt`; calling from the base class
        # will raise NotImplementedError.
        return get_day_of_month(other, self._day_opt)

    def is_on_offset(self, dt) -> bool:
        if self.normalize and not is_normalized(dt):
            return False

        # Default (slow) method for determining if some date is a member of the
        # date range generated by this offset. Subclasses may have this
        # re-implemented in a nicer way.
        a = dt
        b = (dt + self) - self
        return a == b

    # ------------------------------------------------------------------

    # Staticmethod so we can call from Tick.__init__, will be unnecessary
    #  once BaseOffset is a cdef class and is inherited by Tick
    @staticmethod
    def _validate_n(n):
        """
        Require that `n` be an integer.

        Parameters
        ----------
        n : int

        Returns
        -------
        nint : int

        Raises
        ------
        TypeError if `int(n)` raises
        ValueError if n != int(n)
        """
        if util.is_timedelta64_object(n):
            raise TypeError(f'`n` argument must be an integer, got {type(n)}')
        try:
            nint = int(n)
        except (ValueError, TypeError):
            raise TypeError(f'`n` argument must be an integer, got {type(n)}')
        if n != nint:
            raise ValueError(f'`n` argument must be an integer, got {n}')
        return nint

    def __setstate__(self, state):
        """Reconstruct an instance from a pickled state"""
        if isinstance(self, MonthOffset):
            # We can't just override MonthOffset.__setstate__ because of the
            #  combination of MRO resolution and cython not handling
            #  multiple inheritance nicely for cdef classes.
            state.pop("_use_relativedelta", False)
            state.pop("offset", None)
            state.pop("_offset", None)
            state.pop("kwds", {})

        if 'offset' in state:
            # Older (<0.22.0) versions have offset attribute instead of _offset
            if '_offset' in state:  # pragma: no cover
                raise AssertionError('Unexpected key `_offset`')
            state['_offset'] = state.pop('offset')
            state['kwds']['offset'] = state['_offset']

        if '_offset' in state and not isinstance(state['_offset'], timedelta):
            # relativedelta, we need to populate using its kwds
            offset = state['_offset']
            odict = offset.__dict__
            kwds = {key: odict[key] for key in odict if odict[key]}
            state.update(kwds)

        self.n = state.pop("n")
        self.normalize = state.pop("normalize")
        self._cache = state.pop("_cache", {})

        if not len(state):
            # FIXME: kludge because some classes no longer have a __dict__,
            #  so we need to short-circuit before raising on the next line
            return

        self.__dict__.update(state)

        if 'weekmask' in state and 'holidays' in state:
            weekmask = state.pop("weekmask")
            holidays = state.pop("holidays")
            calendar, holidays = _get_calendar(weekmask=weekmask,
                                               holidays=holidays,
                                               calendar=None)
            self.calendar = calendar
            self.holidays = holidays

    def __getstate__(self):
        """Return a pickleable state"""
        state = getattr(self, "__dict__", {}).copy()
        state["n"] = self.n
        state["normalize"] = self.normalize

        # we don't want to actually pickle the calendar object
        # as its a np.busyday; we recreate on deserialization
        if 'calendar' in state:
            del state['calendar']
        try:
            state['kwds'].pop('calendar')
        except KeyError:
            pass

        return state

    @property
    def nanos(self):
        raise ValueError(f"{self} is a non-fixed frequency")

    def onOffset(self, dt) -> bool:
        warnings.warn(
            "onOffset is a deprecated, use is_on_offset instead",
            FutureWarning,
            stacklevel=1,
        )
        return self.is_on_offset(dt)

    def isAnchored(self) -> bool:
        warnings.warn(
            "isAnchored is a deprecated, use is_anchored instead",
            FutureWarning,
            stacklevel=1,
        )
        return self.is_anchored()

    def is_anchored(self) -> bool:
        # TODO: Does this make sense for the general case?  It would help
        # if there were a canonical docstring for what is_anchored means.
        return self.n == 1


cdef class SingleConstructorOffset(BaseOffset):
    @classmethod
    def _from_name(cls, suffix=None):
        # default _from_name calls cls with no args
        if suffix:
            raise ValueError(f"Bad freq suffix {suffix}")
        return cls()


# ---------------------------------------------------------------------
# Tick Offsets

cdef class Tick(SingleConstructorOffset):
    # ensure that reversed-ops with numpy scalars return NotImplemented
    __array_priority__ = 1000
    _adjust_dst = False
    _prefix = "undefined"
    _attributes = frozenset(["n", "normalize"])

    def __init__(self, n=1, normalize=False):
        n = self._validate_n(n)
        self.n = n
        self.normalize = False
        self._cache = {}
        if normalize:
            # GH#21427
            raise ValueError(
                "Tick offset with `normalize=True` are not allowed."
            )

    def _repr_attrs(self) -> str:
        # Since cdef classes have no __dict__, we need to override
        return ""

    @property
    def delta(self):
        from .timedeltas import Timedelta
        return self.n * Timedelta(self._nanos_inc)

    @property
    def nanos(self) -> int64_t:
        return self.n * self._nanos_inc

    def is_on_offset(self, dt) -> bool:
        return True

    def is_anchored(self) -> bool:
        return False

    # --------------------------------------------------------------------
    # Comparison and Arithmetic Methods

    def __eq__(self, other):
        if isinstance(other, str):
            try:
                # GH#23524 if to_offset fails, we are dealing with an
                #  incomparable type so == is False and != is True
                other = to_offset(other)
            except ValueError:
                # e.g. "infer"
                return False
        return self.delta == other

    def __ne__(self, other):
        return not (self == other)

    def __le__(self, other):
        return self.delta.__le__(other)

    def __lt__(self, other):
        return self.delta.__lt__(other)

    def __ge__(self, other):
        return self.delta.__ge__(other)

    def __gt__(self, other):
        return self.delta.__gt__(other)

    def __truediv__(self, other):
        if not isinstance(self, Tick):
            # cython semantics mean the args are sometimes swapped
            result = other.delta.__rtruediv__(self)
        else:
            result = self.delta.__truediv__(other)
        return _wrap_timedelta_result(result)

    def __add__(self, other):
        if not isinstance(self, Tick):
            # cython semantics; this is __radd__
            return other.__add__(self)

        if isinstance(other, Tick):
            if type(self) == type(other):
                return type(self)(self.n + other.n)
            else:
                return delta_to_tick(self.delta + other.delta)
        try:
            return self.apply(other)
        except ApplyTypeError:
            # Includes pd.Period
            return NotImplemented
        except OverflowError as err:
            raise OverflowError(
                f"the add operation between {self} and {other} will overflow"
            ) from err

    def apply(self, other):
        # Timestamp can handle tz and nano sec, thus no need to use apply_wraps
        if isinstance(other, ABCTimestamp):

            # GH#15126
            # in order to avoid a recursive
            # call of __add__ and __radd__ if there is
            # an exception, when we call using the + operator,
            # we directly call the known method
            result = other.__add__(self)
            if result is NotImplemented:
                raise OverflowError
            return result
        elif other is NaT:
            return NaT
        elif is_datetime64_object(other) or PyDate_Check(other):
            # PyDate_Check includes date, datetime
            from pandas import Timestamp
            return Timestamp(other) + self

        if PyDelta_Check(other):
            return other + self.delta
        elif isinstance(other, type(self)):
            # TODO: this is reached in tests that specifically call apply,
            #  but should not be reached "naturally" because __add__ should
            #  catch this case first.
            return type(self)(self.n + other.n)

        raise ApplyTypeError(f"Unhandled type: {type(other).__name__}")

    # --------------------------------------------------------------------
    # Pickle Methods

    def __reduce__(self):
        return (type(self), (self.n,))

    def __setstate__(self, state):
        self.n = state["n"]
        self.normalize = False


cdef class Day(Tick):
    _nanos_inc = 24 * 3600 * 1_000_000_000
    _prefix = "D"


cdef class Hour(Tick):
    _nanos_inc = 3600 * 1_000_000_000
    _prefix = "H"


cdef class Minute(Tick):
    _nanos_inc = 60 * 1_000_000_000
    _prefix = "T"


cdef class Second(Tick):
    _nanos_inc = 1_000_000_000
    _prefix = "S"


cdef class Milli(Tick):
    _nanos_inc = 1_000_000
    _prefix = "L"


cdef class Micro(Tick):
    _nanos_inc = 1000
    _prefix = "U"


cdef class Nano(Tick):
    _nanos_inc = 1
    _prefix = "N"


def delta_to_tick(delta: timedelta) -> Tick:
    if delta.microseconds == 0 and getattr(delta, "nanoseconds", 0) == 0:
        # nanoseconds only for pd.Timedelta
        if delta.seconds == 0:
            return Day(delta.days)
        else:
            seconds = delta.days * 86400 + delta.seconds
            if seconds % 3600 == 0:
                return Hour(seconds / 3600)
            elif seconds % 60 == 0:
                return Minute(seconds / 60)
            else:
                return Second(seconds)
    else:
        nanos = delta_to_nanoseconds(delta)
        if nanos % 1_000_000 == 0:
            return Milli(nanos // 1_000_000)
        elif nanos % 1000 == 0:
            return Micro(nanos // 1000)
        else:  # pragma: no cover
            return Nano(nanos)


# --------------------------------------------------------------------


cdef class BusinessMixin(SingleConstructorOffset):
    """
    Mixin to business types to provide related functions.
    """

    cdef readonly:
        timedelta _offset

    def __init__(self, n=1, normalize=False, offset=timedelta(0)):
        BaseOffset.__init__(self, n, normalize)
        self._offset = offset

    @property
    def offset(self):
        """
        Alias for self._offset.
        """
        # Alias for backward compat
        return self._offset

    def _repr_attrs(self) -> str:
        if self.offset:
            attrs = [f"offset={repr(self.offset)}"]
        else:
            attrs = []
        out = ""
        if attrs:
            out += ": " + ", ".join(attrs)
        return out

    cpdef __setstate__(self, state):
        # We need to use a cdef/cpdef method to set the readonly _offset attribute
        if "_offset" in state:
            self._offset = state.pop("_offset")
        elif "offset" in state:
            self._offset = state.pop("offset")
        BaseOffset.__setstate__(self, state)


cdef class BusinessHourMixin(BusinessMixin):
    _adjust_dst = False

    cdef readonly:
        tuple start, end

    def __init__(
            self, n=1, normalize=False, start="09:00", end="17:00", offset=timedelta(0)
        ):
        BusinessMixin.__init__(self, n, normalize, offset)

        # must be validated here to equality check
        if np.ndim(start) == 0:
            # i.e. not is_list_like
            start = [start]
        if not len(start):
            raise ValueError("Must include at least 1 start time")

        if np.ndim(end) == 0:
            # i.e. not is_list_like
            end = [end]
        if not len(end):
            raise ValueError("Must include at least 1 end time")

        start = np.array([_validate_business_time(x) for x in start])
        end = np.array([_validate_business_time(x) for x in end])

        # Validation of input
        if len(start) != len(end):
            raise ValueError("number of starting time and ending time must be the same")
        num_openings = len(start)

        # sort starting and ending time by starting time
        index = np.argsort(start)

        # convert to tuple so that start and end are hashable
        start = tuple(start[index])
        end = tuple(end[index])

        total_secs = 0
        for i in range(num_openings):
            total_secs += self._get_business_hours_by_sec(start[i], end[i])
            total_secs += self._get_business_hours_by_sec(
                end[i], start[(i + 1) % num_openings]
            )
        if total_secs != 24 * 60 * 60:
            raise ValueError(
                "invalid starting and ending time(s): "
                "opening hours should not touch or overlap with "
                "one another"
            )

        self.start = start
        self.end = end

    def __reduce__(self):
        return type(self), (self.n, self.normalize, self.start, self.end, self.offset)

    def _repr_attrs(self) -> str:
        out = super()._repr_attrs()
        hours = ",".join(
            f'{st.strftime("%H:%M")}-{en.strftime("%H:%M")}'
            for st, en in zip(self.start, self.end)
        )
        attrs = [f"{self._prefix}={hours}"]
        out += ": " + ", ".join(attrs)
        return out

    def _get_business_hours_by_sec(self, start, end):
        """
        Return business hours in a day by seconds.
        """
        # create dummy datetime to calculate business hours in a day
        dtstart = datetime(2014, 4, 1, start.hour, start.minute)
        day = 1 if start < end else 2
        until = datetime(2014, 4, day, end.hour, end.minute)
        return int((until - dtstart).total_seconds())

    def _get_closing_time(self, dt):
        """
        Get the closing time of a business hour interval by its opening time.

        Parameters
        ----------
        dt : datetime
            Opening time of a business hour interval.

        Returns
        -------
        result : datetime
            Corresponding closing time.
        """
        for i, st in enumerate(self.start):
            if st.hour == dt.hour and st.minute == dt.minute:
                return dt + timedelta(
                    seconds=self._get_business_hours_by_sec(st, self.end[i])
                )
        assert False


class CustomMixin:
    """
    Mixin for classes that define and validate calendar, holidays,
    and weekdays attributes.
    """

    def __init__(self, weekmask, holidays, calendar):
        calendar, holidays = _get_calendar(
            weekmask=weekmask, holidays=holidays, calendar=calendar
        )
        # Custom offset instances are identified by the
        # following two attributes. See DateOffset._params()
        # holidays, weekmask

        object.__setattr__(self, "weekmask", weekmask)
        object.__setattr__(self, "holidays", holidays)
        object.__setattr__(self, "calendar", calendar)


class WeekOfMonthMixin(SingleConstructorOffset):
    """
    Mixin for methods common to WeekOfMonth and LastWeekOfMonth.
    """
    def __init__(self, n=1, normalize=False, weekday=0):
        BaseOffset.__init__(self, n, normalize)
        object.__setattr__(self, "weekday", weekday)

        if weekday < 0 or weekday > 6:
            raise ValueError(f"Day must be 0<=day<=6, got {weekday}")

    @apply_wraps
    def apply(self, other):
        compare_day = self._get_offset_day(other)

        months = self.n
        if months > 0 and compare_day > other.day:
            months -= 1
        elif months <= 0 and compare_day < other.day:
            months += 1

        shifted = shift_month(other, months, "start")
        to_day = self._get_offset_day(shifted)
        return shift_day(shifted, to_day - shifted.day)

    def is_on_offset(self, dt) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        return dt.day == self._get_offset_day(dt)

    @property
    def rule_code(self) -> str:
        weekday = int_to_weekday.get(self.weekday, "")
        if self.week == -1:
            # LastWeekOfMonth
            return f"{self._prefix}-{weekday}"
        return f"{self._prefix}-{self.week + 1}{weekday}"


# ----------------------------------------------------------------------
# Year-Based Offset Classes

cdef class YearOffset(SingleConstructorOffset):
    """
    DateOffset that just needs a month.
    """
    _attributes = frozenset(["n", "normalize", "month"])

    # _default_month: int  # FIXME: python annotation here breaks things

    cdef readonly:
        int month

    def __init__(self, n=1, normalize=False, month=None):
        BaseOffset.__init__(self, n, normalize)

        month = month if month is not None else self._default_month
        self.month = month

        if month < 1 or month > 12:
            raise ValueError("Month must go from 1 to 12")

    cpdef __setstate__(self, state):
        self.month = state.pop("month")
        self.n = state.pop("n")
        self.normalize = state.pop("normalize")
        self._cache = {}

    def __reduce__(self):
        return type(self), (self.n, self.normalize, self.month)

    @classmethod
    def _from_name(cls, suffix=None):
        kwargs = {}
        if suffix:
            kwargs["month"] = MONTH_TO_CAL_NUM[suffix]
        return cls(**kwargs)

    @property
    def rule_code(self) -> str:
        month = MONTH_ALIASES[self.month]
        return f"{self._prefix}-{month}"

    def is_on_offset(self, dt) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        return dt.month == self.month and dt.day == self._get_offset_day(dt)

    def _get_offset_day(self, other) -> int:
        # override BaseOffset method to use self.month instead of other.month
        # TODO: there may be a more performant way to do this
        return get_day_of_month(
            other.replace(month=self.month), self._day_opt
        )

    @apply_wraps
    def apply(self, other):
        years = roll_yearday(other, self.n, self.month, self._day_opt)
        months = years * 12 + (self.month - other.month)
        return shift_month(other, months, self._day_opt)

    @apply_index_wraps
    def apply_index(self, dtindex):
        shifted = shift_quarters(
            dtindex.asi8, self.n, self.month, self._day_opt, modby=12
        )
        return type(dtindex)._simple_new(shifted, dtype=dtindex.dtype)


cdef class BYearEnd(YearOffset):
    """
    DateOffset increments between business EOM dates.
    """

    _outputName = "BusinessYearEnd"
    _default_month = 12
    _prefix = "BA"
    _day_opt = "business_end"


cdef class BYearBegin(YearOffset):
    """
    DateOffset increments between business year begin dates.
    """

    _outputName = "BusinessYearBegin"
    _default_month = 1
    _prefix = "BAS"
    _day_opt = "business_start"


cdef class YearEnd(YearOffset):
    """
    DateOffset increments between calendar year ends.
    """

    _default_month = 12
    _prefix = "A"
    _day_opt = "end"


cdef class YearBegin(YearOffset):
    """
    DateOffset increments between calendar year begin dates.
    """

    _default_month = 1
    _prefix = "AS"
    _day_opt = "start"


# ----------------------------------------------------------------------
# Quarter-Based Offset Classes

cdef class QuarterOffset(SingleConstructorOffset):
    _attributes = frozenset(["n", "normalize", "startingMonth"])
    # TODO: Consider combining QuarterOffset and YearOffset __init__ at some
    #       point.  Also apply_index, is_on_offset, rule_code if
    #       startingMonth vs month attr names are resolved

    # FIXME: python annotations here breaks things
    # _default_startingMonth: int
    # _from_name_startingMonth: int

    cdef readonly:
        int startingMonth

    def __init__(self, n=1, normalize=False, startingMonth=None):
        BaseOffset.__init__(self, n, normalize)

        if startingMonth is None:
            startingMonth = self._default_startingMonth
        self.startingMonth = startingMonth

    cpdef __setstate__(self, state):
        self.startingMonth = state.pop("startingMonth")
        self.n = state.pop("n")
        self.normalize = state.pop("normalize")

    def __reduce__(self):
        return type(self), (self.n, self.normalize, self.startingMonth)

    @classmethod
    def _from_name(cls, suffix=None):
        kwargs = {}
        if suffix:
            kwargs["startingMonth"] = MONTH_TO_CAL_NUM[suffix]
        else:
            if cls._from_name_startingMonth is not None:
                kwargs["startingMonth"] = cls._from_name_startingMonth
        return cls(**kwargs)

    @property
    def rule_code(self) -> str:
        month = MONTH_ALIASES[self.startingMonth]
        return f"{self._prefix}-{month}"

    def is_anchored(self) -> bool:
        return self.n == 1 and self.startingMonth is not None

    def is_on_offset(self, dt) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        mod_month = (dt.month - self.startingMonth) % 3
        return mod_month == 0 and dt.day == self._get_offset_day(dt)

    @apply_wraps
    def apply(self, other):
        # months_since: find the calendar quarter containing other.month,
        # e.g. if other.month == 8, the calendar quarter is [Jul, Aug, Sep].
        # Then find the month in that quarter containing an is_on_offset date for
        # self.  `months_since` is the number of months to shift other.month
        # to get to this on-offset month.
        months_since = other.month % 3 - self.startingMonth % 3
        qtrs = roll_qtrday(
            other, self.n, self.startingMonth, day_opt=self._day_opt, modby=3
        )
        months = qtrs * 3 - months_since
        return shift_month(other, months, self._day_opt)

    @apply_index_wraps
    def apply_index(self, dtindex):
        shifted = shift_quarters(
            dtindex.asi8, self.n, self.startingMonth, self._day_opt
        )
        return type(dtindex)._simple_new(shifted, dtype=dtindex.dtype)


cdef class BQuarterEnd(QuarterOffset):
    """
    DateOffset increments between business Quarter dates.

    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/30/2007, 6/29/2007, ...
    """
    _outputName = "BusinessQuarterEnd"
    _default_startingMonth = 3
    _from_name_startingMonth = 12
    _prefix = "BQ"
    _day_opt = "business_end"


# TODO: This is basically the same as BQuarterEnd
cdef class BQuarterBegin(QuarterOffset):
    _outputName = "BusinessQuarterBegin"
    # I suspect this is wrong for *all* of them.
    # TODO: What does the above comment refer to?
    _default_startingMonth = 3
    _from_name_startingMonth = 1
    _prefix = "BQS"
    _day_opt = "business_start"


cdef class QuarterEnd(QuarterOffset):
    """
    DateOffset increments between business Quarter dates.

    startingMonth = 1 corresponds to dates like 1/31/2007, 4/30/2007, ...
    startingMonth = 2 corresponds to dates like 2/28/2007, 5/31/2007, ...
    startingMonth = 3 corresponds to dates like 3/31/2007, 6/30/2007, ...
    """
    _outputName = "QuarterEnd"
    _default_startingMonth = 3
    _prefix = "Q"
    _day_opt = "end"


cdef class QuarterBegin(QuarterOffset):
    _outputName = "QuarterBegin"
    _default_startingMonth = 3
    _from_name_startingMonth = 1
    _prefix = "QS"
    _day_opt = "start"


# ----------------------------------------------------------------------
# Month-Based Offset Classes

cdef class MonthOffset(SingleConstructorOffset):
    def is_on_offset(self, dt) -> bool:
        if self.normalize and not is_normalized(dt):
            return False
        return dt.day == self._get_offset_day(dt)

    @apply_wraps
    def apply(self, other):
        compare_day = self._get_offset_day(other)
        n = roll_convention(other.day, self.n, compare_day)
        return shift_month(other, n, self._day_opt)

    @apply_index_wraps
    def apply_index(self, dtindex):
        shifted = shift_months(dtindex.asi8, self.n, self._day_opt)
        return type(dtindex)._simple_new(shifted, dtype=dtindex.dtype)


cdef class MonthEnd(MonthOffset):
    """
    DateOffset of one month end.
    """
    _prefix = "M"
    _day_opt = "end"


cdef class MonthBegin(MonthOffset):
    """
    DateOffset of one month at beginning.
    """
    _prefix = "MS"
    _day_opt = "start"


cdef class BusinessMonthEnd(MonthOffset):
    """
    DateOffset increments between business EOM dates.
    """
    _prefix = "BM"
    _day_opt = "business_end"


cdef class BusinessMonthBegin(MonthOffset):
    """
    DateOffset of one business month at beginning.
    """
    _prefix = "BMS"
    _day_opt = "business_start"


# ---------------------------------------------------------------------
# Special Offset Classes

cdef class FY5253Mixin(SingleConstructorOffset):
    cdef readonly:
        int startingMonth
        int weekday
        str variation

    def __init__(
        self, n=1, normalize=False, weekday=0, startingMonth=1, variation="nearest"
    ):
        BaseOffset.__init__(self, n, normalize)
        self.startingMonth = startingMonth
        self.weekday = weekday
        self.variation = variation

        if self.n == 0:
            raise ValueError("N cannot be 0")

        if self.variation not in ["nearest", "last"]:
            raise ValueError(f"{self.variation} is not a valid variation")

    def is_anchored(self) -> bool:
        return (
            self.n == 1 and self.startingMonth is not None and self.weekday is not None
        )

    # --------------------------------------------------------------------
    # Name-related methods

    @property
    def rule_code(self) -> str:
        prefix = self._prefix
        suffix = self.get_rule_code_suffix()
        return f"{prefix}-{suffix}"

    def _get_suffix_prefix(self) -> str:
        if self.variation == "nearest":
            return "N"
        else:
            return "L"

    def get_rule_code_suffix(self) -> str:
        prefix = self._get_suffix_prefix()
        month = MONTH_ALIASES[self.startingMonth]
        weekday = int_to_weekday[self.weekday]
        return f"{prefix}-{month}-{weekday}"


# ----------------------------------------------------------------------
# RelativeDelta Arithmetic

def shift_day(other: datetime, days: int) -> datetime:
    """
    Increment the datetime `other` by the given number of days, retaining
    the time-portion of the datetime.  For tz-naive datetimes this is
    equivalent to adding a timedelta.  For tz-aware datetimes it is similar to
    dateutil's relativedelta.__add__, but handles pytz tzinfo objects.

    Parameters
    ----------
    other : datetime or Timestamp
    days : int

    Returns
    -------
    shifted: datetime or Timestamp
    """
    if other.tzinfo is None:
        return other + timedelta(days=days)

    tz = other.tzinfo
    naive = other.replace(tzinfo=None)
    shifted = naive + timedelta(days=days)
    return localize_pydatetime(shifted, tz)


cdef inline int year_add_months(npy_datetimestruct dts, int months) nogil:
    """new year number after shifting npy_datetimestruct number of months"""
    return dts.year + (dts.month + months - 1) // 12


cdef inline int month_add_months(npy_datetimestruct dts, int months) nogil:
    """
    New month number after shifting npy_datetimestruct
    number of months.
    """
    cdef:
        int new_month = (dts.month + months) % 12
    return 12 if new_month == 0 else new_month


@cython.wraparound(False)
@cython.boundscheck(False)
cdef shift_quarters(
    const int64_t[:] dtindex,
    int quarters,
    int q1start_month,
    object day,
    int modby=3,
):
    """
    Given an int64 array representing nanosecond timestamps, shift all elements
    by the specified number of quarters using DateOffset semantics.

    Parameters
    ----------
    dtindex : int64_t[:] timestamps for input dates
    quarters : int number of quarters to shift
    q1start_month : int month in which Q1 begins by convention
    day : {'start', 'end', 'business_start', 'business_end'}
    modby : int (3 for quarters, 12 for years)

    Returns
    -------
    out : ndarray[int64_t]
    """
    cdef:
        Py_ssize_t i
        npy_datetimestruct dts
        int count = len(dtindex)
        int months_to_roll, months_since, n, compare_day
        bint roll_check
        int64_t[:] out = np.empty(count, dtype='int64')

    if day == 'start':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                n = quarters

                months_since = (dts.month - q1start_month) % modby

                # offset semantics - if on the anchor point and going backwards
                # shift to next
                if n <= 0 and (months_since != 0 or
                               (months_since == 0 and dts.day > 1)):
                    n += 1

                dts.year = year_add_months(dts, modby * n - months_since)
                dts.month = month_add_months(dts, modby * n - months_since)
                dts.day = 1

                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'end':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                n = quarters

                months_since = (dts.month - q1start_month) % modby

                if n <= 0 and months_since != 0:
                    # The general case of this condition would be
                    # `months_since != 0 or (months_since == 0 and
                    #    dts.day > get_days_in_month(dts.year, dts.month))`
                    # but the get_days_in_month inequality would never hold.
                    n += 1
                elif n > 0 and (months_since == 0 and
                                dts.day < get_days_in_month(dts.year,
                                                            dts.month)):
                    n -= 1

                dts.year = year_add_months(dts, modby * n - months_since)
                dts.month = month_add_months(dts, modby * n - months_since)
                dts.day = get_days_in_month(dts.year, dts.month)

                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'business_start':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                n = quarters

                months_since = (dts.month - q1start_month) % modby
                # compare_day is only relevant for comparison in the case
                # where months_since == 0.
                compare_day = get_firstbday(dts.year, dts.month)

                if n <= 0 and (months_since != 0 or
                               (months_since == 0 and dts.day > compare_day)):
                    # make sure to roll forward, so negate
                    n += 1
                elif n > 0 and (months_since == 0 and dts.day < compare_day):
                    # pretend to roll back if on same month but
                    # before compare_day
                    n -= 1

                dts.year = year_add_months(dts, modby * n - months_since)
                dts.month = month_add_months(dts, modby * n - months_since)

                dts.day = get_firstbday(dts.year, dts.month)

                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'business_end':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                n = quarters

                months_since = (dts.month - q1start_month) % modby
                # compare_day is only relevant for comparison in the case
                # where months_since == 0.
                compare_day = get_lastbday(dts.year, dts.month)

                if n <= 0 and (months_since != 0 or
                               (months_since == 0 and dts.day > compare_day)):
                    # make sure to roll forward, so negate
                    n += 1
                elif n > 0 and (months_since == 0 and dts.day < compare_day):
                    # pretend to roll back if on same month but
                    # before compare_day
                    n -= 1

                dts.year = year_add_months(dts, modby * n - months_since)
                dts.month = month_add_months(dts, modby * n - months_since)

                dts.day = get_lastbday(dts.year, dts.month)

                out[i] = dtstruct_to_dt64(&dts)

    else:
        raise ValueError("day must be None, 'start', 'end', "
                         "'business_start', or 'business_end'")

    return np.asarray(out)


@cython.wraparound(False)
@cython.boundscheck(False)
def shift_months(const int64_t[:] dtindex, int months, object day=None):
    """
    Given an int64-based datetime index, shift all elements
    specified number of months using DateOffset semantics

    day: {None, 'start', 'end'}
       * None: day of month
       * 'start' 1st day of month
       * 'end' last day of month
    """
    cdef:
        Py_ssize_t i
        npy_datetimestruct dts
        int count = len(dtindex)
        int months_to_roll
        bint roll_check
        int64_t[:] out = np.empty(count, dtype='int64')

    if day is None:
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                dts.year = year_add_months(dts, months)
                dts.month = month_add_months(dts, months)

                dts.day = min(dts.day, get_days_in_month(dts.year, dts.month))
                out[i] = dtstruct_to_dt64(&dts)
    elif day == 'start':
        roll_check = False
        if months <= 0:
            months += 1
            roll_check = True
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months

                # offset semantics - if on the anchor point and going backwards
                # shift to next
                if roll_check and dts.day == 1:
                    months_to_roll -= 1

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)
                dts.day = 1

                out[i] = dtstruct_to_dt64(&dts)
    elif day == 'end':
        roll_check = False
        if months > 0:
            months -= 1
            roll_check = True
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months

                # similar semantics - when adding shift forward by one
                # month if already at an end of month
                if roll_check and dts.day == get_days_in_month(dts.year,
                                                               dts.month):
                    months_to_roll += 1

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)

                dts.day = get_days_in_month(dts.year, dts.month)
                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'business_start':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months
                compare_day = get_firstbday(dts.year, dts.month)

                months_to_roll = roll_convention(dts.day, months_to_roll,
                                                 compare_day)

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)

                dts.day = get_firstbday(dts.year, dts.month)
                out[i] = dtstruct_to_dt64(&dts)

    elif day == 'business_end':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months
                compare_day = get_lastbday(dts.year, dts.month)

                months_to_roll = roll_convention(dts.day, months_to_roll,
                                                 compare_day)

                dts.year = year_add_months(dts, months_to_roll)
                dts.month = month_add_months(dts, months_to_roll)

                dts.day = get_lastbday(dts.year, dts.month)
                out[i] = dtstruct_to_dt64(&dts)

    else:
        raise ValueError("day must be None, 'start', 'end', "
                         "'business_start', or 'business_end'")

    return np.asarray(out)


def shift_month(stamp: datetime, months: int,
                day_opt: object=None) -> datetime:
    """
    Given a datetime (or Timestamp) `stamp`, an integer `months` and an
    option `day_opt`, return a new datetimelike that many months later,
    with day determined by `day_opt` using relativedelta semantics.

    Scalar analogue of shift_months

    Parameters
    ----------
    stamp : datetime or Timestamp
    months : int
    day_opt : None, 'start', 'end', 'business_start', 'business_end', or int
        None: returned datetimelike has the same day as the input, or the
              last day of the month if the new month is too short
        'start': returned datetimelike has day=1
        'end': returned datetimelike has day on the last day of the month
        'business_start': returned datetimelike has day on the first
            business day of the month
        'business_end': returned datetimelike has day on the last
            business day of the month
        int: returned datetimelike has day equal to day_opt

    Returns
    -------
    shifted : datetime or Timestamp (same as input `stamp`)
    """
    cdef:
        int year, month, day
        int days_in_month, dy

    dy = (stamp.month + months) // 12
    month = (stamp.month + months) % 12

    if month == 0:
        month = 12
        dy -= 1
    year = stamp.year + dy

    if day_opt is None:
        days_in_month = get_days_in_month(year, month)
        day = min(stamp.day, days_in_month)
    elif day_opt == 'start':
        day = 1
    elif day_opt == 'end':
        day = get_days_in_month(year, month)
    elif day_opt == 'business_start':
        # first business day of month
        day = get_firstbday(year, month)
    elif day_opt == 'business_end':
        # last business day of month
        day = get_lastbday(year, month)
    elif is_integer_object(day_opt):
        days_in_month = get_days_in_month(year, month)
        day = min(day_opt, days_in_month)
    else:
        raise ValueError(day_opt)
    return stamp.replace(year=year, month=month, day=day)


cdef int get_day_of_month(datetime other, day_opt) except? -1:
    """
    Find the day in `other`'s month that satisfies a DateOffset's is_on_offset
    policy, as described by the `day_opt` argument.

    Parameters
    ----------
    other : datetime or Timestamp
    day_opt : 'start', 'end', 'business_start', 'business_end', or int
        'start': returns 1
        'end': returns last day of the month
        'business_start': returns the first business day of the month
        'business_end': returns the last business day of the month
        int: returns the day in the month indicated by `other`, or the last of
            day the month if the value exceeds in that month's number of days.

    Returns
    -------
    day_of_month : int

    Examples
    -------
    >>> other = datetime(2017, 11, 14)
    >>> get_day_of_month(other, 'start')
    1
    >>> get_day_of_month(other, 'end')
    30

    """
    cdef:
        int days_in_month

    if day_opt == 'start':
        return 1
    elif day_opt == 'end':
        days_in_month = get_days_in_month(other.year, other.month)
        return days_in_month
    elif day_opt == 'business_start':
        # first business day of month
        return get_firstbday(other.year, other.month)
    elif day_opt == 'business_end':
        # last business day of month
        return get_lastbday(other.year, other.month)
    elif is_integer_object(day_opt):
        days_in_month = get_days_in_month(other.year, other.month)
        return min(day_opt, days_in_month)
    elif day_opt is None:
        # Note: unlike `shift_month`, get_day_of_month does not
        # allow day_opt = None
        raise NotImplementedError
    else:
        raise ValueError(day_opt)


cpdef int roll_convention(int other, int n, int compare) nogil:
    """
    Possibly increment or decrement the number of periods to shift
    based on rollforward/rollbackward conventions.

    Parameters
    ----------
    other : int, generally the day component of a datetime
    n : number of periods to increment, before adjusting for rolling
    compare : int, generally the day component of a datetime, in the same
              month as the datetime form which `other` was taken.

    Returns
    -------
    n : int number of periods to increment
    """
    if n > 0 and other < compare:
        n -= 1
    elif n <= 0 and other > compare:
        # as if rolled forward already
        n += 1
    return n


def roll_qtrday(other: datetime, n: int, month: int,
                day_opt: object, modby: int=3) -> int:
    """
    Possibly increment or decrement the number of periods to shift
    based on rollforward/rollbackward conventions.

    Parameters
    ----------
    other : datetime or Timestamp
    n : number of periods to increment, before adjusting for rolling
    month : int reference month giving the first month of the year
    day_opt : 'start', 'end', 'business_start', 'business_end', or int
        The convention to use in finding the day in a given month against
        which to compare for rollforward/rollbackward decisions.
    modby : int 3 for quarters, 12 for years

    Returns
    -------
    n : int number of periods to increment

    See Also
    --------
    get_day_of_month : Find the day in a month provided an offset.
    """
    cdef:
        int months_since
    # TODO: Merge this with roll_yearday by setting modby=12 there?
    #       code de-duplication versus perf hit?
    # TODO: with small adjustments this could be used in shift_quarters
    months_since = other.month % modby - month % modby

    if n > 0:
        if months_since < 0 or (months_since == 0 and
                                other.day < get_day_of_month(other,
                                                             day_opt)):
            # pretend to roll back if on same month but
            # before compare_day
            n -= 1
    else:
        if months_since > 0 or (months_since == 0 and
                                other.day > get_day_of_month(other,
                                                             day_opt)):
            # make sure to roll forward, so negate
            n += 1
    return n


def roll_yearday(other: datetime, n: int, month: int, day_opt: object) -> int:
    """
    Possibly increment or decrement the number of periods to shift
    based on rollforward/rollbackward conventions.

    Parameters
    ----------
    other : datetime or Timestamp
    n : number of periods to increment, before adjusting for rolling
    month : reference month giving the first month of the year
    day_opt : 'start', 'end', 'business_start', 'business_end', or int
        The day of the month to compare against that of `other` when
        incrementing or decrementing the number of periods:

        'start': 1
        'end': last day of the month
        'business_start': first business day of the month
        'business_end': last business day of the month
        int: day in the month indicated by `other`, or the last of day
            the month if the value exceeds in that month's number of days.

    Returns
    -------
    n : int number of periods to increment

    Notes
    -----
    * Mirrors `roll_check` in shift_months

    Examples
    -------
    >>> month = 3
    >>> day_opt = 'start'              # `other` will be compared to March 1
    >>> other = datetime(2017, 2, 10)  # before March 1
    >>> roll_yearday(other, 2, month, day_opt)
    1
    >>> roll_yearday(other, -7, month, day_opt)
    -7
    >>>
    >>> other = Timestamp('2014-03-15', tz='US/Eastern')  # after March 1
    >>> roll_yearday(other, 2, month, day_opt)
    2
    >>> roll_yearday(other, -7, month, day_opt)
    -6

    >>> month = 6
    >>> day_opt = 'end'                # `other` will be compared to June 30
    >>> other = datetime(1999, 6, 29)  # before June 30
    >>> roll_yearday(other, 5, month, day_opt)
    4
    >>> roll_yearday(other, -7, month, day_opt)
    -7
    >>>
    >>> other = Timestamp(2072, 8, 24, 6, 17, 18)  # after June 30
    >>> roll_yearday(other, 5, month, day_opt)
    5
    >>> roll_yearday(other, -7, month, day_opt)
    -6

    """
    # Note: The other.day < ... condition will never hold when day_opt=='start'
    # and the other.day > ... condition will never hold when day_opt=='end'.
    # At some point these extra checks may need to be optimized away.
    # But that point isn't today.
    if n > 0:
        if other.month < month or (other.month == month and
                                   other.day < get_day_of_month(other,
                                                                day_opt)):
            n -= 1
    else:
        if other.month > month or (other.month == month and
                                   other.day > get_day_of_month(other,
                                                                day_opt)):
            n += 1
    return n
