cimport cython
from datetime import (
    timedelta,
    timezone,
)
import zoneinfo
from zoneinfo._zoneinfo import ZoneInfo as _ZoneInfo

from pandas.compat._optional import import_optional_dependency

from cpython.datetime cimport (
    datetime,
    timedelta,
    tzinfo,
)

# dateutil compat

from dateutil.tz import (
    gettz as dateutil_gettz,
    tzfile as _dateutil_tzfile,
    tzlocal as _dateutil_tzlocal,
    tzutc as _dateutil_tzutc,
)
import numpy as np

pytz = import_optional_dependency("pytz", errors="ignore")

cimport numpy as cnp
from numpy cimport int64_t

cnp.import_array()

# ----------------------------------------------------------------------
from pandas._libs.tslibs.util cimport (
    get_nat,
    is_integer_object,
)


cdef int64_t NPY_NAT = get_nat()
cdef tzinfo utc_stdlib = timezone.utc
cdef tzinfo utc_pytz = pytz.UTC if pytz else None
cdef tzinfo utc_dateutil_str = dateutil_gettz("UTC")  # NB: *not* the same as tzutc()

cdef tzinfo utc_zoneinfo = None
cdef type ZoneInfo = zoneinfo.ZoneInfo


# ----------------------------------------------------------------------

cdef bint is_utc_zoneinfo(tzinfo tz):
    # Workaround for cases with missing tzdata
    #  https://github.com/pandas-dev/pandas/pull/46425#discussion_r830633025
    if tz is None:
        return False

    global utc_zoneinfo
    if utc_zoneinfo is None:
        try:
            utc_zoneinfo = zoneinfo.ZoneInfo("UTC")
        except zoneinfo.ZoneInfoNotFoundError:
            return False

    return tz is utc_zoneinfo


cpdef inline bint is_utc(tzinfo tz):
    return (
        tz is utc_stdlib
        or isinstance(tz, _dateutil_tzutc)
        or tz is utc_dateutil_str
        or is_utc_zoneinfo(tz)
        or (utc_pytz is not None and tz is utc_pytz)
    )


cdef bint is_zoneinfo(tzinfo tz):
    return isinstance(tz, ZoneInfo)


cdef bint is_tzlocal(tzinfo tz):
    return isinstance(tz, _dateutil_tzlocal)


cdef bint treat_tz_as_pytz(tzinfo tz):
    return (hasattr(tz, "_utc_transition_times") and
            hasattr(tz, "_transition_info"))


cdef bint treat_tz_as_dateutil(tzinfo tz):
    return hasattr(tz, "_trans_list") and hasattr(tz, "_trans_idx")


# Returns str or tzinfo object
cpdef inline object get_timezone(tzinfo tz):
    """
    We need to do several things here:
    1) Distinguish between pytz and dateutil timezones
    2) Not be over-specific (e.g. US/Eastern with/without DST is same *zone*
       but a different tz object)
    3) Provide something to serialize when we're storing a datetime object
       in pytables.

    We return a string prefaced with dateutil if it's a dateutil tz, else just
    the tz name. It needs to be a string so that we can serialize it with
    UJSON/pytables. maybe_get_tz (below) is the inverse of this process.
    """
    if tz is None:
        raise TypeError("tz argument cannot be None")
    if is_utc(tz):
        return tz
    elif is_zoneinfo(tz):
        return tz.key
    elif treat_tz_as_pytz(tz):
        zone = tz.zone
        if zone is None:
            return tz
        return zone
    elif treat_tz_as_dateutil(tz):
        if ".tar.gz" in tz._filename:
            raise ValueError(
                "Bad tz filename. Dateutil on python 3 on windows has a "
                "bug which causes tzfile._filename to be the same for all "
                "timezone files. Please construct dateutil timezones "
                'implicitly by passing a string like "dateutil/Europe'
                '/London" when you construct your pandas objects instead '
                "of passing a timezone object. See "
                "https://github.com/pandas-dev/pandas/pull/7362")
        return "dateutil/" + tz._filename
    else:
        return tz


cpdef inline tzinfo maybe_get_tz(object tz):
    """
    (Maybe) Construct a timezone object from a string. If tz is a string, use
    it to construct a timezone object. Otherwise, just return tz.
    """
    if isinstance(tz, str):
        if tz == "tzlocal()":
            tz = _dateutil_tzlocal()
        elif tz.startswith("dateutil/"):
            zone = tz[9:]
            tz = dateutil_gettz(zone)
            # On Python 3 on Windows, the filename is not always set correctly.
            if isinstance(tz, _dateutil_tzfile) and ".tar.gz" in tz._filename:
                tz._filename = zone
        elif tz[0] in {"-", "+"}:
            hours = int(tz[0:3])
            minutes = int(tz[0] + tz[4:6])
            tz = timezone(timedelta(hours=hours, minutes=minutes))
        elif tz[0:4] in {"UTC-", "UTC+"}:
            hours = int(tz[3:6])
            minutes = int(tz[3] + tz[7:9])
            tz = timezone(timedelta(hours=hours, minutes=minutes))
        elif tz == "UTC" or tz == "utc":
            tz = utc_stdlib
        else:
            tz = zoneinfo.ZoneInfo(tz)
    elif is_integer_object(tz):
        tz = timezone(timedelta(seconds=tz))
    elif isinstance(tz, tzinfo):
        if treat_tz_as_pytz(tz) and pytz is None:
            # call again for raising proper error
            import_optional_dependency("pytz")
        pass
    elif tz is None:
        pass
    else:
        raise TypeError(type(tz))
    return tz


def _p_tz_cache_key(tz: tzinfo):
    """
    Python interface for cache function to facilitate testing.
    """
    return tz_cache_key(tz)


# Timezone data caches, key is the pytz string or dateutil file name.
dst_cache = {}


cdef object tz_cache_key(tzinfo tz):
    """
    Return the key in the cache for the timezone info object or None
    if unknown.

    The key is currently the tz string for pytz timezones, the filename for
    dateutil timezones.

    Notes
    -----
    This cannot just be the hash of a timezone object. Unfortunately, the
    hashes of two dateutil tz objects which represent the same timezone are
    not equal (even though the tz objects will compare equal and represent
    the same tz file). Also, pytz objects are not always hashable so we use
    str(tz) instead.
    """
    if pytz is not None and isinstance(tz, pytz.tzinfo.BaseTzInfo):
        return tz.zone
    elif isinstance(tz, _dateutil_tzfile):
        if ".tar.gz" in tz._filename:
            raise ValueError("Bad tz filename. Dateutil on python 3 on "
                             "windows has a bug which causes tzfile._filename "
                             "to be the same for all timezone files. Please "
                             "construct dateutil timezones implicitly by "
                             'passing a string like "dateutil/Europe/London" '
                             "when you construct your pandas objects instead "
                             "of passing a timezone object. See "
                             "https://github.com/pandas-dev/pandas/pull/7362")
        return "dateutil" + tz._filename
    elif is_zoneinfo(tz):
        return "zoneinfo/" + tz.key
    else:
        return None


# ----------------------------------------------------------------------
# UTC Offsets


cdef timedelta get_utcoffset(tzinfo tz, datetime obj):
    try:
        return tz._utcoffset
    except AttributeError:
        return tz.utcoffset(obj)


cpdef inline bint is_fixed_offset(tzinfo tz):
    if treat_tz_as_dateutil(tz):
        if len(tz._trans_idx) == 0 and len(tz._trans_list) == 0:
            return 1
        else:
            return 0
    elif treat_tz_as_pytz(tz) and pytz is not None:
        if (len(tz._transition_info) == 0
                and len(tz._utc_transition_times) == 0):
            return 1
        else:
            return 0
    elif is_zoneinfo(tz):
        tz_py = _ZoneInfo(tz.key)
        if tz_py._fixed_offset:
            return 1
        else:
            return 0
    # This also implicitly accepts datetime.timezone objects which are
    # considered fixed
    return 1


cdef object _get_utc_trans_times_from_dateutil_tz(tzinfo tz):
    """
    Transition times in dateutil timezones are stored in local non-dst
    time.  This code converts them to UTC. It's the reverse of the code
    in dateutil.tz.tzfile.__init__.
    """
    new_trans = list(tz._trans_list)
    last_std_offset = 0
    for i, (trans, tti) in enumerate(zip(tz._trans_list, tz._trans_idx, strict=True)):
        if not tti.isdst:
            last_std_offset = tti.offset
        new_trans[i] = trans - last_std_offset
    return new_trans


@cython.wraparound(False)
@cython.boundscheck(False)
cdef int64_t[::1] unbox_utcoffsets(object transinfo):
    cdef:
        Py_ssize_t i
        cnp.npy_intp sz
        int64_t[::1] arr

    sz = len(transinfo)
    arr = cnp.PyArray_EMPTY(1, &sz, cnp.NPY_INT64, 0)

    for i in range(sz):
        arr[i] = int(transinfo[i][0].total_seconds()) * 1_000_000_000

    return arr


# ----------------------------------------------------------------------
# Daylight Savings


cdef tuple _get_zoneinfo_trans_and_deltas(tzinfo tz):
    """
    Get transition times and UTC offsets for a ZoneInfo timezone.

    Uses zoneinfo's Python fallback implementation to get transition data,
    including future transitions generated from POSIX TZ rules.

    Parameters
    ----------
    tz : ZoneInfo

    Returns
    -------
    trans : ndarray[int64_t]
        Nanosecond UTC times of DST transitions.
    deltas : ndarray[int64_t]
        Nanosecond UTC offsets corresponding to DST transitions.
    is_fixed : bint
        True if this is a fixed-offset timezone with no transitions.
    """
    cdef:
        int64_t fixed_offset_seconds, last_hist_ts, max_i8_seconds, start_utc, end_utc
        list trans_utc, deltas_seconds, future_trans
        int year, last_year, std_offset, dst_offset

    tz_py = _ZoneInfo(tz.key)

    if tz_py._fixed_offset:
        fixed_offset_seconds = int(tz_py._tz_after.utcoff.total_seconds())
        trans = np.array([NPY_NAT + 1], dtype=np.int64)
        deltas = np.array([fixed_offset_seconds], dtype="i8") * 1_000_000_000
        return trans, deltas, True

    trans_utc = list(tz_py._trans_utc)
    deltas_seconds = [int(info.utcoff.total_seconds()) for info in tz_py._ttinfos]

    has_future_dst = (
        hasattr(tz_py, "_tz_after")
        and tz_py._tz_after is not None
        and hasattr(tz_py._tz_after, "transitions")
    )
    if has_future_dst:
        # ZoneInfo's _trans_utc only contains historical data (typically up to
        # the late 1990s for many timezones). POSIX TZ rules stored in _tz_after
        # are required to generate transitions for current and future dates.
        # Without this block, DST-observing timezones like Europe/Amsterdam
        # would return incorrect offsets for any date after ~1996.
        tz_after = tz_py._tz_after
        std_offset = int(tz_after.std.utcoff.total_seconds())
        dst_offset = int(tz_after.dst.utcoff.total_seconds())

        if trans_utc:
            last_hist_ts = trans_utc[-1]
            try:
                last_year = datetime.fromtimestamp(last_hist_ts, timezone.utc).year
            except (OSError, OverflowError, ValueError):
                last_year = 1970
        else:
            last_hist_ts = 0
            last_year = 1970

        future_trans = []
        max_i8_seconds = np.iinfo(np.int64).max // 1_000_000_000
        for year in range(last_year, 2263):
            try:
                year_trans = tz_after.transitions(year)
                if not year_trans:
                    break
                start_local, end_local = year_trans
                start_utc = start_local - std_offset
                end_utc = end_local - dst_offset

                # transitions are is currently converted to nanoseconds below, so guard
                # against overflow by only including those that fit in the ns range
                if end_utc <= max_i8_seconds:
                    if start_utc > last_hist_ts:
                        future_trans.append((start_utc, dst_offset))
                    if end_utc > last_hist_ts:
                        future_trans.append((end_utc, std_offset))
            except Exception:
                break

        future_trans.sort()

        for t, d in future_trans:
            trans_utc.append(t)
            deltas_seconds.append(d)

    trans = np.array(trans_utc, dtype="i8") * 1_000_000_000
    trans = np.hstack([np.array([NPY_NAT + 1], dtype=np.int64), trans])

    first_offset_seconds = int(tz_py._tti_before.utcoff.total_seconds())
    deltas = np.array(deltas_seconds, dtype="i8") * 1_000_000_000
    deltas = np.hstack([[first_offset_seconds * 1_000_000_000], deltas])

    return trans, deltas, False


cdef object get_dst_info(tzinfo tz):
    """
    Returns
    -------
    ndarray[int64_t]
        Nanosecond UTC times of DST transitions.
    ndarray[int64_t]
        Nanosecond UTC offsets corresponding to DST transitions.
    str
        Describing the type of tzinfo object.
    """
    cache_key = tz_cache_key(tz)
    if cache_key is None:
        # e.g. pytz.FixedOffset, matplotlib.dates._UTC,
        # psycopg2.tz.FixedOffsetTimezone
        num = int(get_utcoffset(tz, None).total_seconds()) * 1_000_000_000
        # If we have e.g. ZoneInfo here, the get_utcoffset call will return None,
        #  so the total_seconds() call will raise AttributeError.
        return (np.array([NPY_NAT + 1], dtype=np.int64),
                np.array([num], dtype=np.int64),
                "unknown")

    if cache_key not in dst_cache:
        if treat_tz_as_pytz(tz):
            trans = np.array(tz._utc_transition_times, dtype="M8[ns]")
            trans = trans.view("i8")
            if tz._utc_transition_times[0].year == 1:
                trans[0] = NPY_NAT + 1
            deltas = unbox_utcoffsets(tz._transition_info)
            typ = "pytz"

        elif treat_tz_as_dateutil(tz):
            if len(tz._trans_list):
                # get utc trans times
                trans_list = _get_utc_trans_times_from_dateutil_tz(tz)
                trans = np.hstack([
                    np.array([0], dtype="M8[s]"),  # place holder for 1st item
                    np.array(trans_list, dtype="M8[s]")]).astype(
                    "M8[ns]")  # all trans listed
                trans = trans.view("i8")
                trans[0] = NPY_NAT + 1

                # deltas
                deltas = np.array([v.offset for v in (
                    tz._ttinfo_before,) + tz._trans_idx], dtype="i8")
                deltas *= 1_000_000_000
                typ = "dateutil"

            elif is_fixed_offset(tz):
                trans = np.array([NPY_NAT + 1], dtype=np.int64)
                deltas = np.array([tz._ttinfo_std.offset],
                                  dtype="i8") * 1_000_000_000
                typ = "fixed"
            else:
                # 2018-07-12 this is not reached in the tests, and this case
                # is not handled in any of the functions that call
                # get_dst_info.  If this case _were_ hit the calling
                # functions would then hit an IndexError because they assume
                # `deltas` is non-empty.
                # (under the just-deleted code that returned empty arrays)
                raise AssertionError("dateutil tzinfo is not a FixedOffset "
                                     "and has an empty `_trans_list`.", tz)

        elif is_zoneinfo(tz):
            trans, deltas, is_fixed = _get_zoneinfo_trans_and_deltas(tz)
            typ = "fixed" if is_fixed else "zoneinfo"

        else:
            # static tzinfo, we can get here with pytz.StaticTZInfo
            #  which are not caught by treat_tz_as_pytz
            trans = np.array([NPY_NAT + 1], dtype=np.int64)
            num = int(get_utcoffset(tz, None).total_seconds()) * 1_000_000_000
            deltas = np.array([num], dtype=np.int64)
            typ = "static"

        dst_cache[cache_key] = (trans, deltas, typ)

    return dst_cache[cache_key]


def infer_tzinfo(datetime start, datetime end):
    if start is not None and end is not None:
        tz = start.tzinfo
        if not tz_compare(tz, end.tzinfo):
            raise AssertionError(f"Inputs must both have the same timezone, "
                                 f"{tz} != {end.tzinfo}")
    elif start is not None:
        tz = start.tzinfo
    elif end is not None:
        tz = end.tzinfo
    else:
        tz = None
    return tz


cpdef bint tz_compare(tzinfo start, tzinfo end):
    """
    Compare string representations of timezones

    The same timezone can be represented as different instances of
    timezones. For example
    `<DstTzInfo 'Europe/Paris' LMT+0:09:00 STD>` and
    `<DstTzInfo 'Europe/Paris' CET+1:00:00 STD>` are essentially same
    timezones but aren't evaluated such, but the string representation
    for both of these is `'Europe/Paris'`.

    This exists only to add a notion of equality to pytz-style zones
    that is compatible with the notion of equality expected of tzinfo
    subclasses.

    Parameters
    ----------
    start : tzinfo
    end : tzinfo

    Returns:
    -------
    bool
    """
    # GH 18523
    if is_utc(start):
        # GH#38851 consider pytz/dateutil/stdlib UTCs as equivalent
        return is_utc(end)
    elif is_utc(end):
        # Ensure we don't treat tzlocal as equal to UTC when running in UTC
        return False
    elif start is None or end is None:
        return start is None and end is None
    return get_timezone(start) == get_timezone(end)


def tz_standardize(tz: tzinfo) -> tzinfo:
    """
    If the passed tz is a pytz timezone object, "normalize" it to a
    consistent version

    Parameters
    ----------
    tz : tzinfo

    Returns
    -------
    tzinfo

    Examples
    --------
    >>> from datetime import datetime
    >>> from pytz import timezone
    >>> tz = timezone('US/Pacific').normalize(
    ...     datetime(2014, 1, 1, tzinfo=pytz.utc)
    ... ).tzinfo
    >>> tz
    <DstTzInfo 'US/Pacific' PST-1 day, 16:00:00 STD>
    >>> tz_standardize(tz)
    <DstTzInfo 'US/Pacific' LMT-1 day, 16:07:00 STD>

    >>> tz = timezone('US/Pacific')
    >>> tz
    <DstTzInfo 'US/Pacific' LMT-1 day, 16:07:00 STD>
    >>> tz_standardize(tz)
    <DstTzInfo 'US/Pacific' LMT-1 day, 16:07:00 STD>
    """
    if treat_tz_as_pytz(tz) and pytz is not None:
        return pytz.timezone(str(tz))
    return tz
