# -*- coding: utf-8 -*-
# cython: profile=False
# cython: linetrace=False
# distutils: define_macros=CYTHON_TRACE=0
# distutils: define_macros=CYTHON_TRACE_NOGIL=0

cimport cython
from cython cimport Py_ssize_t

# dateutil compat
from dateutil.tz import (
    tzutc as _dateutil_tzutc,
    tzlocal as _dateutil_tzlocal,
    tzfile as _dateutil_tzfile)

import sys
if sys.platform == 'win32' or sys.platform == 'cygwin':
    # equiv pd.compat.is_platform_windows()
    from dateutil.zoneinfo import gettz as dateutil_gettz
else:
    from dateutil.tz import gettz as dateutil_gettz


from pytz.tzinfo import BaseTzInfo as _pytz_BaseTzInfo
import pytz
UTC = pytz.utc


import numpy as np
cimport numpy as np
from numpy cimport ndarray, int64_t
np.import_array()

# ----------------------------------------------------------------------
from util cimport is_string_object, is_integer_object, get_nat

cdef int64_t NPY_NAT = get_nat()

# ----------------------------------------------------------------------

cdef inline bint is_utc(object tz):
    return tz is UTC or isinstance(tz, _dateutil_tzutc)


cdef inline bint is_tzlocal(object tz):
    return isinstance(tz, _dateutil_tzlocal)


cdef inline bint treat_tz_as_pytz(object tz):
    return hasattr(tz, '_utc_transition_times') and hasattr(
        tz, '_transition_info')


cdef inline bint treat_tz_as_dateutil(object tz):
    return hasattr(tz, '_trans_list') and hasattr(tz, '_trans_idx')


cpdef inline object get_timezone(object tz):
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
    if is_utc(tz):
        return 'UTC'
    else:
        if treat_tz_as_dateutil(tz):
            if '.tar.gz' in tz._filename:
                raise ValueError(
                    'Bad tz filename. Dateutil on python 3 on windows has a '
                    'bug which causes tzfile._filename to be the same for all '
                    'timezone files. Please construct dateutil timezones '
                    'implicitly by passing a string like "dateutil/Europe'
                    '/London" when you construct your pandas objects instead '
                    'of passing a timezone object. See '
                    'https://github.com/pandas-dev/pandas/pull/7362')
            return 'dateutil/' + tz._filename
        else:
            # tz is a pytz timezone or unknown.
            try:
                zone = tz.zone
                if zone is None:
                    return tz
                return zone
            except AttributeError:
                return tz


cpdef inline object maybe_get_tz(object tz):
    """
    (Maybe) Construct a timezone object from a string. If tz is a string, use
    it to construct a timezone object. Otherwise, just return tz.
    """
    if is_string_object(tz):
        if tz == 'tzlocal()':
            tz = _dateutil_tzlocal()
        elif tz.startswith('dateutil/'):
            zone = tz[9:]
            tz = dateutil_gettz(zone)
            # On Python 3 on Windows, the filename is not always set correctly.
            if isinstance(tz, _dateutil_tzfile) and '.tar.gz' in tz._filename:
                tz._filename = zone
        else:
            tz = pytz.timezone(tz)
    elif is_integer_object(tz):
        tz = pytz.FixedOffset(tz / 60)
    return tz


def _p_tz_cache_key(tz):
    """ Python interface for cache function to facilitate testing."""
    return tz_cache_key(tz)


# Timezone data caches, key is the pytz string or dateutil file name.
dst_cache = {}


cdef inline object tz_cache_key(object tz):
    """
    Return the key in the cache for the timezone info object or None
    if unknown.

    The key is currently the tz string for pytz timezones, the filename for
    dateutil timezones.

    Notes
    =====
    This cannot just be the hash of a timezone object. Unfortunately, the
    hashes of two dateutil tz objects which represent the same timezone are
    not equal (even though the tz objects will compare equal and represent
    the same tz file). Also, pytz objects are not always hashable so we use
    str(tz) instead.
    """
    if isinstance(tz, _pytz_BaseTzInfo):
        return tz.zone
    elif isinstance(tz, _dateutil_tzfile):
        if '.tar.gz' in tz._filename:
            raise ValueError('Bad tz filename. Dateutil on python 3 on '
                             'windows has a bug which causes tzfile._filename '
                             'to be the same for all timezone files. Please '
                             'construct dateutil timezones implicitly by '
                             'passing a string like "dateutil/Europe/London" '
                             'when you construct your pandas objects instead '
                             'of passing a timezone object. See '
                             'https://github.com/pandas-dev/pandas/pull/7362')
        return 'dateutil' + tz._filename
    else:
        return None


# ----------------------------------------------------------------------
# UTC Offsets


cpdef get_utcoffset(tzinfo, obj):
    try:
        return tzinfo._utcoffset
    except AttributeError:
        return tzinfo.utcoffset(obj)


cdef inline bint is_fixed_offset(object tz):
    if treat_tz_as_dateutil(tz):
        if len(tz._trans_idx) == 0 and len(tz._trans_list) == 0:
            return 1
        else:
            return 0
    elif treat_tz_as_pytz(tz):
        if (len(tz._transition_info) == 0
                and len(tz._utc_transition_times) == 0):
            return 1
        else:
            return 0
    return 1


cdef object get_utc_trans_times_from_dateutil_tz(object tz):
    """
    Transition times in dateutil timezones are stored in local non-dst
    time.  This code converts them to UTC. It's the reverse of the code
    in dateutil.tz.tzfile.__init__.
    """
    new_trans = list(tz._trans_list)
    last_std_offset = 0
    for i, (trans, tti) in enumerate(zip(tz._trans_list, tz._trans_idx)):
        if not tti.isdst:
            last_std_offset = tti.offset
        new_trans[i] = trans - last_std_offset
    return new_trans


cpdef ndarray unbox_utcoffsets(object transinfo):
    cdef:
        Py_ssize_t i, sz
        ndarray[int64_t] arr

    sz = len(transinfo)
    arr = np.empty(sz, dtype='i8')

    for i in range(sz):
        arr[i] = int(transinfo[i][0].total_seconds()) * 1000000000

    return arr


# ----------------------------------------------------------------------
# Daylight Savings


cdef object get_dst_info(object tz):
    """
    return a tuple of :
      (UTC times of DST transitions,
       UTC offsets in microseconds corresponding to DST transitions,
       string of type of transitions)

    """
    cache_key = tz_cache_key(tz)
    if cache_key is None:
        num = int(get_utcoffset(tz, None).total_seconds()) * 1000000000
        return (np.array([NPY_NAT + 1], dtype=np.int64),
                np.array([num], dtype=np.int64),
                None)

    if cache_key not in dst_cache:
        if treat_tz_as_pytz(tz):
            trans = np.array(tz._utc_transition_times, dtype='M8[ns]')
            trans = trans.view('i8')
            try:
                if tz._utc_transition_times[0].year == 1:
                    trans[0] = NPY_NAT + 1
            except Exception:
                pass
            deltas = unbox_utcoffsets(tz._transition_info)
            typ = 'pytz'

        elif treat_tz_as_dateutil(tz):
            if len(tz._trans_list):
                # get utc trans times
                trans_list = get_utc_trans_times_from_dateutil_tz(tz)
                trans = np.hstack([
                    np.array([0], dtype='M8[s]'),  # place holder for 1st item
                    np.array(trans_list, dtype='M8[s]')]).astype(
                    'M8[ns]')  # all trans listed
                trans = trans.view('i8')
                trans[0] = NPY_NAT + 1

                # deltas
                deltas = np.array([v.offset for v in (
                    tz._ttinfo_before,) + tz._trans_idx], dtype='i8')
                deltas *= 1000000000
                typ = 'dateutil'

            elif is_fixed_offset(tz):
                trans = np.array([NPY_NAT + 1], dtype=np.int64)
                deltas = np.array([tz._ttinfo_std.offset],
                                  dtype='i8') * 1000000000
                typ = 'fixed'
            else:
                trans = np.array([], dtype='M8[ns]')
                deltas = np.array([], dtype='i8')
                typ = None

        else:
            # static tzinfo
            trans = np.array([NPY_NAT + 1], dtype=np.int64)
            num = int(get_utcoffset(tz, None).total_seconds()) * 1000000000
            deltas = np.array([num], dtype=np.int64)
            typ = 'static'

        dst_cache[cache_key] = (trans, deltas, typ)

    return dst_cache[cache_key]


def infer_tzinfo(start, end):
    if start is not None and end is not None:
        tz = start.tzinfo
        if not (get_timezone(tz) == get_timezone(end.tzinfo)):
            msg = 'Inputs must both have the same timezone, {tz1} != {tz2}'
            raise AssertionError(msg.format(tz1=tz, tz2=end.tzinfo))
    elif start is not None:
        tz = start.tzinfo
    elif end is not None:
        tz = end.tzinfo
    else:
        tz = None
    return tz
