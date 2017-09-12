# -*- coding: utf-8 -*-
# cython: profile=False

# dateutil compat
from dateutil.tz import (
	tzutc as _dateutil_tzutc,
	tzlocal as _dateutil_tzlocal)

import pytz
UTC = pytz.utc


cdef inline bint _is_utc(object tz):
    return tz is UTC or isinstance(tz, _dateutil_tzutc)


cdef inline bint _is_tzlocal(object tz):
    return isinstance(tz, _dateutil_tzlocal)


cdef inline bint _treat_tz_as_pytz(object tz):
    return hasattr(tz, '_utc_transition_times') and hasattr(
        tz, '_transition_info')


cdef inline bint _treat_tz_as_dateutil(object tz):
    return hasattr(tz, '_trans_list') and hasattr(tz, '_trans_idx')


cdef inline object _get_zone(object tz):
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
    if _is_utc(tz):
        return 'UTC'
    else:
        if _treat_tz_as_dateutil(tz):
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


def get_timezone(tz):
    return _get_zone(tz)

#----------------------------------------------------------------------
# UTC Offsets

cpdef _get_utcoffset(tzinfo, obj):
    try:
        return tzinfo._utcoffset
    except AttributeError:
        return tzinfo.utcoffset(obj)
