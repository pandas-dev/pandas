# -*- coding: utf-8 -*-
# cython: profile=False

cimport cython
from cython cimport Py_ssize_t

import numpy as np
cimport numpy as np
from numpy cimport int64_t, ndarray
np.import_array()

import pytz

from cpython.datetime cimport datetime

from np_datetime cimport (check_dts_bounds,
                          pandas_datetimestruct,
                          dt64_to_dtstruct, dtstruct_to_dt64)

cimport util

from timezones cimport (
    is_utc, is_tzlocal, is_fixed_offset,
    treat_tz_as_dateutil, treat_tz_as_pytz,
    get_utcoffset, get_dst_info, get_timezone)

# ----------------------------------------------------------------------
# Constants
cdef int64_t NPY_NAT = util.get_nat()

cdef int64_t DAY_NS = 86400000000000LL

UTC = pytz.UTC


# ----------------------------------------------------------------------
# _TSObject Conversion

# lightweight C object to hold datetime & int64 pair
cdef class _TSObject:
    # cdef:
    #    pandas_datetimestruct dts      # pandas_datetimestruct
    #    int64_t value               # numpy dt64
    #    object tzinfo

    property value:
        def __get__(self):
            return self.value


cdef inline void _localize_tso(_TSObject obj, object tz):
    """
    Take a TSObject in UTC and localizes to timezone tz.
    """
    cdef:
        ndarray[int64_t] trans, deltas
        Py_ssize_t delta, posn

    if is_utc(tz):
        obj.tzinfo = tz
    elif is_tzlocal(tz):
        dt64_to_dtstruct(obj.value, &obj.dts)
        dt = datetime(obj.dts.year, obj.dts.month, obj.dts.day, obj.dts.hour,
                      obj.dts.min, obj.dts.sec, obj.dts.us, tz)
        delta = int(get_utcoffset(tz, dt).total_seconds()) * 1000000000
        if obj.value != NPY_NAT:
            dt64_to_dtstruct(obj.value + delta, &obj.dts)
        else:
            dt64_to_dtstruct(obj.value, &obj.dts)
        obj.tzinfo = tz
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = get_dst_info(tz)

        pos = trans.searchsorted(obj.value, side='right') - 1

        # static/pytz/dateutil specific code
        if is_fixed_offset(tz):
            # statictzinfo
            if len(deltas) > 0 and obj.value != NPY_NAT:
                dt64_to_dtstruct(obj.value + deltas[0], &obj.dts)
            else:
                dt64_to_dtstruct(obj.value, &obj.dts)
            obj.tzinfo = tz
        elif treat_tz_as_pytz(tz):
            inf = tz._transition_info[pos]
            if obj.value != NPY_NAT:
                dt64_to_dtstruct(obj.value + deltas[pos], &obj.dts)
            else:
                dt64_to_dtstruct(obj.value, &obj.dts)
            obj.tzinfo = tz._tzinfos[inf]
        elif treat_tz_as_dateutil(tz):
            if obj.value != NPY_NAT:
                dt64_to_dtstruct(obj.value + deltas[pos], &obj.dts)
            else:
                dt64_to_dtstruct(obj.value, &obj.dts)
            obj.tzinfo = tz
        else:
            obj.tzinfo = tz


# ----------------------------------------------------------------------
# Localization / Timezone Conversion


cpdef int64_t tz_convert_single(int64_t val, object tz1, object tz2):
    """
    Convert the val (in i8) from timezone1 to timezone2

    This is a single timezone versoin of tz_convert

    Parameters
    ----------
    val : int64
    tz1 : string / timezone object
    tz2 : string / timezone object

    Returns
    -------
    int64 converted

    """

    cdef:
        ndarray[int64_t] trans, deltas
        Py_ssize_t pos
        int64_t v, offset, utc_date
        pandas_datetimestruct dts

    if val == NPY_NAT:
        return val

    # Convert to UTC
    if is_tzlocal(tz1):
        dt64_to_dtstruct(val, &dts)
        dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                      dts.min, dts.sec, dts.us, tz1)
        delta = int(get_utcoffset(tz1, dt).total_seconds()) * 1000000000
        utc_date = val - delta
    elif get_timezone(tz1) != 'UTC':
        trans, deltas, typ = get_dst_info(tz1)
        pos = trans.searchsorted(val, side='right') - 1
        if pos < 0:
            raise ValueError('First time before start of DST info')
        offset = deltas[pos]
        utc_date = val - offset
    else:
        utc_date = val

    if get_timezone(tz2) == 'UTC':
        return utc_date
    if is_tzlocal(tz2):
        dt64_to_dtstruct(val, &dts)
        dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                      dts.min, dts.sec, dts.us, tz2)
        delta = int(get_utcoffset(tz2, dt).total_seconds()) * 1000000000
        return utc_date + delta

    # Convert UTC to other timezone
    trans, deltas, typ = get_dst_info(tz2)

    pos = trans.searchsorted(utc_date, side='right') - 1
    if pos < 0:
        raise ValueError('First time before start of DST info')

    offset = deltas[pos]
    return utc_date + offset


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_convert(ndarray[int64_t] vals, object tz1, object tz2):
    """
    Convert the values (in i8) from timezone1 to timezone2

    Parameters
    ----------
    vals : int64 ndarray
    tz1 : string / timezone object
    tz2 : string / timezone object

    Returns
    -------
    int64 ndarray of converted
    """

    cdef:
        ndarray[int64_t] utc_dates, tt, result, trans, deltas
        Py_ssize_t i, j, pos, n = len(vals)
        ndarray[Py_ssize_t] posn
        int64_t v, offset, delta
        pandas_datetimestruct dts

    if len(vals) == 0:
        return np.array([], dtype=np.int64)

    # Convert to UTC
    if get_timezone(tz1) != 'UTC':
        utc_dates = np.empty(n, dtype=np.int64)
        if is_tzlocal(tz1):
            for i in range(n):
                v = vals[i]
                if v == NPY_NAT:
                    utc_dates[i] = NPY_NAT
                else:
                    dt64_to_dtstruct(v, &dts)
                    dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                                  dts.min, dts.sec, dts.us, tz1)
                    delta = (int(get_utcoffset(tz1, dt).total_seconds())
                             * 1000000000)
                    utc_dates[i] = v - delta
        else:
            trans, deltas, typ = get_dst_info(tz1)

            # all-NaT
            tt = vals[vals != NPY_NAT]
            if not len(tt):
                return vals

            posn = trans.searchsorted(tt, side='right')
            j = 0
            for i in range(n):
                v = vals[i]
                if v == NPY_NAT:
                    utc_dates[i] = NPY_NAT
                else:
                    pos = posn[j] - 1
                    j = j + 1
                    if pos < 0:
                        raise ValueError('First time before start of DST info')
                    offset = deltas[pos]
                    utc_dates[i] = v - offset
    else:
        utc_dates = vals

    if get_timezone(tz2) == 'UTC':
        return utc_dates

    result = np.zeros(n, dtype=np.int64)
    if is_tzlocal(tz2):
        for i in range(n):
            v = utc_dates[i]
            if v == NPY_NAT:
                result[i] = NPY_NAT
            else:
                dt64_to_dtstruct(v, &dts)
                dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                              dts.min, dts.sec, dts.us, tz2)
                delta = (int(get_utcoffset(tz2, dt).total_seconds())
                             * 1000000000)
                result[i] = v + delta
        return result

    # Convert UTC to other timezone
    trans, deltas, typ = get_dst_info(tz2)

    # use first non-NaT element
    # if all-NaT, return all-NaT
    if (result == NPY_NAT).all():
        return result

    # if all NaT, return all NaT
    tt = utc_dates[utc_dates!=NPY_NAT]
    if not len(tt):
        return utc_dates

    posn = trans.searchsorted(tt, side='right')

    j = 0
    for i in range(n):
        v = utc_dates[i]
        if vals[i] == NPY_NAT:
            result[i] = vals[i]
        else:
            pos = posn[j] - 1
            j = j + 1
            if pos < 0:
                raise ValueError('First time before start of DST info')
            offset = deltas[pos]
            result[i] = v + offset
    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_localize_to_utc(ndarray[int64_t] vals, object tz, object ambiguous=None,
                       object errors='raise'):
    """
    Localize tzinfo-naive i8 to given time zone (using pytz). If
    there are ambiguities in the values, raise AmbiguousTimeError.

    Returns
    -------
    localized : DatetimeIndex
    """
    cdef:
        ndarray[int64_t] trans, deltas, idx_shifted
        ndarray ambiguous_array
        Py_ssize_t i, idx, pos, ntrans, n = len(vals)
        int64_t *tdata
        int64_t v, left, right
        ndarray[int64_t] result, result_a, result_b, dst_hours
        pandas_datetimestruct dts
        bint infer_dst = False, is_dst = False, fill = False
        bint is_coerce = errors == 'coerce', is_raise = errors == 'raise'

    # Vectorized version of DstTzInfo.localize

    assert is_coerce or is_raise

    if tz == UTC or tz is None:
        return vals

    result = np.empty(n, dtype=np.int64)

    if is_tzlocal(tz):
        for i in range(n):
            v = vals[i]
            dt64_to_dtstruct(v, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(get_utcoffset(tz, dt).total_seconds()) * 1000000000
            result[i] = v - delta
        return result

    if util.is_string_object(ambiguous):
        if ambiguous == 'infer':
            infer_dst = True
        elif ambiguous == 'NaT':
            fill = True
    elif isinstance(ambiguous, bool):
        is_dst = True
        if ambiguous:
            ambiguous_array = np.ones(len(vals), dtype=bool)
        else:
            ambiguous_array = np.zeros(len(vals), dtype=bool)
    elif hasattr(ambiguous, '__iter__'):
        is_dst = True
        if len(ambiguous) != len(vals):
            raise ValueError("Length of ambiguous bool-array must be "
                             "the same size as vals")
        ambiguous_array = np.asarray(ambiguous)

    trans, deltas, typ = get_dst_info(tz)

    tdata = <int64_t*> trans.data
    ntrans = len(trans)

    result_a = np.empty(n, dtype=np.int64)
    result_b = np.empty(n, dtype=np.int64)
    result_a.fill(NPY_NAT)
    result_b.fill(NPY_NAT)

    # left side
    idx_shifted = (np.maximum(0, trans.searchsorted(
        vals - DAY_NS, side='right') - 1)).astype(np.int64)

    for i in range(n):
        v = vals[i] - deltas[idx_shifted[i]]
        pos = bisect_right_i8(tdata, v, ntrans) - 1

        # timestamp falls to the left side of the DST transition
        if v + deltas[pos] == vals[i]:
            result_a[i] = v

    # right side
    idx_shifted = (np.maximum(0, trans.searchsorted(
        vals + DAY_NS, side='right') - 1)).astype(np.int64)

    for i in range(n):
        v = vals[i] - deltas[idx_shifted[i]]
        pos = bisect_right_i8(tdata, v, ntrans) - 1

        # timestamp falls to the right side of the DST transition
        if v + deltas[pos] == vals[i]:
            result_b[i] = v

    if infer_dst:
        dst_hours = np.empty(n, dtype=np.int64)
        dst_hours.fill(NPY_NAT)

        # Get the ambiguous hours (given the above, these are the hours
        # where result_a != result_b and neither of them are NAT)
        both_nat = np.logical_and(result_a != NPY_NAT, result_b != NPY_NAT)
        both_eq = result_a == result_b
        trans_idx = np.squeeze(np.nonzero(np.logical_and(both_nat, ~both_eq)))
        if trans_idx.size == 1:
            stamp = _render_tstamp(vals[trans_idx])
            raise pytz.AmbiguousTimeError(
                "Cannot infer dst time from %s as there "
                "are no repeated times" % stamp)
        # Split the array into contiguous chunks (where the difference between
        # indices is 1).  These are effectively dst transitions in different
        # years which is useful for checking that there is not an ambiguous
        # transition in an individual year.
        if trans_idx.size > 0:
            one_diff = np.where(np.diff(trans_idx) != 1)[0] +1
            trans_grp = np.array_split(trans_idx, one_diff)

            # Iterate through each day, if there are no hours where the
            # delta is negative (indicates a repeat of hour) the switch
            # cannot be inferred
            for grp in trans_grp:

                delta = np.diff(result_a[grp])
                if grp.size == 1 or np.all(delta > 0):
                    stamp = _render_tstamp(vals[grp[0]])
                    raise pytz.AmbiguousTimeError(stamp)

                # Find the index for the switch and pull from a for dst and b
                # for standard
                switch_idx = (delta <= 0).nonzero()[0]
                if switch_idx.size > 1:
                    raise pytz.AmbiguousTimeError(
                        "There are %i dst switches when "
                        "there should only be 1." % switch_idx.size)
                switch_idx = switch_idx[0] + 1 # Pull the only index and adjust
                a_idx = grp[:switch_idx]
                b_idx = grp[switch_idx:]
                dst_hours[grp] = np.hstack((result_a[a_idx], result_b[b_idx]))

    for i in range(n):
        left = result_a[i]
        right = result_b[i]
        if vals[i] == NPY_NAT:
            result[i] = vals[i]
        elif left != NPY_NAT and right != NPY_NAT:
            if left == right:
                result[i] = left
            else:
                if infer_dst and dst_hours[i] != NPY_NAT:
                    result[i] = dst_hours[i]
                elif is_dst:
                    if ambiguous_array[i]:
                        result[i] = left
                    else:
                        result[i] = right
                elif fill:
                    result[i] = NPY_NAT
                else:
                    stamp = _render_tstamp(vals[i])
                    raise pytz.AmbiguousTimeError(
                        "Cannot infer dst time from %r, try using the "
                        "'ambiguous' argument" % stamp)
        elif left != NPY_NAT:
            result[i] = left
        elif right != NPY_NAT:
            result[i] = right
        else:
            if is_coerce:
                result[i] = NPY_NAT
            else:
                stamp = _render_tstamp(vals[i])
                raise pytz.NonExistentTimeError(stamp)

    return result


cdef inline bisect_right_i8(int64_t *data, int64_t val, Py_ssize_t n):
    cdef Py_ssize_t pivot, left = 0, right = n

    assert n >= 1

    # edge cases
    if val > data[n - 1]:
        return n

    if val < data[0]:
        return 0

    while left < right:
        pivot = left + (right - left) // 2

        if data[pivot] <= val:
            left = pivot + 1
        else:
            right = pivot

    return left


cdef inline str _render_tstamp(int64_t val):
    """ Helper function to render exception messages"""
    from pandas._libs.tslib import Timestamp
    return str(Timestamp(val))
