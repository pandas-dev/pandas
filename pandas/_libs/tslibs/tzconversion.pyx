"""
timezone conversion
"""
import cython
from cython import Py_ssize_t

from cpython.datetime cimport (
    PyDateTime_IMPORT, PyDelta_Check, datetime, tzinfo)
PyDateTime_IMPORT

import pytz
from dateutil.tz import tzutc

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, int64_t, uint8_t, intp_t
cnp.import_array()

from pandas._libs.tslibs.ccalendar import DAY_SECONDS, HOUR_SECONDS
from pandas._libs.tslibs.nattype cimport NPY_NAT
from pandas._libs.tslibs.np_datetime cimport (
    npy_datetimestruct, dt64_to_dtstruct)
from pandas._libs.tslibs.timedeltas cimport delta_to_nanoseconds
from pandas._libs.tslibs.timezones cimport (
    get_dst_info, is_tzlocal, is_utc, get_timezone, get_utcoffset)


# TODO: cdef scalar version to call from convert_str_to_tsobject
@cython.boundscheck(False)
@cython.wraparound(False)
def tz_localize_to_utc(ndarray[int64_t] vals, object tz, object ambiguous=None,
                       object nonexistent=None):
    """
    Localize tzinfo-naive i8 to given time zone (using pytz). If
    there are ambiguities in the values, raise AmbiguousTimeError.

    Parameters
    ----------
    vals : ndarray[int64_t]
    tz : tzinfo or None
    ambiguous : str, bool, or arraylike
        When clocks moved backward due to DST, ambiguous times may arise.
        For example in Central European Time (UTC+01), when going from 03:00
        DST to 02:00 non-DST, 02:30:00 local time occurs both at 00:30:00 UTC
        and at 01:30:00 UTC. In such a situation, the `ambiguous` parameter
        dictates how ambiguous times should be handled.

        - 'infer' will attempt to infer fall dst-transition hours based on
          order
        - bool-ndarray where True signifies a DST time, False signifies a
          non-DST time (note that this flag is only applicable for ambiguous
          times, but the array must have the same length as vals)
        - bool if True, treat all vals as DST. If False, treat them as non-DST
        - 'NaT' will return NaT where there are ambiguous times

    nonexistent : {None, "NaT", "shift_forward", "shift_backward", "raise", \
timedelta-like}
        How to handle non-existent times when converting wall times to UTC

        .. versionadded:: 0.24.0

    Returns
    -------
    localized : ndarray[int64_t]
    """
    cdef:
        int64_t[:] deltas, idx_shifted, idx_shifted_left, idx_shifted_right
        ndarray[uint8_t, cast=True] ambiguous_array, both_nat, both_eq
        Py_ssize_t i, idx, pos, ntrans, n = len(vals)
        Py_ssize_t delta_idx_offset, delta_idx, pos_left, pos_right
        int64_t *tdata
        int64_t v, left, right, val, v_left, v_right, new_local, remaining_mins
        int64_t first_delta
        int64_t HOURS_NS = HOUR_SECONDS * 1000000000, shift_delta = 0
        ndarray[int64_t] trans, result, result_a, result_b, dst_hours, delta
        ndarray trans_idx, grp, a_idx, b_idx, one_diff
        npy_datetimestruct dts
        bint infer_dst = False, is_dst = False, fill = False
        bint shift_forward = False, shift_backward = False
        bint fill_nonexist = False
        list trans_grp
        str stamp

    # Vectorized version of DstTzInfo.localize
    if is_utc(tz) or tz is None:
        return vals

    result = np.empty(n, dtype=np.int64)

    if is_tzlocal(tz):
        for i in range(n):
            v = vals[i]
            if v == NPY_NAT:
                result[i] = NPY_NAT
            else:
                result[i] = _tz_convert_tzlocal_utc(v, tz, to_utc=True)
        return result

    # silence false-positive compiler warning
    ambiguous_array = np.empty(0, dtype=bool)
    if isinstance(ambiguous, str):
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
        ambiguous_array = np.asarray(ambiguous, dtype=bool)

    if nonexistent == 'NaT':
        fill_nonexist = True
    elif nonexistent == 'shift_forward':
        shift_forward = True
    elif nonexistent == 'shift_backward':
        shift_backward = True
    elif PyDelta_Check(nonexistent):
        shift_delta = delta_to_nanoseconds(nonexistent)
    elif nonexistent not in ('raise', None):
        msg = ("nonexistent must be one of {'NaT', 'raise', 'shift_forward', "
               "shift_backwards} or a timedelta object")
        raise ValueError(msg)

    trans, deltas, _ = get_dst_info(tz)

    tdata = <int64_t*>cnp.PyArray_DATA(trans)
    ntrans = len(trans)

    # Determine whether each date lies left of the DST transition (store in
    # result_a) or right of the DST transition (store in result_b)
    result_a = np.empty(n, dtype=np.int64)
    result_b = np.empty(n, dtype=np.int64)
    result_a[:] = NPY_NAT
    result_b[:] = NPY_NAT

    idx_shifted_left = (np.maximum(0, trans.searchsorted(
        vals - DAY_SECONDS * 1000000000, side='right') - 1)).astype(np.int64)

    idx_shifted_right = (np.maximum(0, trans.searchsorted(
        vals + DAY_SECONDS * 1000000000, side='right') - 1)).astype(np.int64)

    for i in range(n):
        val = vals[i]
        v_left = val - deltas[idx_shifted_left[i]]
        pos_left = bisect_right_i8(tdata, v_left, ntrans) - 1
        # timestamp falls to the left side of the DST transition
        if v_left + deltas[pos_left] == val:
            result_a[i] = v_left

        v_right = val - deltas[idx_shifted_right[i]]
        pos_right = bisect_right_i8(tdata, v_right, ntrans) - 1
        # timestamp falls to the right side of the DST transition
        if v_right + deltas[pos_right] == val:
            result_b[i] = v_right

    # silence false-positive compiler warning
    dst_hours = np.empty(0, dtype=np.int64)
    if infer_dst:
        dst_hours = np.empty(n, dtype=np.int64)
        dst_hours[:] = NPY_NAT

        # Get the ambiguous hours (given the above, these are the hours
        # where result_a != result_b and neither of them are NAT)
        both_nat = np.logical_and(result_a != NPY_NAT, result_b != NPY_NAT)
        both_eq = result_a == result_b
        trans_idx = np.squeeze(np.nonzero(np.logical_and(both_nat, ~both_eq)))
        if trans_idx.size == 1:
            stamp = _render_tstamp(vals[trans_idx])
            raise pytz.AmbiguousTimeError(
                f"Cannot infer dst time from {stamp} as there "
                f"are no repeated times")
        # Split the array into contiguous chunks (where the difference between
        # indices is 1).  These are effectively dst transitions in different
        # years which is useful for checking that there is not an ambiguous
        # transition in an individual year.
        if trans_idx.size > 0:
            one_diff = np.where(np.diff(trans_idx) != 1)[0] + 1
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
                        f"There are {switch_idx.size} dst switches when "
                        f"there should only be 1.")
                switch_idx = switch_idx[0] + 1
                # Pull the only index and adjust
                a_idx = grp[:switch_idx]
                b_idx = grp[switch_idx:]
                dst_hours[grp] = np.hstack((result_a[a_idx], result_b[b_idx]))

    for i in range(n):
        val = vals[i]
        left = result_a[i]
        right = result_b[i]
        if val == NPY_NAT:
            result[i] = val
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
                    stamp = _render_tstamp(val)
                    raise pytz.AmbiguousTimeError(
                        f"Cannot infer dst time from {stamp}, try using the "
                        f"'ambiguous' argument")
        elif left != NPY_NAT:
            result[i] = left
        elif right != NPY_NAT:
            result[i] = right
        else:
            # Handle nonexistent times
            if shift_forward or shift_backward or shift_delta != 0:
                # Shift the nonexistent time to the closest existing time
                remaining_mins = val % HOURS_NS
                if shift_delta != 0:
                    # Validate that we don't relocalize on another nonexistent
                    # time
                    if -1 < shift_delta + remaining_mins < HOURS_NS:
                        raise ValueError(
                            f"The provided timedelta will relocalize on a "
                            f"nonexistent time: {nonexistent}"
                        )
                    new_local = val + shift_delta
                elif shift_forward:
                    new_local = val + (HOURS_NS - remaining_mins)
                else:
                    # Subtract 1 since the beginning hour is _inclusive_ of
                    # nonexistent times
                    new_local = val - remaining_mins - 1
                delta_idx = trans.searchsorted(new_local, side='right')
                # Shift the delta_idx by if the UTC offset of
                # the target tz is greater than 0 and we're moving forward
                # or vice versa
                first_delta = deltas[0]
                if (shift_forward or shift_delta > 0) and first_delta > 0:
                    delta_idx_offset = 1
                elif (shift_backward or shift_delta < 0) and first_delta < 0:
                    delta_idx_offset = 1
                else:
                    delta_idx_offset = 0
                delta_idx = delta_idx - delta_idx_offset
                result[i] = new_local - deltas[delta_idx]
            elif fill_nonexist:
                result[i] = NPY_NAT
            else:
                stamp = _render_tstamp(val)
                raise pytz.NonExistentTimeError(stamp)

    return result


cdef inline Py_ssize_t bisect_right_i8(int64_t *data,
                                       int64_t val, Py_ssize_t n):
    cdef:
        Py_ssize_t pivot, left = 0, right = n

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
    from pandas._libs.tslibs.timestamps import Timestamp
    return str(Timestamp(val))


# ----------------------------------------------------------------------
# Timezone Conversion

cdef int64_t tz_convert_utc_to_tzlocal(int64_t utc_val, tzinfo tz):
    """
    Parameters
    ----------
    utc_val : int64_t
    tz : tzinfo

    Returns
    -------
    local_val : int64_t
    """
    return _tz_convert_tzlocal_utc(utc_val, tz, to_utc=False)


cpdef int64_t tz_convert_single(int64_t val, object tz1, object tz2):
    """
    Convert the val (in i8) from timezone1 to timezone2

    This is a single timezone version of tz_convert

    Parameters
    ----------
    val : int64
    tz1 : string / timezone object
    tz2 : string / timezone object

    Returns
    -------
    converted: int64
    """
    cdef:
        int64_t[:] deltas
        Py_ssize_t pos
        int64_t v, offset, utc_date
        npy_datetimestruct dts
        int64_t arr[1]

    # See GH#17734 We should always be converting either from UTC or to UTC
    assert is_utc(tz1) or is_utc(tz2)

    if val == NPY_NAT:
        return val

    # Convert to UTC
    if is_tzlocal(tz1):
        utc_date = _tz_convert_tzlocal_utc(val, tz1, to_utc=True)
    elif not is_utc(get_timezone(tz1)):
        arr[0] = val
        utc_date = _tz_convert_dst(arr, tz1, to_utc=True)[0]
    else:
        utc_date = val

    if is_utc(get_timezone(tz2)):
        return utc_date
    elif is_tzlocal(tz2):
        return _tz_convert_tzlocal_utc(utc_date, tz2, to_utc=False)
    else:
        # Convert UTC to other timezone
        arr[0] = utc_date
        # Note: at least with cython 0.28.3, doing a lookup `[0]` in the next
        # line is sensitive to the declared return type of _tz_convert_dst;
        # if it is declared as returning ndarray[int64_t], a compile-time error
        # is raised.
        return _tz_convert_dst(arr, tz2, to_utc=False)[0]


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_convert(int64_t[:] vals, object tz1, object tz2):
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
        int64_t[:] utc_dates, converted

    if len(vals) == 0:
        return np.array([], dtype=np.int64)

    # Convert to UTC
    utc_dates = _tz_convert_one_way(vals, tz1, to_utc=True)
    converted = _tz_convert_one_way(utc_dates, tz2, to_utc=False)
    return np.array(converted, dtype=np.int64)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int64_t[:] _tz_convert_one_way(int64_t[:] vals, object tz, bint to_utc):
    """
    Convert the given values (in i8) either to UTC or from UTC.

    Parameters
    ----------
    vals : int64 ndarray
    tz1 : string / timezone object
    to_utc : bint

    Returns
    -------
    converted : ndarray[int64_t]
    """
    cdef:
        int64_t[:] converted, result
        Py_ssize_t i, n = len(vals)
        int64_t val

    if not is_utc(get_timezone(tz)):
        converted = np.empty(n, dtype=np.int64)
        if is_tzlocal(tz):
            for i in range(n):
                val = vals[i]
                if val == NPY_NAT:
                    converted[i] = NPY_NAT
                else:
                    converted[i] = _tz_convert_tzlocal_utc(val, tz, to_utc)
        else:
            converted = _tz_convert_dst(vals, tz, to_utc)
    else:
        converted = vals

    return converted


cdef inline int64_t _tzlocal_get_offset_components(int64_t val, tzinfo tz,
                                                   bint to_utc,
                                                   bint *fold=NULL):
    """
    Calculate offset in nanoseconds needed to convert the i8 representation of
    a datetime from a tzlocal timezone to UTC, or vice-versa.

    Parameters
    ----------
    val : int64_t
    tz : tzinfo
    to_utc : bint
        True if converting tzlocal _to_ UTC, False if going the other direction
    fold : bint*, default NULL
        pointer to fold: whether datetime ends up in a fold or not
        after adjustment

    Returns
    -------
    delta : int64_t

    Notes
    -----
    Sets fold by pointer
    """
    cdef:
        npy_datetimestruct dts
        datetime dt
        int64_t delta

    dt64_to_dtstruct(val, &dts)
    dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                  dts.min, dts.sec, dts.us)
    # get_utcoffset (tz.utcoffset under the hood) only makes sense if datetime
    # is _wall time_, so if val is a UTC timestamp convert to wall time
    if not to_utc:
        dt = dt.replace(tzinfo=tzutc())
        dt = dt.astimezone(tz)

    if fold is not NULL:
        fold[0] = dt.fold

    return int(get_utcoffset(tz, dt).total_seconds()) * 1000000000


cdef int64_t _tz_convert_tzlocal_utc(int64_t val, tzinfo tz, bint to_utc=True):
    """
    Convert the i8 representation of a datetime from a tzlocal timezone to
    UTC, or vice-versa.

    Private, not intended for use outside of tslibs.conversion

    Parameters
    ----------
    val : int64_t
    tz : tzinfo
    to_utc : bint
        True if converting tzlocal _to_ UTC, False if going the other direction

    Returns
    -------
    result : int64_t
    """
    cdef int64_t delta

    delta = _tzlocal_get_offset_components(val, tz, to_utc, NULL)

    if to_utc:
        return val - delta
    else:
        return val + delta


cdef int64_t _tz_convert_tzlocal_fromutc(int64_t val, tzinfo tz, bint *fold):
    """
    Convert the i8 representation of a datetime from UTC to local timezone,
    set fold by pointer

    Private, not intended for use outside of tslibs.conversion

    Parameters
    ----------
    val : int64_t
    tz : tzinfo
    fold : bint*
        pointer to fold: whether datetime ends up in a fold or not
        after adjustment

    Returns
    -------
    result : int64_t

    Notes
    -----
    Sets fold by pointer
    """
    cdef int64_t delta

    delta = _tzlocal_get_offset_components(val, tz, False, fold)

    return val + delta


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int64_t[:] _tz_convert_dst(
    const int64_t[:] values, tzinfo tz, bint to_utc=True,
):
    """
    tz_convert for non-UTC non-tzlocal cases where we have to check
    DST transitions pointwise.

    Parameters
    ----------
    values : ndarray[int64_t]
    tz : tzinfo
    to_utc : bool
        True if converting _to_ UTC, False if converting _from_ utc

    Returns
    -------
    result : ndarray[int64_t]
    """
    cdef:
        Py_ssize_t n = len(values)
        Py_ssize_t i
        intp_t[:] pos
        int64_t[:] result = np.empty(n, dtype=np.int64)
        ndarray[int64_t] trans
        int64_t[:] deltas
        int64_t v
        bint tz_is_local

    tz_is_local = is_tzlocal(tz)

    if not tz_is_local:
        # get_dst_info cannot extract offsets from tzlocal because its
        # dependent on a datetime
        trans, deltas, _ = get_dst_info(tz)
        if not to_utc:
            # We add `offset` below instead of subtracting it
            deltas = -1 * np.array(deltas, dtype='i8')

        # Previously, this search was done pointwise to try and benefit
        # from getting to skip searches for iNaTs. However, it seems call
        # overhead dominates the search time so doing it once in bulk
        # is substantially faster (GH#24603)
        pos = trans.searchsorted(values, side='right') - 1

    for i in range(n):
        v = values[i]
        if v == NPY_NAT:
            result[i] = v
        elif tz_is_local:
            result[i] = _tz_convert_tzlocal_utc(v, tz, to_utc=to_utc)
        else:
            if pos[i] < 0:
                raise ValueError('First time before start of DST info')
            result[i] = v - deltas[pos[i]]

    return result
