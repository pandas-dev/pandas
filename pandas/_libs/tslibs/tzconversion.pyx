"""
timezone conversion
"""
import cython
from cython import Py_ssize_t

from cpython.datetime cimport (
    PyDelta_Check,
    datetime,
    datetime_new,
    import_datetime,
    timedelta,
    tzinfo,
)

import_datetime()

import numpy as np
import pytz

cimport numpy as cnp
from numpy cimport (
    int64_t,
    intp_t,
    ndarray,
    uint8_t,
)

cnp.import_array()

from pandas._libs.tslibs.ccalendar cimport (
    DAY_NANOS,
    HOUR_NANOS,
)
from pandas._libs.tslibs.nattype cimport NPY_NAT
from pandas._libs.tslibs.np_datetime cimport (
    dt64_to_dtstruct,
    npy_datetimestruct,
)
from pandas._libs.tslibs.timezones cimport (
    get_dst_info,
    is_fixed_offset,
    is_tzlocal,
    is_utc,
    utc_pytz,
)


cdef int64_t tz_localize_to_utc_single(
    int64_t val, tzinfo tz, object ambiguous=None, object nonexistent=None,
) except? -1:
    """See tz_localize_to_utc.__doc__"""
    cdef:
        int64_t delta
        int64_t[::1] deltas

    if val == NPY_NAT:
        return val

    elif is_utc(tz) or tz is None:
        return val

    elif is_tzlocal(tz):
        return val - _tz_localize_using_tzinfo_api(val, tz, to_utc=True)

    elif is_fixed_offset(tz):
        _, deltas, _ = get_dst_info(tz)
        delta = deltas[0]
        return val - delta

    else:
        return tz_localize_to_utc(
            np.array([val], dtype="i8"),
            tz,
            ambiguous=ambiguous,
            nonexistent=nonexistent,
        )[0]


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_localize_to_utc(ndarray[int64_t] vals, tzinfo tz, object ambiguous=None,
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

    Returns
    -------
    localized : ndarray[int64_t]
    """
    cdef:
        const int64_t[::1] deltas
        ndarray[uint8_t, cast=True] ambiguous_array
        Py_ssize_t i, isl, isr, idx, pos, ntrans, n = vals.shape[0]
        Py_ssize_t delta_idx_offset, delta_idx, pos_left, pos_right
        int64_t *tdata
        int64_t v, left, right, val, v_left, v_right, new_local, remaining_mins
        int64_t first_delta, delta
        int64_t shift_delta = 0
        ndarray[int64_t] trans, result_a, result_b, dst_hours
        int64_t[::1] result
        npy_datetimestruct dts
        bint infer_dst = False, is_dst = False, fill = False
        bint shift_forward = False, shift_backward = False
        bint fill_nonexist = False
        str stamp

    # Vectorized version of DstTzInfo.localize
    if is_utc(tz) or tz is None:
        return vals.copy()

    result = np.empty(n, dtype=np.int64)

    if is_tzlocal(tz):
        for i in range(n):
            v = vals[i]
            if v == NPY_NAT:
                result[i] = NPY_NAT
            else:
                result[i] = v - _tz_localize_using_tzinfo_api(v, tz, to_utc=True)
        return result.base  # to return underlying ndarray

    elif is_fixed_offset(tz):
        _, deltas, _ = get_dst_info(tz)
        delta = deltas[0]
        for i in range(n):
            v = vals[i]
            if v == NPY_NAT:
                result[i] = NPY_NAT
            else:
                result[i] = v - delta
        return result.base  # to return underlying ndarray

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
        from .timedeltas import delta_to_nanoseconds
        shift_delta = delta_to_nanoseconds(nonexistent)
    elif nonexistent not in ('raise', None):
        msg = ("nonexistent must be one of {'NaT', 'raise', 'shift_forward', "
               "shift_backwards} or a timedelta object")
        raise ValueError(msg)

    trans, deltas, _ = get_dst_info(tz)

    tdata = <int64_t*>cnp.PyArray_DATA(trans)
    ntrans = trans.shape[0]

    # Determine whether each date lies left of the DST transition (store in
    # result_a) or right of the DST transition (store in result_b)
    result_a = np.empty(n, dtype=np.int64)
    result_b = np.empty(n, dtype=np.int64)
    result_a[:] = NPY_NAT
    result_b[:] = NPY_NAT

    for i in range(n):
        # This loops resembles the "Find the two best possibilities" block
        #  in pytz's DstTZInfo.localize method.
        val = vals[i]
        if val == NPY_NAT:
            continue

        # TODO: be careful of overflow in val-DAY_NANOS
        isl = bisect_right_i8(tdata, val - DAY_NANOS, ntrans) - 1
        if isl < 0:
            isl = 0

        v_left = val - deltas[isl]
        pos_left = bisect_right_i8(tdata, v_left, ntrans) - 1
        # timestamp falls to the left side of the DST transition
        if v_left + deltas[pos_left] == val:
            result_a[i] = v_left

        # TODO: be careful of overflow in val+DAY_NANOS
        isr = bisect_right_i8(tdata, val + DAY_NANOS, ntrans) - 1
        if isr < 0:
            isr = 0

        v_right = val - deltas[isr]
        pos_right = bisect_right_i8(tdata, v_right, ntrans) - 1
        # timestamp falls to the right side of the DST transition
        if v_right + deltas[pos_right] == val:
            result_b[i] = v_right

    # silence false-positive compiler warning
    dst_hours = np.empty(0, dtype=np.int64)
    if infer_dst:
        dst_hours = _get_dst_hours(vals, result_a, result_b)

    # Pre-compute delta_idx_offset that will be used if we go down non-existent
    #  paths.
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
                        "'ambiguous' argument"
                    )
        elif left != NPY_NAT:
            result[i] = left
        elif right != NPY_NAT:
            result[i] = right
        else:
            # Handle nonexistent times
            if shift_forward or shift_backward or shift_delta != 0:
                # Shift the nonexistent time to the closest existing time
                remaining_mins = val % HOUR_NANOS
                if shift_delta != 0:
                    # Validate that we don't relocalize on another nonexistent
                    # time
                    if -1 < shift_delta + remaining_mins < HOUR_NANOS:
                        raise ValueError(
                            "The provided timedelta will relocalize on a "
                            f"nonexistent time: {nonexistent}"
                        )
                    new_local = val + shift_delta
                elif shift_forward:
                    new_local = val + (HOUR_NANOS - remaining_mins)
                else:
                    # Subtract 1 since the beginning hour is _inclusive_ of
                    # nonexistent times
                    new_local = val - remaining_mins - 1

                delta_idx = bisect_right_i8(tdata, new_local, ntrans)

                delta_idx = delta_idx - delta_idx_offset
                result[i] = new_local - deltas[delta_idx]
            elif fill_nonexist:
                result[i] = NPY_NAT
            else:
                stamp = _render_tstamp(val)
                raise pytz.NonExistentTimeError(stamp)

    return result.base  # .base to get underlying ndarray


cdef inline Py_ssize_t bisect_right_i8(int64_t *data,
                                       int64_t val, Py_ssize_t n):
    # Caller is responsible for checking n > 0
    # This looks very similar to local_search_right in the ndarray.searchsorted
    #  implementation.
    cdef:
        Py_ssize_t pivot, left = 0, right = n

    # edge cases
    if val > data[n - 1]:
        return n

    # Caller is responsible for ensuring 'val >= data[0]'. This is
    #  ensured by the fact that 'data' comes from get_dst_info where data[0]
    #  is *always* NPY_NAT+1. If that ever changes, we will need to restore
    #  the following disabled check.
    # if val < data[0]:
    #    return 0

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


cdef ndarray[int64_t] _get_dst_hours(
    # vals only needed here to potential render an exception message
    const int64_t[:] vals,
    ndarray[int64_t] result_a,
    ndarray[int64_t] result_b,
):
    cdef:
        Py_ssize_t i, n = vals.shape[0]
        ndarray[uint8_t, cast=True] mismatch
        ndarray[int64_t] delta, dst_hours
        ndarray[intp_t] switch_idxs, trans_idx, grp, a_idx, b_idx, one_diff
        list trans_grp
        intp_t switch_idx
        int64_t left, right

    dst_hours = np.empty(n, dtype=np.int64)
    dst_hours[:] = NPY_NAT

    mismatch = np.zeros(n, dtype=bool)

    for i in range(n):
        left = result_a[i]
        right = result_b[i]

        # Get the ambiguous hours (given the above, these are the hours
        # where result_a != result_b and neither of them are NAT)
        if left != right and left != NPY_NAT and right != NPY_NAT:
            mismatch[i] = 1

    trans_idx = mismatch.nonzero()[0]

    if trans_idx.size == 1:
        stamp = _render_tstamp(vals[trans_idx[0]])
        raise pytz.AmbiguousTimeError(
            f"Cannot infer dst time from {stamp} as there "
            "are no repeated times"
        )

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
            switch_idxs = (delta <= 0).nonzero()[0]
            if switch_idxs.size > 1:
                raise pytz.AmbiguousTimeError(
                    f"There are {switch_idxs.size} dst switches when "
                    "there should only be 1."
                )

            switch_idx = switch_idxs[0] + 1
            # Pull the only index and adjust
            a_idx = grp[:switch_idx]
            b_idx = grp[switch_idx:]
            dst_hours[grp] = np.hstack((result_a[a_idx], result_b[b_idx]))

    return dst_hours


# ----------------------------------------------------------------------
# Timezone Conversion

cdef int64_t localize_tzinfo_api(
    int64_t utc_val, tzinfo tz, bint* fold=NULL
) except? -1:
    """
    Parameters
    ----------
    utc_val : int64_t
    tz : tzinfo
    fold : bint*
        pointer to fold: whether datetime ends up in a fold or not
        after adjustment

    Returns
    -------
    delta : int64_t
        Value to add when converting from utc.
    """
    return _tz_localize_using_tzinfo_api(utc_val, tz, to_utc=False, fold=fold)


def py_tz_convert_from_utc_single(int64_t utc_val, tzinfo tz):
    # The 'bint* fold=NULL' in tz_convert_from_utc_single means we cannot
    #  make it cdef, so this is version exposed for testing from python.
    return tz_convert_from_utc_single(utc_val, tz)


cdef int64_t tz_convert_from_utc_single(
    int64_t utc_val,
    tzinfo tz,
    bint* fold=NULL,
    Py_ssize_t* outpos=NULL,
) except? -1:
    """
    Convert the val (in i8) from UTC to tz

    This is a single value version of tz_convert_from_utc.

    Parameters
    ----------
    utc_val : int64
    tz : tzinfo
    fold : bint*, default NULL
    outpos : Py_ssize_t*, default NULL

    Returns
    -------
    converted: int64
    """
    cdef:
        int64_t delta
        int64_t[::1] deltas
        ndarray[int64_t, ndim=1] trans
        int64_t* tdata
        intp_t pos

    if utc_val == NPY_NAT:
        return utc_val

    if is_utc(tz):
        return utc_val
    elif is_tzlocal(tz):
        return utc_val + _tz_localize_using_tzinfo_api(utc_val, tz, to_utc=False)
    else:
        trans, deltas, typ = get_dst_info(tz)
        tdata = <int64_t*>cnp.PyArray_DATA(trans)

        if typ == "dateutil":
            pos = bisect_right_i8(tdata, utc_val, trans.shape[0]) - 1

            if fold is not NULL:
                fold[0] = infer_dateutil_fold(utc_val, trans, deltas, pos)
            return utc_val + deltas[pos]

        elif typ == "pytz":
            pos = bisect_right_i8(tdata, utc_val, trans.shape[0]) - 1

            # We need to get 'pos' back to the caller so it can pick the
            #  correct "standardized" tzinfo objecg.
            if outpos is not NULL:
                outpos[0] = pos
            return utc_val + deltas[pos]

        else:
            # All other cases have len(deltas) == 1. As of 2018-07-17
            #  (and 2022-03-07), all test cases that get here have
            #  is_fixed_offset(tz).
            return utc_val + deltas[0]


def tz_convert_from_utc(const int64_t[:] vals, tzinfo tz):
    """
    Convert the values (in i8) from UTC to tz

    Parameters
    ----------
    vals : int64 ndarray
    tz : tzinfo

    Returns
    -------
    int64 ndarray of converted
    """
    cdef:
        const int64_t[:] converted

    if vals.shape[0] == 0:
        return np.array([], dtype=np.int64)

    converted = _tz_convert_from_utc(vals, tz)
    return np.asarray(converted, dtype=np.int64)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef const int64_t[:] _tz_convert_from_utc(const int64_t[:] stamps, tzinfo tz):
    """
    Convert the given values (in i8) either to UTC or from UTC.

    Parameters
    ----------
    stamps : int64 ndarray
    tz : tzinfo

    Returns
    -------
    converted : ndarray[int64_t]
    """
    cdef:
        Py_ssize_t i, ntrans = -1, n = stamps.shape[0]
        ndarray[int64_t] trans
        int64_t[::1] deltas
        int64_t* tdata = NULL
        intp_t pos
        int64_t utc_val, local_val, delta = NPY_NAT
        bint use_utc = False, use_tzlocal = False, use_fixed = False
        str typ

        int64_t[::1] result

    if is_utc(tz) or tz is None:
        # Much faster than going through the "standard" pattern below
        return stamps.copy()

    if is_utc(tz) or tz is None:
        use_utc = True
    elif is_tzlocal(tz):
        use_tzlocal = True
    else:
        trans, deltas, typ = get_dst_info(tz)
        ntrans = trans.shape[0]
        if typ not in ["pytz", "dateutil"]:
            # static/fixed; in this case we know that len(delta) == 1
            use_fixed = True
            delta = deltas[0]
        else:
            tdata = <int64_t*>cnp.PyArray_DATA(trans)

    result = np.empty(n, dtype=np.int64)

    for i in range(n):
        utc_val = stamps[i]
        if utc_val == NPY_NAT:
            result[i] = NPY_NAT
            continue

        # The pattern used in vectorized.pyx checks for use_utc here,
        #  but we handle that case above.
        if use_tzlocal:
            local_val = utc_val + _tz_localize_using_tzinfo_api(utc_val, tz, to_utc=False)
        elif use_fixed:
            local_val = utc_val + delta
        else:
            pos = bisect_right_i8(tdata, utc_val, ntrans) - 1
            local_val = utc_val + deltas[pos]

        result[i] = local_val

    return result


# OSError may be thrown by tzlocal on windows at or close to 1970-01-01
#  see https://github.com/pandas-dev/pandas/pull/37591#issuecomment-720628241
cdef int64_t _tz_localize_using_tzinfo_api(
    int64_t val, tzinfo tz, bint to_utc=True, bint* fold=NULL
) except? -1:
    """
    Convert the i8 representation of a datetime from a general-case timezone to
    UTC, or vice-versa using the datetime/tzinfo API.

    Private, not intended for use outside of tslibs.tzconversion.

    Parameters
    ----------
    val : int64_t
    tz : tzinfo
    to_utc : bint
        True if converting _to_ UTC, False if going the other direction.
    fold : bint*, default NULL
        pointer to fold: whether datetime ends up in a fold or not
        after adjustment.
        Only passed with to_utc=False.

    Returns
    -------
    delta : int64_t
        Value to add when converting from utc, subtract when converting to utc.

    Notes
    -----
    Sets fold by pointer
    """
    cdef:
        npy_datetimestruct dts
        datetime dt
        int64_t delta
        timedelta td

    dt64_to_dtstruct(val, &dts)

    # datetime_new is cython-optimized constructor
    if not to_utc:
        # tz.utcoffset only makes sense if datetime
        # is _wall time_, so if val is a UTC timestamp convert to wall time
        dt = datetime_new(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, utc_pytz)
        dt = dt.astimezone(tz)

        if fold is not NULL:
            # NB: fold is only passed with to_utc=False
            fold[0] = dt.fold
    else:
        dt = datetime_new(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, None)

    td = tz.utcoffset(dt)
    delta = int(td.total_seconds() * 1_000_000_000)
    return delta


# NB: relies on dateutil internals, subject to change.
cdef bint infer_dateutil_fold(
    int64_t value,
    const int64_t[::1] trans,
    const int64_t[::1] deltas,
    intp_t pos,
):
    """
    Infer _TSObject fold property from value by assuming 0 and then setting
    to 1 if necessary.

    Parameters
    ----------
    value : int64_t
    trans : ndarray[int64_t]
        ndarray of offset transition points in nanoseconds since epoch.
    deltas : int64_t[:]
        array of offsets corresponding to transition points in trans.
    pos : intp_t
        Position of the last transition point before taking fold into account.

    Returns
    -------
    bint
        Due to daylight saving time, one wall clock time can occur twice
        when shifting from summer to winter time; fold describes whether the
        datetime-like corresponds  to the first (0) or the second time (1)
        the wall clock hits the ambiguous time

    References
    ----------
    .. [1] "PEP 495 - Local Time Disambiguation"
           https://www.python.org/dev/peps/pep-0495/#the-fold-attribute
    """
    cdef:
        bint fold = 0

    if pos > 0:
        fold_delta = deltas[pos - 1] - deltas[pos]
        if value - fold_delta < trans[pos]:
            fold = 1

    return fold
