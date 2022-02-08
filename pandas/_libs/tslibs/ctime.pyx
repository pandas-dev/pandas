"""
Cython implementation of (parts of) the standard library time module.
"""

from cpython.exc cimport PyErr_SetFromErrno
from libc.stdint cimport int64_t


cdef extern from "Python.h":
    ctypedef int64_t _PyTime_t
    _PyTime_t _PyTime_GetSystemClock() nogil
    double _PyTime_AsSecondsDouble(_PyTime_t t) nogil

from libc.time cimport (
    localtime as libc_localtime,
    time_t,
    tm,
)


def pytime():
    """
    python-exposed for testing
    """
    return time()


def pylocaltime():
    """
    python-exposed for testing
    """
    lt = localtime()
    # https://github.com/pandas-dev/pandas/pull/45864#issuecomment-1033021599
    return {
        "tm_year": lt.tm_year,
        "tm_mon": lt.tm_mon,
        "tm_mday": lt.tm_mday,
        "tm_hour": lt.tm_hour,
        "tm_min": lt.tm_min,
        "tm_sec": lt.tm_sec,
        "tm_wday": lt.tm_wday,
        "tm_yday": lt.tm_yday,
        "tm_isdst": lt.tm_isdst,
    }


cdef inline double time() nogil:
    cdef:
        _PyTime_t tic

    tic = _PyTime_GetSystemClock()
    return _PyTime_AsSecondsDouble(tic)


cdef inline int _raise_from_errno() except -1 with gil:
    PyErr_SetFromErrno(RuntimeError)
    return <int>-1  # Let the C compiler know that this function always raises.


cdef inline tm localtime() nogil except *:
    """
    Analogue to the stdlib time.localtime.  The returned struct
    has some entries that the stdlib version does not: tm_gmtoff, tm_zone
    """
    cdef:
        time_t tic = <time_t>time()
        tm* result

    result = libc_localtime(&tic)
    if result is NULL:
        _raise_from_errno()
    # Fix 0-based date values (and the 1900-based year).
    # See tmtotuple() in https://github.com/python/cpython/blob/master/Modules/timemodule.c
    result.tm_year += 1900
    result.tm_mon += 1
    result.tm_wday = (result.tm_wday + 6) % 7
    result.tm_yday += 1
    return result[0]
