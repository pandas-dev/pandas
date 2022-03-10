"""
We define base classes that will be inherited by Timestamp, Timedelta, etc
in order to allow for fast isinstance checks without circular dependency issues.

This is analogous to core.dtypes.generic.
"""

from cpython.datetime cimport (
    datetime,
    time,
)


cdef class ABCTimestamp(datetime):
    pass


cdef:
    # Use strftime to get the locale-specific version of am and pm
    # we will use these when formatting datetime as string
    bytes AM_LOCAL = time(1).strftime("%p").encode("utf-8")
    bytes PM_LOCAL = time(13).strftime("%p").encode("utf-8")


def get_local_ampm():
    """Return the utf-8 encoded am and pm strings used in the current locale

    Note that the considered locale is the one active at module import time.
    """
    global AM_LOCAL, PM_LOCAL
    return AM_LOCAL, PM_LOCAL
