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


def get_local_ampm():
    """Return the am and pm strings used in the current locale.

    We will use these when formatting datetime as string using string templates, which
    is faster than strftime when executed on arrays. The strftime directive where we
    need these is `%p`.

    Note that the considered locale is the one active when this function is called
    (not at cython compilation time).
    """
    return time(1).strftime("%p"), time(13).strftime("%p")
