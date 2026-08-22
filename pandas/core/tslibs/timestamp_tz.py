"""TimestampTZ - a thin Timestamp subclass that requires timezone-aware datetimes.

This module provides TimestampTZ, a small subclass of pandas.Timestamp that
validates that the resulting Timestamp is timezone-aware. It keeps construction
and timezone parsing/validation delegated to pandas.Timestamp and only enforces
that the produced Timestamp has a tz. This avoids duplicating timezone parsing
logic and keeps the class lightweight.

Examples
--------
>>> from pandas.core.tslibs.timestamp_tz import TimestampTZ
>>> TimestampTZ("2020-01-01T00:00:00+00:00")
TimestampTZ('2020-01-01 00:00:00+00:00', tz='UTC')

You can also pass a tz keyword to localize naive inputs:
>>> TimestampTZ("2020-01-01T00:00:00", tz="UTC")
TimestampTZ('2020-01-01 00:00:00+00:00', tz='UTC')

Naive datetimes without an explicit tz will raise:
>>> TimestampTZ("2020-01-01T00:00:00")
Traceback (most recent call last):
    ...
ValueError: TimestampTZ requires a timezone-aware Timestamp

Notes
-----
- This class intentionally does minimal work: it delegates parsing and
  localization to pandas.Timestamp and only tests for presence of `tz`.
- Because pandas.Timestamp is an extension/Cython-backed type, some
  operations on TimestampTZ may return base Timestamp instances. If preserving
  the subclass across all operations is required, consider implementing a
  wrapper type instead.

"""

from pandas import Timestamp

__all__ = ["TimestampTZ"]


class TimestampTZ(Timestamp):
    """A Timestamp subclass that requires timezone-aware datetimes.

    Parameters
    ----------
    *args, **kwargs
        Same as pandas.Timestamp. A `tz` keyword may be provided to localize
        naive inputs. After construction, the resulting Timestamp must have a
        timezone; otherwise a ValueError is raised.

    Examples
    --------
    >>> TimestampTZ("2020-01-01T00:00:00+00:00")
    TimestampTZ('2020-01-01 00:00:00+00:00', tz='UTC')

    >>> TimestampTZ("2020-01-01T00:00:00", tz="UTC")
    TimestampTZ('2020-01-01 00:00:00+00:00', tz='UTC')

    >>> TimestampTZ("2020-01-01T00:00:00")
    Traceback (most recent call last):
        ...
    ValueError: TimestampTZ requires a timezone-aware Timestamp
    """

    def __new__(cls, *args, tz=None, **kwargs):
        # Delegate construction and timezone parsing/localization to pandas.Timestamp
        # Timestamp accepts a `tz` kw and will localize/convert as appropriate.
        ts = super().__new__(cls, *args, tz=tz, **kwargs)

        # Ensure the resulting Timestamp is timezone-aware
        if ts.tz is None:
            raise ValueError("TimestampTZ requires a timezone-aware Timestamp")

        return ts
