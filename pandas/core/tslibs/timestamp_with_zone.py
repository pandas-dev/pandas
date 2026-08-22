"""TimestampWithZone - a Timestamp subclass that enforces timezone-aware datetimes.

This module provides TimestampWithZone, a subclass of pandas.Timestamp that
validates that the resulting Timestamp is timezone-aware. It keeps construction
and timezone parsing/validation delegated to pandas.Timestamp and only enforces
that the produced Timestamp has a tz. This avoids duplicating timezone parsing
logic and keeps the class lightweight.

Examples
--------
>>> from pandas.core.tslibs.timestamp_with_zone import TimestampWithZone
>>> TimestampWithZone("2020-01-01T00:00:00+00:00")
TimestampWithZone('2020-01-01 00:00:00+00:00', tz='UTC')

You can also pass a tz keyword to localize naive inputs:
>>> TimestampWithZone("2020-01-01T00:00:00", tz="UTC")
TimestampWithZone('2020-01-01 00:00:00+00:00', tz='UTC')

Naive datetimes without an explicit tz will raise:
>>> TimestampWithZone("2020-01-01T00:00:00")
Traceback (most recent call last):
    ...
ValueError: TimestampWithZone requires a timezone-aware Timestamp

Notes
-----
- This class intentionally does minimal work: it delegates parsing and
  localization to pandas.Timestamp and only tests for presence of `tz`.
- Because pandas.Timestamp is an extension/Cython-backed type, some
  operations on TimestampWithZone may return base Timestamp instances. If preserving
  the subclass across all operations is required, consider implementing a
  wrapper type instead.

"""

from pandas import Timestamp

__all__ = ["TimestampWithZone"]


class TimestampWithZone(Timestamp):
    """A Timestamp subclass that enforces timezone-aware datetimes.

    This class ensures that any Timestamp created is always timezone-aware,
    raising a ValueError if a naive datetime is provided without an explicit tz.

    Parameters
    ----------
    *args, **kwargs
        Same as pandas.Timestamp. A `tz` keyword may be provided to localize
        naive inputs. After construction, the resulting Timestamp must have a
        timezone; otherwise a ValueError is raised.

    Raises
    ------
    ValueError
        If the resulting Timestamp does not have a timezone (tz is None).

    Examples
    --------
    Create from an ISO string with timezone:

    >>> TimestampWithZone("2020-01-01T00:00:00+00:00")  # doctest: +SKIP
    TimestampWithZone('2020-01-01 00:00:00+00:00', tz='UTC')

    Localize a naive datetime by providing tz keyword:

    >>> TimestampWithZone("2020-01-01T00:00:00", tz="UTC")  # doctest: +SKIP
    TimestampWithZone('2020-01-01 00:00:00+00:00', tz='UTC')

    Creating from a naive datetime string without tz keyword raises:

    >>> TimestampWithZone("2020-01-01T00:00:00")  # doctest: +SKIP
    Traceback (most recent call last):
        ...
    ValueError: TimestampWithZone requires a timezone-aware Timestamp
    """

    def __new__(cls, *args, tz=None, **kwargs):
        """Create a new TimestampWithZone instance.

        Parameters
        ----------
        *args
            Positional arguments passed to pandas.Timestamp.
        tz : str, pytz.timezone, dateutil.tz.tzfile, or None, optional
            Time zone for the Timestamp. If provided, will localize naive inputs.
        **kwargs
            Keyword arguments passed to pandas.Timestamp.

        Returns
        -------
        TimestampWithZone
            A timezone-aware Timestamp instance.

        Raises
        ------
        ValueError
            If the resulting Timestamp is not timezone-aware.
        """
        # Delegate construction and timezone parsing/localization to pandas.Timestamp
        # Timestamp accepts a `tz` kw and will localize/convert as appropriate.
        ts = super().__new__(cls, *args, tz=tz, **kwargs)

        # Ensure the resulting Timestamp is timezone-aware
        if ts.tz is None:
            raise ValueError("TimestampWithZone requires a timezone-aware Timestamp")

        return ts

    def __repr__(self) -> str:
        """Return string representation of the Timestamp."""
        return (
            f"{self.__class__.__name__}('{self.isoformat()}', tz='{self.tz}')"
        )
