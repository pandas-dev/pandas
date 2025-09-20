"""
Provides the default implementation of :class:`ArrowFactory <arrow.factory.ArrowFactory>`
methods for use as a module API.

"""

from datetime import date, datetime
from datetime import tzinfo as dt_tzinfo
from time import struct_time
from typing import Any, List, Optional, Tuple, Type, Union, overload

from arrow.arrow import TZ_EXPR, Arrow
from arrow.constants import DEFAULT_LOCALE
from arrow.factory import ArrowFactory

# internal default factory.
_factory = ArrowFactory()

# TODO: Use Positional Only Argument (https://www.python.org/dev/peps/pep-0570/)
#  after Python 3.7 deprecation


@overload
def get(
    *,
    locale: str = DEFAULT_LOCALE,
    tzinfo: Optional[TZ_EXPR] = None,
    normalize_whitespace: bool = False,
) -> Arrow:
    ...  # pragma: no cover


@overload
def get(
    *args: int,
    locale: str = DEFAULT_LOCALE,
    tzinfo: Optional[TZ_EXPR] = None,
    normalize_whitespace: bool = False,
) -> Arrow:
    ...  # pragma: no cover


@overload
def get(
    __obj: Union[
        Arrow,
        datetime,
        date,
        struct_time,
        dt_tzinfo,
        int,
        float,
        str,
        Tuple[int, int, int],
    ],
    *,
    locale: str = DEFAULT_LOCALE,
    tzinfo: Optional[TZ_EXPR] = None,
    normalize_whitespace: bool = False,
) -> Arrow:
    ...  # pragma: no cover


@overload
def get(
    __arg1: Union[datetime, date],
    __arg2: TZ_EXPR,
    *,
    locale: str = DEFAULT_LOCALE,
    tzinfo: Optional[TZ_EXPR] = None,
    normalize_whitespace: bool = False,
) -> Arrow:
    ...  # pragma: no cover


@overload
def get(
    __arg1: str,
    __arg2: Union[str, List[str]],
    *,
    locale: str = DEFAULT_LOCALE,
    tzinfo: Optional[TZ_EXPR] = None,
    normalize_whitespace: bool = False,
) -> Arrow:
    ...  # pragma: no cover


def get(*args: Any, **kwargs: Any) -> Arrow:
    """Calls the default :class:`ArrowFactory <arrow.factory.ArrowFactory>` ``get`` method."""

    return _factory.get(*args, **kwargs)


get.__doc__ = _factory.get.__doc__


def utcnow() -> Arrow:
    """Calls the default :class:`ArrowFactory <arrow.factory.ArrowFactory>` ``utcnow`` method."""

    return _factory.utcnow()


utcnow.__doc__ = _factory.utcnow.__doc__


def now(tz: Optional[TZ_EXPR] = None) -> Arrow:
    """Calls the default :class:`ArrowFactory <arrow.factory.ArrowFactory>` ``now`` method."""

    return _factory.now(tz)


now.__doc__ = _factory.now.__doc__


def factory(type: Type[Arrow]) -> ArrowFactory:
    """Returns an :class:`.ArrowFactory` for the specified :class:`Arrow <arrow.arrow.Arrow>`
    or derived type.

    :param type: the type, :class:`Arrow <arrow.arrow.Arrow>` or derived.

    """

    return ArrowFactory(type)


__all__ = ["get", "utcnow", "now", "factory"]
