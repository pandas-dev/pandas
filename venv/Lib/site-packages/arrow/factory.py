"""
Implements the :class:`ArrowFactory <arrow.factory.ArrowFactory>` class,
providing factory methods for common :class:`Arrow <arrow.arrow.Arrow>`
construction scenarios.

"""

import calendar
from datetime import date, datetime, timezone
from datetime import tzinfo as dt_tzinfo
from decimal import Decimal
from time import struct_time
from typing import Any, List, Optional, Tuple, Type, Union, overload

from arrow import parser
from arrow.arrow import TZ_EXPR, Arrow
from arrow.constants import DEFAULT_LOCALE
from arrow.util import is_timestamp, iso_to_gregorian


class ArrowFactory:
    """A factory for generating :class:`Arrow <arrow.arrow.Arrow>` objects.

    :param type: (optional) the :class:`Arrow <arrow.arrow.Arrow>`-based class to construct from.
        Defaults to :class:`Arrow <arrow.arrow.Arrow>`.

    """

    type: Type[Arrow]

    def __init__(self, type: Type[Arrow] = Arrow) -> None:
        self.type = type

    @overload
    def get(
        self,
        *,
        locale: str = DEFAULT_LOCALE,
        tzinfo: Optional[TZ_EXPR] = None,
        normalize_whitespace: bool = False,
    ) -> Arrow: ...  # pragma: no cover

    @overload
    def get(
        self,
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
    ) -> Arrow: ...  # pragma: no cover

    @overload
    def get(
        self,
        __arg1: Union[datetime, date],
        __arg2: TZ_EXPR,
        *,
        locale: str = DEFAULT_LOCALE,
        tzinfo: Optional[TZ_EXPR] = None,
        normalize_whitespace: bool = False,
    ) -> Arrow: ...  # pragma: no cover

    @overload
    def get(
        self,
        __arg1: str,
        __arg2: Union[str, List[str]],
        *,
        locale: str = DEFAULT_LOCALE,
        tzinfo: Optional[TZ_EXPR] = None,
        normalize_whitespace: bool = False,
    ) -> Arrow: ...  # pragma: no cover

    def get(self, *args: Any, **kwargs: Any) -> Arrow:
        """Returns an :class:`Arrow <arrow.arrow.Arrow>` object based on flexible inputs.

        :param locale: (optional) a ``str`` specifying a locale for the parser. Defaults to 'en-us'.
        :param tzinfo: (optional) a :ref:`timezone expression <tz-expr>` or tzinfo object.
            Replaces the timezone unless using an input form that is explicitly UTC or specifies
            the timezone in a positional argument. Defaults to UTC.
        :param normalize_whitespace: (optional) a ``bool`` specifying whether or not to normalize
            redundant whitespace (spaces, tabs, and newlines) in a datetime string before parsing.
            Defaults to false.

        Usage::

            >>> import arrow

        **No inputs** to get current UTC time::

            >>> arrow.get()
            <Arrow [2013-05-08T05:51:43.316458+00:00]>

        **One** :class:`Arrow <arrow.arrow.Arrow>` object, to get a copy.

            >>> arw = arrow.utcnow()
            >>> arrow.get(arw)
            <Arrow [2013-10-23T15:21:54.354846+00:00]>

        **One** ``float`` or ``int``, convertible to a floating-point timestamp, to get
        that timestamp in UTC::

            >>> arrow.get(1367992474.293378)
            <Arrow [2013-05-08T05:54:34.293378+00:00]>

            >>> arrow.get(1367992474)
            <Arrow [2013-05-08T05:54:34+00:00]>

        **One** ISO 8601-formatted ``str``, to parse it::

            >>> arrow.get('2013-09-29T01:26:43.830580')
            <Arrow [2013-09-29T01:26:43.830580+00:00]>

        **One** ISO 8601-formatted ``str``, in basic format, to parse it::

            >>> arrow.get('20160413T133656.456289')
            <Arrow [2016-04-13T13:36:56.456289+00:00]>

        **One** ``tzinfo``, to get the current time **converted** to that timezone::

            >>> arrow.get(tz.tzlocal())
            <Arrow [2013-05-07T22:57:28.484717-07:00]>

        **One** naive ``datetime``, to get that datetime in UTC::

            >>> arrow.get(datetime(2013, 5, 5))
            <Arrow [2013-05-05T00:00:00+00:00]>

        **One** aware ``datetime``, to get that datetime::

            >>> arrow.get(datetime(2013, 5, 5, tzinfo=tz.tzlocal()))
            <Arrow [2013-05-05T00:00:00-07:00]>

        **One** naive ``date``, to get that date in UTC::

            >>> arrow.get(date(2013, 5, 5))
            <Arrow [2013-05-05T00:00:00+00:00]>

        **One** time.struct time::

            >>> arrow.get(gmtime(0))
            <Arrow [1970-01-01T00:00:00+00:00]>

        **One** iso calendar ``tuple``, to get that week date in UTC::

            >>> arrow.get((2013, 18, 7))
            <Arrow [2013-05-05T00:00:00+00:00]>

        **Two** arguments, a naive or aware ``datetime``, and a replacement
        :ref:`timezone expression <tz-expr>`::

            >>> arrow.get(datetime(2013, 5, 5), 'US/Pacific')
            <Arrow [2013-05-05T00:00:00-07:00]>

        **Two** arguments, a naive ``date``, and a replacement
        :ref:`timezone expression <tz-expr>`::

            >>> arrow.get(date(2013, 5, 5), 'US/Pacific')
            <Arrow [2013-05-05T00:00:00-07:00]>

        **Two** arguments, both ``str``, to parse the first according to the format of the second::

            >>> arrow.get('2013-05-05 12:30:45 America/Chicago', 'YYYY-MM-DD HH:mm:ss ZZZ')
            <Arrow [2013-05-05T12:30:45-05:00]>

        **Two** arguments, first a ``str`` to parse and second a ``list`` of formats to try::

            >>> arrow.get('2013-05-05 12:30:45', ['MM/DD/YYYY', 'YYYY-MM-DD HH:mm:ss'])
            <Arrow [2013-05-05T12:30:45+00:00]>

        **Three or more** arguments, as for the direct constructor of an ``Arrow`` object::

            >>> arrow.get(2013, 5, 5, 12, 30, 45)
            <Arrow [2013-05-05T12:30:45+00:00]>

        """

        arg_count = len(args)
        locale = kwargs.pop("locale", DEFAULT_LOCALE)
        tz = kwargs.get("tzinfo", None)
        normalize_whitespace = kwargs.pop("normalize_whitespace", False)

        # if kwargs given, send to constructor unless only tzinfo provided
        if len(kwargs) > 1:
            arg_count = 3

        # tzinfo kwarg is not provided
        if len(kwargs) == 1 and tz is None:
            arg_count = 3

        # () -> now, @ tzinfo or utc
        if arg_count == 0:
            if isinstance(tz, str):
                tz = parser.TzinfoParser.parse(tz)
                return self.type.now(tzinfo=tz)

            if isinstance(tz, dt_tzinfo):
                return self.type.now(tzinfo=tz)

            return self.type.utcnow()

        if arg_count == 1:
            arg = args[0]
            if isinstance(arg, Decimal):
                arg = float(arg)

            # (None) -> raises an exception
            if arg is None:
                raise TypeError("Cannot parse argument of type None.")

            # try (int, float) -> from timestamp @ tzinfo
            elif not isinstance(arg, str) and is_timestamp(arg):
                if tz is None:
                    # set to UTC by default
                    tz = timezone.utc
                return self.type.fromtimestamp(arg, tzinfo=tz)

            # (Arrow) -> from the object's datetime @ tzinfo
            elif isinstance(arg, Arrow):
                return self.type.fromdatetime(arg.datetime, tzinfo=tz)

            # (datetime) -> from datetime @ tzinfo
            elif isinstance(arg, datetime):
                return self.type.fromdatetime(arg, tzinfo=tz)

            # (date) -> from date @ tzinfo
            elif isinstance(arg, date):
                return self.type.fromdate(arg, tzinfo=tz)

            # (tzinfo) -> now @ tzinfo
            elif isinstance(arg, dt_tzinfo):
                return self.type.now(tzinfo=arg)

            # (str) -> parse @ tzinfo
            elif isinstance(arg, str):
                dt = parser.DateTimeParser(locale).parse_iso(arg, normalize_whitespace)
                return self.type.fromdatetime(dt, tzinfo=tz)

            # (struct_time) -> from struct_time
            elif isinstance(arg, struct_time):
                return self.type.utcfromtimestamp(calendar.timegm(arg))

            # (iso calendar) -> convert then from date @ tzinfo
            elif isinstance(arg, tuple) and len(arg) == 3:
                d = iso_to_gregorian(*arg)
                return self.type.fromdate(d, tzinfo=tz)

            else:
                raise TypeError(f"Cannot parse single argument of type {type(arg)!r}.")

        elif arg_count == 2:
            arg_1, arg_2 = args[0], args[1]

            if isinstance(arg_1, datetime):
                # (datetime, tzinfo/str) -> fromdatetime @ tzinfo
                if isinstance(arg_2, (dt_tzinfo, str)):
                    return self.type.fromdatetime(arg_1, tzinfo=arg_2)
                else:
                    raise TypeError(
                        f"Cannot parse two arguments of types 'datetime', {type(arg_2)!r}."
                    )

            elif isinstance(arg_1, date):
                # (date, tzinfo/str) -> fromdate @ tzinfo
                if isinstance(arg_2, (dt_tzinfo, str)):
                    return self.type.fromdate(arg_1, tzinfo=arg_2)
                else:
                    raise TypeError(
                        f"Cannot parse two arguments of types 'date', {type(arg_2)!r}."
                    )

            # (str, format) -> parse @ tzinfo
            elif isinstance(arg_1, str) and isinstance(arg_2, (str, list)):
                dt = parser.DateTimeParser(locale).parse(
                    args[0], args[1], normalize_whitespace
                )
                return self.type.fromdatetime(dt, tzinfo=tz)

            else:
                raise TypeError(
                    f"Cannot parse two arguments of types {type(arg_1)!r} and {type(arg_2)!r}."
                )

        # 3+ args -> datetime-like via constructor
        else:
            return self.type(*args, **kwargs)

    def utcnow(self) -> Arrow:
        """Returns an :class:`Arrow <arrow.arrow.Arrow>` object, representing "now" in UTC time.

        Usage::

            >>> import arrow
            >>> arrow.utcnow()
            <Arrow [2013-05-08T05:19:07.018993+00:00]>
        """

        return self.type.utcnow()

    def now(self, tz: Optional[TZ_EXPR] = None) -> Arrow:
        """Returns an :class:`Arrow <arrow.arrow.Arrow>` object, representing "now" in the given
        timezone.

        :param tz: (optional) A :ref:`timezone expression <tz-expr>`.  Defaults to local time.

        Usage::

            >>> import arrow
            >>> arrow.now()
            <Arrow [2013-05-07T22:19:11.363410-07:00]>

            >>> arrow.now('US/Pacific')
            <Arrow [2013-05-07T22:19:15.251821-07:00]>

            >>> arrow.now('+02:00')
            <Arrow [2013-05-08T07:19:25.618646+02:00]>

            >>> arrow.now('local')
            <Arrow [2013-05-07T22:19:39.130059-07:00]>
        """

        if tz is None:
            tz = datetime.now().astimezone().tzinfo
        elif not isinstance(tz, dt_tzinfo):
            tz = parser.TzinfoParser.parse(tz)

        return self.type.now(tz)
