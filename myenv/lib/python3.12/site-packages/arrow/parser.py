"""Provides the :class:`Arrow <arrow.parser.DateTimeParser>` class, a better way to parse datetime strings."""

import re
import sys
from datetime import datetime, timedelta
from datetime import tzinfo as dt_tzinfo
from functools import lru_cache
from typing import (
    Any,
    ClassVar,
    Dict,
    Iterable,
    List,
    Match,
    Optional,
    Pattern,
    SupportsFloat,
    SupportsInt,
    Tuple,
    Union,
    cast,
    overload,
)

from dateutil import tz

from arrow import locales
from arrow.constants import DEFAULT_LOCALE
from arrow.util import next_weekday, normalize_timestamp

if sys.version_info < (3, 8):  # pragma: no cover
    from typing_extensions import Literal, TypedDict
else:
    from typing import Literal, TypedDict  # pragma: no cover


class ParserError(ValueError):
    pass


# Allows for ParserErrors to be propagated from _build_datetime()
# when day_of_year errors occur.
# Before this, the ParserErrors were caught by the try/except in
# _parse_multiformat() and the appropriate error message was not
# transmitted to the user.
class ParserMatchError(ParserError):
    pass


_WEEKDATE_ELEMENT = Union[str, bytes, SupportsInt, bytearray]

_FORMAT_TYPE = Literal[
    "YYYY",
    "YY",
    "MM",
    "M",
    "DDDD",
    "DDD",
    "DD",
    "D",
    "HH",
    "H",
    "hh",
    "h",
    "mm",
    "m",
    "ss",
    "s",
    "X",
    "x",
    "ZZZ",
    "ZZ",
    "Z",
    "S",
    "W",
    "MMMM",
    "MMM",
    "Do",
    "dddd",
    "ddd",
    "d",
    "a",
    "A",
]


class _Parts(TypedDict, total=False):
    year: int
    month: int
    day_of_year: int
    day: int
    hour: int
    minute: int
    second: int
    microsecond: int
    timestamp: float
    expanded_timestamp: int
    tzinfo: dt_tzinfo
    am_pm: Literal["am", "pm"]
    day_of_week: int
    weekdate: Tuple[_WEEKDATE_ELEMENT, _WEEKDATE_ELEMENT, Optional[_WEEKDATE_ELEMENT]]


class DateTimeParser:
    _FORMAT_RE: ClassVar[Pattern[str]] = re.compile(
        r"(YYY?Y?|MM?M?M?|Do|DD?D?D?|d?d?d?d|HH?|hh?|mm?|ss?|S+|ZZ?Z?|a|A|x|X|W)"
    )
    _ESCAPE_RE: ClassVar[Pattern[str]] = re.compile(r"\[[^\[\]]*\]")

    _ONE_OR_TWO_DIGIT_RE: ClassVar[Pattern[str]] = re.compile(r"\d{1,2}")
    _ONE_OR_TWO_OR_THREE_DIGIT_RE: ClassVar[Pattern[str]] = re.compile(r"\d{1,3}")
    _ONE_OR_MORE_DIGIT_RE: ClassVar[Pattern[str]] = re.compile(r"\d+")
    _TWO_DIGIT_RE: ClassVar[Pattern[str]] = re.compile(r"\d{2}")
    _THREE_DIGIT_RE: ClassVar[Pattern[str]] = re.compile(r"\d{3}")
    _FOUR_DIGIT_RE: ClassVar[Pattern[str]] = re.compile(r"\d{4}")
    _TZ_Z_RE: ClassVar[Pattern[str]] = re.compile(r"([\+\-])(\d{2})(?:(\d{2}))?|Z")
    _TZ_ZZ_RE: ClassVar[Pattern[str]] = re.compile(r"([\+\-])(\d{2})(?:\:(\d{2}))?|Z")
    _TZ_NAME_RE: ClassVar[Pattern[str]] = re.compile(r"\w[\w+\-/]+")
    # NOTE: timestamps cannot be parsed from natural language strings (by removing the ^...$) because it will
    # break cases like "15 Jul 2000" and a format list (see issue #447)
    _TIMESTAMP_RE: ClassVar[Pattern[str]] = re.compile(r"^\-?\d+\.?\d+$")
    _TIMESTAMP_EXPANDED_RE: ClassVar[Pattern[str]] = re.compile(r"^\-?\d+$")
    _TIME_RE: ClassVar[Pattern[str]] = re.compile(
        r"^(\d{2})(?:\:?(\d{2}))?(?:\:?(\d{2}))?(?:([\.\,])(\d+))?$"
    )
    _WEEK_DATE_RE: ClassVar[Pattern[str]] = re.compile(
        r"(?P<year>\d{4})[\-]?W(?P<week>\d{2})[\-]?(?P<day>\d)?"
    )

    _BASE_INPUT_RE_MAP: ClassVar[Dict[_FORMAT_TYPE, Pattern[str]]] = {
        "YYYY": _FOUR_DIGIT_RE,
        "YY": _TWO_DIGIT_RE,
        "MM": _TWO_DIGIT_RE,
        "M": _ONE_OR_TWO_DIGIT_RE,
        "DDDD": _THREE_DIGIT_RE,
        "DDD": _ONE_OR_TWO_OR_THREE_DIGIT_RE,
        "DD": _TWO_DIGIT_RE,
        "D": _ONE_OR_TWO_DIGIT_RE,
        "HH": _TWO_DIGIT_RE,
        "H": _ONE_OR_TWO_DIGIT_RE,
        "hh": _TWO_DIGIT_RE,
        "h": _ONE_OR_TWO_DIGIT_RE,
        "mm": _TWO_DIGIT_RE,
        "m": _ONE_OR_TWO_DIGIT_RE,
        "ss": _TWO_DIGIT_RE,
        "s": _ONE_OR_TWO_DIGIT_RE,
        "X": _TIMESTAMP_RE,
        "x": _TIMESTAMP_EXPANDED_RE,
        "ZZZ": _TZ_NAME_RE,
        "ZZ": _TZ_ZZ_RE,
        "Z": _TZ_Z_RE,
        "S": _ONE_OR_MORE_DIGIT_RE,
        "W": _WEEK_DATE_RE,
    }

    SEPARATORS: ClassVar[List[str]] = ["-", "/", "."]

    locale: locales.Locale
    _input_re_map: Dict[_FORMAT_TYPE, Pattern[str]]

    def __init__(self, locale: str = DEFAULT_LOCALE, cache_size: int = 0) -> None:
        self.locale = locales.get_locale(locale)
        self._input_re_map = self._BASE_INPUT_RE_MAP.copy()
        self._input_re_map.update(
            {
                "MMMM": self._generate_choice_re(
                    self.locale.month_names[1:], re.IGNORECASE
                ),
                "MMM": self._generate_choice_re(
                    self.locale.month_abbreviations[1:], re.IGNORECASE
                ),
                "Do": re.compile(self.locale.ordinal_day_re),
                "dddd": self._generate_choice_re(
                    self.locale.day_names[1:], re.IGNORECASE
                ),
                "ddd": self._generate_choice_re(
                    self.locale.day_abbreviations[1:], re.IGNORECASE
                ),
                "d": re.compile(r"[1-7]"),
                "a": self._generate_choice_re(
                    (self.locale.meridians["am"], self.locale.meridians["pm"])
                ),
                # note: 'A' token accepts both 'am/pm' and 'AM/PM' formats to
                # ensure backwards compatibility of this token
                "A": self._generate_choice_re(self.locale.meridians.values()),
            }
        )
        if cache_size > 0:
            self._generate_pattern_re = lru_cache(maxsize=cache_size)(  # type: ignore
                self._generate_pattern_re
            )

    # TODO: since we support more than ISO 8601, we should rename this function
    # IDEA: break into multiple functions
    def parse_iso(
        self, datetime_string: str, normalize_whitespace: bool = False
    ) -> datetime:
        if normalize_whitespace:
            datetime_string = re.sub(r"\s+", " ", datetime_string.strip())

        has_space_divider = " " in datetime_string
        has_t_divider = "T" in datetime_string

        num_spaces = datetime_string.count(" ")
        if has_space_divider and num_spaces != 1 or has_t_divider and num_spaces > 0:
            raise ParserError(
                f"Expected an ISO 8601-like string, but was given {datetime_string!r}. "
                "Try passing in a format string to resolve this."
            )

        has_time = has_space_divider or has_t_divider
        has_tz = False

        # date formats (ISO 8601 and others) to test against
        # NOTE: YYYYMM is omitted to avoid confusion with YYMMDD (no longer part of ISO 8601, but is still often used)
        formats = [
            "YYYY-MM-DD",
            "YYYY-M-DD",
            "YYYY-M-D",
            "YYYY/MM/DD",
            "YYYY/M/DD",
            "YYYY/M/D",
            "YYYY.MM.DD",
            "YYYY.M.DD",
            "YYYY.M.D",
            "YYYYMMDD",
            "YYYY-DDDD",
            "YYYYDDDD",
            "YYYY-MM",
            "YYYY/MM",
            "YYYY.MM",
            "YYYY",
            "W",
        ]

        if has_time:
            if has_space_divider:
                date_string, time_string = datetime_string.split(" ", 1)
            else:
                date_string, time_string = datetime_string.split("T", 1)

            time_parts = re.split(
                r"[\+\-Z]", time_string, maxsplit=1, flags=re.IGNORECASE
            )

            time_components: Optional[Match[str]] = self._TIME_RE.match(time_parts[0])

            if time_components is None:
                raise ParserError(
                    "Invalid time component provided. "
                    "Please specify a format or provide a valid time component in the basic or extended ISO 8601 time format."
                )

            (
                hours,
                minutes,
                seconds,
                subseconds_sep,
                subseconds,
            ) = time_components.groups()

            has_tz = len(time_parts) == 2
            has_minutes = minutes is not None
            has_seconds = seconds is not None
            has_subseconds = subseconds is not None

            is_basic_time_format = ":" not in time_parts[0]
            tz_format = "Z"

            # use 'ZZ' token instead since tz offset is present in non-basic format
            if has_tz and ":" in time_parts[1]:
                tz_format = "ZZ"

            time_sep = "" if is_basic_time_format else ":"

            if has_subseconds:
                time_string = "HH{time_sep}mm{time_sep}ss{subseconds_sep}S".format(
                    time_sep=time_sep, subseconds_sep=subseconds_sep
                )
            elif has_seconds:
                time_string = "HH{time_sep}mm{time_sep}ss".format(time_sep=time_sep)
            elif has_minutes:
                time_string = f"HH{time_sep}mm"
            else:
                time_string = "HH"

            if has_space_divider:
                formats = [f"{f} {time_string}" for f in formats]
            else:
                formats = [f"{f}T{time_string}" for f in formats]

        if has_time and has_tz:
            # Add "Z" or "ZZ" to the format strings to indicate to
            # _parse_token() that a timezone needs to be parsed
            formats = [f"{f}{tz_format}" for f in formats]

        return self._parse_multiformat(datetime_string, formats)

    def parse(
        self,
        datetime_string: str,
        fmt: Union[List[str], str],
        normalize_whitespace: bool = False,
    ) -> datetime:
        if normalize_whitespace:
            datetime_string = re.sub(r"\s+", " ", datetime_string)

        if isinstance(fmt, list):
            return self._parse_multiformat(datetime_string, fmt)

        try:
            fmt_tokens: List[_FORMAT_TYPE]
            fmt_pattern_re: Pattern[str]
            fmt_tokens, fmt_pattern_re = self._generate_pattern_re(fmt)
        except re.error as e:
            raise ParserMatchError(
                f"Failed to generate regular expression pattern: {e}."
            )

        match = fmt_pattern_re.search(datetime_string)

        if match is None:
            raise ParserMatchError(
                f"Failed to match {fmt!r} when parsing {datetime_string!r}."
            )

        parts: _Parts = {}
        for token in fmt_tokens:
            value: Union[Tuple[str, str, str], str]
            if token == "Do":
                value = match.group("value")
            elif token == "W":
                value = (match.group("year"), match.group("week"), match.group("day"))
            else:
                value = match.group(token)

            if value is None:
                raise ParserMatchError(
                    f"Unable to find a match group for the specified token {token!r}."
                )

            self._parse_token(token, value, parts)  # type: ignore[arg-type]

        return self._build_datetime(parts)

    def _generate_pattern_re(self, fmt: str) -> Tuple[List[_FORMAT_TYPE], Pattern[str]]:
        # fmt is a string of tokens like 'YYYY-MM-DD'
        # we construct a new string by replacing each
        # token by its pattern:
        # 'YYYY-MM-DD' -> '(?P<YYYY>\d{4})-(?P<MM>\d{2})-(?P<DD>\d{2})'
        tokens: List[_FORMAT_TYPE] = []
        offset = 0

        # Escape all special RegEx chars
        escaped_fmt = re.escape(fmt)

        # Extract the bracketed expressions to be reinserted later.
        escaped_fmt = re.sub(self._ESCAPE_RE, "#", escaped_fmt)

        # Any number of S is the same as one.
        # TODO: allow users to specify the number of digits to parse
        escaped_fmt = re.sub(r"S+", "S", escaped_fmt)

        escaped_data = re.findall(self._ESCAPE_RE, fmt)

        fmt_pattern = escaped_fmt

        for m in self._FORMAT_RE.finditer(escaped_fmt):
            token: _FORMAT_TYPE = cast(_FORMAT_TYPE, m.group(0))
            try:
                input_re = self._input_re_map[token]
            except KeyError:
                raise ParserError(f"Unrecognized token {token!r}.")
            input_pattern = f"(?P<{token}>{input_re.pattern})"
            tokens.append(token)
            # a pattern doesn't have the same length as the token
            # it replaces! We keep the difference in the offset variable.
            # This works because the string is scanned left-to-right and matches
            # are returned in the order found by finditer.
            fmt_pattern = (
                fmt_pattern[: m.start() + offset]
                + input_pattern
                + fmt_pattern[m.end() + offset :]
            )
            offset += len(input_pattern) - (m.end() - m.start())

        final_fmt_pattern = ""
        split_fmt = fmt_pattern.split(r"\#")

        # Due to the way Python splits, 'split_fmt' will always be longer
        for i in range(len(split_fmt)):
            final_fmt_pattern += split_fmt[i]
            if i < len(escaped_data):
                final_fmt_pattern += escaped_data[i][1:-1]

        # Wrap final_fmt_pattern in a custom word boundary to strictly
        # match the formatting pattern and filter out date and time formats
        # that include junk such as: blah1998-09-12 blah, blah 1998-09-12blah,
        # blah1998-09-12blah. The custom word boundary matches every character
        # that is not a whitespace character to allow for searching for a date
        # and time string in a natural language sentence. Therefore, searching
        # for a string of the form YYYY-MM-DD in "blah 1998-09-12 blah" will
        # work properly.
        # Certain punctuation before or after the target pattern such as
        # "1998-09-12," is permitted. For the full list of valid punctuation,
        # see the documentation.

        starting_word_boundary = (
            r"(?<!\S\S)"  # Don't have two consecutive non-whitespace characters. This ensures that we allow cases
            # like .11.25.2019 but not 1.11.25.2019 (for pattern MM.DD.YYYY)
            r"(?<![^\,\.\;\:\?\!\"\'\`\[\]\{\}\(\)<>\s])"  # This is the list of punctuation that is ok before the
            # pattern (i.e. "It can't not be these characters before the pattern")
            r"(\b|^)"
            # The \b is to block cases like 1201912 but allow 201912 for pattern YYYYMM. The ^ was necessary to allow a
            # negative number through i.e. before epoch numbers
        )
        ending_word_boundary = (
            r"(?=[\,\.\;\:\?\!\"\'\`\[\]\{\}\(\)\<\>]?"  # Positive lookahead stating that these punctuation marks
            # can appear after the pattern at most 1 time
            r"(?!\S))"  # Don't allow any non-whitespace character after the punctuation
        )
        bounded_fmt_pattern = r"{}{}{}".format(
            starting_word_boundary, final_fmt_pattern, ending_word_boundary
        )

        return tokens, re.compile(bounded_fmt_pattern, flags=re.IGNORECASE)

    @overload
    def _parse_token(
        self,
        token: Literal[
            "YYYY",
            "YY",
            "MM",
            "M",
            "DDDD",
            "DDD",
            "DD",
            "D",
            "Do",
            "HH",
            "hh",
            "h",
            "H",
            "mm",
            "m",
            "ss",
            "s",
            "x",
        ],
        value: Union[str, bytes, SupportsInt, bytearray],
        parts: _Parts,
    ) -> None:
        ...  # pragma: no cover

    @overload
    def _parse_token(
        self,
        token: Literal["X"],
        value: Union[str, bytes, SupportsFloat, bytearray],
        parts: _Parts,
    ) -> None:
        ...  # pragma: no cover

    @overload
    def _parse_token(
        self,
        token: Literal["MMMM", "MMM", "dddd", "ddd", "S"],
        value: Union[str, bytes, bytearray],
        parts: _Parts,
    ) -> None:
        ...  # pragma: no cover

    @overload
    def _parse_token(
        self,
        token: Literal["a", "A", "ZZZ", "ZZ", "Z"],
        value: Union[str, bytes],
        parts: _Parts,
    ) -> None:
        ...  # pragma: no cover

    @overload
    def _parse_token(
        self,
        token: Literal["W"],
        value: Tuple[_WEEKDATE_ELEMENT, _WEEKDATE_ELEMENT, Optional[_WEEKDATE_ELEMENT]],
        parts: _Parts,
    ) -> None:
        ...  # pragma: no cover

    def _parse_token(
        self,
        token: Any,
        value: Any,
        parts: _Parts,
    ) -> None:
        if token == "YYYY":
            parts["year"] = int(value)

        elif token == "YY":
            value = int(value)
            parts["year"] = 1900 + value if value > 68 else 2000 + value

        elif token in ["MMMM", "MMM"]:
            # FIXME: month_number() is nullable
            parts["month"] = self.locale.month_number(value.lower())  # type: ignore[typeddict-item]

        elif token in ["MM", "M"]:
            parts["month"] = int(value)

        elif token in ["DDDD", "DDD"]:
            parts["day_of_year"] = int(value)

        elif token in ["DD", "D"]:
            parts["day"] = int(value)

        elif token == "Do":
            parts["day"] = int(value)

        elif token == "dddd":
            # locale day names are 1-indexed
            day_of_week = [x.lower() for x in self.locale.day_names].index(
                value.lower()
            )
            parts["day_of_week"] = day_of_week - 1

        elif token == "ddd":
            # locale day abbreviations are 1-indexed
            day_of_week = [x.lower() for x in self.locale.day_abbreviations].index(
                value.lower()
            )
            parts["day_of_week"] = day_of_week - 1

        elif token.upper() in ["HH", "H"]:
            parts["hour"] = int(value)

        elif token in ["mm", "m"]:
            parts["minute"] = int(value)

        elif token in ["ss", "s"]:
            parts["second"] = int(value)

        elif token == "S":
            # We have the *most significant* digits of an arbitrary-precision integer.
            # We want the six most significant digits as an integer, rounded.
            # IDEA: add nanosecond support somehow? Need datetime support for it first.
            value = value.ljust(7, "0")

            # floating-point (IEEE-754) defaults to half-to-even rounding
            seventh_digit = int(value[6])
            if seventh_digit == 5:
                rounding = int(value[5]) % 2
            elif seventh_digit > 5:
                rounding = 1
            else:
                rounding = 0

            parts["microsecond"] = int(value[:6]) + rounding

        elif token == "X":
            parts["timestamp"] = float(value)

        elif token == "x":
            parts["expanded_timestamp"] = int(value)

        elif token in ["ZZZ", "ZZ", "Z"]:
            parts["tzinfo"] = TzinfoParser.parse(value)

        elif token in ["a", "A"]:
            if value in (self.locale.meridians["am"], self.locale.meridians["AM"]):
                parts["am_pm"] = "am"
                if "hour" in parts and not 0 <= parts["hour"] <= 12:
                    raise ParserMatchError(
                        f"Hour token value must be between 0 and 12 inclusive for token {token!r}."
                    )
            elif value in (self.locale.meridians["pm"], self.locale.meridians["PM"]):
                parts["am_pm"] = "pm"
        elif token == "W":
            parts["weekdate"] = value

    @staticmethod
    def _build_datetime(parts: _Parts) -> datetime:
        weekdate = parts.get("weekdate")

        if weekdate is not None:
            year, week = int(weekdate[0]), int(weekdate[1])

            if weekdate[2] is not None:
                _day = int(weekdate[2])
            else:
                # day not given, default to 1
                _day = 1

            date_string = f"{year}-{week}-{_day}"

            #  tokens for ISO 8601 weekdates
            dt = datetime.strptime(date_string, "%G-%V-%u")

            parts["year"] = dt.year
            parts["month"] = dt.month
            parts["day"] = dt.day

        timestamp = parts.get("timestamp")

        if timestamp is not None:
            return datetime.fromtimestamp(timestamp, tz=tz.tzutc())

        expanded_timestamp = parts.get("expanded_timestamp")

        if expanded_timestamp is not None:
            return datetime.fromtimestamp(
                normalize_timestamp(expanded_timestamp),
                tz=tz.tzutc(),
            )

        day_of_year = parts.get("day_of_year")

        if day_of_year is not None:
            _year = parts.get("year")
            month = parts.get("month")
            if _year is None:
                raise ParserError(
                    "Year component is required with the DDD and DDDD tokens."
                )

            if month is not None:
                raise ParserError(
                    "Month component is not allowed with the DDD and DDDD tokens."
                )

            date_string = f"{_year}-{day_of_year}"
            try:
                dt = datetime.strptime(date_string, "%Y-%j")
            except ValueError:
                raise ParserError(
                    f"The provided day of year {day_of_year!r} is invalid."
                )

            parts["year"] = dt.year
            parts["month"] = dt.month
            parts["day"] = dt.day

        day_of_week: Optional[int] = parts.get("day_of_week")
        day = parts.get("day")

        # If day is passed, ignore day of week
        if day_of_week is not None and day is None:
            year = parts.get("year", 1970)
            month = parts.get("month", 1)
            day = 1

            # dddd => first day of week after epoch
            # dddd YYYY => first day of week in specified year
            # dddd MM YYYY => first day of week in specified year and month
            # dddd MM => first day after epoch in specified month
            next_weekday_dt = next_weekday(datetime(year, month, day), day_of_week)
            parts["year"] = next_weekday_dt.year
            parts["month"] = next_weekday_dt.month
            parts["day"] = next_weekday_dt.day

        am_pm = parts.get("am_pm")
        hour = parts.get("hour", 0)

        if am_pm == "pm" and hour < 12:
            hour += 12
        elif am_pm == "am" and hour == 12:
            hour = 0

        # Support for midnight at the end of day
        if hour == 24:
            if parts.get("minute", 0) != 0:
                raise ParserError("Midnight at the end of day must not contain minutes")
            if parts.get("second", 0) != 0:
                raise ParserError("Midnight at the end of day must not contain seconds")
            if parts.get("microsecond", 0) != 0:
                raise ParserError(
                    "Midnight at the end of day must not contain microseconds"
                )
            hour = 0
            day_increment = 1
        else:
            day_increment = 0

        # account for rounding up to 1000000
        microsecond = parts.get("microsecond", 0)
        if microsecond == 1000000:
            microsecond = 0
            second_increment = 1
        else:
            second_increment = 0

        increment = timedelta(days=day_increment, seconds=second_increment)

        return (
            datetime(
                year=parts.get("year", 1),
                month=parts.get("month", 1),
                day=parts.get("day", 1),
                hour=hour,
                minute=parts.get("minute", 0),
                second=parts.get("second", 0),
                microsecond=microsecond,
                tzinfo=parts.get("tzinfo"),
            )
            + increment
        )

    def _parse_multiformat(self, string: str, formats: Iterable[str]) -> datetime:
        _datetime: Optional[datetime] = None

        for fmt in formats:
            try:
                _datetime = self.parse(string, fmt)
                break
            except ParserMatchError:
                pass

        if _datetime is None:
            supported_formats = ", ".join(formats)
            raise ParserError(
                f"Could not match input {string!r} to any of the following formats: {supported_formats}."
            )

        return _datetime

    # generates a capture group of choices separated by an OR operator
    @staticmethod
    def _generate_choice_re(
        choices: Iterable[str], flags: Union[int, re.RegexFlag] = 0
    ) -> Pattern[str]:
        return re.compile(r"({})".format("|".join(choices)), flags=flags)


class TzinfoParser:
    _TZINFO_RE: ClassVar[Pattern[str]] = re.compile(
        r"^(?:\(UTC)*([\+\-])?(\d{2})(?:\:?(\d{2}))?"
    )

    @classmethod
    def parse(cls, tzinfo_string: str) -> dt_tzinfo:
        tzinfo: Optional[dt_tzinfo] = None

        if tzinfo_string == "local":
            tzinfo = tz.tzlocal()

        elif tzinfo_string in ["utc", "UTC", "Z"]:
            tzinfo = tz.tzutc()

        else:
            iso_match = cls._TZINFO_RE.match(tzinfo_string)

            if iso_match:
                sign: Optional[str]
                hours: str
                minutes: Union[str, int, None]
                sign, hours, minutes = iso_match.groups()
                seconds = int(hours) * 3600 + int(minutes or 0) * 60

                if sign == "-":
                    seconds *= -1

                tzinfo = tz.tzoffset(None, seconds)

            else:
                tzinfo = tz.gettz(tzinfo_string)

        if tzinfo is None:
            raise ParserError(f"Could not parse timezone expression {tzinfo_string!r}.")

        return tzinfo
