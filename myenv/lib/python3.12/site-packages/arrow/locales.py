"""Provides internationalization for arrow in over 60 languages and dialects."""

import sys
from math import trunc
from typing import (
    Any,
    ClassVar,
    Dict,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
    cast,
)

if sys.version_info < (3, 8):  # pragma: no cover
    from typing_extensions import Literal
else:
    from typing import Literal  # pragma: no cover

TimeFrameLiteral = Literal[
    "now",
    "second",
    "seconds",
    "minute",
    "minutes",
    "hour",
    "hours",
    "day",
    "days",
    "week",
    "weeks",
    "month",
    "months",
    "quarter",
    "quarters",
    "year",
    "years",
]

_TimeFrameElements = Union[
    str, Sequence[str], Mapping[str, str], Mapping[str, Sequence[str]]
]

_locale_map: Dict[str, Type["Locale"]] = {}


def get_locale(name: str) -> "Locale":
    """Returns an appropriate :class:`Locale <arrow.locales.Locale>`
    corresponding to an input locale name.

    :param name: the name of the locale.

    """

    normalized_locale_name = name.lower().replace("_", "-")
    locale_cls = _locale_map.get(normalized_locale_name)

    if locale_cls is None:
        raise ValueError(f"Unsupported locale {normalized_locale_name!r}.")

    return locale_cls()


def get_locale_by_class_name(name: str) -> "Locale":
    """Returns an appropriate :class:`Locale <arrow.locales.Locale>`
    corresponding to an locale class name.

    :param name: the name of the locale class.

    """
    locale_cls: Optional[Type[Locale]] = globals().get(name)

    if locale_cls is None:
        raise ValueError(f"Unsupported locale {name!r}.")

    return locale_cls()


class Locale:
    """Represents locale-specific data and functionality."""

    names: ClassVar[List[str]] = []

    timeframes: ClassVar[Mapping[TimeFrameLiteral, _TimeFrameElements]] = {
        "now": "",
        "second": "",
        "seconds": "",
        "minute": "",
        "minutes": "",
        "hour": "",
        "hours": "",
        "day": "",
        "days": "",
        "week": "",
        "weeks": "",
        "month": "",
        "months": "",
        "quarter": "",
        "quarters": "",
        "year": "",
        "years": "",
    }

    meridians: ClassVar[Dict[str, str]] = {"am": "", "pm": "", "AM": "", "PM": ""}

    past: ClassVar[str]
    future: ClassVar[str]
    and_word: ClassVar[Optional[str]] = None

    month_names: ClassVar[List[str]] = []
    month_abbreviations: ClassVar[List[str]] = []

    day_names: ClassVar[List[str]] = []
    day_abbreviations: ClassVar[List[str]] = []

    ordinal_day_re: ClassVar[str] = r"(\d+)"

    _month_name_to_ordinal: Optional[Dict[str, int]]

    def __init_subclass__(cls, **kwargs: Any) -> None:
        for locale_name in cls.names:
            if locale_name in _locale_map:
                raise LookupError(f"Duplicated locale name: {locale_name}")

            _locale_map[locale_name.lower().replace("_", "-")] = cls

    def __init__(self) -> None:
        self._month_name_to_ordinal = None

    def describe(
        self,
        timeframe: TimeFrameLiteral,
        delta: Union[float, int] = 0,
        only_distance: bool = False,
    ) -> str:
        """Describes a delta within a timeframe in plain language.

        :param timeframe: a string representing a timeframe.
        :param delta: a quantity representing a delta in a timeframe.
        :param only_distance: return only distance eg: "11 seconds" without "in" or "ago" keywords
        """

        humanized = self._format_timeframe(timeframe, trunc(delta))
        if not only_distance:
            humanized = self._format_relative(humanized, timeframe, delta)

        return humanized

    def describe_multi(
        self,
        timeframes: Sequence[Tuple[TimeFrameLiteral, Union[int, float]]],
        only_distance: bool = False,
    ) -> str:
        """Describes a delta within multiple timeframes in plain language.

        :param timeframes: a list of string, quantity pairs each representing a timeframe and delta.
        :param only_distance: return only distance eg: "2 hours and 11 seconds" without "in" or "ago" keywords
        """

        parts = [
            self._format_timeframe(timeframe, trunc(delta))
            for timeframe, delta in timeframes
        ]
        if self.and_word:
            parts.insert(-1, self.and_word)
        humanized = " ".join(parts)

        if not only_distance:
            # Needed to determine the correct relative string to use
            timeframe_value = 0

            for _, unit_value in timeframes:
                if trunc(unit_value) != 0:
                    timeframe_value = trunc(unit_value)
                    break

            # Note it doesn't matter the timeframe unit we use on the call, only the value
            humanized = self._format_relative(humanized, "seconds", timeframe_value)

        return humanized

    def day_name(self, day: int) -> str:
        """Returns the day name for a specified day of the week.

        :param day: the ``int`` day of the week (1-7).

        """

        return self.day_names[day]

    def day_abbreviation(self, day: int) -> str:
        """Returns the day abbreviation for a specified day of the week.

        :param day: the ``int`` day of the week (1-7).

        """

        return self.day_abbreviations[day]

    def month_name(self, month: int) -> str:
        """Returns the month name for a specified month of the year.

        :param month: the ``int`` month of the year (1-12).

        """

        return self.month_names[month]

    def month_abbreviation(self, month: int) -> str:
        """Returns the month abbreviation for a specified month of the year.

        :param month: the ``int`` month of the year (1-12).

        """

        return self.month_abbreviations[month]

    def month_number(self, name: str) -> Optional[int]:
        """Returns the month number for a month specified by name or abbreviation.

        :param name: the month name or abbreviation.

        """

        if self._month_name_to_ordinal is None:
            self._month_name_to_ordinal = self._name_to_ordinal(self.month_names)
            self._month_name_to_ordinal.update(
                self._name_to_ordinal(self.month_abbreviations)
            )

        return self._month_name_to_ordinal.get(name)

    def year_full(self, year: int) -> str:
        """Returns the year for specific locale if available

        :param year: the ``int`` year (4-digit)
        """
        return f"{year:04d}"

    def year_abbreviation(self, year: int) -> str:
        """Returns the year for specific locale if available

        :param year: the ``int`` year (4-digit)
        """
        return f"{year:04d}"[2:]

    def meridian(self, hour: int, token: Any) -> Optional[str]:
        """Returns the meridian indicator for a specified hour and format token.

        :param hour: the ``int`` hour of the day.
        :param token: the format token.
        """

        if token == "a":
            return self.meridians["am"] if hour < 12 else self.meridians["pm"]
        if token == "A":
            return self.meridians["AM"] if hour < 12 else self.meridians["PM"]
        return None

    def ordinal_number(self, n: int) -> str:
        """Returns the ordinal format of a given integer

        :param n: an integer
        """
        return self._ordinal_number(n)

    def _ordinal_number(self, n: int) -> str:
        return f"{n}"

    def _name_to_ordinal(self, lst: Sequence[str]) -> Dict[str, int]:
        return {elem.lower(): i for i, elem in enumerate(lst[1:], 1)}

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        # TODO: remove cast
        return cast(str, self.timeframes[timeframe]).format(trunc(abs(delta)))

    def _format_relative(
        self,
        humanized: str,
        timeframe: TimeFrameLiteral,
        delta: Union[float, int],
    ) -> str:
        if timeframe == "now":
            return humanized

        direction = self.past if delta < 0 else self.future

        return direction.format(humanized)


class EnglishLocale(Locale):
    names = [
        "en",
        "en-us",
        "en-gb",
        "en-au",
        "en-be",
        "en-jp",
        "en-za",
        "en-ca",
        "en-ph",
    ]

    past = "{0} ago"
    future = "in {0}"
    and_word = "and"

    timeframes = {
        "now": "just now",
        "second": "a second",
        "seconds": "{0} seconds",
        "minute": "a minute",
        "minutes": "{0} minutes",
        "hour": "an hour",
        "hours": "{0} hours",
        "day": "a day",
        "days": "{0} days",
        "week": "a week",
        "weeks": "{0} weeks",
        "month": "a month",
        "months": "{0} months",
        "quarter": "a quarter",
        "quarters": "{0} quarters",
        "year": "a year",
        "years": "{0} years",
    }

    meridians = {"am": "am", "pm": "pm", "AM": "AM", "PM": "PM"}

    month_names = [
        "",
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ]
    month_abbreviations = [
        "",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]

    day_names = [
        "",
        "Monday",
        "Tuesday",
        "Wednesday",
        "Thursday",
        "Friday",
        "Saturday",
        "Sunday",
    ]
    day_abbreviations = ["", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]

    ordinal_day_re = r"((?P<value>[2-3]?1(?=st)|[2-3]?2(?=nd)|[2-3]?3(?=rd)|[1-3]?[04-9](?=th)|1[1-3](?=th))(st|nd|rd|th))"

    def _ordinal_number(self, n: int) -> str:
        if n % 100 not in (11, 12, 13):
            remainder = abs(n) % 10
            if remainder == 1:
                return f"{n}st"
            elif remainder == 2:
                return f"{n}nd"
            elif remainder == 3:
                return f"{n}rd"
        return f"{n}th"

    def describe(
        self,
        timeframe: TimeFrameLiteral,
        delta: Union[int, float] = 0,
        only_distance: bool = False,
    ) -> str:
        """Describes a delta within a timeframe in plain language.

        :param timeframe: a string representing a timeframe.
        :param delta: a quantity representing a delta in a timeframe.
        :param only_distance: return only distance eg: "11 seconds" without "in" or "ago" keywords
        """

        humanized = super().describe(timeframe, delta, only_distance)
        if only_distance and timeframe == "now":
            humanized = "instantly"

        return humanized


class ItalianLocale(Locale):
    names = ["it", "it-it"]
    past = "{0} fa"
    future = "tra {0}"
    and_word = "e"

    timeframes = {
        "now": "adesso",
        "second": "un secondo",
        "seconds": "{0} qualche secondo",
        "minute": "un minuto",
        "minutes": "{0} minuti",
        "hour": "un'ora",
        "hours": "{0} ore",
        "day": "un giorno",
        "days": "{0} giorni",
        "week": "una settimana",
        "weeks": "{0} settimane",
        "month": "un mese",
        "months": "{0} mesi",
        "year": "un anno",
        "years": "{0} anni",
    }

    month_names = [
        "",
        "gennaio",
        "febbraio",
        "marzo",
        "aprile",
        "maggio",
        "giugno",
        "luglio",
        "agosto",
        "settembre",
        "ottobre",
        "novembre",
        "dicembre",
    ]
    month_abbreviations = [
        "",
        "gen",
        "feb",
        "mar",
        "apr",
        "mag",
        "giu",
        "lug",
        "ago",
        "set",
        "ott",
        "nov",
        "dic",
    ]

    day_names = [
        "",
        "lunedì",
        "martedì",
        "mercoledì",
        "giovedì",
        "venerdì",
        "sabato",
        "domenica",
    ]
    day_abbreviations = ["", "lun", "mar", "mer", "gio", "ven", "sab", "dom"]

    ordinal_day_re = r"((?P<value>[1-3]?[0-9](?=[ºª]))[ºª])"

    def _ordinal_number(self, n: int) -> str:
        return f"{n}º"


class SpanishLocale(Locale):
    names = ["es", "es-es"]
    past = "hace {0}"
    future = "en {0}"
    and_word = "y"

    timeframes = {
        "now": "ahora",
        "second": "un segundo",
        "seconds": "{0} segundos",
        "minute": "un minuto",
        "minutes": "{0} minutos",
        "hour": "una hora",
        "hours": "{0} horas",
        "day": "un día",
        "days": "{0} días",
        "week": "una semana",
        "weeks": "{0} semanas",
        "month": "un mes",
        "months": "{0} meses",
        "year": "un año",
        "years": "{0} años",
    }

    meridians = {"am": "am", "pm": "pm", "AM": "AM", "PM": "PM"}

    month_names = [
        "",
        "enero",
        "febrero",
        "marzo",
        "abril",
        "mayo",
        "junio",
        "julio",
        "agosto",
        "septiembre",
        "octubre",
        "noviembre",
        "diciembre",
    ]
    month_abbreviations = [
        "",
        "ene",
        "feb",
        "mar",
        "abr",
        "may",
        "jun",
        "jul",
        "ago",
        "sep",
        "oct",
        "nov",
        "dic",
    ]

    day_names = [
        "",
        "lunes",
        "martes",
        "miércoles",
        "jueves",
        "viernes",
        "sábado",
        "domingo",
    ]
    day_abbreviations = ["", "lun", "mar", "mie", "jue", "vie", "sab", "dom"]

    ordinal_day_re = r"((?P<value>[1-3]?[0-9](?=[ºª]))[ºª])"

    def _ordinal_number(self, n: int) -> str:
        return f"{n}º"


class FrenchBaseLocale(Locale):
    past = "il y a {0}"
    future = "dans {0}"
    and_word = "et"

    timeframes = {
        "now": "maintenant",
        "second": "une seconde",
        "seconds": "{0} secondes",
        "minute": "une minute",
        "minutes": "{0} minutes",
        "hour": "une heure",
        "hours": "{0} heures",
        "day": "un jour",
        "days": "{0} jours",
        "week": "une semaine",
        "weeks": "{0} semaines",
        "month": "un mois",
        "months": "{0} mois",
        "year": "un an",
        "years": "{0} ans",
    }

    month_names = [
        "",
        "janvier",
        "février",
        "mars",
        "avril",
        "mai",
        "juin",
        "juillet",
        "août",
        "septembre",
        "octobre",
        "novembre",
        "décembre",
    ]

    day_names = [
        "",
        "lundi",
        "mardi",
        "mercredi",
        "jeudi",
        "vendredi",
        "samedi",
        "dimanche",
    ]
    day_abbreviations = ["", "lun", "mar", "mer", "jeu", "ven", "sam", "dim"]

    ordinal_day_re = (
        r"((?P<value>\b1(?=er\b)|[1-3]?[02-9](?=e\b)|[1-3]1(?=e\b))(er|e)\b)"
    )

    def _ordinal_number(self, n: int) -> str:
        if abs(n) == 1:
            return f"{n}er"
        return f"{n}e"


class FrenchLocale(FrenchBaseLocale, Locale):
    names = ["fr", "fr-fr"]

    month_abbreviations = [
        "",
        "janv",
        "févr",
        "mars",
        "avr",
        "mai",
        "juin",
        "juil",
        "août",
        "sept",
        "oct",
        "nov",
        "déc",
    ]


class FrenchCanadianLocale(FrenchBaseLocale, Locale):
    names = ["fr-ca"]

    month_abbreviations = [
        "",
        "janv",
        "févr",
        "mars",
        "avr",
        "mai",
        "juin",
        "juill",
        "août",
        "sept",
        "oct",
        "nov",
        "déc",
    ]


class GreekLocale(Locale):
    names = ["el", "el-gr"]

    past = "{0} πριν"
    future = "σε {0}"
    and_word = "και"

    timeframes = {
        "now": "τώρα",
        "second": "ένα δεύτερο",
        "seconds": "{0} δευτερόλεπτα",
        "minute": "ένα λεπτό",
        "minutes": "{0} λεπτά",
        "hour": "μία ώρα",
        "hours": "{0} ώρες",
        "day": "μία μέρα",
        "days": "{0} μέρες",
        "week": "μία εβδομάδα",
        "weeks": "{0} εβδομάδες",
        "month": "ένα μήνα",
        "months": "{0} μήνες",
        "year": "ένα χρόνο",
        "years": "{0} χρόνια",
    }

    month_names = [
        "",
        "Ιανουαρίου",
        "Φεβρουαρίου",
        "Μαρτίου",
        "Απριλίου",
        "Μαΐου",
        "Ιουνίου",
        "Ιουλίου",
        "Αυγούστου",
        "Σεπτεμβρίου",
        "Οκτωβρίου",
        "Νοεμβρίου",
        "Δεκεμβρίου",
    ]
    month_abbreviations = [
        "",
        "Ιαν",
        "Φεβ",
        "Μαρ",
        "Απρ",
        "Μαϊ",
        "Ιον",
        "Ιολ",
        "Αυγ",
        "Σεπ",
        "Οκτ",
        "Νοε",
        "Δεκ",
    ]

    day_names = [
        "",
        "Δευτέρα",
        "Τρίτη",
        "Τετάρτη",
        "Πέμπτη",
        "Παρασκευή",
        "Σάββατο",
        "Κυριακή",
    ]
    day_abbreviations = ["", "Δευ", "Τρι", "Τετ", "Πεμ", "Παρ", "Σαβ", "Κυρ"]


class JapaneseLocale(Locale):
    names = ["ja", "ja-jp"]

    past = "{0}前"
    future = "{0}後"
    and_word = ""

    timeframes = {
        "now": "現在",
        "second": "1秒",
        "seconds": "{0}秒",
        "minute": "1分",
        "minutes": "{0}分",
        "hour": "1時間",
        "hours": "{0}時間",
        "day": "1日",
        "days": "{0}日",
        "week": "1週間",
        "weeks": "{0}週間",
        "month": "1ヶ月",
        "months": "{0}ヶ月",
        "year": "1年",
        "years": "{0}年",
    }

    month_names = [
        "",
        "1月",
        "2月",
        "3月",
        "4月",
        "5月",
        "6月",
        "7月",
        "8月",
        "9月",
        "10月",
        "11月",
        "12月",
    ]
    month_abbreviations = [
        "",
        " 1",
        " 2",
        " 3",
        " 4",
        " 5",
        " 6",
        " 7",
        " 8",
        " 9",
        "10",
        "11",
        "12",
    ]

    day_names = ["", "月曜日", "火曜日", "水曜日", "木曜日", "金曜日", "土曜日", "日曜日"]
    day_abbreviations = ["", "月", "火", "水", "木", "金", "土", "日"]


class SwedishLocale(Locale):
    names = ["sv", "sv-se"]

    past = "för {0} sen"
    future = "om {0}"
    and_word = "och"

    timeframes = {
        "now": "just nu",
        "second": "en sekund",
        "seconds": "{0} sekunder",
        "minute": "en minut",
        "minutes": "{0} minuter",
        "hour": "en timme",
        "hours": "{0} timmar",
        "day": "en dag",
        "days": "{0} dagar",
        "week": "en vecka",
        "weeks": "{0} veckor",
        "month": "en månad",
        "months": "{0} månader",
        "year": "ett år",
        "years": "{0} år",
    }

    month_names = [
        "",
        "januari",
        "februari",
        "mars",
        "april",
        "maj",
        "juni",
        "juli",
        "augusti",
        "september",
        "oktober",
        "november",
        "december",
    ]
    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mar",
        "apr",
        "maj",
        "jun",
        "jul",
        "aug",
        "sep",
        "okt",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "måndag",
        "tisdag",
        "onsdag",
        "torsdag",
        "fredag",
        "lördag",
        "söndag",
    ]
    day_abbreviations = ["", "mån", "tis", "ons", "tor", "fre", "lör", "sön"]


class FinnishLocale(Locale):
    names = ["fi", "fi-fi"]

    # The finnish grammar is very complex, and its hard to convert
    # 1-to-1 to something like English.

    past = "{0} sitten"
    future = "{0} kuluttua"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "juuri nyt",
        "second": {"past": "sekunti", "future": "sekunnin"},
        "seconds": {"past": "{0} sekuntia", "future": "{0} sekunnin"},
        "minute": {"past": "minuutti", "future": "minuutin"},
        "minutes": {"past": "{0} minuuttia", "future": "{0} minuutin"},
        "hour": {"past": "tunti", "future": "tunnin"},
        "hours": {"past": "{0} tuntia", "future": "{0} tunnin"},
        "day": {"past": "päivä", "future": "päivän"},
        "days": {"past": "{0} päivää", "future": "{0} päivän"},
        "week": {"past": "viikko", "future": "viikon"},
        "weeks": {"past": "{0} viikkoa", "future": "{0} viikon"},
        "month": {"past": "kuukausi", "future": "kuukauden"},
        "months": {"past": "{0} kuukautta", "future": "{0} kuukauden"},
        "year": {"past": "vuosi", "future": "vuoden"},
        "years": {"past": "{0} vuotta", "future": "{0} vuoden"},
    }

    # Months and days are lowercase in Finnish
    month_names = [
        "",
        "tammikuu",
        "helmikuu",
        "maaliskuu",
        "huhtikuu",
        "toukokuu",
        "kesäkuu",
        "heinäkuu",
        "elokuu",
        "syyskuu",
        "lokakuu",
        "marraskuu",
        "joulukuu",
    ]

    month_abbreviations = [
        "",
        "tammi",
        "helmi",
        "maalis",
        "huhti",
        "touko",
        "kesä",
        "heinä",
        "elo",
        "syys",
        "loka",
        "marras",
        "joulu",
    ]

    day_names = [
        "",
        "maanantai",
        "tiistai",
        "keskiviikko",
        "torstai",
        "perjantai",
        "lauantai",
        "sunnuntai",
    ]

    day_abbreviations = ["", "ma", "ti", "ke", "to", "pe", "la", "su"]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]

        if isinstance(form, Mapping):
            if delta < 0:
                form = form["past"]
            else:
                form = form["future"]

        return form.format(abs(delta))

    def _ordinal_number(self, n: int) -> str:
        return f"{n}."


class ChineseCNLocale(Locale):
    names = ["zh", "zh-cn"]

    past = "{0}前"
    future = "{0}后"

    timeframes = {
        "now": "刚才",
        "second": "1秒",
        "seconds": "{0}秒",
        "minute": "1分钟",
        "minutes": "{0}分钟",
        "hour": "1小时",
        "hours": "{0}小时",
        "day": "1天",
        "days": "{0}天",
        "week": "1周",
        "weeks": "{0}周",
        "month": "1个月",
        "months": "{0}个月",
        "year": "1年",
        "years": "{0}年",
    }

    month_names = [
        "",
        "一月",
        "二月",
        "三月",
        "四月",
        "五月",
        "六月",
        "七月",
        "八月",
        "九月",
        "十月",
        "十一月",
        "十二月",
    ]
    month_abbreviations = [
        "",
        " 1",
        " 2",
        " 3",
        " 4",
        " 5",
        " 6",
        " 7",
        " 8",
        " 9",
        "10",
        "11",
        "12",
    ]

    day_names = ["", "星期一", "星期二", "星期三", "星期四", "星期五", "星期六", "星期日"]
    day_abbreviations = ["", "一", "二", "三", "四", "五", "六", "日"]


class ChineseTWLocale(Locale):
    names = ["zh-tw"]

    past = "{0}前"
    future = "{0}後"
    and_word = "和"

    timeframes = {
        "now": "剛才",
        "second": "1秒",
        "seconds": "{0}秒",
        "minute": "1分鐘",
        "minutes": "{0}分鐘",
        "hour": "1小時",
        "hours": "{0}小時",
        "day": "1天",
        "days": "{0}天",
        "week": "1週",
        "weeks": "{0}週",
        "month": "1個月",
        "months": "{0}個月",
        "year": "1年",
        "years": "{0}年",
    }

    month_names = [
        "",
        "1月",
        "2月",
        "3月",
        "4月",
        "5月",
        "6月",
        "7月",
        "8月",
        "9月",
        "10月",
        "11月",
        "12月",
    ]
    month_abbreviations = [
        "",
        " 1",
        " 2",
        " 3",
        " 4",
        " 5",
        " 6",
        " 7",
        " 8",
        " 9",
        "10",
        "11",
        "12",
    ]

    day_names = ["", "週一", "週二", "週三", "週四", "週五", "週六", "週日"]
    day_abbreviations = ["", "一", "二", "三", "四", "五", "六", "日"]


class HongKongLocale(Locale):
    names = ["zh-hk"]

    past = "{0}前"
    future = "{0}後"

    timeframes = {
        "now": "剛才",
        "second": "1秒",
        "seconds": "{0}秒",
        "minute": "1分鐘",
        "minutes": "{0}分鐘",
        "hour": "1小時",
        "hours": "{0}小時",
        "day": "1天",
        "days": "{0}天",
        "week": "1星期",
        "weeks": "{0}星期",
        "month": "1個月",
        "months": "{0}個月",
        "year": "1年",
        "years": "{0}年",
    }

    month_names = [
        "",
        "1月",
        "2月",
        "3月",
        "4月",
        "5月",
        "6月",
        "7月",
        "8月",
        "9月",
        "10月",
        "11月",
        "12月",
    ]
    month_abbreviations = [
        "",
        " 1",
        " 2",
        " 3",
        " 4",
        " 5",
        " 6",
        " 7",
        " 8",
        " 9",
        "10",
        "11",
        "12",
    ]

    day_names = ["", "星期一", "星期二", "星期三", "星期四", "星期五", "星期六", "星期日"]
    day_abbreviations = ["", "一", "二", "三", "四", "五", "六", "日"]


class KoreanLocale(Locale):
    names = ["ko", "ko-kr"]

    past = "{0} 전"
    future = "{0} 후"

    timeframes = {
        "now": "지금",
        "second": "1초",
        "seconds": "{0}초",
        "minute": "1분",
        "minutes": "{0}분",
        "hour": "한시간",
        "hours": "{0}시간",
        "day": "하루",
        "days": "{0}일",
        "week": "1주",
        "weeks": "{0}주",
        "month": "한달",
        "months": "{0}개월",
        "year": "1년",
        "years": "{0}년",
    }

    special_dayframes = {
        -3: "그끄제",
        -2: "그제",
        -1: "어제",
        1: "내일",
        2: "모레",
        3: "글피",
        4: "그글피",
    }

    special_yearframes = {-2: "제작년", -1: "작년", 1: "내년", 2: "내후년"}

    month_names = [
        "",
        "1월",
        "2월",
        "3월",
        "4월",
        "5월",
        "6월",
        "7월",
        "8월",
        "9월",
        "10월",
        "11월",
        "12월",
    ]
    month_abbreviations = [
        "",
        " 1",
        " 2",
        " 3",
        " 4",
        " 5",
        " 6",
        " 7",
        " 8",
        " 9",
        "10",
        "11",
        "12",
    ]

    day_names = ["", "월요일", "화요일", "수요일", "목요일", "금요일", "토요일", "일요일"]
    day_abbreviations = ["", "월", "화", "수", "목", "금", "토", "일"]

    def _ordinal_number(self, n: int) -> str:
        ordinals = ["0", "첫", "두", "세", "네", "다섯", "여섯", "일곱", "여덟", "아홉", "열"]
        if n < len(ordinals):
            return f"{ordinals[n]}번째"
        return f"{n}번째"

    def _format_relative(
        self,
        humanized: str,
        timeframe: TimeFrameLiteral,
        delta: Union[float, int],
    ) -> str:
        if timeframe in ("day", "days"):
            special = self.special_dayframes.get(int(delta))
            if special:
                return special
        elif timeframe in ("year", "years"):
            special = self.special_yearframes.get(int(delta))
            if special:
                return special

        return super()._format_relative(humanized, timeframe, delta)


# derived locale types & implementations.
class DutchLocale(Locale):
    names = ["nl", "nl-nl"]

    past = "{0} geleden"
    future = "over {0}"

    timeframes = {
        "now": "nu",
        "second": "een seconde",
        "seconds": "{0} seconden",
        "minute": "een minuut",
        "minutes": "{0} minuten",
        "hour": "een uur",
        "hours": "{0} uur",
        "day": "een dag",
        "days": "{0} dagen",
        "week": "een week",
        "weeks": "{0} weken",
        "month": "een maand",
        "months": "{0} maanden",
        "year": "een jaar",
        "years": "{0} jaar",
    }

    # In Dutch names of months and days are not starting with a capital letter
    # like in the English language.
    month_names = [
        "",
        "januari",
        "februari",
        "maart",
        "april",
        "mei",
        "juni",
        "juli",
        "augustus",
        "september",
        "oktober",
        "november",
        "december",
    ]
    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mrt",
        "apr",
        "mei",
        "jun",
        "jul",
        "aug",
        "sep",
        "okt",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "maandag",
        "dinsdag",
        "woensdag",
        "donderdag",
        "vrijdag",
        "zaterdag",
        "zondag",
    ]
    day_abbreviations = ["", "ma", "di", "wo", "do", "vr", "za", "zo"]


class SlavicBaseLocale(Locale):
    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]
        delta = abs(delta)

        if isinstance(form, Mapping):
            if delta % 10 == 1 and delta % 100 != 11:
                form = form["singular"]
            elif 2 <= delta % 10 <= 4 and (delta % 100 < 10 or delta % 100 >= 20):
                form = form["dual"]
            else:
                form = form["plural"]

        return form.format(delta)


class BelarusianLocale(SlavicBaseLocale):
    names = ["be", "be-by"]

    past = "{0} таму"
    future = "праз {0}"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "зараз",
        "second": "секунду",
        "seconds": "{0} некалькі секунд",
        "minute": "хвіліну",
        "minutes": {
            "singular": "{0} хвіліну",
            "dual": "{0} хвіліны",
            "plural": "{0} хвілін",
        },
        "hour": "гадзіну",
        "hours": {
            "singular": "{0} гадзіну",
            "dual": "{0} гадзіны",
            "plural": "{0} гадзін",
        },
        "day": "дзень",
        "days": {"singular": "{0} дзень", "dual": "{0} дні", "plural": "{0} дзён"},
        "month": "месяц",
        "months": {
            "singular": "{0} месяц",
            "dual": "{0} месяцы",
            "plural": "{0} месяцаў",
        },
        "year": "год",
        "years": {"singular": "{0} год", "dual": "{0} гады", "plural": "{0} гадоў"},
    }

    month_names = [
        "",
        "студзеня",
        "лютага",
        "сакавіка",
        "красавіка",
        "траўня",
        "чэрвеня",
        "ліпеня",
        "жніўня",
        "верасня",
        "кастрычніка",
        "лістапада",
        "снежня",
    ]
    month_abbreviations = [
        "",
        "студ",
        "лют",
        "сак",
        "крас",
        "трав",
        "чэрв",
        "ліп",
        "жнів",
        "вер",
        "каст",
        "ліст",
        "снеж",
    ]

    day_names = [
        "",
        "панядзелак",
        "аўторак",
        "серада",
        "чацвер",
        "пятніца",
        "субота",
        "нядзеля",
    ]
    day_abbreviations = ["", "пн", "ат", "ср", "чц", "пт", "сб", "нд"]


class PolishLocale(SlavicBaseLocale):
    names = ["pl", "pl-pl"]

    past = "{0} temu"
    future = "za {0}"

    # The nouns should be in genitive case (Polish: "dopełniacz")
    # in order to correctly form `past` & `future` expressions.
    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "teraz",
        "second": "sekundę",
        "seconds": {
            "singular": "{0} sekund",
            "dual": "{0} sekundy",
            "plural": "{0} sekund",
        },
        "minute": "minutę",
        "minutes": {
            "singular": "{0} minut",
            "dual": "{0} minuty",
            "plural": "{0} minut",
        },
        "hour": "godzinę",
        "hours": {
            "singular": "{0} godzin",
            "dual": "{0} godziny",
            "plural": "{0} godzin",
        },
        "day": "dzień",
        "days": "{0} dni",
        "week": "tydzień",
        "weeks": {
            "singular": "{0} tygodni",
            "dual": "{0} tygodnie",
            "plural": "{0} tygodni",
        },
        "month": "miesiąc",
        "months": {
            "singular": "{0} miesięcy",
            "dual": "{0} miesiące",
            "plural": "{0} miesięcy",
        },
        "year": "rok",
        "years": {"singular": "{0} lat", "dual": "{0} lata", "plural": "{0} lat"},
    }

    month_names = [
        "",
        "styczeń",
        "luty",
        "marzec",
        "kwiecień",
        "maj",
        "czerwiec",
        "lipiec",
        "sierpień",
        "wrzesień",
        "październik",
        "listopad",
        "grudzień",
    ]
    month_abbreviations = [
        "",
        "sty",
        "lut",
        "mar",
        "kwi",
        "maj",
        "cze",
        "lip",
        "sie",
        "wrz",
        "paź",
        "lis",
        "gru",
    ]

    day_names = [
        "",
        "poniedziałek",
        "wtorek",
        "środa",
        "czwartek",
        "piątek",
        "sobota",
        "niedziela",
    ]
    day_abbreviations = ["", "Pn", "Wt", "Śr", "Czw", "Pt", "So", "Nd"]


class RussianLocale(SlavicBaseLocale):
    names = ["ru", "ru-ru"]

    past = "{0} назад"
    future = "через {0}"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "сейчас",
        "second": "секунда",
        "seconds": {
            "singular": "{0} секунду",
            "dual": "{0} секунды",
            "plural": "{0} секунд",
        },
        "minute": "минуту",
        "minutes": {
            "singular": "{0} минуту",
            "dual": "{0} минуты",
            "plural": "{0} минут",
        },
        "hour": "час",
        "hours": {"singular": "{0} час", "dual": "{0} часа", "plural": "{0} часов"},
        "day": "день",
        "days": {"singular": "{0} день", "dual": "{0} дня", "plural": "{0} дней"},
        "week": "неделю",
        "weeks": {
            "singular": "{0} неделю",
            "dual": "{0} недели",
            "plural": "{0} недель",
        },
        "month": "месяц",
        "months": {
            "singular": "{0} месяц",
            "dual": "{0} месяца",
            "plural": "{0} месяцев",
        },
        "quarter": "квартал",
        "quarters": {
            "singular": "{0} квартал",
            "dual": "{0} квартала",
            "plural": "{0} кварталов",
        },
        "year": "год",
        "years": {"singular": "{0} год", "dual": "{0} года", "plural": "{0} лет"},
    }

    month_names = [
        "",
        "января",
        "февраля",
        "марта",
        "апреля",
        "мая",
        "июня",
        "июля",
        "августа",
        "сентября",
        "октября",
        "ноября",
        "декабря",
    ]
    month_abbreviations = [
        "",
        "янв",
        "фев",
        "мар",
        "апр",
        "май",
        "июн",
        "июл",
        "авг",
        "сен",
        "окт",
        "ноя",
        "дек",
    ]

    day_names = [
        "",
        "понедельник",
        "вторник",
        "среда",
        "четверг",
        "пятница",
        "суббота",
        "воскресенье",
    ]
    day_abbreviations = ["", "пн", "вт", "ср", "чт", "пт", "сб", "вс"]


class AfrikaansLocale(Locale):
    names = ["af", "af-nl"]

    past = "{0} gelede"
    future = "in {0}"

    timeframes = {
        "now": "nou",
        "second": "n sekonde",
        "seconds": "{0} sekondes",
        "minute": "minuut",
        "minutes": "{0} minute",
        "hour": "uur",
        "hours": "{0} ure",
        "day": "een dag",
        "days": "{0} dae",
        "month": "een maand",
        "months": "{0} maande",
        "year": "een jaar",
        "years": "{0} jaar",
    }

    month_names = [
        "",
        "Januarie",
        "Februarie",
        "Maart",
        "April",
        "Mei",
        "Junie",
        "Julie",
        "Augustus",
        "September",
        "Oktober",
        "November",
        "Desember",
    ]
    month_abbreviations = [
        "",
        "Jan",
        "Feb",
        "Mrt",
        "Apr",
        "Mei",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Okt",
        "Nov",
        "Des",
    ]

    day_names = [
        "",
        "Maandag",
        "Dinsdag",
        "Woensdag",
        "Donderdag",
        "Vrydag",
        "Saterdag",
        "Sondag",
    ]
    day_abbreviations = ["", "Ma", "Di", "Wo", "Do", "Vr", "Za", "So"]


class BulgarianLocale(SlavicBaseLocale):
    names = ["bg", "bg-bg"]

    past = "{0} назад"
    future = "напред {0}"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "сега",
        "second": "секунда",
        "seconds": "{0} няколко секунди",
        "minute": "минута",
        "minutes": {
            "singular": "{0} минута",
            "dual": "{0} минути",
            "plural": "{0} минути",
        },
        "hour": "час",
        "hours": {"singular": "{0} час", "dual": "{0} часа", "plural": "{0} часа"},
        "day": "ден",
        "days": {"singular": "{0} ден", "dual": "{0} дни", "plural": "{0} дни"},
        "month": "месец",
        "months": {
            "singular": "{0} месец",
            "dual": "{0} месеца",
            "plural": "{0} месеца",
        },
        "year": "година",
        "years": {
            "singular": "{0} година",
            "dual": "{0} години",
            "plural": "{0} години",
        },
    }

    month_names = [
        "",
        "януари",
        "февруари",
        "март",
        "април",
        "май",
        "юни",
        "юли",
        "август",
        "септември",
        "октомври",
        "ноември",
        "декември",
    ]
    month_abbreviations = [
        "",
        "ян",
        "февр",
        "март",
        "апр",
        "май",
        "юни",
        "юли",
        "авг",
        "септ",
        "окт",
        "ноем",
        "дек",
    ]

    day_names = [
        "",
        "понеделник",
        "вторник",
        "сряда",
        "четвъртък",
        "петък",
        "събота",
        "неделя",
    ]
    day_abbreviations = ["", "пон", "вт", "ср", "четв", "пет", "съб", "нед"]


class UkrainianLocale(SlavicBaseLocale):
    names = ["ua", "uk", "uk-ua"]

    past = "{0} тому"
    future = "за {0}"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "зараз",
        "second": "секунда",
        "seconds": "{0} кілька секунд",
        "minute": "хвилину",
        "minutes": {
            "singular": "{0} хвилину",
            "dual": "{0} хвилини",
            "plural": "{0} хвилин",
        },
        "hour": "годину",
        "hours": {
            "singular": "{0} годину",
            "dual": "{0} години",
            "plural": "{0} годин",
        },
        "day": "день",
        "days": {"singular": "{0} день", "dual": "{0} дні", "plural": "{0} днів"},
        "month": "місяць",
        "months": {
            "singular": "{0} місяць",
            "dual": "{0} місяці",
            "plural": "{0} місяців",
        },
        "year": "рік",
        "years": {"singular": "{0} рік", "dual": "{0} роки", "plural": "{0} років"},
    }

    month_names = [
        "",
        "січня",
        "лютого",
        "березня",
        "квітня",
        "травня",
        "червня",
        "липня",
        "серпня",
        "вересня",
        "жовтня",
        "листопада",
        "грудня",
    ]
    month_abbreviations = [
        "",
        "січ",
        "лют",
        "бер",
        "квіт",
        "трав",
        "черв",
        "лип",
        "серп",
        "вер",
        "жовт",
        "лист",
        "груд",
    ]

    day_names = [
        "",
        "понеділок",
        "вівторок",
        "середа",
        "четвер",
        "п’ятниця",
        "субота",
        "неділя",
    ]
    day_abbreviations = ["", "пн", "вт", "ср", "чт", "пт", "сб", "нд"]


class MacedonianLocale(SlavicBaseLocale):
    names = ["mk", "mk-mk"]

    past = "пред {0}"
    future = "за {0}"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "сега",
        "second": "една секунда",
        "seconds": {
            "singular": "{0} секунда",
            "dual": "{0} секунди",
            "plural": "{0} секунди",
        },
        "minute": "една минута",
        "minutes": {
            "singular": "{0} минута",
            "dual": "{0} минути",
            "plural": "{0} минути",
        },
        "hour": "еден саат",
        "hours": {"singular": "{0} саат", "dual": "{0} саати", "plural": "{0} саати"},
        "day": "еден ден",
        "days": {"singular": "{0} ден", "dual": "{0} дена", "plural": "{0} дена"},
        "week": "една недела",
        "weeks": {
            "singular": "{0} недела",
            "dual": "{0} недели",
            "plural": "{0} недели",
        },
        "month": "еден месец",
        "months": {
            "singular": "{0} месец",
            "dual": "{0} месеци",
            "plural": "{0} месеци",
        },
        "year": "една година",
        "years": {
            "singular": "{0} година",
            "dual": "{0} години",
            "plural": "{0} години",
        },
    }

    meridians = {"am": "дп", "pm": "пп", "AM": "претпладне", "PM": "попладне"}

    month_names = [
        "",
        "Јануари",
        "Февруари",
        "Март",
        "Април",
        "Мај",
        "Јуни",
        "Јули",
        "Август",
        "Септември",
        "Октомври",
        "Ноември",
        "Декември",
    ]
    month_abbreviations = [
        "",
        "Јан",
        "Фев",
        "Мар",
        "Апр",
        "Мај",
        "Јун",
        "Јул",
        "Авг",
        "Септ",
        "Окт",
        "Ноем",
        "Декем",
    ]

    day_names = [
        "",
        "Понеделник",
        "Вторник",
        "Среда",
        "Четврток",
        "Петок",
        "Сабота",
        "Недела",
    ]
    day_abbreviations = [
        "",
        "Пон",
        "Вт",
        "Сре",
        "Чет",
        "Пет",
        "Саб",
        "Нед",
    ]


class GermanBaseLocale(Locale):
    past = "vor {0}"
    future = "in {0}"
    and_word = "und"

    timeframes: ClassVar[Dict[TimeFrameLiteral, str]] = {
        "now": "gerade eben",
        "second": "einer Sekunde",
        "seconds": "{0} Sekunden",
        "minute": "einer Minute",
        "minutes": "{0} Minuten",
        "hour": "einer Stunde",
        "hours": "{0} Stunden",
        "day": "einem Tag",
        "days": "{0} Tagen",
        "week": "einer Woche",
        "weeks": "{0} Wochen",
        "month": "einem Monat",
        "months": "{0} Monaten",
        "year": "einem Jahr",
        "years": "{0} Jahren",
    }

    timeframes_only_distance = timeframes.copy()
    timeframes_only_distance["second"] = "eine Sekunde"
    timeframes_only_distance["minute"] = "eine Minute"
    timeframes_only_distance["hour"] = "eine Stunde"
    timeframes_only_distance["day"] = "ein Tag"
    timeframes_only_distance["days"] = "{0} Tage"
    timeframes_only_distance["week"] = "eine Woche"
    timeframes_only_distance["month"] = "ein Monat"
    timeframes_only_distance["months"] = "{0} Monate"
    timeframes_only_distance["year"] = "ein Jahr"
    timeframes_only_distance["years"] = "{0} Jahre"

    month_names = [
        "",
        "Januar",
        "Februar",
        "März",
        "April",
        "Mai",
        "Juni",
        "Juli",
        "August",
        "September",
        "Oktober",
        "November",
        "Dezember",
    ]

    month_abbreviations = [
        "",
        "Jan",
        "Feb",
        "Mär",
        "Apr",
        "Mai",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Okt",
        "Nov",
        "Dez",
    ]

    day_names = [
        "",
        "Montag",
        "Dienstag",
        "Mittwoch",
        "Donnerstag",
        "Freitag",
        "Samstag",
        "Sonntag",
    ]

    day_abbreviations = ["", "Mo", "Di", "Mi", "Do", "Fr", "Sa", "So"]

    def _ordinal_number(self, n: int) -> str:
        return f"{n}."

    def describe(
        self,
        timeframe: TimeFrameLiteral,
        delta: Union[int, float] = 0,
        only_distance: bool = False,
    ) -> str:
        """Describes a delta within a timeframe in plain language.

        :param timeframe: a string representing a timeframe.
        :param delta: a quantity representing a delta in a timeframe.
        :param only_distance: return only distance eg: "11 seconds" without "in" or "ago" keywords
        """

        if not only_distance:
            return super().describe(timeframe, delta, only_distance)

        # German uses a different case without 'in' or 'ago'
        humanized: str = self.timeframes_only_distance[timeframe].format(
            trunc(abs(delta))
        )

        return humanized


class GermanLocale(GermanBaseLocale, Locale):
    names = ["de", "de-de"]


class SwissLocale(GermanBaseLocale, Locale):
    names = ["de-ch"]


class AustrianLocale(GermanBaseLocale, Locale):
    names = ["de-at"]

    month_names = [
        "",
        "Jänner",
        "Februar",
        "März",
        "April",
        "Mai",
        "Juni",
        "Juli",
        "August",
        "September",
        "Oktober",
        "November",
        "Dezember",
    ]


class NorwegianLocale(Locale):
    names = ["nb", "nb-no"]

    past = "for {0} siden"
    future = "om {0}"

    timeframes = {
        "now": "nå nettopp",
        "second": "ett sekund",
        "seconds": "{0} sekunder",
        "minute": "ett minutt",
        "minutes": "{0} minutter",
        "hour": "en time",
        "hours": "{0} timer",
        "day": "en dag",
        "days": "{0} dager",
        "week": "en uke",
        "weeks": "{0} uker",
        "month": "en måned",
        "months": "{0} måneder",
        "year": "ett år",
        "years": "{0} år",
    }

    month_names = [
        "",
        "januar",
        "februar",
        "mars",
        "april",
        "mai",
        "juni",
        "juli",
        "august",
        "september",
        "oktober",
        "november",
        "desember",
    ]
    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mar",
        "apr",
        "mai",
        "jun",
        "jul",
        "aug",
        "sep",
        "okt",
        "nov",
        "des",
    ]

    day_names = [
        "",
        "mandag",
        "tirsdag",
        "onsdag",
        "torsdag",
        "fredag",
        "lørdag",
        "søndag",
    ]
    day_abbreviations = ["", "ma", "ti", "on", "to", "fr", "lø", "sø"]

    def _ordinal_number(self, n: int) -> str:
        return f"{n}."


class NewNorwegianLocale(Locale):
    names = ["nn", "nn-no"]

    past = "for {0} sidan"
    future = "om {0}"

    timeframes = {
        "now": "no nettopp",
        "second": "eitt sekund",
        "seconds": "{0} sekund",
        "minute": "eitt minutt",
        "minutes": "{0} minutt",
        "hour": "ein time",
        "hours": "{0} timar",
        "day": "ein dag",
        "days": "{0} dagar",
        "week": "ei veke",
        "weeks": "{0} veker",
        "month": "ein månad",
        "months": "{0} månader",
        "year": "eitt år",
        "years": "{0} år",
    }

    month_names = [
        "",
        "januar",
        "februar",
        "mars",
        "april",
        "mai",
        "juni",
        "juli",
        "august",
        "september",
        "oktober",
        "november",
        "desember",
    ]
    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mar",
        "apr",
        "mai",
        "jun",
        "jul",
        "aug",
        "sep",
        "okt",
        "nov",
        "des",
    ]

    day_names = [
        "",
        "måndag",
        "tysdag",
        "onsdag",
        "torsdag",
        "fredag",
        "laurdag",
        "sundag",
    ]
    day_abbreviations = ["", "må", "ty", "on", "to", "fr", "la", "su"]

    def _ordinal_number(self, n: int) -> str:
        return f"{n}."


class PortugueseLocale(Locale):
    names = ["pt", "pt-pt"]

    past = "há {0}"
    future = "em {0}"
    and_word = "e"

    timeframes = {
        "now": "agora",
        "second": "um segundo",
        "seconds": "{0} segundos",
        "minute": "um minuto",
        "minutes": "{0} minutos",
        "hour": "uma hora",
        "hours": "{0} horas",
        "day": "um dia",
        "days": "{0} dias",
        "week": "uma semana",
        "weeks": "{0} semanas",
        "month": "um mês",
        "months": "{0} meses",
        "year": "um ano",
        "years": "{0} anos",
    }

    month_names = [
        "",
        "Janeiro",
        "Fevereiro",
        "Março",
        "Abril",
        "Maio",
        "Junho",
        "Julho",
        "Agosto",
        "Setembro",
        "Outubro",
        "Novembro",
        "Dezembro",
    ]
    month_abbreviations = [
        "",
        "Jan",
        "Fev",
        "Mar",
        "Abr",
        "Mai",
        "Jun",
        "Jul",
        "Ago",
        "Set",
        "Out",
        "Nov",
        "Dez",
    ]

    day_names = [
        "",
        "Segunda-feira",
        "Terça-feira",
        "Quarta-feira",
        "Quinta-feira",
        "Sexta-feira",
        "Sábado",
        "Domingo",
    ]
    day_abbreviations = ["", "Seg", "Ter", "Qua", "Qui", "Sex", "Sab", "Dom"]


class BrazilianPortugueseLocale(PortugueseLocale):
    names = ["pt-br"]

    past = "faz {0}"


class TagalogLocale(Locale):
    names = ["tl", "tl-ph"]

    past = "nakaraang {0}"
    future = "{0} mula ngayon"

    timeframes = {
        "now": "ngayon lang",
        "second": "isang segundo",
        "seconds": "{0} segundo",
        "minute": "isang minuto",
        "minutes": "{0} minuto",
        "hour": "isang oras",
        "hours": "{0} oras",
        "day": "isang araw",
        "days": "{0} araw",
        "week": "isang linggo",
        "weeks": "{0} linggo",
        "month": "isang buwan",
        "months": "{0} buwan",
        "year": "isang taon",
        "years": "{0} taon",
    }

    month_names = [
        "",
        "Enero",
        "Pebrero",
        "Marso",
        "Abril",
        "Mayo",
        "Hunyo",
        "Hulyo",
        "Agosto",
        "Setyembre",
        "Oktubre",
        "Nobyembre",
        "Disyembre",
    ]
    month_abbreviations = [
        "",
        "Ene",
        "Peb",
        "Mar",
        "Abr",
        "May",
        "Hun",
        "Hul",
        "Ago",
        "Set",
        "Okt",
        "Nob",
        "Dis",
    ]

    day_names = [
        "",
        "Lunes",
        "Martes",
        "Miyerkules",
        "Huwebes",
        "Biyernes",
        "Sabado",
        "Linggo",
    ]
    day_abbreviations = ["", "Lun", "Mar", "Miy", "Huw", "Biy", "Sab", "Lin"]

    meridians = {"am": "nu", "pm": "nh", "AM": "ng umaga", "PM": "ng hapon"}

    def _ordinal_number(self, n: int) -> str:
        return f"ika-{n}"


class VietnameseLocale(Locale):
    names = ["vi", "vi-vn"]

    past = "{0} trước"
    future = "{0} nữa"

    timeframes = {
        "now": "hiện tại",
        "second": "một giây",
        "seconds": "{0} giây",
        "minute": "một phút",
        "minutes": "{0} phút",
        "hour": "một giờ",
        "hours": "{0} giờ",
        "day": "một ngày",
        "days": "{0} ngày",
        "week": "một tuần",
        "weeks": "{0} tuần",
        "month": "một tháng",
        "months": "{0} tháng",
        "year": "một năm",
        "years": "{0} năm",
    }

    month_names = [
        "",
        "Tháng Một",
        "Tháng Hai",
        "Tháng Ba",
        "Tháng Tư",
        "Tháng Năm",
        "Tháng Sáu",
        "Tháng Bảy",
        "Tháng Tám",
        "Tháng Chín",
        "Tháng Mười",
        "Tháng Mười Một",
        "Tháng Mười Hai",
    ]
    month_abbreviations = [
        "",
        "Tháng 1",
        "Tháng 2",
        "Tháng 3",
        "Tháng 4",
        "Tháng 5",
        "Tháng 6",
        "Tháng 7",
        "Tháng 8",
        "Tháng 9",
        "Tháng 10",
        "Tháng 11",
        "Tháng 12",
    ]

    day_names = [
        "",
        "Thứ Hai",
        "Thứ Ba",
        "Thứ Tư",
        "Thứ Năm",
        "Thứ Sáu",
        "Thứ Bảy",
        "Chủ Nhật",
    ]
    day_abbreviations = ["", "Thứ 2", "Thứ 3", "Thứ 4", "Thứ 5", "Thứ 6", "Thứ 7", "CN"]


class TurkishLocale(Locale):
    names = ["tr", "tr-tr"]

    past = "{0} önce"
    future = "{0} sonra"
    and_word = "ve"

    timeframes = {
        "now": "şimdi",
        "second": "bir saniye",
        "seconds": "{0} saniye",
        "minute": "bir dakika",
        "minutes": "{0} dakika",
        "hour": "bir saat",
        "hours": "{0} saat",
        "day": "bir gün",
        "days": "{0} gün",
        "week": "bir hafta",
        "weeks": "{0} hafta",
        "month": "bir ay",
        "months": "{0} ay",
        "year": "bir yıl",
        "years": "{0} yıl",
    }

    meridians = {"am": "öö", "pm": "ös", "AM": "ÖÖ", "PM": "ÖS"}

    month_names = [
        "",
        "Ocak",
        "Şubat",
        "Mart",
        "Nisan",
        "Mayıs",
        "Haziran",
        "Temmuz",
        "Ağustos",
        "Eylül",
        "Ekim",
        "Kasım",
        "Aralık",
    ]
    month_abbreviations = [
        "",
        "Oca",
        "Şub",
        "Mar",
        "Nis",
        "May",
        "Haz",
        "Tem",
        "Ağu",
        "Eyl",
        "Eki",
        "Kas",
        "Ara",
    ]

    day_names = [
        "",
        "Pazartesi",
        "Salı",
        "Çarşamba",
        "Perşembe",
        "Cuma",
        "Cumartesi",
        "Pazar",
    ]
    day_abbreviations = ["", "Pzt", "Sal", "Çar", "Per", "Cum", "Cmt", "Paz"]


class AzerbaijaniLocale(Locale):
    names = ["az", "az-az"]

    past = "{0} əvvəl"
    future = "{0} sonra"

    timeframes = {
        "now": "indi",
        "second": "bir saniyə",
        "seconds": "{0} saniyə",
        "minute": "bir dəqiqə",
        "minutes": "{0} dəqiqə",
        "hour": "bir saat",
        "hours": "{0} saat",
        "day": "bir gün",
        "days": "{0} gün",
        "week": "bir həftə",
        "weeks": "{0} həftə",
        "month": "bir ay",
        "months": "{0} ay",
        "year": "bir il",
        "years": "{0} il",
    }

    month_names = [
        "",
        "Yanvar",
        "Fevral",
        "Mart",
        "Aprel",
        "May",
        "İyun",
        "İyul",
        "Avqust",
        "Sentyabr",
        "Oktyabr",
        "Noyabr",
        "Dekabr",
    ]
    month_abbreviations = [
        "",
        "Yan",
        "Fev",
        "Mar",
        "Apr",
        "May",
        "İyn",
        "İyl",
        "Avq",
        "Sen",
        "Okt",
        "Noy",
        "Dek",
    ]

    day_names = [
        "",
        "Bazar ertəsi",
        "Çərşənbə axşamı",
        "Çərşənbə",
        "Cümə axşamı",
        "Cümə",
        "Şənbə",
        "Bazar",
    ]
    day_abbreviations = ["", "Ber", "Çax", "Çər", "Cax", "Cüm", "Şnb", "Bzr"]


class ArabicLocale(Locale):
    names = [
        "ar",
        "ar-ae",
        "ar-bh",
        "ar-dj",
        "ar-eg",
        "ar-eh",
        "ar-er",
        "ar-km",
        "ar-kw",
        "ar-ly",
        "ar-om",
        "ar-qa",
        "ar-sa",
        "ar-sd",
        "ar-so",
        "ar-ss",
        "ar-td",
        "ar-ye",
    ]

    past = "منذ {0}"
    future = "خلال {0}"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "الآن",
        "second": "ثانية",
        "seconds": {"2": "ثانيتين", "ten": "{0} ثوان", "higher": "{0} ثانية"},
        "minute": "دقيقة",
        "minutes": {"2": "دقيقتين", "ten": "{0} دقائق", "higher": "{0} دقيقة"},
        "hour": "ساعة",
        "hours": {"2": "ساعتين", "ten": "{0} ساعات", "higher": "{0} ساعة"},
        "day": "يوم",
        "days": {"2": "يومين", "ten": "{0} أيام", "higher": "{0} يوم"},
        "week": "اسبوع",
        "weeks": {"2": "اسبوعين", "ten": "{0} أسابيع", "higher": "{0} اسبوع"},
        "month": "شهر",
        "months": {"2": "شهرين", "ten": "{0} أشهر", "higher": "{0} شهر"},
        "year": "سنة",
        "years": {"2": "سنتين", "ten": "{0} سنوات", "higher": "{0} سنة"},
    }

    month_names = [
        "",
        "يناير",
        "فبراير",
        "مارس",
        "أبريل",
        "مايو",
        "يونيو",
        "يوليو",
        "أغسطس",
        "سبتمبر",
        "أكتوبر",
        "نوفمبر",
        "ديسمبر",
    ]
    month_abbreviations = [
        "",
        "يناير",
        "فبراير",
        "مارس",
        "أبريل",
        "مايو",
        "يونيو",
        "يوليو",
        "أغسطس",
        "سبتمبر",
        "أكتوبر",
        "نوفمبر",
        "ديسمبر",
    ]

    day_names = [
        "",
        "الإثنين",
        "الثلاثاء",
        "الأربعاء",
        "الخميس",
        "الجمعة",
        "السبت",
        "الأحد",
    ]
    day_abbreviations = ["", "إثنين", "ثلاثاء", "أربعاء", "خميس", "جمعة", "سبت", "أحد"]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]
        delta = abs(delta)
        if isinstance(form, Mapping):
            if delta == 2:
                form = form["2"]
            elif 2 < delta <= 10:
                form = form["ten"]
            else:
                form = form["higher"]

        return form.format(delta)


class LevantArabicLocale(ArabicLocale):
    names = ["ar-iq", "ar-jo", "ar-lb", "ar-ps", "ar-sy"]
    month_names = [
        "",
        "كانون الثاني",
        "شباط",
        "آذار",
        "نيسان",
        "أيار",
        "حزيران",
        "تموز",
        "آب",
        "أيلول",
        "تشرين الأول",
        "تشرين الثاني",
        "كانون الأول",
    ]
    month_abbreviations = [
        "",
        "كانون الثاني",
        "شباط",
        "آذار",
        "نيسان",
        "أيار",
        "حزيران",
        "تموز",
        "آب",
        "أيلول",
        "تشرين الأول",
        "تشرين الثاني",
        "كانون الأول",
    ]


class AlgeriaTunisiaArabicLocale(ArabicLocale):
    names = ["ar-tn", "ar-dz"]
    month_names = [
        "",
        "جانفي",
        "فيفري",
        "مارس",
        "أفريل",
        "ماي",
        "جوان",
        "جويلية",
        "أوت",
        "سبتمبر",
        "أكتوبر",
        "نوفمبر",
        "ديسمبر",
    ]
    month_abbreviations = [
        "",
        "جانفي",
        "فيفري",
        "مارس",
        "أفريل",
        "ماي",
        "جوان",
        "جويلية",
        "أوت",
        "سبتمبر",
        "أكتوبر",
        "نوفمبر",
        "ديسمبر",
    ]


class MauritaniaArabicLocale(ArabicLocale):
    names = ["ar-mr"]
    month_names = [
        "",
        "يناير",
        "فبراير",
        "مارس",
        "إبريل",
        "مايو",
        "يونيو",
        "يوليو",
        "أغشت",
        "شتمبر",
        "أكتوبر",
        "نوفمبر",
        "دجمبر",
    ]
    month_abbreviations = [
        "",
        "يناير",
        "فبراير",
        "مارس",
        "إبريل",
        "مايو",
        "يونيو",
        "يوليو",
        "أغشت",
        "شتمبر",
        "أكتوبر",
        "نوفمبر",
        "دجمبر",
    ]


class MoroccoArabicLocale(ArabicLocale):
    names = ["ar-ma"]
    month_names = [
        "",
        "يناير",
        "فبراير",
        "مارس",
        "أبريل",
        "ماي",
        "يونيو",
        "يوليوز",
        "غشت",
        "شتنبر",
        "أكتوبر",
        "نونبر",
        "دجنبر",
    ]
    month_abbreviations = [
        "",
        "يناير",
        "فبراير",
        "مارس",
        "أبريل",
        "ماي",
        "يونيو",
        "يوليوز",
        "غشت",
        "شتنبر",
        "أكتوبر",
        "نونبر",
        "دجنبر",
    ]


class IcelandicLocale(Locale):
    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]

        if isinstance(form, Mapping):
            if delta < 0:
                form = form["past"]
            elif delta > 0:
                form = form["future"]
            else:
                raise ValueError(
                    "Icelandic Locale does not support units with a delta of zero. "
                    "Please consider making a contribution to fix this issue."
                )
                # FIXME: handle when delta is 0

        return form.format(abs(delta))

    names = ["is", "is-is"]

    past = "fyrir {0} síðan"
    future = "eftir {0}"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "rétt í þessu",
        "second": {"past": "sekúndu", "future": "sekúndu"},
        "seconds": {"past": "{0} nokkrum sekúndum", "future": "nokkrar sekúndur"},
        "minute": {"past": "einni mínútu", "future": "eina mínútu"},
        "minutes": {"past": "{0} mínútum", "future": "{0} mínútur"},
        "hour": {"past": "einum tíma", "future": "einn tíma"},
        "hours": {"past": "{0} tímum", "future": "{0} tíma"},
        "day": {"past": "einum degi", "future": "einn dag"},
        "days": {"past": "{0} dögum", "future": "{0} daga"},
        "month": {"past": "einum mánuði", "future": "einn mánuð"},
        "months": {"past": "{0} mánuðum", "future": "{0} mánuði"},
        "year": {"past": "einu ári", "future": "eitt ár"},
        "years": {"past": "{0} árum", "future": "{0} ár"},
    }

    meridians = {"am": "f.h.", "pm": "e.h.", "AM": "f.h.", "PM": "e.h."}

    month_names = [
        "",
        "janúar",
        "febrúar",
        "mars",
        "apríl",
        "maí",
        "júní",
        "júlí",
        "ágúst",
        "september",
        "október",
        "nóvember",
        "desember",
    ]
    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mar",
        "apr",
        "maí",
        "jún",
        "júl",
        "ágú",
        "sep",
        "okt",
        "nóv",
        "des",
    ]

    day_names = [
        "",
        "mánudagur",
        "þriðjudagur",
        "miðvikudagur",
        "fimmtudagur",
        "föstudagur",
        "laugardagur",
        "sunnudagur",
    ]
    day_abbreviations = ["", "mán", "þri", "mið", "fim", "fös", "lau", "sun"]


class DanishLocale(Locale):
    names = ["da", "da-dk"]

    past = "for {0} siden"
    future = "om {0}"
    and_word = "og"

    timeframes = {
        "now": "lige nu",
        "second": "et sekund",
        "seconds": "{0} sekunder",
        "minute": "et minut",
        "minutes": "{0} minutter",
        "hour": "en time",
        "hours": "{0} timer",
        "day": "en dag",
        "days": "{0} dage",
        "week": "en uge",
        "weeks": "{0} uger",
        "month": "en måned",
        "months": "{0} måneder",
        "year": "et år",
        "years": "{0} år",
    }

    month_names = [
        "",
        "januar",
        "februar",
        "marts",
        "april",
        "maj",
        "juni",
        "juli",
        "august",
        "september",
        "oktober",
        "november",
        "december",
    ]
    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mar",
        "apr",
        "maj",
        "jun",
        "jul",
        "aug",
        "sep",
        "okt",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "mandag",
        "tirsdag",
        "onsdag",
        "torsdag",
        "fredag",
        "lørdag",
        "søndag",
    ]
    day_abbreviations = ["", "man", "tir", "ons", "tor", "fre", "lør", "søn"]

    def _ordinal_number(self, n: int) -> str:
        return f"{n}."


class MalayalamLocale(Locale):
    names = ["ml"]

    past = "{0} മുമ്പ്"
    future = "{0} ശേഷം"

    timeframes = {
        "now": "ഇപ്പോൾ",
        "second": "ഒരു നിമിഷം",
        "seconds": "{0} സെക്കന്റ്‌",
        "minute": "ഒരു മിനിറ്റ്",
        "minutes": "{0} മിനിറ്റ്",
        "hour": "ഒരു മണിക്കൂർ",
        "hours": "{0} മണിക്കൂർ",
        "day": "ഒരു ദിവസം ",
        "days": "{0} ദിവസം ",
        "month": "ഒരു മാസം ",
        "months": "{0} മാസം ",
        "year": "ഒരു വർഷം ",
        "years": "{0} വർഷം ",
    }

    meridians = {
        "am": "രാവിലെ",
        "pm": "ഉച്ചക്ക് ശേഷം",
        "AM": "രാവിലെ",
        "PM": "ഉച്ചക്ക് ശേഷം",
    }

    month_names = [
        "",
        "ജനുവരി",
        "ഫെബ്രുവരി",
        "മാർച്ച്‌",
        "ഏപ്രിൽ ",
        "മെയ്‌ ",
        "ജൂണ്‍",
        "ജൂലൈ",
        "ഓഗസ്റ്റ്‌",
        "സെപ്റ്റംബർ",
        "ഒക്ടോബർ",
        "നവംബർ",
        "ഡിസംബർ",
    ]
    month_abbreviations = [
        "",
        "ജനു",
        "ഫെബ് ",
        "മാർ",
        "ഏപ്രിൽ",
        "മേയ്",
        "ജൂണ്‍",
        "ജൂലൈ",
        "ഓഗസ്റ",
        "സെപ്റ്റ",
        "ഒക്ടോ",
        "നവം",
        "ഡിസം",
    ]

    day_names = ["", "തിങ്കള്‍", "ചൊവ്വ", "ബുധന്‍", "വ്യാഴം", "വെള്ളി", "ശനി", "ഞായര്‍"]
    day_abbreviations = [
        "",
        "തിങ്കള്‍",
        "ചൊവ്വ",
        "ബുധന്‍",
        "വ്യാഴം",
        "വെള്ളി",
        "ശനി",
        "ഞായര്‍",
    ]


class HindiLocale(Locale):
    names = ["hi", "hi-in"]

    past = "{0} पहले"
    future = "{0} बाद"

    timeframes = {
        "now": "अभी",
        "second": "एक पल",
        "seconds": "{0} सेकंड्",
        "minute": "एक मिनट ",
        "minutes": "{0} मिनट ",
        "hour": "एक घंटा",
        "hours": "{0} घंटे",
        "day": "एक दिन",
        "days": "{0} दिन",
        "month": "एक माह ",
        "months": "{0} महीने ",
        "year": "एक वर्ष ",
        "years": "{0} साल ",
    }

    meridians = {"am": "सुबह", "pm": "शाम", "AM": "सुबह", "PM": "शाम"}

    month_names = [
        "",
        "जनवरी",
        "फरवरी",
        "मार्च",
        "अप्रैल ",
        "मई",
        "जून",
        "जुलाई",
        "अगस्त",
        "सितंबर",
        "अक्टूबर",
        "नवंबर",
        "दिसंबर",
    ]
    month_abbreviations = [
        "",
        "जन",
        "फ़र",
        "मार्च",
        "अप्रै",
        "मई",
        "जून",
        "जुलाई",
        "आग",
        "सित",
        "अकत",
        "नवे",
        "दिस",
    ]

    day_names = [
        "",
        "सोमवार",
        "मंगलवार",
        "बुधवार",
        "गुरुवार",
        "शुक्रवार",
        "शनिवार",
        "रविवार",
    ]
    day_abbreviations = ["", "सोम", "मंगल", "बुध", "गुरुवार", "शुक्र", "शनि", "रवि"]


class CzechLocale(Locale):
    names = ["cs", "cs-cz"]

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "Teď",
        "second": {"past": "vteřina", "future": "vteřina"},
        "seconds": {
            "zero": "vteřina",
            "past": "{0} sekundami",
            "future-singular": "{0} sekundy",
            "future-paucal": "{0} sekund",
        },
        "minute": {"past": "minutou", "future": "minutu"},
        "minutes": {
            "zero": "{0} minut",
            "past": "{0} minutami",
            "future-singular": "{0} minuty",
            "future-paucal": "{0} minut",
        },
        "hour": {"past": "hodinou", "future": "hodinu"},
        "hours": {
            "zero": "{0} hodin",
            "past": "{0} hodinami",
            "future-singular": "{0} hodiny",
            "future-paucal": "{0} hodin",
        },
        "day": {"past": "dnem", "future": "den"},
        "days": {
            "zero": "{0} dnů",
            "past": "{0} dny",
            "future-singular": "{0} dny",
            "future-paucal": "{0} dnů",
        },
        "week": {"past": "týdnem", "future": "týden"},
        "weeks": {
            "zero": "{0} týdnů",
            "past": "{0} týdny",
            "future-singular": "{0} týdny",
            "future-paucal": "{0} týdnů",
        },
        "month": {"past": "měsícem", "future": "měsíc"},
        "months": {
            "zero": "{0} měsíců",
            "past": "{0} měsíci",
            "future-singular": "{0} měsíce",
            "future-paucal": "{0} měsíců",
        },
        "year": {"past": "rokem", "future": "rok"},
        "years": {
            "zero": "{0} let",
            "past": "{0} lety",
            "future-singular": "{0} roky",
            "future-paucal": "{0} let",
        },
    }

    past = "Před {0}"
    future = "Za {0}"

    month_names = [
        "",
        "leden",
        "únor",
        "březen",
        "duben",
        "květen",
        "červen",
        "červenec",
        "srpen",
        "září",
        "říjen",
        "listopad",
        "prosinec",
    ]
    month_abbreviations = [
        "",
        "led",
        "úno",
        "bře",
        "dub",
        "kvě",
        "čvn",
        "čvc",
        "srp",
        "zář",
        "říj",
        "lis",
        "pro",
    ]

    day_names = [
        "",
        "pondělí",
        "úterý",
        "středa",
        "čtvrtek",
        "pátek",
        "sobota",
        "neděle",
    ]
    day_abbreviations = ["", "po", "út", "st", "čt", "pá", "so", "ne"]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        """Czech aware time frame format function, takes into account
        the differences between past and future forms."""
        abs_delta = abs(delta)
        form = self.timeframes[timeframe]

        if isinstance(form, str):
            return form.format(abs_delta)

        if delta == 0:
            key = "zero"  # And *never* use 0 in the singular!
        elif delta < 0:
            key = "past"
        else:
            # Needed since both regular future and future-singular and future-paucal cases
            if "future-singular" not in form:
                key = "future"
            elif 2 <= abs_delta % 10 <= 4 and (
                abs_delta % 100 < 10 or abs_delta % 100 >= 20
            ):
                key = "future-singular"
            else:
                key = "future-paucal"

        form: str = form[key]
        return form.format(abs_delta)


class SlovakLocale(Locale):
    names = ["sk", "sk-sk"]

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "Teraz",
        "second": {"past": "sekundou", "future": "sekundu"},
        "seconds": {
            "zero": "{0} sekúnd",
            "past": "{0} sekundami",
            "future-singular": "{0} sekundy",
            "future-paucal": "{0} sekúnd",
        },
        "minute": {"past": "minútou", "future": "minútu"},
        "minutes": {
            "zero": "{0} minút",
            "past": "{0} minútami",
            "future-singular": "{0} minúty",
            "future-paucal": "{0} minút",
        },
        "hour": {"past": "hodinou", "future": "hodinu"},
        "hours": {
            "zero": "{0} hodín",
            "past": "{0} hodinami",
            "future-singular": "{0} hodiny",
            "future-paucal": "{0} hodín",
        },
        "day": {"past": "dňom", "future": "deň"},
        "days": {
            "zero": "{0} dní",
            "past": "{0} dňami",
            "future-singular": "{0} dni",
            "future-paucal": "{0} dní",
        },
        "week": {"past": "týždňom", "future": "týždeň"},
        "weeks": {
            "zero": "{0} týždňov",
            "past": "{0} týždňami",
            "future-singular": "{0} týždne",
            "future-paucal": "{0} týždňov",
        },
        "month": {"past": "mesiacom", "future": "mesiac"},
        "months": {
            "zero": "{0} mesiacov",
            "past": "{0} mesiacmi",
            "future-singular": "{0} mesiace",
            "future-paucal": "{0} mesiacov",
        },
        "year": {"past": "rokom", "future": "rok"},
        "years": {
            "zero": "{0} rokov",
            "past": "{0} rokmi",
            "future-singular": "{0} roky",
            "future-paucal": "{0} rokov",
        },
    }

    past = "Pred {0}"
    future = "O {0}"
    and_word = "a"

    month_names = [
        "",
        "január",
        "február",
        "marec",
        "apríl",
        "máj",
        "jún",
        "júl",
        "august",
        "september",
        "október",
        "november",
        "december",
    ]
    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mar",
        "apr",
        "máj",
        "jún",
        "júl",
        "aug",
        "sep",
        "okt",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "pondelok",
        "utorok",
        "streda",
        "štvrtok",
        "piatok",
        "sobota",
        "nedeľa",
    ]
    day_abbreviations = ["", "po", "ut", "st", "št", "pi", "so", "ne"]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        """Slovak aware time frame format function, takes into account
        the differences between past and future forms."""
        abs_delta = abs(delta)
        form = self.timeframes[timeframe]

        if isinstance(form, str):
            return form.format(abs_delta)

        if delta == 0:
            key = "zero"  # And *never* use 0 in the singular!
        elif delta < 0:
            key = "past"
        else:
            if "future-singular" not in form:
                key = "future"
            elif 2 <= abs_delta % 10 <= 4 and (
                abs_delta % 100 < 10 or abs_delta % 100 >= 20
            ):
                key = "future-singular"
            else:
                key = "future-paucal"

        form: str = form[key]
        return form.format(abs_delta)


class FarsiLocale(Locale):
    names = ["fa", "fa-ir"]

    past = "{0} قبل"
    future = "در {0}"

    timeframes = {
        "now": "اکنون",
        "second": "یک لحظه",
        "seconds": "{0} ثانیه",
        "minute": "یک دقیقه",
        "minutes": "{0} دقیقه",
        "hour": "یک ساعت",
        "hours": "{0} ساعت",
        "day": "یک روز",
        "days": "{0} روز",
        "month": "یک ماه",
        "months": "{0} ماه",
        "year": "یک سال",
        "years": "{0} سال",
    }

    meridians = {
        "am": "قبل از ظهر",
        "pm": "بعد از ظهر",
        "AM": "قبل از ظهر",
        "PM": "بعد از ظهر",
    }

    month_names = [
        "",
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ]
    month_abbreviations = [
        "",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]

    day_names = [
        "",
        "دو شنبه",
        "سه شنبه",
        "چهارشنبه",
        "پنجشنبه",
        "جمعه",
        "شنبه",
        "یکشنبه",
    ]
    day_abbreviations = ["", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]


class HebrewLocale(Locale):
    names = ["he", "he-il"]

    past = "לפני {0}"
    future = "בעוד {0}"
    and_word = "ו"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "הרגע",
        "second": "שנייה",
        "seconds": "{0} שניות",
        "minute": "דקה",
        "minutes": "{0} דקות",
        "hour": "שעה",
        "hours": {"2": "שעתיים", "ten": "{0} שעות", "higher": "{0} שעות"},
        "day": "יום",
        "days": {"2": "יומיים", "ten": "{0} ימים", "higher": "{0} יום"},
        "week": "שבוע",
        "weeks": {"2": "שבועיים", "ten": "{0} שבועות", "higher": "{0} שבועות"},
        "month": "חודש",
        "months": {"2": "חודשיים", "ten": "{0} חודשים", "higher": "{0} חודשים"},
        "year": "שנה",
        "years": {"2": "שנתיים", "ten": "{0} שנים", "higher": "{0} שנה"},
    }

    meridians = {
        "am": 'לפנ"צ',
        "pm": 'אחר"צ',
        "AM": "לפני הצהריים",
        "PM": "אחרי הצהריים",
    }

    month_names = [
        "",
        "ינואר",
        "פברואר",
        "מרץ",
        "אפריל",
        "מאי",
        "יוני",
        "יולי",
        "אוגוסט",
        "ספטמבר",
        "אוקטובר",
        "נובמבר",
        "דצמבר",
    ]
    month_abbreviations = [
        "",
        "ינו׳",
        "פבר׳",
        "מרץ",
        "אפר׳",
        "מאי",
        "יוני",
        "יולי",
        "אוג׳",
        "ספט׳",
        "אוק׳",
        "נוב׳",
        "דצמ׳",
    ]

    day_names = ["", "שני", "שלישי", "רביעי", "חמישי", "שישי", "שבת", "ראשון"]
    day_abbreviations = ["", "ב׳", "ג׳", "ד׳", "ה׳", "ו׳", "ש׳", "א׳"]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]
        delta = abs(delta)
        if isinstance(form, Mapping):
            if delta == 2:
                form = form["2"]
            elif delta == 0 or 2 < delta <= 10:
                form = form["ten"]
            else:
                form = form["higher"]

        return form.format(delta)

    def describe_multi(
        self,
        timeframes: Sequence[Tuple[TimeFrameLiteral, Union[int, float]]],
        only_distance: bool = False,
    ) -> str:
        """Describes a delta within multiple timeframes in plain language.
        In Hebrew, the and word behaves a bit differently.

        :param timeframes: a list of string, quantity pairs each representing a timeframe and delta.
        :param only_distance: return only distance eg: "2 hours and 11 seconds" without "in" or "ago" keywords
        """

        humanized = ""
        for index, (timeframe, delta) in enumerate(timeframes):
            last_humanized = self._format_timeframe(timeframe, trunc(delta))
            if index == 0:
                humanized = last_humanized
            elif index == len(timeframes) - 1:  # Must have at least 2 items
                humanized += " " + self.and_word
                if last_humanized[0].isdecimal():
                    humanized += "־"
                humanized += last_humanized
            else:  # Don't add for the last one
                humanized += ", " + last_humanized

        if not only_distance:
            humanized = self._format_relative(humanized, timeframe, trunc(delta))

        return humanized


class MarathiLocale(Locale):
    names = ["mr"]

    past = "{0} आधी"
    future = "{0} नंतर"

    timeframes = {
        "now": "सद्य",
        "second": "एक सेकंद",
        "seconds": "{0} सेकंद",
        "minute": "एक मिनिट ",
        "minutes": "{0} मिनिट ",
        "hour": "एक तास",
        "hours": "{0} तास",
        "day": "एक दिवस",
        "days": "{0} दिवस",
        "month": "एक महिना ",
        "months": "{0} महिने ",
        "year": "एक वर्ष ",
        "years": "{0} वर्ष ",
    }

    meridians = {"am": "सकाळ", "pm": "संध्याकाळ", "AM": "सकाळ", "PM": "संध्याकाळ"}

    month_names = [
        "",
        "जानेवारी",
        "फेब्रुवारी",
        "मार्च",
        "एप्रिल",
        "मे",
        "जून",
        "जुलै",
        "अॉगस्ट",
        "सप्टेंबर",
        "अॉक्टोबर",
        "नोव्हेंबर",
        "डिसेंबर",
    ]
    month_abbreviations = [
        "",
        "जान",
        "फेब्रु",
        "मार्च",
        "एप्रि",
        "मे",
        "जून",
        "जुलै",
        "अॉग",
        "सप्टें",
        "अॉक्टो",
        "नोव्हें",
        "डिसें",
    ]

    day_names = [
        "",
        "सोमवार",
        "मंगळवार",
        "बुधवार",
        "गुरुवार",
        "शुक्रवार",
        "शनिवार",
        "रविवार",
    ]
    day_abbreviations = ["", "सोम", "मंगळ", "बुध", "गुरु", "शुक्र", "शनि", "रवि"]


class CatalanLocale(Locale):
    names = ["ca", "ca-es", "ca-ad", "ca-fr", "ca-it"]
    past = "Fa {0}"
    future = "En {0}"
    and_word = "i"

    timeframes = {
        "now": "Ara mateix",
        "second": "un segon",
        "seconds": "{0} segons",
        "minute": "un minut",
        "minutes": "{0} minuts",
        "hour": "una hora",
        "hours": "{0} hores",
        "day": "un dia",
        "days": "{0} dies",
        "month": "un mes",
        "months": "{0} mesos",
        "year": "un any",
        "years": "{0} anys",
    }

    month_names = [
        "",
        "gener",
        "febrer",
        "març",
        "abril",
        "maig",
        "juny",
        "juliol",
        "agost",
        "setembre",
        "octubre",
        "novembre",
        "desembre",
    ]
    month_abbreviations = [
        "",
        "gen.",
        "febr.",
        "març",
        "abr.",
        "maig",
        "juny",
        "jul.",
        "ag.",
        "set.",
        "oct.",
        "nov.",
        "des.",
    ]
    day_names = [
        "",
        "dilluns",
        "dimarts",
        "dimecres",
        "dijous",
        "divendres",
        "dissabte",
        "diumenge",
    ]
    day_abbreviations = [
        "",
        "dl.",
        "dt.",
        "dc.",
        "dj.",
        "dv.",
        "ds.",
        "dg.",
    ]


class BasqueLocale(Locale):
    names = ["eu", "eu-eu"]
    past = "duela {0}"
    future = "{0}"  # I don't know what's the right phrase in Basque for the future.

    timeframes = {
        "now": "Orain",
        "second": "segundo bat",
        "seconds": "{0} segundu",
        "minute": "minutu bat",
        "minutes": "{0} minutu",
        "hour": "ordu bat",
        "hours": "{0} ordu",
        "day": "egun bat",
        "days": "{0} egun",
        "month": "hilabete bat",
        "months": "{0} hilabet",
        "year": "urte bat",
        "years": "{0} urte",
    }

    month_names = [
        "",
        "urtarrilak",
        "otsailak",
        "martxoak",
        "apirilak",
        "maiatzak",
        "ekainak",
        "uztailak",
        "abuztuak",
        "irailak",
        "urriak",
        "azaroak",
        "abenduak",
    ]
    month_abbreviations = [
        "",
        "urt",
        "ots",
        "mar",
        "api",
        "mai",
        "eka",
        "uzt",
        "abu",
        "ira",
        "urr",
        "aza",
        "abe",
    ]
    day_names = [
        "",
        "astelehena",
        "asteartea",
        "asteazkena",
        "osteguna",
        "ostirala",
        "larunbata",
        "igandea",
    ]
    day_abbreviations = ["", "al", "ar", "az", "og", "ol", "lr", "ig"]


class HungarianLocale(Locale):
    names = ["hu", "hu-hu"]

    past = "{0} ezelőtt"
    future = "{0} múlva"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "éppen most",
        "second": {"past": "egy második", "future": "egy második"},
        "seconds": {"past": "{0} másodpercekkel", "future": "{0} pár másodperc"},
        "minute": {"past": "egy perccel", "future": "egy perc"},
        "minutes": {"past": "{0} perccel", "future": "{0} perc"},
        "hour": {"past": "egy órával", "future": "egy óra"},
        "hours": {"past": "{0} órával", "future": "{0} óra"},
        "day": {"past": "egy nappal", "future": "egy nap"},
        "days": {"past": "{0} nappal", "future": "{0} nap"},
        "week": {"past": "egy héttel", "future": "egy hét"},
        "weeks": {"past": "{0} héttel", "future": "{0} hét"},
        "month": {"past": "egy hónappal", "future": "egy hónap"},
        "months": {"past": "{0} hónappal", "future": "{0} hónap"},
        "year": {"past": "egy évvel", "future": "egy év"},
        "years": {"past": "{0} évvel", "future": "{0} év"},
    }

    month_names = [
        "",
        "január",
        "február",
        "március",
        "április",
        "május",
        "június",
        "július",
        "augusztus",
        "szeptember",
        "október",
        "november",
        "december",
    ]
    month_abbreviations = [
        "",
        "jan",
        "febr",
        "márc",
        "ápr",
        "máj",
        "jún",
        "júl",
        "aug",
        "szept",
        "okt",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "hétfő",
        "kedd",
        "szerda",
        "csütörtök",
        "péntek",
        "szombat",
        "vasárnap",
    ]
    day_abbreviations = ["", "hét", "kedd", "szer", "csüt", "pént", "szom", "vas"]

    meridians = {"am": "de", "pm": "du", "AM": "DE", "PM": "DU"}

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]

        if isinstance(form, Mapping):
            if delta > 0:
                form = form["future"]
            else:
                form = form["past"]

        return form.format(abs(delta))


class EsperantoLocale(Locale):
    names = ["eo", "eo-xx"]
    past = "antaŭ {0}"
    future = "post {0}"

    timeframes = {
        "now": "nun",
        "second": "sekundo",
        "seconds": "{0} kelkaj sekundoj",
        "minute": "unu minuto",
        "minutes": "{0} minutoj",
        "hour": "un horo",
        "hours": "{0} horoj",
        "day": "unu tago",
        "days": "{0} tagoj",
        "month": "unu monato",
        "months": "{0} monatoj",
        "year": "unu jaro",
        "years": "{0} jaroj",
    }

    month_names = [
        "",
        "januaro",
        "februaro",
        "marto",
        "aprilo",
        "majo",
        "junio",
        "julio",
        "aŭgusto",
        "septembro",
        "oktobro",
        "novembro",
        "decembro",
    ]
    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mar",
        "apr",
        "maj",
        "jun",
        "jul",
        "aŭg",
        "sep",
        "okt",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "lundo",
        "mardo",
        "merkredo",
        "ĵaŭdo",
        "vendredo",
        "sabato",
        "dimanĉo",
    ]
    day_abbreviations = ["", "lun", "mar", "mer", "ĵaŭ", "ven", "sab", "dim"]

    meridians = {"am": "atm", "pm": "ptm", "AM": "ATM", "PM": "PTM"}

    ordinal_day_re = r"((?P<value>[1-3]?[0-9](?=a))a)"

    def _ordinal_number(self, n: int) -> str:
        return f"{n}a"


class ThaiLocale(Locale):
    names = ["th", "th-th"]

    past = "{0} ที่ผ่านมา"
    future = "ในอีก {0}"

    timeframes = {
        "now": "ขณะนี้",
        "second": "วินาที",
        "seconds": "{0} ไม่กี่วินาที",
        "minute": "1 นาที",
        "minutes": "{0} นาที",
        "hour": "1 ชั่วโมง",
        "hours": "{0} ชั่วโมง",
        "day": "1 วัน",
        "days": "{0} วัน",
        "month": "1 เดือน",
        "months": "{0} เดือน",
        "year": "1 ปี",
        "years": "{0} ปี",
    }

    month_names = [
        "",
        "มกราคม",
        "กุมภาพันธ์",
        "มีนาคม",
        "เมษายน",
        "พฤษภาคม",
        "มิถุนายน",
        "กรกฎาคม",
        "สิงหาคม",
        "กันยายน",
        "ตุลาคม",
        "พฤศจิกายน",
        "ธันวาคม",
    ]
    month_abbreviations = [
        "",
        "ม.ค.",
        "ก.พ.",
        "มี.ค.",
        "เม.ย.",
        "พ.ค.",
        "มิ.ย.",
        "ก.ค.",
        "ส.ค.",
        "ก.ย.",
        "ต.ค.",
        "พ.ย.",
        "ธ.ค.",
    ]

    day_names = ["", "จันทร์", "อังคาร", "พุธ", "พฤหัสบดี", "ศุกร์", "เสาร์", "อาทิตย์"]
    day_abbreviations = ["", "จ", "อ", "พ", "พฤ", "ศ", "ส", "อา"]

    meridians = {"am": "am", "pm": "pm", "AM": "AM", "PM": "PM"}

    BE_OFFSET = 543

    def year_full(self, year: int) -> str:
        """Thai always use Buddhist Era (BE) which is CE + 543"""
        year += self.BE_OFFSET
        return f"{year:04d}"

    def year_abbreviation(self, year: int) -> str:
        """Thai always use Buddhist Era (BE) which is CE + 543"""
        year += self.BE_OFFSET
        return f"{year:04d}"[2:]

    def _format_relative(
        self,
        humanized: str,
        timeframe: TimeFrameLiteral,
        delta: Union[float, int],
    ) -> str:
        """Thai normally doesn't have any space between words"""
        if timeframe == "now":
            return humanized

        direction = self.past if delta < 0 else self.future
        relative_string = direction.format(humanized)

        if timeframe == "seconds":
            relative_string = relative_string.replace(" ", "")

        return relative_string


class LaotianLocale(Locale):
    names = ["lo", "lo-la"]

    past = "{0} ກ່ອນຫນ້ານີ້"
    future = "ໃນ {0}"

    timeframes = {
        "now": "ດຽວນີ້",
        "second": "ວິນາທີ",
        "seconds": "{0} ວິນາທີ",
        "minute": "ນາທີ",
        "minutes": "{0} ນາທີ",
        "hour": "ຊົ່ວໂມງ",
        "hours": "{0} ຊົ່ວໂມງ",
        "day": "ມື້",
        "days": "{0} ມື້",
        "week": "ອາທິດ",
        "weeks": "{0} ອາທິດ",
        "month": "ເດືອນ",
        "months": "{0} ເດືອນ",
        "year": "ປີ",
        "years": "{0} ປີ",
    }

    month_names = [
        "",
        "ມັງກອນ",  # mangkon
        "ກຸມພາ",  # kumpha
        "ມີນາ",  # mina
        "ເມສາ",  # mesa
        "ພຶດສະພາ",  # phudsapha
        "ມິຖຸນາ",  # mithuna
        "ກໍລະກົດ",  # kolakod
        "ສິງຫາ",  # singha
        "ກັນຍາ",  # knaia
        "ຕຸລາ",  # tula
        "ພະຈິກ",  # phachik
        "ທັນວາ",  # thanuaa
    ]
    month_abbreviations = [
        "",
        "ມັງກອນ",
        "ກຸມພາ",
        "ມີນາ",
        "ເມສາ",
        "ພຶດສະພາ",
        "ມິຖຸນາ",
        "ກໍລະກົດ",
        "ສິງຫາ",
        "ກັນຍາ",
        "ຕຸລາ",
        "ພະຈິກ",
        "ທັນວາ",
    ]

    day_names = [
        "",
        "ວັນຈັນ",  # vanchan
        "ວັນອັງຄານ",  # vnoangkhan
        "ວັນພຸດ",  # vanphud
        "ວັນພະຫັດ",  # vanphahad
        "ວັນ​ສຸກ",  # vansuk
        "ວັນເສົາ",  # vansao
        "ວັນອາທິດ",  # vnoathid
    ]
    day_abbreviations = [
        "",
        "ວັນຈັນ",
        "ວັນອັງຄານ",
        "ວັນພຸດ",
        "ວັນພະຫັດ",
        "ວັນ​ສຸກ",
        "ວັນເສົາ",
        "ວັນອາທິດ",
    ]

    BE_OFFSET = 543

    def year_full(self, year: int) -> str:
        """Lao always use Buddhist Era (BE) which is CE + 543"""
        year += self.BE_OFFSET
        return f"{year:04d}"

    def year_abbreviation(self, year: int) -> str:
        """Lao always use Buddhist Era (BE) which is CE + 543"""
        year += self.BE_OFFSET
        return f"{year:04d}"[2:]

    def _format_relative(
        self,
        humanized: str,
        timeframe: TimeFrameLiteral,
        delta: Union[float, int],
    ) -> str:
        """Lao normally doesn't have any space between words"""
        if timeframe == "now":
            return humanized

        direction = self.past if delta < 0 else self.future
        relative_string = direction.format(humanized)

        if timeframe == "seconds":
            relative_string = relative_string.replace(" ", "")

        return relative_string


class BengaliLocale(Locale):
    names = ["bn", "bn-bd", "bn-in"]

    past = "{0} আগে"
    future = "{0} পরে"

    timeframes = {
        "now": "এখন",
        "second": "একটি দ্বিতীয়",
        "seconds": "{0} সেকেন্ড",
        "minute": "এক মিনিট",
        "minutes": "{0} মিনিট",
        "hour": "এক ঘণ্টা",
        "hours": "{0} ঘণ্টা",
        "day": "এক দিন",
        "days": "{0} দিন",
        "month": "এক মাস",
        "months": "{0} মাস ",
        "year": "এক বছর",
        "years": "{0} বছর",
    }

    meridians = {"am": "সকাল", "pm": "বিকাল", "AM": "সকাল", "PM": "বিকাল"}

    month_names = [
        "",
        "জানুয়ারি",
        "ফেব্রুয়ারি",
        "মার্চ",
        "এপ্রিল",
        "মে",
        "জুন",
        "জুলাই",
        "আগস্ট",
        "সেপ্টেম্বর",
        "অক্টোবর",
        "নভেম্বর",
        "ডিসেম্বর",
    ]
    month_abbreviations = [
        "",
        "জানু",
        "ফেব",
        "মার্চ",
        "এপ্রি",
        "মে",
        "জুন",
        "জুল",
        "অগা",
        "সেপ্ট",
        "অক্টো",
        "নভে",
        "ডিসে",
    ]

    day_names = [
        "",
        "সোমবার",
        "মঙ্গলবার",
        "বুধবার",
        "বৃহস্পতিবার",
        "শুক্রবার",
        "শনিবার",
        "রবিবার",
    ]
    day_abbreviations = ["", "সোম", "মঙ্গল", "বুধ", "বৃহঃ", "শুক্র", "শনি", "রবি"]

    def _ordinal_number(self, n: int) -> str:
        if n > 10 or n == 0:
            return f"{n}তম"
        if n in [1, 5, 7, 8, 9, 10]:
            return f"{n}ম"
        if n in [2, 3]:
            return f"{n}য়"
        if n == 4:
            return f"{n}র্থ"
        if n == 6:
            return f"{n}ষ্ঠ"
        return ""


class RomanshLocale(Locale):
    names = ["rm", "rm-ch"]

    past = "avant {0}"
    future = "en {0}"

    timeframes = {
        "now": "en quest mument",
        "second": "in secunda",
        "seconds": "{0} secundas",
        "minute": "ina minuta",
        "minutes": "{0} minutas",
        "hour": "in'ura",
        "hours": "{0} ura",
        "day": "in di",
        "days": "{0} dis",
        "week": "in'emna",
        "weeks": "{0} emnas",
        "month": "in mais",
        "months": "{0} mais",
        "year": "in onn",
        "years": "{0} onns",
    }

    month_names = [
        "",
        "schaner",
        "favrer",
        "mars",
        "avrigl",
        "matg",
        "zercladur",
        "fanadur",
        "avust",
        "settember",
        "october",
        "november",
        "december",
    ]

    month_abbreviations = [
        "",
        "schan",
        "fav",
        "mars",
        "avr",
        "matg",
        "zer",
        "fan",
        "avu",
        "set",
        "oct",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "glindesdi",
        "mardi",
        "mesemna",
        "gievgia",
        "venderdi",
        "sonda",
        "dumengia",
    ]

    day_abbreviations = ["", "gli", "ma", "me", "gie", "ve", "so", "du"]


class RomanianLocale(Locale):
    names = ["ro", "ro-ro"]

    past = "{0} în urmă"
    future = "peste {0}"
    and_word = "și"

    timeframes = {
        "now": "acum",
        "second": "o secunda",
        "seconds": "{0} câteva secunde",
        "minute": "un minut",
        "minutes": "{0} minute",
        "hour": "o oră",
        "hours": "{0} ore",
        "day": "o zi",
        "days": "{0} zile",
        "month": "o lună",
        "months": "{0} luni",
        "year": "un an",
        "years": "{0} ani",
    }

    month_names = [
        "",
        "ianuarie",
        "februarie",
        "martie",
        "aprilie",
        "mai",
        "iunie",
        "iulie",
        "august",
        "septembrie",
        "octombrie",
        "noiembrie",
        "decembrie",
    ]
    month_abbreviations = [
        "",
        "ian",
        "febr",
        "mart",
        "apr",
        "mai",
        "iun",
        "iul",
        "aug",
        "sept",
        "oct",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "luni",
        "marți",
        "miercuri",
        "joi",
        "vineri",
        "sâmbătă",
        "duminică",
    ]
    day_abbreviations = ["", "Lun", "Mar", "Mie", "Joi", "Vin", "Sâm", "Dum"]


class SlovenianLocale(Locale):
    names = ["sl", "sl-si"]

    past = "pred {0}"
    future = "čez {0}"
    and_word = "in"

    timeframes = {
        "now": "zdaj",
        "second": "sekundo",
        "seconds": "{0} sekund",
        "minute": "minuta",
        "minutes": "{0} minutami",
        "hour": "uro",
        "hours": "{0} ur",
        "day": "dan",
        "days": "{0} dni",
        "month": "mesec",
        "months": "{0} mesecev",
        "year": "leto",
        "years": "{0} let",
    }

    meridians = {"am": "", "pm": "", "AM": "", "PM": ""}

    month_names = [
        "",
        "Januar",
        "Februar",
        "Marec",
        "April",
        "Maj",
        "Junij",
        "Julij",
        "Avgust",
        "September",
        "Oktober",
        "November",
        "December",
    ]

    month_abbreviations = [
        "",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "Maj",
        "Jun",
        "Jul",
        "Avg",
        "Sep",
        "Okt",
        "Nov",
        "Dec",
    ]

    day_names = [
        "",
        "Ponedeljek",
        "Torek",
        "Sreda",
        "Četrtek",
        "Petek",
        "Sobota",
        "Nedelja",
    ]

    day_abbreviations = ["", "Pon", "Tor", "Sre", "Čet", "Pet", "Sob", "Ned"]


class IndonesianLocale(Locale):
    names = ["id", "id-id"]

    past = "{0} yang lalu"
    future = "dalam {0}"
    and_word = "dan"

    timeframes = {
        "now": "baru saja",
        "second": "1 sebentar",
        "seconds": "{0} detik",
        "minute": "1 menit",
        "minutes": "{0} menit",
        "hour": "1 jam",
        "hours": "{0} jam",
        "day": "1 hari",
        "days": "{0} hari",
        "week": "1 minggu",
        "weeks": "{0} minggu",
        "month": "1 bulan",
        "months": "{0} bulan",
        "quarter": "1 kuartal",
        "quarters": "{0} kuartal",
        "year": "1 tahun",
        "years": "{0} tahun",
    }

    meridians = {"am": "", "pm": "", "AM": "", "PM": ""}

    month_names = [
        "",
        "Januari",
        "Februari",
        "Maret",
        "April",
        "Mei",
        "Juni",
        "Juli",
        "Agustus",
        "September",
        "Oktober",
        "November",
        "Desember",
    ]

    month_abbreviations = [
        "",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "Mei",
        "Jun",
        "Jul",
        "Ags",
        "Sept",
        "Okt",
        "Nov",
        "Des",
    ]

    day_names = ["", "Senin", "Selasa", "Rabu", "Kamis", "Jumat", "Sabtu", "Minggu"]

    day_abbreviations = [
        "",
        "Senin",
        "Selasa",
        "Rabu",
        "Kamis",
        "Jumat",
        "Sabtu",
        "Minggu",
    ]


class NepaliLocale(Locale):
    names = ["ne", "ne-np"]

    past = "{0} पहिले"
    future = "{0} पछी"

    timeframes = {
        "now": "अहिले",
        "second": "एक सेकेन्ड",
        "seconds": "{0} सेकण्ड",
        "minute": "मिनेट",
        "minutes": "{0} मिनेट",
        "hour": "एक घण्टा",
        "hours": "{0} घण्टा",
        "day": "एक दिन",
        "days": "{0} दिन",
        "month": "एक महिना",
        "months": "{0} महिना",
        "year": "एक बर्ष",
        "years": "{0} बर्ष",
    }

    meridians = {"am": "पूर्वाह्न", "pm": "अपरान्ह", "AM": "पूर्वाह्न", "PM": "अपरान्ह"}

    month_names = [
        "",
        "जनवरी",
        "फेब्रुअरी",
        "मार्च",
        "एप्रील",
        "मे",
        "जुन",
        "जुलाई",
        "अगष्ट",
        "सेप्टेम्बर",
        "अक्टोबर",
        "नोवेम्बर",
        "डिसेम्बर",
    ]
    month_abbreviations = [
        "",
        "जन",
        "फेब",
        "मार्च",
        "एप्रील",
        "मे",
        "जुन",
        "जुलाई",
        "अग",
        "सेप",
        "अक्ट",
        "नोव",
        "डिस",
    ]

    day_names = [
        "",
        "सोमवार",
        "मंगलवार",
        "बुधवार",
        "बिहिवार",
        "शुक्रवार",
        "शनिवार",
        "आइतवार",
    ]

    day_abbreviations = ["", "सोम", "मंगल", "बुध", "बिहि", "शुक्र", "शनि", "आइत"]


class EstonianLocale(Locale):
    names = ["ee", "et"]

    past = "{0} tagasi"
    future = "{0} pärast"
    and_word = "ja"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Mapping[str, str]]] = {
        "now": {"past": "just nüüd", "future": "just nüüd"},
        "second": {"past": "üks sekund", "future": "ühe sekundi"},
        "seconds": {"past": "{0} sekundit", "future": "{0} sekundi"},
        "minute": {"past": "üks minut", "future": "ühe minuti"},
        "minutes": {"past": "{0} minutit", "future": "{0} minuti"},
        "hour": {"past": "tund aega", "future": "tunni aja"},
        "hours": {"past": "{0} tundi", "future": "{0} tunni"},
        "day": {"past": "üks päev", "future": "ühe päeva"},
        "days": {"past": "{0} päeva", "future": "{0} päeva"},
        "month": {"past": "üks kuu", "future": "ühe kuu"},
        "months": {"past": "{0} kuud", "future": "{0} kuu"},
        "year": {"past": "üks aasta", "future": "ühe aasta"},
        "years": {"past": "{0} aastat", "future": "{0} aasta"},
    }

    month_names = [
        "",
        "Jaanuar",
        "Veebruar",
        "Märts",
        "Aprill",
        "Mai",
        "Juuni",
        "Juuli",
        "August",
        "September",
        "Oktoober",
        "November",
        "Detsember",
    ]
    month_abbreviations = [
        "",
        "Jan",
        "Veb",
        "Mär",
        "Apr",
        "Mai",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Okt",
        "Nov",
        "Dets",
    ]

    day_names = [
        "",
        "Esmaspäev",
        "Teisipäev",
        "Kolmapäev",
        "Neljapäev",
        "Reede",
        "Laupäev",
        "Pühapäev",
    ]
    day_abbreviations = ["", "Esm", "Teis", "Kolm", "Nelj", "Re", "Lau", "Püh"]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]
        if delta > 0:
            _form = form["future"]
        else:
            _form = form["past"]
        return _form.format(abs(delta))


class LatvianLocale(Locale):
    names = ["lv", "lv-lv"]

    past = "pirms {0}"
    future = "pēc {0}"
    and_word = "un"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "tagad",
        "second": "sekundes",
        "seconds": "{0} sekundēm",
        "minute": "minūtes",
        "minutes": "{0} minūtēm",
        "hour": "stundas",
        "hours": "{0} stundām",
        "day": "dienas",
        "days": "{0} dienām",
        "week": "nedēļas",
        "weeks": "{0} nedēļām",
        "month": "mēneša",
        "months": "{0} mēnešiem",
        "year": "gada",
        "years": "{0} gadiem",
    }

    month_names = [
        "",
        "janvāris",
        "februāris",
        "marts",
        "aprīlis",
        "maijs",
        "jūnijs",
        "jūlijs",
        "augusts",
        "septembris",
        "oktobris",
        "novembris",
        "decembris",
    ]

    month_abbreviations = [
        "",
        "jan",
        "feb",
        "marts",
        "apr",
        "maijs",
        "jūnijs",
        "jūlijs",
        "aug",
        "sept",
        "okt",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "pirmdiena",
        "otrdiena",
        "trešdiena",
        "ceturtdiena",
        "piektdiena",
        "sestdiena",
        "svētdiena",
    ]

    day_abbreviations = [
        "",
        "pi",
        "ot",
        "tr",
        "ce",
        "pi",
        "se",
        "sv",
    ]


class SwahiliLocale(Locale):
    names = [
        "sw",
        "sw-ke",
        "sw-tz",
    ]

    past = "{0} iliyopita"
    future = "muda wa {0}"
    and_word = "na"

    timeframes = {
        "now": "sasa hivi",
        "second": "sekunde",
        "seconds": "sekunde {0}",
        "minute": "dakika moja",
        "minutes": "dakika {0}",
        "hour": "saa moja",
        "hours": "saa {0}",
        "day": "siku moja",
        "days": "siku {0}",
        "week": "wiki moja",
        "weeks": "wiki {0}",
        "month": "mwezi moja",
        "months": "miezi {0}",
        "year": "mwaka moja",
        "years": "miaka {0}",
    }

    meridians = {"am": "asu", "pm": "mch", "AM": "ASU", "PM": "MCH"}

    month_names = [
        "",
        "Januari",
        "Februari",
        "Machi",
        "Aprili",
        "Mei",
        "Juni",
        "Julai",
        "Agosti",
        "Septemba",
        "Oktoba",
        "Novemba",
        "Desemba",
    ]
    month_abbreviations = [
        "",
        "Jan",
        "Feb",
        "Mac",
        "Apr",
        "Mei",
        "Jun",
        "Jul",
        "Ago",
        "Sep",
        "Okt",
        "Nov",
        "Des",
    ]

    day_names = [
        "",
        "Jumatatu",
        "Jumanne",
        "Jumatano",
        "Alhamisi",
        "Ijumaa",
        "Jumamosi",
        "Jumapili",
    ]
    day_abbreviations = [
        "",
        "Jumatatu",
        "Jumanne",
        "Jumatano",
        "Alhamisi",
        "Ijumaa",
        "Jumamosi",
        "Jumapili",
    ]


class CroatianLocale(Locale):
    names = ["hr", "hr-hr"]

    past = "prije {0}"
    future = "za {0}"
    and_word = "i"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "upravo sad",
        "second": "sekundu",
        "seconds": {"double": "{0} sekunde", "higher": "{0} sekundi"},
        "minute": "minutu",
        "minutes": {"double": "{0} minute", "higher": "{0} minuta"},
        "hour": "sat",
        "hours": {"double": "{0} sata", "higher": "{0} sati"},
        "day": "jedan dan",
        "days": {"double": "{0} dana", "higher": "{0} dana"},
        "week": "tjedan",
        "weeks": {"double": "{0} tjedna", "higher": "{0} tjedana"},
        "month": "mjesec",
        "months": {"double": "{0} mjeseca", "higher": "{0} mjeseci"},
        "year": "godinu",
        "years": {"double": "{0} godine", "higher": "{0} godina"},
    }

    month_names = [
        "",
        "siječanj",
        "veljača",
        "ožujak",
        "travanj",
        "svibanj",
        "lipanj",
        "srpanj",
        "kolovoz",
        "rujan",
        "listopad",
        "studeni",
        "prosinac",
    ]

    month_abbreviations = [
        "",
        "siječ",
        "velj",
        "ožuj",
        "trav",
        "svib",
        "lip",
        "srp",
        "kol",
        "ruj",
        "list",
        "stud",
        "pros",
    ]

    day_names = [
        "",
        "ponedjeljak",
        "utorak",
        "srijeda",
        "četvrtak",
        "petak",
        "subota",
        "nedjelja",
    ]

    day_abbreviations = [
        "",
        "po",
        "ut",
        "sr",
        "če",
        "pe",
        "su",
        "ne",
    ]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]
        delta = abs(delta)
        if isinstance(form, Mapping):
            if 1 < delta <= 4:
                form = form["double"]
            else:
                form = form["higher"]

        return form.format(delta)


class LatinLocale(Locale):
    names = ["la", "la-va"]

    past = "ante {0}"
    future = "in {0}"
    and_word = "et"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "nunc",
        "second": "secundum",
        "seconds": "{0} secundis",
        "minute": "minutam",
        "minutes": "{0} minutis",
        "hour": "horam",
        "hours": "{0} horas",
        "day": "diem",
        "days": "{0} dies",
        "week": "hebdomadem",
        "weeks": "{0} hebdomades",
        "month": "mensem",
        "months": "{0} mensis",
        "year": "annum",
        "years": "{0} annos",
    }

    month_names = [
        "",
        "Ianuarius",
        "Februarius",
        "Martius",
        "Aprilis",
        "Maius",
        "Iunius",
        "Iulius",
        "Augustus",
        "September",
        "October",
        "November",
        "December",
    ]

    month_abbreviations = [
        "",
        "Ian",
        "Febr",
        "Mart",
        "Apr",
        "Mai",
        "Iun",
        "Iul",
        "Aug",
        "Sept",
        "Oct",
        "Nov",
        "Dec",
    ]

    day_names = [
        "",
        "dies Lunae",
        "dies Martis",
        "dies Mercurii",
        "dies Iovis",
        "dies Veneris",
        "dies Saturni",
        "dies Solis",
    ]

    day_abbreviations = [
        "",
        "dies Lunae",
        "dies Martis",
        "dies Mercurii",
        "dies Iovis",
        "dies Veneris",
        "dies Saturni",
        "dies Solis",
    ]


class LithuanianLocale(Locale):
    names = ["lt", "lt-lt"]

    past = "prieš {0}"
    future = "po {0}"
    and_word = "ir"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "dabar",
        "second": "sekundės",
        "seconds": "{0} sekundžių",
        "minute": "minutės",
        "minutes": "{0} minučių",
        "hour": "valandos",
        "hours": "{0} valandų",
        "day": "dieną",
        "days": "{0} dienų",
        "week": "savaitės",
        "weeks": "{0} savaičių",
        "month": "mėnesio",
        "months": "{0} mėnesių",
        "year": "metų",
        "years": "{0} metų",
    }

    month_names = [
        "",
        "sausis",
        "vasaris",
        "kovas",
        "balandis",
        "gegužė",
        "birželis",
        "liepa",
        "rugpjūtis",
        "rugsėjis",
        "spalis",
        "lapkritis",
        "gruodis",
    ]

    month_abbreviations = [
        "",
        "saus",
        "vas",
        "kovas",
        "bal",
        "geg",
        "birž",
        "liepa",
        "rugp",
        "rugs",
        "spalis",
        "lapkr",
        "gr",
    ]

    day_names = [
        "",
        "pirmadienis",
        "antradienis",
        "trečiadienis",
        "ketvirtadienis",
        "penktadienis",
        "šeštadienis",
        "sekmadienis",
    ]

    day_abbreviations = [
        "",
        "pi",
        "an",
        "tr",
        "ke",
        "pe",
        "še",
        "se",
    ]


class MalayLocale(Locale):
    names = ["ms", "ms-my", "ms-bn"]

    past = "{0} yang lalu"
    future = "dalam {0}"
    and_word = "dan"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "sekarang",
        "second": "saat",
        "seconds": "{0} saat",
        "minute": "minit",
        "minutes": "{0} minit",
        "hour": "jam",
        "hours": "{0} jam",
        "day": "hari",
        "days": "{0} hari",
        "week": "minggu",
        "weeks": "{0} minggu",
        "month": "bulan",
        "months": "{0} bulan",
        "year": "tahun",
        "years": "{0} tahun",
    }

    month_names = [
        "",
        "Januari",
        "Februari",
        "Mac",
        "April",
        "Mei",
        "Jun",
        "Julai",
        "Ogos",
        "September",
        "Oktober",
        "November",
        "Disember",
    ]

    month_abbreviations = [
        "",
        "Jan.",
        "Feb.",
        "Mac",
        "Apr.",
        "Mei",
        "Jun",
        "Julai",
        "Og.",
        "Sept.",
        "Okt.",
        "Nov.",
        "Dis.",
    ]

    day_names = [
        "",
        "Isnin",
        "Selasa",
        "Rabu",
        "Khamis",
        "Jumaat",
        "Sabtu",
        "Ahad",
    ]

    day_abbreviations = [
        "",
        "Isnin",
        "Selasa",
        "Rabu",
        "Khamis",
        "Jumaat",
        "Sabtu",
        "Ahad",
    ]


class MalteseLocale(Locale):
    names = ["mt", "mt-mt"]

    past = "{0} ilu"
    future = "fi {0}"
    and_word = "u"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "issa",
        "second": "sekonda",
        "seconds": "{0} sekondi",
        "minute": "minuta",
        "minutes": "{0} minuti",
        "hour": "siegħa",
        "hours": {"dual": "{0} sagħtejn", "plural": "{0} sigħat"},
        "day": "jum",
        "days": {"dual": "{0} jumejn", "plural": "{0} ijiem"},
        "week": "ġimgħa",
        "weeks": {"dual": "{0} ġimagħtejn", "plural": "{0} ġimgħat"},
        "month": "xahar",
        "months": {"dual": "{0} xahrejn", "plural": "{0} xhur"},
        "year": "sena",
        "years": {"dual": "{0} sentejn", "plural": "{0} snin"},
    }

    month_names = [
        "",
        "Jannar",
        "Frar",
        "Marzu",
        "April",
        "Mejju",
        "Ġunju",
        "Lulju",
        "Awwissu",
        "Settembru",
        "Ottubru",
        "Novembru",
        "Diċembru",
    ]

    month_abbreviations = [
        "",
        "Jan",
        "Fr",
        "Mar",
        "Apr",
        "Mejju",
        "Ġun",
        "Lul",
        "Aw",
        "Sett",
        "Ott",
        "Nov",
        "Diċ",
    ]

    day_names = [
        "",
        "It-Tnejn",
        "It-Tlieta",
        "L-Erbgħa",
        "Il-Ħamis",
        "Il-Ġimgħa",
        "Is-Sibt",
        "Il-Ħadd",
    ]

    day_abbreviations = [
        "",
        "T",
        "TL",
        "E",
        "Ħ",
        "Ġ",
        "S",
        "Ħ",
    ]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]
        delta = abs(delta)
        if isinstance(form, Mapping):
            if delta == 2:
                form = form["dual"]
            else:
                form = form["plural"]

        return form.format(delta)


class SamiLocale(Locale):
    names = ["se", "se-fi", "se-no", "se-se"]

    past = "{0} dassái"
    future = "{0} "  # NOTE: couldn't find preposition for Sami here, none needed?

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "dál",
        "second": "sekunda",
        "seconds": "{0} sekundda",
        "minute": "minuhta",
        "minutes": "{0} minuhta",
        "hour": "diimmu",
        "hours": "{0} diimmu",
        "day": "beaivvi",
        "days": "{0} beaivvi",
        "week": "vahku",
        "weeks": "{0} vahku",
        "month": "mánu",
        "months": "{0} mánu",
        "year": "jagi",
        "years": "{0} jagi",
    }

    month_names = [
        "",
        "Ođđajagimánnu",
        "Guovvamánnu",
        "Njukčamánnu",
        "Cuoŋománnu",
        "Miessemánnu",
        "Geassemánnu",
        "Suoidnemánnu",
        "Borgemánnu",
        "Čakčamánnu",
        "Golggotmánnu",
        "Skábmamánnu",
        "Juovlamánnu",
    ]

    month_abbreviations = [
        "",
        "Ođđajagimánnu",
        "Guovvamánnu",
        "Njukčamánnu",
        "Cuoŋománnu",
        "Miessemánnu",
        "Geassemánnu",
        "Suoidnemánnu",
        "Borgemánnu",
        "Čakčamánnu",
        "Golggotmánnu",
        "Skábmamánnu",
        "Juovlamánnu",
    ]

    day_names = [
        "",
        "Mánnodat",
        "Disdat",
        "Gaskavahkku",
        "Duorastat",
        "Bearjadat",
        "Lávvordat",
        "Sotnabeaivi",
    ]

    day_abbreviations = [
        "",
        "Mánnodat",
        "Disdat",
        "Gaskavahkku",
        "Duorastat",
        "Bearjadat",
        "Lávvordat",
        "Sotnabeaivi",
    ]


class OdiaLocale(Locale):
    names = ["or", "or-in"]

    past = "{0} ପୂର୍ବେ"
    future = "{0} ପରେ"

    timeframes = {
        "now": "ବର୍ତ୍ତମାନ",
        "second": "ଏକ ସେକେଣ୍ଡ",
        "seconds": "{0} ସେକେଣ୍ଡ",
        "minute": "ଏକ ମିନଟ",
        "minutes": "{0} ମିନଟ",
        "hour": "ଏକ ଘଣ୍ଟା",
        "hours": "{0} ଘଣ୍ଟା",
        "day": "ଏକ ଦିନ",
        "days": "{0} ଦିନ",
        "month": "ଏକ ମାସ",
        "months": "{0} ମାସ ",
        "year": "ଏକ ବର୍ଷ",
        "years": "{0} ବର୍ଷ",
    }

    meridians = {"am": "ପୂର୍ବାହ୍ନ", "pm": "ଅପରାହ୍ନ", "AM": "ପୂର୍ବାହ୍ନ", "PM": "ଅପରାହ୍ନ"}

    month_names = [
        "",
        "ଜାନୁଆରୀ",
        "ଫେବୃଆରୀ",
        "ମାର୍ଚ୍ଚ୍",
        "ଅପ୍ରେଲ",
        "ମଇ",
        "ଜୁନ୍",
        "ଜୁଲାଇ",
        "ଅଗଷ୍ଟ",
        "ସେପ୍ଟେମ୍ବର",
        "ଅକ୍ଟୋବର୍",
        "ନଭେମ୍ବର୍",
        "ଡିସେମ୍ବର୍",
    ]
    month_abbreviations = [
        "",
        "ଜାନୁ",
        "ଫେବୃ",
        "ମାର୍ଚ୍ଚ୍",
        "ଅପ୍ରେ",
        "ମଇ",
        "ଜୁନ୍",
        "ଜୁଲା",
        "ଅଗ",
        "ସେପ୍ଟେ",
        "ଅକ୍ଟୋ",
        "ନଭେ",
        "ଡିସେ",
    ]

    day_names = [
        "",
        "ସୋମବାର",
        "ମଙ୍ଗଳବାର",
        "ବୁଧବାର",
        "ଗୁରୁବାର",
        "ଶୁକ୍ରବାର",
        "ଶନିବାର",
        "ରବିବାର",
    ]
    day_abbreviations = [
        "",
        "ସୋମ",
        "ମଙ୍ଗଳ",
        "ବୁଧ",
        "ଗୁରୁ",
        "ଶୁକ୍ର",
        "ଶନି",
        "ରବି",
    ]

    def _ordinal_number(self, n: int) -> str:
        if n > 10 or n == 0:
            return f"{n}ତମ"
        if n in [1, 5, 7, 8, 9, 10]:
            return f"{n}ମ"
        if n in [2, 3]:
            return f"{n}ୟ"
        if n == 4:
            return f"{n}ର୍ଥ"
        if n == 6:
            return f"{n}ଷ୍ଠ"
        return ""


class SerbianLocale(Locale):
    names = ["sr", "sr-rs", "sr-sp"]

    past = "pre {0}"
    future = "za {0}"
    and_word = "i"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[str, Mapping[str, str]]]] = {
        "now": "sada",
        "second": "sekundu",
        "seconds": {"double": "{0} sekunde", "higher": "{0} sekundi"},
        "minute": "minutu",
        "minutes": {"double": "{0} minute", "higher": "{0} minuta"},
        "hour": "sat",
        "hours": {"double": "{0} sata", "higher": "{0} sati"},
        "day": "dan",
        "days": {"double": "{0} dana", "higher": "{0} dana"},
        "week": "nedelju",
        "weeks": {"double": "{0} nedelje", "higher": "{0} nedelja"},
        "month": "mesec",
        "months": {"double": "{0} meseca", "higher": "{0} meseci"},
        "year": "godinu",
        "years": {"double": "{0} godine", "higher": "{0} godina"},
    }

    month_names = [
        "",
        "januar",  # јануар
        "februar",  # фебруар
        "mart",  # март
        "april",  # април
        "maj",  # мај
        "jun",  # јун
        "jul",  # јул
        "avgust",  # август
        "septembar",  # септембар
        "oktobar",  # октобар
        "novembar",  # новембар
        "decembar",  # децембар
    ]

    month_abbreviations = [
        "",
        "jan",
        "feb",
        "mar",
        "apr",
        "maj",
        "jun",
        "jul",
        "avg",
        "sep",
        "okt",
        "nov",
        "dec",
    ]

    day_names = [
        "",
        "ponedeljak",  # понедељак
        "utorak",  # уторак
        "sreda",  # среда
        "četvrtak",  # четвртак
        "petak",  # петак
        "subota",  # субота
        "nedelja",  # недеља
    ]

    day_abbreviations = [
        "",
        "po",  # по
        "ut",  # ут
        "sr",  # ср
        "če",  # че
        "pe",  # пе
        "su",  # су
        "ne",  # не
    ]

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        form = self.timeframes[timeframe]
        delta = abs(delta)
        if isinstance(form, Mapping):
            if 1 < delta <= 4:
                form = form["double"]
            else:
                form = form["higher"]

        return form.format(delta)


class LuxembourgishLocale(Locale):
    names = ["lb", "lb-lu"]

    past = "virun {0}"
    future = "an {0}"
    and_word = "an"

    timeframes: ClassVar[Dict[TimeFrameLiteral, str]] = {
        "now": "just elo",
        "second": "enger Sekonn",
        "seconds": "{0} Sekonnen",
        "minute": "enger Minutt",
        "minutes": "{0} Minutten",
        "hour": "enger Stonn",
        "hours": "{0} Stonnen",
        "day": "engem Dag",
        "days": "{0} Deeg",
        "week": "enger Woch",
        "weeks": "{0} Wochen",
        "month": "engem Mount",
        "months": "{0} Méint",
        "year": "engem Joer",
        "years": "{0} Jahren",
    }

    timeframes_only_distance = timeframes.copy()
    timeframes_only_distance["second"] = "eng Sekonn"
    timeframes_only_distance["minute"] = "eng Minutt"
    timeframes_only_distance["hour"] = "eng Stonn"
    timeframes_only_distance["day"] = "een Dag"
    timeframes_only_distance["days"] = "{0} Deeg"
    timeframes_only_distance["week"] = "eng Woch"
    timeframes_only_distance["month"] = "ee Mount"
    timeframes_only_distance["months"] = "{0} Méint"
    timeframes_only_distance["year"] = "ee Joer"
    timeframes_only_distance["years"] = "{0} Joer"

    month_names = [
        "",
        "Januar",
        "Februar",
        "Mäerz",
        "Abrëll",
        "Mee",
        "Juni",
        "Juli",
        "August",
        "September",
        "Oktouber",
        "November",
        "Dezember",
    ]

    month_abbreviations = [
        "",
        "Jan",
        "Feb",
        "Mäe",
        "Abr",
        "Mee",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Okt",
        "Nov",
        "Dez",
    ]

    day_names = [
        "",
        "Méindeg",
        "Dënschdeg",
        "Mëttwoch",
        "Donneschdeg",
        "Freideg",
        "Samschdeg",
        "Sonndeg",
    ]

    day_abbreviations = ["", "Méi", "Dën", "Mët", "Don", "Fre", "Sam", "Son"]

    def _ordinal_number(self, n: int) -> str:
        return f"{n}."

    def describe(
        self,
        timeframe: TimeFrameLiteral,
        delta: Union[int, float] = 0,
        only_distance: bool = False,
    ) -> str:
        if not only_distance:
            return super().describe(timeframe, delta, only_distance)

        # Luxembourgish uses a different case without 'in' or 'ago'
        humanized: str = self.timeframes_only_distance[timeframe].format(
            trunc(abs(delta))
        )

        return humanized


class ZuluLocale(Locale):
    names = ["zu", "zu-za"]

    past = "{0} edlule"
    future = "{0} "
    and_word = "futhi"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[Mapping[str, str], str]]] = {
        "now": "manje",
        "second": {"past": "umzuzwana", "future": "ngomzuzwana"},
        "seconds": {"past": "{0} imizuzwana", "future": "{0} ngemizuzwana"},
        "minute": {"past": "umzuzu", "future": "ngomzuzu"},
        "minutes": {"past": "{0} imizuzu", "future": "{0} ngemizuzu"},
        "hour": {"past": "ihora", "future": "ngehora"},
        "hours": {"past": "{0} amahora", "future": "{0} emahoreni"},
        "day": {"past": "usuku", "future": "ngosuku"},
        "days": {"past": "{0} izinsuku", "future": "{0} ezinsukwini"},
        "week": {"past": "isonto", "future": "ngesonto"},
        "weeks": {"past": "{0} amasonto", "future": "{0} emasontweni"},
        "month": {"past": "inyanga", "future": "ngenyanga"},
        "months": {"past": "{0} izinyanga", "future": "{0} ezinyangeni"},
        "year": {"past": "unyaka", "future": "ngonyak"},
        "years": {"past": "{0} iminyaka", "future": "{0} eminyakeni"},
    }

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        """Zulu aware time frame format function, takes into account
        the differences between past and future forms."""
        abs_delta = abs(delta)
        form = self.timeframes[timeframe]

        if isinstance(form, str):
            return form.format(abs_delta)

        if delta > 0:
            key = "future"
        else:
            key = "past"
        form = form[key]

        return form.format(abs_delta)

    month_names = [
        "",
        "uMasingane",
        "uNhlolanja",
        "uNdasa",
        "UMbasa",
        "UNhlaba",
        "UNhlangulana",
        "uNtulikazi",
        "UNcwaba",
        "uMandulo",
        "uMfumfu",
        "uLwezi",
        "uZibandlela",
    ]

    month_abbreviations = [
        "",
        "uMasingane",
        "uNhlolanja",
        "uNdasa",
        "UMbasa",
        "UNhlaba",
        "UNhlangulana",
        "uNtulikazi",
        "UNcwaba",
        "uMandulo",
        "uMfumfu",
        "uLwezi",
        "uZibandlela",
    ]

    day_names = [
        "",
        "uMsombuluko",
        "uLwesibili",
        "uLwesithathu",
        "uLwesine",
        "uLwesihlanu",
        "uMgqibelo",
        "iSonto",
    ]

    day_abbreviations = [
        "",
        "uMsombuluko",
        "uLwesibili",
        "uLwesithathu",
        "uLwesine",
        "uLwesihlanu",
        "uMgqibelo",
        "iSonto",
    ]


class TamilLocale(Locale):
    names = ["ta", "ta-in", "ta-lk"]

    past = "{0} நேரத்திற்கு முன்பு"
    future = "இல் {0}"

    timeframes = {
        "now": "இப்போது",
        "second": "ஒரு இரண்டாவது",
        "seconds": "{0} விநாடிகள்",
        "minute": "ஒரு நிமிடம்",
        "minutes": "{0} நிமிடங்கள்",
        "hour": "ஒரு மணி",
        "hours": "{0} மணிநேரம்",
        "day": "ஒரு நாள்",
        "days": "{0} நாட்கள்",
        "week": "ஒரு வாரம்",
        "weeks": "{0} வாரங்கள்",
        "month": "ஒரு மாதம்",
        "months": "{0} மாதங்கள்",
        "year": "ஒரு ஆண்டு",
        "years": "{0} ஆண்டுகள்",
    }

    month_names = [
        "",
        "சித்திரை",
        "வைகாசி",
        "ஆனி",
        "ஆடி",
        "ஆவணி",
        "புரட்டாசி",
        "ஐப்பசி",
        "கார்த்திகை",
        "மார்கழி",
        "தை",
        "மாசி",
        "பங்குனி",
    ]

    month_abbreviations = [
        "",
        "ஜன",
        "பிப்",
        "மார்",
        "ஏப்",
        "மே",
        "ஜூன்",
        "ஜூலை",
        "ஆக",
        "செப்",
        "அக்",
        "நவ",
        "டிச",
    ]

    day_names = [
        "",
        "திங்கட்கிழமை",
        "செவ்வாய்க்கிழமை",
        "புதன்கிழமை",
        "வியாழக்கிழமை",
        "வெள்ளிக்கிழமை",
        "சனிக்கிழமை",
        "ஞாயிற்றுக்கிழமை",
    ]

    day_abbreviations = [
        "",
        "திங்கட்",
        "செவ்வாய்",
        "புதன்",
        "வியாழன்",
        "வெள்ளி",
        "சனி",
        "ஞாயிறு",
    ]

    def _ordinal_number(self, n: int) -> str:
        if n == 1:
            return f"{n}வது"
        elif n >= 0:
            return f"{n}ஆம்"
        else:
            return ""


class AlbanianLocale(Locale):
    names = ["sq", "sq-al"]

    past = "{0} më parë"
    future = "në {0}"
    and_word = "dhe"

    timeframes = {
        "now": "tani",
        "second": "sekondë",
        "seconds": "{0} sekonda",
        "minute": "minutë",
        "minutes": "{0} minuta",
        "hour": "orë",
        "hours": "{0} orë",
        "day": "ditë",
        "days": "{0} ditë",
        "week": "javë",
        "weeks": "{0} javë",
        "month": "muaj",
        "months": "{0} muaj",
        "year": "vit",
        "years": "{0} vjet",
    }

    month_names = [
        "",
        "janar",
        "shkurt",
        "mars",
        "prill",
        "maj",
        "qershor",
        "korrik",
        "gusht",
        "shtator",
        "tetor",
        "nëntor",
        "dhjetor",
    ]

    month_abbreviations = [
        "",
        "jan",
        "shk",
        "mar",
        "pri",
        "maj",
        "qer",
        "korr",
        "gush",
        "sht",
        "tet",
        "nën",
        "dhj",
    ]

    day_names = [
        "",
        "e hënë",
        "e martë",
        "e mërkurë",
        "e enjte",
        "e premte",
        "e shtunë",
        "e diel",
    ]

    day_abbreviations = [
        "",
        "hën",
        "mar",
        "mër",
        "enj",
        "pre",
        "sht",
        "die",
    ]


class GeorgianLocale(Locale):
    names = ["ka", "ka-ge"]

    past = "{0} წინ"  # ts’in
    future = "{0} შემდეგ"  # shemdeg
    and_word = "და"  # da

    timeframes = {
        "now": "ახლა",  # akhla
        # When a cardinal qualifies a noun, it stands in the singular
        "second": "წამის",  # ts’amis
        "seconds": "{0} წამის",
        "minute": "წუთის",  # ts’utis
        "minutes": "{0} წუთის",
        "hour": "საათის",  # saatis
        "hours": "{0} საათის",
        "day": "დღის",  # dghis
        "days": "{0} დღის",
        "week": "კვირის",  # k’viris
        "weeks": "{0} კვირის",
        "month": "თვის",  # tvis
        "months": "{0} თვის",
        "year": "წლის",  # ts’lis
        "years": "{0} წლის",
    }

    month_names = [
        # modern month names
        "",
        "იანვარი",  # Ianvari
        "თებერვალი",  # Tebervali
        "მარტი",  # Mart'i
        "აპრილი",  # Ap'rili
        "მაისი",  # Maisi
        "ივნისი",  # Ivnisi
        "ივლისი",  # Ivlisi
        "აგვისტო",  # Agvist'o
        "სექტემბერი",  # Sekt'emberi
        "ოქტომბერი",  # Okt'omberi
        "ნოემბერი",  # Noemberi
        "დეკემბერი",  # Dek'emberi
    ]

    month_abbreviations = [
        # no abbr. found yet
        "",
        "იანვარი",  # Ianvari
        "თებერვალი",  # Tebervali
        "მარტი",  # Mart'i
        "აპრილი",  # Ap'rili
        "მაისი",  # Maisi
        "ივნისი",  # Ivnisi
        "ივლისი",  # Ivlisi
        "აგვისტო",  # Agvist'o
        "სექტემბერი",  # Sekt'emberi
        "ოქტომბერი",  # Okt'omberi
        "ნოემბერი",  # Noemberi
        "დეკემბერი",  # Dek'emberi
    ]

    day_names = [
        "",
        "ორშაბათი",  # orshabati
        "სამშაბათი",  # samshabati
        "ოთხშაბათი",  # otkhshabati
        "ხუთშაბათი",  # khutshabati
        "პარასკევი",  # p’arask’evi
        "შაბათი",  # shabati
        # "k’vira" also serves as week; to avoid confusion "k’vira-dge" can be used for Sunday
        "კვირა",  # k’vira
    ]

    day_abbreviations = [
        "",
        "ორშაბათი",  # orshabati
        "სამშაბათი",  # samshabati
        "ოთხშაბათი",  # otkhshabati
        "ხუთშაბათი",  # khutshabati
        "პარასკევი",  # p’arask’evi
        "შაბათი",  # shabati
        "კვირა",  # k’vira
    ]


class SinhalaLocale(Locale):
    names = ["si", "si-lk"]

    past = "{0}ට පෙර"
    future = "{0}"
    and_word = "සහ"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[Mapping[str, str], str]]] = {
        "now": "දැන්",
        "second": {
            "past": "තත්පරයක",
            "future": "තත්පරයකින්",
        },  # ක් is the article
        "seconds": {
            "past": "තත්පර {0} ක",
            "future": "තත්පර {0} කින්",
        },
        "minute": {
            "past": "විනාඩියක",
            "future": "විනාඩියකින්",
        },
        "minutes": {
            "past": "විනාඩි {0} ක",
            "future": "මිනිත්තු {0} කින්",
        },
        "hour": {"past": "පැයක", "future": "පැයකින්"},
        "hours": {
            "past": "පැය {0} ක",
            "future": "පැය {0} කින්",
        },
        "day": {"past": "දිනක", "future": "දිනකට"},
        "days": {
            "past": "දින {0} ක",
            "future": "දින {0} කින්",
        },
        "week": {"past": "සතියක", "future": "සතියකින්"},
        "weeks": {
            "past": "සති {0} ක",
            "future": "සති {0} කින්",
        },
        "month": {"past": "මාසයක", "future": "එය මාසය තුළ"},
        "months": {
            "past": "මාස {0} ක",
            "future": "මාස {0} කින්",
        },
        "year": {"past": "වසරක", "future": "වසරක් තුළ"},
        "years": {
            "past": "අවුරුදු {0} ක",
            "future": "අවුරුදු {0} තුළ",
        },
    }
    # Sinhala: the general format to describe timeframe is different from past and future,
    # so we do not copy the original timeframes dictionary
    timeframes_only_distance = {}
    timeframes_only_distance["second"] = "තත්පරයක්"
    timeframes_only_distance["seconds"] = "තත්පර {0}"
    timeframes_only_distance["minute"] = "මිනිත්තුවක්"
    timeframes_only_distance["minutes"] = "විනාඩි {0}"
    timeframes_only_distance["hour"] = "පැයක්"
    timeframes_only_distance["hours"] = "පැය {0}"
    timeframes_only_distance["day"] = "දවසක්"
    timeframes_only_distance["days"] = "දවස් {0}"
    timeframes_only_distance["week"] = "සතියක්"
    timeframes_only_distance["weeks"] = "සති {0}"
    timeframes_only_distance["month"] = "මාසයක්"
    timeframes_only_distance["months"] = "මාස {0}"
    timeframes_only_distance["year"] = "අවුරුද්දක්"
    timeframes_only_distance["years"] = "අවුරුදු {0}"

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        """
        Sinhala awares time frame format function, takes into account
        the differences between general, past, and future forms (three different suffixes).
        """
        abs_delta = abs(delta)
        form = self.timeframes[timeframe]

        if isinstance(form, str):
            return form.format(abs_delta)

        if delta > 0:
            key = "future"
        else:
            key = "past"
        form = form[key]

        return form.format(abs_delta)

    def describe(
        self,
        timeframe: TimeFrameLiteral,
        delta: Union[float, int] = 1,  # key is always future when only_distance=False
        only_distance: bool = False,
    ) -> str:
        """Describes a delta within a timeframe in plain language.

        :param timeframe: a string representing a timeframe.
        :param delta: a quantity representing a delta in a timeframe.
        :param only_distance: return only distance eg: "11 seconds" without "in" or "ago" keywords
        """

        if not only_distance:
            return super().describe(timeframe, delta, only_distance)
        # Sinhala uses a different case without 'in' or 'ago'
        humanized = self.timeframes_only_distance[timeframe].format(trunc(abs(delta)))

        return humanized

    month_names = [
        "",
        "ජනවාරි",
        "පෙබරවාරි",
        "මාර්තු",
        "අප්‍රේල්",
        "මැයි",
        "ජූනි",
        "ජූලි",
        "අගෝස්තු",
        "සැප්තැම්බර්",
        "ඔක්තෝබර්",
        "නොවැම්බර්",
        "දෙසැම්බර්",
    ]

    month_abbreviations = [
        "",
        "ජන",
        "පෙබ",
        "මාර්",
        "අප්‍රේ",
        "මැයි",
        "ජුනි",
        "ජූලි",
        "අගෝ",
        "සැප්",
        "ඔක්",
        "නොවැ",
        "දෙසැ",
    ]

    day_names = [
        "",
        "සදුදා",
        "අඟහරැවදා",
        "බදාදා",
        "බ්‍රහස්‍පතින්‍දා",
        "සිකුරාදා",
        "සෙනසුරාදා",
        "ඉරිදා",
    ]

    day_abbreviations = [
        "",
        "සදුද",
        "බදා",
        "බදා",
        "සිකු",
        "සෙන",
        "අ",
        "ඉරිදා",
    ]


class UrduLocale(Locale):
    names = ["ur", "ur-pk"]

    past = "پہلے {0}"
    future = "میں {0}"
    and_word = "اور"

    timeframes = {
        "now": "ابھی",
        "second": "ایک سیکنڈ",
        "seconds": "{0} سیکنڈ",
        "minute": "ایک منٹ",
        "minutes": "{0} منٹ",
        "hour": "ایک گھنٹے",
        "hours": "{0} گھنٹے",
        "day": "ایک دن",
        "days": "{0} دن",
        "week": "ایک ہفتے",
        "weeks": "{0} ہفتے",
        "month": "ایک مہینہ",
        "months": "{0} ماہ",
        "year": "ایک سال",
        "years": "{0} سال",
    }

    month_names = [
        "",
        "جنوری",
        "فروری",
        "مارچ",
        "اپریل",
        "مئی",
        "جون",
        "جولائی",
        "اگست",
        "ستمبر",
        "اکتوبر",
        "نومبر",
        "دسمبر",
    ]

    month_abbreviations = [
        "",
        "جنوری",
        "فروری",
        "مارچ",
        "اپریل",
        "مئی",
        "جون",
        "جولائی",
        "اگست",
        "ستمبر",
        "اکتوبر",
        "نومبر",
        "دسمبر",
    ]

    day_names = [
        "",
        "سوموار",
        "منگل",
        "بدھ",
        "جمعرات",
        "جمعہ",
        "ہفتہ",
        "اتوار",
    ]

    day_abbreviations = [
        "",
        "سوموار",
        "منگل",
        "بدھ",
        "جمعرات",
        "جمعہ",
        "ہفتہ",
        "اتوار",
    ]


class KazakhLocale(Locale):
    names = ["kk", "kk-kz"]

    past = "{0} бұрын"
    future = "{0} кейін"
    timeframes = {
        "now": "қазір",
        "second": "бір секунд",
        "seconds": "{0} секунд",
        "minute": "бір минут",
        "minutes": "{0} минут",
        "hour": "бір сағат",
        "hours": "{0} сағат",
        "day": "бір күн",
        "days": "{0} күн",
        "week": "бір апта",
        "weeks": "{0} апта",
        "month": "бір ай",
        "months": "{0} ай",
        "year": "бір жыл",
        "years": "{0} жыл",
    }

    month_names = [
        "",
        "Қаңтар",
        "Ақпан",
        "Наурыз",
        "Сәуір",
        "Мамыр",
        "Маусым",
        "Шілде",
        "Тамыз",
        "Қыркүйек",
        "Қазан",
        "Қараша",
        "Желтоқсан",
    ]
    month_abbreviations = [
        "",
        "Қан",
        "Ақп",
        "Нау",
        "Сәу",
        "Мам",
        "Мау",
        "Шіл",
        "Там",
        "Қыр",
        "Қаз",
        "Қар",
        "Жел",
    ]

    day_names = [
        "",
        "Дүйсембі",
        "Сейсенбі",
        "Сәрсенбі",
        "Бейсенбі",
        "Жұма",
        "Сенбі",
        "Жексенбі",
    ]
    day_abbreviations = ["", "Дс", "Сс", "Ср", "Бс", "Жм", "Сб", "Жс"]


class AmharicLocale(Locale):
    names = ["am", "am-et"]

    past = "{0} በፊት"
    future = "{0} ውስጥ"
    and_word = "እና"

    timeframes: ClassVar[Mapping[TimeFrameLiteral, Union[Mapping[str, str], str]]] = {
        "now": "አሁን",
        "second": {
            "past": "ከአንድ ሰከንድ",
            "future": "በአንድ ሰከንድ",
        },
        "seconds": {
            "past": "ከ {0} ሰከንድ",
            "future": "በ {0} ሰከንድ",
        },
        "minute": {
            "past": "ከአንድ ደቂቃ",
            "future": "በአንድ ደቂቃ",
        },
        "minutes": {
            "past": "ከ {0} ደቂቃዎች",
            "future": "በ {0} ደቂቃዎች",
        },
        "hour": {
            "past": "ከአንድ ሰዓት",
            "future": "በአንድ ሰዓት",
        },
        "hours": {
            "past": "ከ {0} ሰዓታት",
            "future": "በ {0} ሰከንድ",
        },
        "day": {
            "past": "ከአንድ ቀን",
            "future": "በአንድ ቀን",
        },
        "days": {
            "past": "ከ {0} ቀናት",
            "future": "በ {0} ቀናት",
        },
        "week": {
            "past": "ከአንድ ሳምንት",
            "future": "በአንድ ሳምንት",
        },
        "weeks": {
            "past": "ከ {0} ሳምንታት",
            "future": "በ {0} ሳምንታት",
        },
        "month": {
            "past": "ከአንድ ወር",
            "future": "በአንድ ወር",
        },
        "months": {
            "past": "ከ {0} ወር",
            "future": "በ {0} ወራት",
        },
        "year": {
            "past": "ከአንድ አመት",
            "future": "በአንድ አመት",
        },
        "years": {
            "past": "ከ {0} ዓመታት",
            "future": "በ {0} ዓመታት",
        },
    }
    # Amharic: the general format to describe timeframe is different from past and future,
    # so we do not copy the original timeframes dictionary
    timeframes_only_distance = {
        "second": "አንድ ሰከንድ",
        "seconds": "{0} ሰከንድ",
        "minute": "አንድ ደቂቃ",
        "minutes": "{0} ደቂቃዎች",
        "hour": "አንድ ሰዓት",
        "hours": "{0} ሰዓት",
        "day": "አንድ ቀን",
        "days": "{0} ቀናት",
        "week": "አንድ ሳምንት",
        "weeks": "{0} ሳምንት",
        "month": "አንድ ወር",
        "months": "{0} ወራት",
        "year": "አንድ አመት",
        "years": "{0} ዓመታት",
    }

    month_names = [
        "",
        "ጃንዩወሪ",
        "ፌብሩወሪ",
        "ማርች",
        "ኤፕሪል",
        "ሜይ",
        "ጁን",
        "ጁላይ",
        "ኦገስት",
        "ሴፕቴምበር",
        "ኦክቶበር",
        "ኖቬምበር",
        "ዲሴምበር",
    ]

    month_abbreviations = [
        "",
        "ጃንዩ",
        "ፌብሩ",
        "ማርች",
        "ኤፕሪ",
        "ሜይ",
        "ጁን",
        "ጁላይ",
        "ኦገስ",
        "ሴፕቴ",
        "ኦክቶ",
        "ኖቬም",
        "ዲሴም",
    ]

    day_names = [
        "",
        "ሰኞ",
        "ማክሰኞ",
        "ረቡዕ",
        "ሐሙስ",
        "ዓርብ",
        "ቅዳሜ",
        "እሑድ",
    ]
    day_abbreviations = ["", "እ", "ሰ", "ማ", "ረ", "ሐ", "ዓ", "ቅ"]

    def _ordinal_number(self, n: int) -> str:
        return f"{n}ኛ"

    def _format_timeframe(self, timeframe: TimeFrameLiteral, delta: int) -> str:
        """
        Amharic awares time frame format function, takes into account
        the differences between general, past, and future forms (three different suffixes).
        """
        abs_delta = abs(delta)
        form = self.timeframes[timeframe]

        if isinstance(form, str):
            return form.format(abs_delta)

        if delta > 0:
            key = "future"
        else:
            key = "past"
        form = form[key]

        return form.format(abs_delta)

    def describe(
        self,
        timeframe: TimeFrameLiteral,
        delta: Union[float, int] = 1,  # key is always future when only_distance=False
        only_distance: bool = False,
    ) -> str:
        """Describes a delta within a timeframe in plain language.

        :param timeframe: a string representing a timeframe.
        :param delta: a quantity representing a delta in a timeframe.
        :param only_distance: return only distance eg: "11 seconds" without "in" or "ago" keywords
        """

        if not only_distance:
            return super().describe(timeframe, delta, only_distance)
        humanized = self.timeframes_only_distance[timeframe].format(trunc(abs(delta)))

        return humanized


class ArmenianLocale(Locale):
    names = ["hy", "hy-am"]
    past = "{0} առաջ"
    future = "{0}ից"
    and_word = "Եվ"  # Yev

    timeframes = {
        "now": "հիմա",
        "second": "վայրկյան",
        "seconds": "{0} վայրկյան",
        "minute": "րոպե",
        "minutes": "{0} րոպե",
        "hour": "ժամ",
        "hours": "{0} ժամ",
        "day": "օր",
        "days": "{0} օր",
        "month": "ամիս",
        "months": "{0} ամիս",
        "year": "տարին",
        "years": "{0} տարին",
        "week": "շաբաթ",
        "weeks": "{0} շաբաթ",
    }

    meridians = {
        "am": "Ամ",
        "pm": "պ.մ.",
        "AM": "Ամ",
        "PM": "պ.մ.",
    }

    month_names = [
        "",
        "հունվար",
        "փետրվար",
        "մարտ",
        "ապրիլ",
        "մայիս",
        "հունիս",
        "հուլիս",
        "օգոստոս",
        "սեպտեմբեր",
        "հոկտեմբեր",
        "նոյեմբեր",
        "դեկտեմբեր",
    ]

    month_abbreviations = [
        "",
        "հունվար",
        "փետրվար",
        "մարտ",
        "ապրիլ",
        "մայիս",
        "հունիս",
        "հուլիս",
        "օգոստոս",
        "սեպտեմբեր",
        "հոկտեմբեր",
        "նոյեմբեր",
        "դեկտեմբեր",
    ]

    day_names = [
        "",
        "երկուշաբթի",
        "երեքշաբթի",
        "չորեքշաբթի",
        "հինգշաբթի",
        "ուրբաթ",
        "շաբաթ",
        "կիրակի",
    ]

    day_abbreviations = [
        "",
        "երկ.",
        "երեք.",
        "չորեք.",
        "հինգ.",
        "ուրբ.",
        "շաբ.",
        "կիր.",
    ]


class UzbekLocale(Locale):
    names = ["uz", "uz-uz"]
    past = "{0}dan avval"
    future = "{0}dan keyin"
    timeframes = {
        "now": "hozir",
        "second": "bir soniya",
        "seconds": "{0} soniya",
        "minute": "bir daqiqa",
        "minutes": "{0} daqiqa",
        "hour": "bir soat",
        "hours": "{0} soat",
        "day": "bir kun",
        "days": "{0} kun",
        "week": "bir hafta",
        "weeks": "{0} hafta",
        "month": "bir oy",
        "months": "{0} oy",
        "year": "bir yil",
        "years": "{0} yil",
    }

    month_names = [
        "",
        "Yanvar",
        "Fevral",
        "Mart",
        "Aprel",
        "May",
        "Iyun",
        "Iyul",
        "Avgust",
        "Sentyabr",
        "Oktyabr",
        "Noyabr",
        "Dekabr",
    ]

    month_abbreviations = [
        "",
        "Yan",
        "Fev",
        "Mar",
        "Apr",
        "May",
        "Iyn",
        "Iyl",
        "Avg",
        "Sen",
        "Okt",
        "Noy",
        "Dek",
    ]

    day_names = [
        "",
        "Dushanba",
        "Seshanba",
        "Chorshanba",
        "Payshanba",
        "Juma",
        "Shanba",
        "Yakshanba",
    ]

    day_abbreviations = ["", "Dush", "Sesh", "Chor", "Pay", "Jum", "Shan", "Yak"]
