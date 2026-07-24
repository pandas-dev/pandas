import datetime
import re
from _typeshed import Incomplete
from collections import OrderedDict
from collections.abc import Iterable, Sequence, Set as AbstractSet
from typing import Any, Final, TypedDict, type_check_only

from dateparser.conf import Settings
from dateparser.date import DateData, DateDataParser, _DetectLanguagesFunction
from dateparser.languages.loader import LocaleDataLoader
from dateparser.languages.locale import Locale
from dateparser.search.text_detection import FullTextLanguageDetector

@type_check_only
class _SearchDates(TypedDict):
    Language: str | None
    Dates: list[tuple[str, datetime.datetime]] | None

RELATIVE_REG: Final[re.Pattern[str]]

def date_is_relative(translation: str) -> bool: ...

class _ExactLanguageSearch:
    loader: LocaleDataLoader
    language: Locale | None
    def __init__(self, loader: LocaleDataLoader) -> None: ...
    def get_current_language(self, shortname: str) -> None: ...
    def search(self, shortname: str, text: str, settings: Settings | None) -> tuple[list[str], list[str]]: ...
    @staticmethod
    def set_relative_base(
        substring: str, already_parsed: list[tuple[DateData, bool]]
    ) -> tuple[str, datetime.datetime | None]: ...
    def choose_best_split(
        self, possible_parsed_splits: list[Incomplete], possible_substrings_splits: list[Incomplete]
    ) -> tuple[Incomplete, Incomplete]: ...
    def split_by(self, item: str, original: str, splitter: str) -> list[list[list[str]]]: ...
    def split_if_not_parsed(self, item: str, original: str) -> list[list[list[str]]]: ...
    def parse_item(
        self,
        parser: DateDataParser,
        item: str,
        translated_item: str,
        parsed: list[tuple[DateData, bool]],
        need_relative_base: bool | None,
    ) -> tuple[DateData, bool]: ...
    def parse_found_objects(
        self,
        parser: DateDataParser,
        to_parse: Iterable[str],
        original: Sequence[str],
        translated: Sequence[str],
        settings: Settings,
    ) -> tuple[list[tuple[DateData, bool]], list[str]]: ...
    def search_parse(self, shortname: str, text: str, settings: Settings) -> list[tuple[str, datetime.datetime]]: ...

class DateSearchWithDetection:
    loader: LocaleDataLoader
    available_language_map: OrderedDict[str, Locale]
    search: _ExactLanguageSearch
    language_detector: FullTextLanguageDetector
    def __init__(self) -> None: ...
    def detect_language(
        self,
        text: str,
        languages: list[str] | tuple[str, ...] | AbstractSet[str] | None,
        settings: Settings | dict[str, Any] | None = None,
        detect_languages_function: _DetectLanguagesFunction | None = None,
    ) -> str | None: ...
    def search_dates(
        self,
        text: str,
        languages: list[str] | tuple[str, ...] | AbstractSet[str] | None = None,
        settings: Settings | dict[str, Any] | None = None,
        detect_languages_function: _DetectLanguagesFunction | None = None,
    ) -> _SearchDates: ...
    def preprocess_text(self, text: str, languages: Iterable[str] | None) -> str: ...
