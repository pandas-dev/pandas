from collections.abc import Iterator

from dateparser.conf import Settings
from dateparser.languages.locale import Locale

class BaseLanguageDetector:
    languages: list[Locale]
    def __init__(self, languages: list[Locale]) -> None: ...
    def iterate_applicable_languages(
        self, date_string: str, modify: bool = False, settings: Settings | None = None
    ) -> Iterator[Locale]: ...

class AutoDetectLanguage(BaseLanguageDetector):
    language_pool: list[Locale]
    allow_redetection: bool
    def __init__(self, languages: list[Locale], allow_redetection: bool = False) -> None: ...
    def iterate_applicable_languages(
        self, date_string: str, modify: bool = False, settings: Settings | None = None
    ) -> Iterator[Locale]: ...

class ExactLanguages(BaseLanguageDetector):
    def __init__(self, languages: list[Locale]) -> None: ...
    def iterate_applicable_languages(
        self, date_string: str, modify: bool = False, settings: Settings | None = None
    ) -> Iterator[Locale]: ...
