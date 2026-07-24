from dateparser.conf import Settings
from dateparser.languages.locale import Locale
from dateparser.search.detection import BaseLanguageDetector

class FullTextLanguageDetector(BaseLanguageDetector):
    language_unique_chars: list[set[str]]
    language_chars: list[set[str]]
    def __init__(self, languages: list[Locale]) -> None: ...
    def get_unique_characters(self, settings: Settings) -> None: ...
    def character_check(self, date_string: str, settings: Settings) -> None: ...
