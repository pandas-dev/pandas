import re
from collections import OrderedDict
from collections.abc import Iterator
from typing import Final

from .locale import Locale

LOCALE_SPLIT_PATTERN: Final[re.Pattern[str]]

class LocaleDataLoader:
    def get_locale_map(
        self,
        languages: list[str] | None = None,
        locales: list[str] | None = None,
        region: str | None = None,
        use_given_order: bool = False,
        allow_conflicting_locales: bool = False,
    ) -> OrderedDict[str, Locale]: ...
    def get_locales(
        self,
        languages: list[str] | None = None,
        locales: list[str] | None = None,
        region: str | None = None,
        use_given_order: bool = False,
        allow_conflicting_locales: bool = False,
    ) -> Iterator[Locale]: ...
    def get_locale(self, shortname: str) -> Locale: ...

default_loader: LocaleDataLoader
