from collections.abc import Set as AbstractSet
from datetime import datetime
from typing import Any, Literal, overload

from dateparser.conf import Settings

from ..date import _DetectLanguagesFunction

@overload
def search_dates(
    text: str,
    languages: list[str] | tuple[str, ...] | AbstractSet[str] | None,
    settings: Settings | dict[str, Any] | None,
    add_detected_language: Literal[True],
    detect_languages_function: _DetectLanguagesFunction | None = None,
) -> list[tuple[str, datetime, str]] | None: ...
@overload
def search_dates(
    text: str,
    languages: list[str] | tuple[str, ...] | AbstractSet[str] | None = None,
    settings: Settings | dict[str, Any] | None = None,
    add_detected_language: Literal[False] = False,
    detect_languages_function: _DetectLanguagesFunction | None = None,
) -> list[tuple[str, datetime]] | None: ...
