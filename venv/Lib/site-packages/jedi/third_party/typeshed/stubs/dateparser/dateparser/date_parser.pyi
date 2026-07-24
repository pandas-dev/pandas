from collections.abc import Callable
from datetime import datetime, tzinfo

from dateparser.conf import Settings

class DateParser:
    def parse(
        self,
        date_string: str,
        parse_method: Callable[[str, Settings, tzinfo | None], tuple[datetime, str]],
        settings: Settings | None = None,
    ) -> tuple[datetime, str]: ...

date_parser: DateParser
