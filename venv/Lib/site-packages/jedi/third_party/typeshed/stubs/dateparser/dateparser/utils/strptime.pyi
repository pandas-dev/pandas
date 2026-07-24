import re
from datetime import datetime
from time import struct_time
from typing import Final, Protocol, type_check_only

TIME_MATCHER: Final[re.Pattern[str]]
MS_SEARCHER: Final[re.Pattern[str]]

@type_check_only
class _strptime_time(Protocol):
    def __call__(self, data_string: str, format: str = "%a %b %d %H:%M:%S %Y") -> struct_time: ...

def patch_strptime() -> _strptime_time: ...
def strptime(date_string: str, format: str) -> datetime: ...
