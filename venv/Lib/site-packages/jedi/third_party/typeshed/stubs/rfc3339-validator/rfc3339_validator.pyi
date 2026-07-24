import re
from typing import Final

__version__: Final[str]
__author__: Final[str]
__email__: Final[str]
RFC3339_REGEX_FLAGS: Final[int]
RFC3339_REGEX: Final[re.Pattern[str]]

def validate_rfc3339(date_string: str) -> bool: ...
