from re import Pattern
from typing import Final

BWS: Final = "[ \t]{0,}?"
CHUNK_EXT: Final[str]
CHUNK_EXT_NAME: Final = r"[!#$%&'*+\-.^_`|~0-9A-Za-z]{1,}"
CHUNK_EXT_RE: Final[Pattern[str]]
CHUNK_EXT_VAL: Final[str]
DIGIT: Final = "[0-9]"
FIELD_CONTENT: Final[str]
FIELD_VALUE: Final[str]
FIELD_VCHAR: Final = r"[\x21-\x7e\x80-\xff]"
HEADER_FIELD_RE: Final[Pattern[str]]
HEXDIG: Final = r"[0-9a-fA-F]"
OBS_TEXT: Final = r"\x80-\xff"
ONLY_DIGIT_RE: Final[Pattern[str]]
ONLY_HEXDIG_RE: Final[Pattern[str]]
OWS: Final = "[ \t]{0,}?"
QDTEXT: Final = "[\t !#-[\\]-~\\x80-\\xff]"
QUOTED_PAIR: Final = "\\\\([\t \\x21-\\x7e\\x80-\\xff])"
QUOTED_PAIR_RE: Final[Pattern[str]]
QUOTED_STRING: Final[str]
QUOTED_STRING_RE: Final[Pattern[str]]
RWS: Final = "[ \t]{1,}?"
TCHAR: Final = r"[!#$%&'*+\-.^_`|~0-9A-Za-z]"
TOKEN: Final = r"[!#$%&'*+\-.^_`|~0-9A-Za-z]{1,}"
VCHAR: Final = r"\x21-\x7e"
WS: Final = "[ \t]"
