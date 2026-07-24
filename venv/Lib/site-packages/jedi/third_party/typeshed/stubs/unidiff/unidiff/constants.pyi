from re import Pattern
from typing import Final

RE_SOURCE_FILENAME: Pattern[str]
RE_TARGET_FILENAME: Pattern[str]
RE_DIFF_GIT_HEADER: Pattern[str]
RE_DIFF_GIT_HEADER_URI_LIKE: Pattern[str]
RE_DIFF_GIT_HEADER_NO_PREFIX: Pattern[str]
RE_DIFF_GIT_DELETED_FILE: Pattern[str]
RE_DIFF_GIT_NEW_FILE: Pattern[str]
RE_HUNK_HEADER: Pattern[str]
RE_HUNK_BODY_LINE: Pattern[str]
RE_HUNK_EMPTY_BODY_LINE: Pattern[str]
RE_NO_NEWLINE_MARKER: Pattern[str]
RE_BINARY_DIFF: Pattern[str]

DEFAULT_ENCODING: Final = "UTF-8"
DEV_NULL: Final = "/dev/null"
LINE_TYPE_ADDED: Final = "+"
LINE_TYPE_REMOVED: Final = "-"
LINE_TYPE_CONTEXT: Final = " "
LINE_TYPE_EMPTY: Final = ""
LINE_TYPE_NO_NEWLINE: Final = "\\"
LINE_VALUE_NO_NEWLINE: Final = " No newline at end of file"
