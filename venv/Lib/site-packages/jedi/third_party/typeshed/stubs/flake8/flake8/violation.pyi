from logging import Logger
from typing import NamedTuple

LOG: Logger

class Violation(NamedTuple):
    code: str
    filename: str
    line_number: int
    column_number: int
    text: str
    physical_line: str | None
    def is_inline_ignored(self, disable_noqa: bool) -> bool: ...
