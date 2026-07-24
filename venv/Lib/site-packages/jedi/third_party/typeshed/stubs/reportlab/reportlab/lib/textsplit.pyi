import re
from _typeshed import Incomplete, ReadableBuffer
from collections.abc import Sequence
from typing import Final

__version__: Final[str]
CANNOT_START_LINE: Final[Sequence[str]]
ALL_CANNOT_START: Final[str]
CANNOT_END_LINE: Final[Sequence[str]]
ALL_CANNOT_END: Final[str]

def is_multi_byte(ch: str | bytes | bytearray) -> bool: ...
def getCharWidths(word: str, fontName: str, fontSize: float) -> list[float]: ...
def wordSplit(word, maxWidths, fontName, fontSize, encoding: str = "utf8") -> list[list[Incomplete]]: ...
def dumbSplit(word, widths, maxWidths) -> list[list[Incomplete]]: ...
def kinsokuShoriSplit(word, widths, availWidth) -> None: ...

rx: re.Pattern[str]

def cjkwrap(text: ReadableBuffer, width: float, encoding: str = "utf8"): ...
