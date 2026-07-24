from typing import Final

openers: Final[str]
closers: Final[str]
delimiters: Final[str]
closing_delimiters: Final[str]
quote_pairs: dict[str, str]

def match_chars(c1: str, c2: str) -> bool: ...
