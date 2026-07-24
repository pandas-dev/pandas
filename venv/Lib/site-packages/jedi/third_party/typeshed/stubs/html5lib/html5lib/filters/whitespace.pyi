import re
from _typeshed import Incomplete

from . import base

spaceCharacters: str
SPACES_REGEX: re.Pattern[str]

class Filter(base.Filter[dict[str, Incomplete]]):
    spacePreserveElements: frozenset[str]

def collapse_spaces(text: str) -> str: ...
