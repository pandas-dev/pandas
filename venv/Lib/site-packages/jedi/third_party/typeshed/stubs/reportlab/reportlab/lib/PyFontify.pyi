import re
from _typeshed import Incomplete
from typing import Final

__version__: Final[str]

def replace(src, sep, rep): ...

keywordsList: list[str]
commentPat: str
pat: str
quotePat: str
tripleQuotePat: str
nonKeyPat: str
keyPat: str
matchPat: str
matchRE: re.Pattern[str]
idKeyPat: str
idRE: re.Pattern[str]

def fontify(pytext, searchfrom: int = 0, searchto=None) -> list[tuple[str, int, int, Incomplete]]: ...
def test(path) -> None: ...
