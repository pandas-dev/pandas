from typing import Final, overload

from . import FixedBody, Observer
from ._libastro import _DateInitType

db: Final[str]
stars: dict[str, FixedBody]

@overload
def star(name: str, observer: Observer, /) -> FixedBody: ...
@overload
def star(name: str, when: _DateInitType = ..., epoch: _DateInitType = ...) -> FixedBody: ...

STAR_NUMBER_NAME: Final[dict[int, str]]
STAR_NAME_NUMBER: Final[dict[str, int]]
