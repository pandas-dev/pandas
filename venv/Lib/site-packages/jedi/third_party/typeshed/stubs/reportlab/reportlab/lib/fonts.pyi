from typing import Final, Literal
from typing_extensions import TypeAlias

_BoolInt: TypeAlias = Literal[0, 1]

__version__: Final[str]

def ps2tt(psfn: str) -> tuple[str, _BoolInt, _BoolInt]: ...
def tt2ps(fn: str, b: _BoolInt, i: _BoolInt) -> str: ...
def addMapping(face: str, bold: _BoolInt, italic: _BoolInt, psname: str) -> None: ...
