from typing import Final, TypeVar

_T = TypeVar("_T")

__author__: Final[str]
__copyright__: Final[str]
__contributors__: Final[list[str]]
__version__: Final[str]
__license__: Final[str]

escape_range: list[tuple[int, int]]

def encode(c: str) -> str: ...
def iri2uri(uri: _T) -> _T: ...
