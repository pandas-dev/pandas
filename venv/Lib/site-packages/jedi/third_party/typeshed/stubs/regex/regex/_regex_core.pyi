import enum
from collections.abc import Callable
from typing import Any, AnyStr, Generic
from typing_extensions import TypeAlias

from ._main import Pattern

__all__ = [
    "A",
    "ASCII",
    "B",
    "BESTMATCH",
    "D",
    "DEBUG",
    "E",
    "ENHANCEMATCH",
    "F",
    "FULLCASE",
    "I",
    "IGNORECASE",
    "L",
    "LOCALE",
    "M",
    "MULTILINE",
    "P",
    "POSIX",
    "R",
    "REVERSE",
    "S",
    "DOTALL",
    "T",
    "TEMPLATE",
    "U",
    "UNICODE",
    "V0",
    "VERSION0",
    "V1",
    "VERSION1",
    "W",
    "WORD",
    "X",
    "VERBOSE",
    "error",
    "Scanner",
    "RegexFlag",
]

class error(Exception):
    def __init__(self, message: str, pattern: AnyStr | None = None, pos: int | None = None) -> None: ...

class RegexFlag(enum.IntFlag):
    A = 0x80
    ASCII = A
    B = 0x1000
    BESTMATCH = B
    D = 0x200
    DEBUG = D
    E = 0x8000
    ENHANCEMATCH = E
    F = 0x4000
    FULLCASE = F
    I = 0x2
    IGNORECASE = I
    L = 0x4
    LOCALE = L
    M = 0x8
    MULTILINE = M
    P = 0x10000
    POSIX = P
    R = 0x400
    REVERSE = R
    T = 0x1
    TEMPLATE = T
    S = 0x10
    DOTALL = S
    U = 0x20
    UNICODE = U
    V0 = 0x2000
    VERSION0 = V0
    V1 = 0x100
    VERSION1 = V1
    W = 0x800
    WORD = W
    X = 0x40
    VERBOSE = X

ASCII = RegexFlag.ASCII
BESTMATCH = RegexFlag.BESTMATCH
DEBUG = RegexFlag.DEBUG
ENHANCEMATCH = RegexFlag.ENHANCEMATCH
FULLCASE = RegexFlag.FULLCASE
IGNORECASE = RegexFlag.IGNORECASE
LOCALE = RegexFlag.LOCALE
MULTILINE = RegexFlag.MULTILINE
POSIX = RegexFlag.POSIX
REVERSE = RegexFlag.REVERSE
TEMPLATE = RegexFlag.TEMPLATE
DOTALL = RegexFlag.DOTALL
UNICODE = RegexFlag.UNICODE
VERBOSE = RegexFlag.VERBOSE
VERSION0 = RegexFlag.VERSION0
VERSION1 = RegexFlag.VERSION1
WORD = RegexFlag.WORD
A = RegexFlag.A
B = RegexFlag.B
D = RegexFlag.D
E = RegexFlag.E
F = RegexFlag.F
I = RegexFlag.I
L = RegexFlag.L
M = RegexFlag.M
P = RegexFlag.P
R = RegexFlag.R
S = RegexFlag.S
U = RegexFlag.U
V0 = RegexFlag.V0
V1 = RegexFlag.V1
W = RegexFlag.W
X = RegexFlag.X
T = RegexFlag.T

DEFAULT_VERSION = VERSION1

_Lexicon: TypeAlias = list[tuple[AnyStr, Callable[[Scanner[AnyStr], AnyStr], Any]]]

class Scanner(Generic[AnyStr]):
    lexicon: _Lexicon[AnyStr]
    scanner: Pattern[AnyStr]

    def __init__(self, lexicon: _Lexicon[AnyStr], flags: int = 0) -> None: ...
    def scan(self, string: AnyStr) -> tuple[list[Any], AnyStr]: ...
