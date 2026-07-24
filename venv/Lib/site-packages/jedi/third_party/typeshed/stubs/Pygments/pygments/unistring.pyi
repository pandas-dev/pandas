from typing import Final, Literal, TypeAlias

_Cats: TypeAlias = Literal[
    "Cc",
    "Cf",
    "Cn",
    "Co",
    "Cs",
    "Ll",
    "Lm",
    "Lo",
    "Lt",
    "Lu",
    "Mc",
    "Me",
    "Mn",
    "Nd",
    "Nl",
    "No",
    "Pc",
    "Pd",
    "Pe",
    "Pf",
    "Pi",
    "Po",
    "Ps",
    "Sc",
    "Sk",
    "Sm",
    "So",
    "Zl",
    "Zp",
    "Zs",
]

Cc: Final[str]
Cf: Final[str]
Cn: Final[str]
Co: Final[str]
Cs: Final[str]
Ll: Final[str]
Lm: Final[str]
Lo: Final[str]
Lt: Final[str]
Lu: Final[str]
Mc: Final[str]
Me: Final[str]
Mn: Final[str]
Nd: Final[str]
Nl: Final[str]
No: Final[str]
Pc: Final[str]
Pd: Final[str]
Pe: Final[str]
Pf: Final[str]
Pi: Final[str]
Po: Final[str]
Ps: Final[str]
Sc: Final[str]
Sk: Final[str]
Sm: Final[str]
So: Final[str]
Zl: Final[str]
Zp: Final[str]
Zs: Final[str]
xid_continue: Final[str]
xid_start: Final[str]
cats: Final[list[_Cats]]

def combine(*args: _Cats) -> str: ...
def allexcept(*args: _Cats) -> str: ...
