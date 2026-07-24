from typing import Final, Protocol, overload, type_check_only

@type_check_only
class _SupportsGet(Protocol):
    @overload
    def get(self, name: str, /) -> str | None: ...
    @overload
    def get(self, name: str, default: str, /) -> str: ...

NO_COLORS: Final[int]
ANSI_8_COLORS: Final[int]
ANSI_16_COLORS: Final[int]
ANSI_256_COLORS: Final[int]
TRUE_COLORS: Final[int]

def detect_color_support(env: _SupportsGet) -> int: ...
