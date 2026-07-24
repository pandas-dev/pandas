from collections.abc import Callable
from re import Match
from typing import Any, Literal

class JSONDecoder:
    encoding: str
    object_hook: Callable[[dict[Any, Any]], Any] | None
    parse_float: Callable[[str], Any] | None
    parse_int: Callable[[str], Any] | None
    parse_constant: Callable[[str], Any] | None
    strict: bool
    object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None
    memo: dict[Any, Any]
    scan_once: Callable[[str, int], tuple[bool, int]]

    def __init__(
        self,
        encoding: str | None = None,
        object_hook: Callable[[dict[Any, Any]], Any] | None = None,
        parse_float: Callable[[str], Any] | None = None,
        parse_int: Callable[[str], Any] | None = None,
        parse_constant: Callable[[str], Any] | None = None,
        strict: bool = True,
        object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = None,
        allow_nan: bool = False,
    ) -> None: ...
    def decode(self, s: str, _w: Callable[[str, int], Match[str]] = ..., _PY3: Literal[True] = True) -> Any: ...
    def raw_decode(
        self, s: str, idx: int = 0, _w: Callable[[str, int], Match[str]] = ..., _PY3: Literal[True] = True
    ) -> tuple[Any, int]: ...

__all__ = ["JSONDecoder"]
