from collections.abc import Callable

from simplejson.decoder import JSONDecoder

class JSONDecodeError(ValueError):
    msg: str
    doc: str
    pos: int
    end: int | None
    lineno: int
    colno: int
    endlineno: int | None
    endcolno: int | None

def make_scanner(context: JSONDecoder) -> Callable[[str, int], tuple[bool, int]]: ...

__all__ = ["make_scanner", "JSONDecodeError"]
