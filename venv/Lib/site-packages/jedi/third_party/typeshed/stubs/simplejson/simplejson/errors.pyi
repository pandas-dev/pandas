__all__ = ["JSONDecodeError"]

def linecol(doc: str, pos: int) -> tuple[int, int]: ...
def errmsg(msg: str, doc: str, pos: int, end: int | None = None) -> str: ...

class JSONDecodeError(ValueError):
    msg: str
    doc: str
    pos: int
    end: int | None
    lineno: int
    colno: int
    endlineno: int | None
    endcolno: int | None
    def __init__(self, msg: str, doc: str, pos: int, end: int | None = None) -> None: ...
    def __reduce__(self) -> tuple[JSONDecodeError, tuple[str, str, int, int | None]]: ...
