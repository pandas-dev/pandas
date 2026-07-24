class UWSGIParseException(Exception): ...

class InvalidUWSGIHeader(UWSGIParseException):
    msg: str
    code: int

    def __init__(self, msg: str = "") -> None: ...

class UnsupportedModifier(UWSGIParseException):
    modifier: int
    code: int

    def __init__(self, modifier: int) -> None: ...

class ForbiddenUWSGIRequest(UWSGIParseException):
    host: str
    code: int

    def __init__(self, host: str) -> None: ...
