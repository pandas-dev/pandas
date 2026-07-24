debug: bool

class Error(Exception):
    msg: str
    def __init__(self, msg: str) -> None: ...

class ValidationError(Error):
    msg: str
    file: str
    def __init__(self, msg: str, file: str) -> None: ...

class ParsingError(Error):
    msg: str
    file: str
    def __init__(self, msg: str, file: str) -> None: ...

class NoKeyError(Error):
    key: str
    group: str
    file: str
    def __init__(self, key: str, group: str, file: str) -> None: ...

class DuplicateKeyError(Error):
    key: str
    group: str
    file: str
    def __init__(self, key: str, group: str, file: str) -> None: ...

class NoGroupError(Error):
    group: str
    file: str
    def __init__(self, group: str, file: str) -> None: ...

class DuplicateGroupError(Error):
    group: str
    file: str
    def __init__(self, group: str, file: str) -> None: ...

class NoThemeError(Error):
    theme: str
    def __init__(self, theme: str) -> None: ...
