from email.message import Message
from typing import IO, Mapping, Optional, Tuple, Union
from urllib.response import addinfourl

# Stubs for urllib.error

class URLError(IOError):
    reason: Union[str, BaseException]
    def __init__(self, reason: Union[str, BaseException], filename: Optional[str] = ...) -> None: ...

class HTTPError(URLError, addinfourl):
    code: int
    def __init__(self, url: str, code: int, msg: str, hdrs: Mapping[str, str], fp: Optional[IO[bytes]]) -> None: ...

class ContentTooShortError(URLError):
    content: Tuple[str, Message]
    def __init__(self, message: str, content: Tuple[str, Message]) -> None: ...
