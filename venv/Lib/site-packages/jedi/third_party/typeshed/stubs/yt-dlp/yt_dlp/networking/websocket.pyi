import abc
from typing import Any

from .common import RequestHandler, Response

class WebSocketResponse(Response):
    # Both raise NotImplementedError.
    def send(self, message: bytes | str) -> Any: ...
    def recv(self) -> Any: ...

class WebSocketRequestHandler(RequestHandler, abc.ABC, metaclass=abc.ABCMeta): ...
