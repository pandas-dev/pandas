from wsproto.frame_protocol import CloseReason

class SimpleWebsocketError(RuntimeError): ...

class ConnectionError(SimpleWebsocketError):
    status_code: int | None
    def __init__(self, status_code: int | None = None) -> None: ...

class ConnectionClosed(SimpleWebsocketError):
    reason: CloseReason
    message: str | None
    def __init__(self, reason: CloseReason = ..., message: str | None = None) -> None: ...
