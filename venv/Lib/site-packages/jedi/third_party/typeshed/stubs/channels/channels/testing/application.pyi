from typing import Any

from asgiref.testing import ApplicationCommunicator as BaseApplicationCommunicator

def no_op() -> None: ...

class ApplicationCommunicator(BaseApplicationCommunicator):
    # ASGI messages are dictionaries with a "type" key and protocol-specific fields.
    # Dictionary values can be strings, bytes, lists, or other types depending on the protocol:
    # - HTTP: {"type": "http.request", "body": b"request data", "headers": [...], ...}
    # - WebSocket: {"type": "websocket.receive", "bytes": b"binary data"} or {"text": "string"}
    # - Custom protocols: Application-specific message dictionaries
    async def send_input(self, message: dict[str, Any]) -> None: ...
    async def receive_output(self, timeout: float = 1) -> dict[str, Any]: ...

    # The following methods are not present in the original source code,
    # but are commonly used in practice. Since the base package doesn't
    # provide type hints for them, they are added here to improve type correctness.
    async def receive_nothing(self, timeout: float = 0.1, interval: float = 0.01) -> bool: ...
    async def wait(self, timeout: float = 1) -> None: ...
    def stop(self, exceptions: bool = True) -> None: ...
