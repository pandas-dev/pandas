from abc import abstractmethod
from logging import Logger
from typing import Generic, TypeVar

from .amqp_object import AMQPObject, Method as AMQPMethod
from .spec import BasicProperties

_M = TypeVar("_M", bound=AMQPMethod)

LOGGER: Logger

class Frame(AMQPObject):
    frame_type: int
    channel_number: int
    def __init__(self, frame_type: int, channel_number: int) -> None: ...
    @abstractmethod
    def marshal(self) -> bytes: ...

class Method(Frame, Generic[_M]):
    method: _M
    def __init__(self, channel_number: int, method: _M) -> None: ...
    def marshal(self) -> bytes: ...

class Header(Frame):
    body_size: int
    properties: BasicProperties
    def __init__(self, channel_number: int, body_size: int, props: BasicProperties) -> None: ...
    def marshal(self) -> bytes: ...

class Body(Frame):
    fragment: bytes
    def __init__(self, channel_number: int, fragment: bytes) -> None: ...
    def marshal(self) -> bytes: ...

class Heartbeat(Frame):
    def __init__(self) -> None: ...
    def marshal(self) -> bytes: ...

class ProtocolHeader(AMQPObject):
    frame_type: int
    major: int
    minor: int
    revision: int
    def __init__(self, major: int | None = None, minor: int | None = None, revision: int | None = None) -> None: ...
    def marshal(self) -> bytes: ...

def decode_frame(data_in: bytes) -> tuple[int, Frame | None]: ...
