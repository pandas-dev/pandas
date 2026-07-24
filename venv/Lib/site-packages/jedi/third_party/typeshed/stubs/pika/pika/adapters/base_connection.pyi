import abc
from _typeshed import Incomplete
from collections.abc import Callable
from logging import Logger
from typing_extensions import Self

from ..adapters.utils import nbio_interface
from ..connection import Connection

LOGGER: Logger

class BaseConnection(Connection, metaclass=abc.ABCMeta):
    def __init__(
        self,
        parameters,
        on_open_callback: Callable[[Self], object] | None,
        on_open_error_callback: Callable[[Self, BaseException], object] | None,
        on_close_callback: Callable[[Self, BaseException], object] | None,
        nbio,
        internal_connection_workflow: bool,
    ) -> None: ...
    @classmethod
    @abc.abstractmethod
    def create_connection(cls, connection_configs, on_done, custom_ioloop=None, workflow=None): ...
    @property
    def ioloop(self): ...

class _StreamingProtocolShim(nbio_interface.AbstractStreamProtocol):
    connection_made: Incomplete
    connection_lost: Incomplete
    eof_received: Incomplete
    data_received: Incomplete
    conn: Incomplete
    def __init__(self, conn) -> None: ...
    def __getattr__(self, attr: str): ...
