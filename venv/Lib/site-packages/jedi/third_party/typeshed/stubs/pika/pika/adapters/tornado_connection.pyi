from _typeshed import Incomplete
from logging import Logger

from pika.adapters import base_connection

LOGGER: Logger

class TornadoConnection(base_connection.BaseConnection):
    def __init__(
        self,
        parameters: Incomplete | None = ...,
        on_open_callback: Incomplete | None = ...,
        on_open_error_callback: Incomplete | None = ...,
        on_close_callback: Incomplete | None = ...,
        custom_ioloop: Incomplete | None = ...,
        internal_connection_workflow: bool = ...,
    ) -> None: ...
    @classmethod
    def create_connection(
        cls, connection_configs, on_done, custom_ioloop: Incomplete | None = ..., workflow: Incomplete | None = ...
    ): ...
