from _typeshed import Incomplete
from collections.abc import Callable
from logging import Logger
from typing import Any, Final
from typing_extensions import Self

from .callback import CallbackManager
from .connection import Connection
from .data import _ArgumentMapping
from .exchange_type import ExchangeType
from .frame import Body, Header, Method
from .spec import Basic, BasicProperties, Confirm, Exchange, Queue, Tx

LOGGER: Logger
MAX_CHANNELS: Final[int]

class Channel:
    CLOSED: Final = 0
    OPENING: Final = 1
    OPEN: Final = 2
    CLOSING: Final = 3

    channel_number: int
    callbacks: CallbackManager
    connection: Connection
    flow_active: bool

    def __init__(self, connection: Connection, channel_number: int, on_open_callback: Callable[[Self], object]) -> None: ...
    def __int__(self) -> int: ...
    def add_callback(self, callback, replies, one_shot: bool = True) -> None: ...
    def add_on_cancel_callback(self, callback) -> None: ...
    def add_on_close_callback(self, callback) -> None: ...
    def add_on_flow_callback(self, callback) -> None: ...
    def add_on_return_callback(self, callback) -> None: ...
    def basic_ack(self, delivery_tag: int = 0, multiple: bool = False) -> None: ...
    def basic_cancel(
        self, consumer_tag: str = "", callback: Callable[[Method[Basic.CancelOk]], object] | None = None
    ) -> None: ...
    def basic_consume(
        self,
        queue: str,
        on_message_callback: Callable[[Channel, Basic.Deliver, BasicProperties, bytes], object],
        auto_ack: bool = False,
        exclusive: bool = False,
        consumer_tag: str | None = None,
        arguments: _ArgumentMapping | None = None,
        callback: Callable[[Method[Basic.ConsumeOk]], object] | None = None,
    ) -> str: ...
    def basic_get(
        self, queue: str, callback: Callable[[Channel, Basic.GetOk, BasicProperties, bytes], object], auto_ack: bool = False
    ) -> None: ...
    def basic_nack(self, delivery_tag: int = 0, multiple: bool = False, requeue: bool = True) -> None: ...
    def basic_publish(
        self,
        exchange: str,
        routing_key: str,
        body: str | bytes,
        properties: BasicProperties | None = None,
        mandatory: bool = False,
    ) -> None: ...
    def basic_qos(
        self,
        prefetch_size: int = 0,
        prefetch_count: int = 0,
        global_qos: bool = False,
        callback: Callable[[Method[Basic.QosOk]], object] | None = None,
    ) -> None: ...
    def basic_reject(self, delivery_tag: int = 0, requeue: bool = True) -> None: ...
    def basic_recover(
        self, requeue: bool = False, callback: Callable[[Method[Basic.RecoverOk]], object] | None = None
    ) -> None: ...
    def close(self, reply_code: int = 0, reply_text: str = "Normal shutdown") -> None: ...
    def confirm_delivery(
        self,
        ack_nack_callback: Callable[[Method[Basic.Ack | Basic.Nack]], object],
        callback: Callable[[Method[Confirm.SelectOk]], object] | None = None,
    ) -> None: ...
    @property
    def consumer_tags(self) -> list[str]: ...
    def exchange_bind(
        self,
        destination: str,
        source: str,
        routing_key: str = "",
        arguments: _ArgumentMapping | None = None,
        callback: Callable[[Method[Exchange.BindOk]], object] | None = None,
    ) -> None: ...
    def exchange_declare(
        self,
        exchange: str,
        exchange_type: ExchangeType | str = ...,
        passive: bool = False,
        durable: bool = False,
        auto_delete: bool = False,
        internal: bool = False,
        arguments: _ArgumentMapping | None = None,
        callback: Callable[[Method[Exchange.DeclareOk]], object] | None = None,
    ) -> None: ...
    def exchange_delete(
        self,
        exchange: str | None = None,
        if_unused: bool = False,
        callback: Callable[[Method[Exchange.DeleteOk]], object] | None = None,
    ) -> None: ...
    def exchange_unbind(
        self,
        destination: str | None = None,
        source: str | None = None,
        routing_key: str = "",
        arguments: _ArgumentMapping | None = None,
        callback: Callable[[Method[Exchange.UnbindOk]], object] | None = None,
    ) -> None: ...
    def flow(self, active: bool, callback: Callable[[bool], object] | None = None) -> None: ...
    @property
    def is_closed(self) -> bool: ...
    @property
    def is_closing(self) -> bool: ...
    @property
    def is_open(self) -> bool: ...
    @property
    def is_opening(self) -> bool: ...
    def open(self) -> None: ...
    def queue_bind(
        self,
        queue: str,
        exchange: str,
        routing_key: str | None = None,
        arguments: _ArgumentMapping | None = None,
        callback: Callable[[Method[Queue.BindOk]], object] | None = None,
    ) -> None: ...
    def queue_declare(
        self,
        queue: str,
        passive: bool = False,
        durable: bool = False,
        exclusive: bool = False,
        auto_delete: bool = False,
        arguments: _ArgumentMapping | None = None,
        callback: Callable[[Method[Queue.DeclareOk]], object] | None = None,
    ) -> None: ...
    def queue_delete(
        self,
        queue: str,
        if_unused: bool = False,
        if_empty: bool = False,
        callback: Callable[[Method[Queue.DeleteOk]], object] | None = None,
    ) -> None: ...
    def queue_purge(self, queue: str, callback: Callable[[Method[Queue.PurgeOk]], object] | None = None) -> None: ...
    def queue_unbind(
        self,
        queue: str,
        exchange: str | None = None,
        routing_key: str | None = None,
        arguments: _ArgumentMapping | None = None,
        callback: Callable[[Method[Queue.UnbindOk]], object] | None = None,
    ): ...
    def tx_commit(self, callback: Callable[[Method[Tx.CommitOk]], object] | None = None) -> None: ...
    def tx_rollback(self, callback: Callable[[Method[Tx.RollbackOk]], object] | None = None) -> None: ...
    def tx_select(self, callback: Callable[[Method[Tx.SelectOk]], object] | None = None) -> None: ...

class ContentFrameAssembler:
    def __init__(self) -> None: ...
    def process(self, frame_value: Method[Any] | Header | Body) -> tuple[Incomplete, Incomplete, bytes] | None: ...
