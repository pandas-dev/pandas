from _typeshed import Incomplete
from collections.abc import Callable
from logging import Logger
from typing import Literal

from pika import amqp_object, frame

LOGGER: Logger

def name_or_value(value: amqp_object.AMQPObject | frame.Frame | int | str) -> str: ...
def sanitize_prefix(function): ...
def check_for_prefix_and_key(function): ...

class CallbackManager:
    CALLS: str
    ARGUMENTS: str
    DUPLICATE_WARNING: str
    CALLBACK: str
    ONE_SHOT: str
    ONLY_CALLER: str
    def __init__(self) -> None: ...
    def add(
        self,
        prefix: str | int,
        key: str | object,
        callback: Callable[[Incomplete], Incomplete],
        one_shot: bool = True,
        only_caller: object | None = None,
        arguments=None,
    ) -> tuple[str | int, str | object]: ...
    def clear(self) -> None: ...
    def cleanup(self, prefix: str | int) -> bool: ...
    def pending(self, prefix: str | int, key: str | object) -> int | None: ...
    def process(self, prefix: str | int, key: str | object, caller, *args, **keywords) -> bool: ...
    def remove(
        self,
        prefix: str | int,
        key: str | object,
        callback_value: Callable[[Incomplete], Incomplete] | None = None,
        arguments=None,
    ) -> Literal[True]: ...
    def remove_all(self, prefix: str | int, key: str | object) -> None: ...
