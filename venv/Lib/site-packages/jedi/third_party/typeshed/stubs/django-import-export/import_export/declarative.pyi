import _typeshed
from logging import Logger
from typing import Any

logger: Logger

class DeclarativeMetaclass(type):
    def __new__(cls: type[_typeshed.Self], name: str, bases: tuple[type[Any], ...], attrs: dict[str, Any]) -> _typeshed.Self: ...

class ModelDeclarativeMetaclass(DeclarativeMetaclass):
    def __new__(cls: type[_typeshed.Self], name: str, bases: tuple[type[Any], ...], attrs: dict[str, Any]) -> _typeshed.Self: ...
