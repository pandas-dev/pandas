import abc
from dataclasses import dataclass
from functools import cached_property
from typing import Final

@dataclass(frozen=True)
class JsRuntimeInfo:
    name: str
    path: str
    version: str
    version_tuple: tuple[int, ...]
    supported: bool = True

class JsRuntime(abc.ABC):
    def __init__(self, path: str | None = None) -> None: ...
    @cached_property
    @abc.abstractmethod
    def info(self) -> JsRuntimeInfo | None: ...

class DenoJsRuntime(JsRuntime):
    MIN_SUPPORTED_VERSION: Final[tuple[int, int, int]]
    @cached_property
    def info(self) -> JsRuntimeInfo | None: ...

class BunJsRuntime(JsRuntime):
    MIN_SUPPORTED_VERSION: Final[tuple[int, int, int]]
    @cached_property
    def info(self) -> JsRuntimeInfo | None: ...

class NodeJsRuntime(JsRuntime):
    MIN_SUPPORTED_VERSION: Final[tuple[int, int, int]]
    @cached_property
    def info(self) -> JsRuntimeInfo | None: ...

class QuickJsRuntime(JsRuntime):
    MIN_SUPPORTED_VERSION: Final[tuple[int, int, int]]
    @cached_property
    def info(self) -> JsRuntimeInfo | None: ...
