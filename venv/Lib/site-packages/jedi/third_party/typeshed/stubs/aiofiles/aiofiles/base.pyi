from asyncio.events import AbstractEventLoop
from collections.abc import Awaitable, Callable, Generator
from concurrent.futures import Executor
from contextlib import AbstractAsyncContextManager
from types import TracebackType
from typing import Any, BinaryIO, Generic, TextIO, TypeVar
from typing_extensions import Self

_T = TypeVar("_T")
_V_co = TypeVar("_V_co", covariant=True)

def wrap(func: Callable[..., _T]) -> Callable[..., Awaitable[_T]]: ...

class AsyncBase(Generic[_T]):
    def __init__(self, file: TextIO | BinaryIO | None, loop: AbstractEventLoop | None, executor: Executor | None) -> None: ...
    def __aiter__(self) -> Self: ...
    async def __anext__(self) -> _T: ...

class AsyncIndirectBase(AsyncBase[_T]):
    def __init__(
        self, name: str, loop: AbstractEventLoop | None, executor: Executor | None, indirect: Callable[[], TextIO | BinaryIO]
    ) -> None: ...

class AiofilesContextManager(Awaitable[_V_co], AbstractAsyncContextManager[_V_co]):
    __slots__ = ("_coro", "_obj")
    def __init__(self, coro: Awaitable[_V_co]) -> None: ...
    def __await__(self) -> Generator[Any, Any, _V_co]: ...
    async def __aenter__(self) -> _V_co: ...
    async def __aexit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...
