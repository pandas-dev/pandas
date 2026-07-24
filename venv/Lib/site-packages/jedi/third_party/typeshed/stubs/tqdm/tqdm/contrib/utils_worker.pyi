from _typeshed import Incomplete
from collections import deque
from collections.abc import Callable
from concurrent.futures import Future, ThreadPoolExecutor
from typing import TypeVar
from typing_extensions import ParamSpec

__all__ = ["MonoWorker"]

_P = ParamSpec("_P")
_R = TypeVar("_R")

class MonoWorker:
    pool: ThreadPoolExecutor
    futures: deque[Future[Incomplete]]
    def submit(self, func: Callable[_P, _R], *args: _P.args, **kwargs: _P.kwargs) -> Future[_R]: ...
