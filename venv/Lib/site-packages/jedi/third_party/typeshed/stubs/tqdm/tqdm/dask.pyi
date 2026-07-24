from _typeshed import Incomplete
from collections.abc import Callable
from typing import ClassVar, type_check_only
from typing_extensions import Self

__all__ = ["TqdmCallback"]

# dask.callbacks.Callback
@type_check_only
class _Callback:
    active: ClassVar[set[tuple[Callable[..., Incomplete] | None, ...]]]
    def __init__(
        self,
        start: Incomplete | None,
        start_state: Incomplete | None,
        pretask: Incomplete | None,
        posttask: Incomplete | None,
        finish: Incomplete | None,
    ) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, *args: object) -> None: ...
    def register(self) -> None: ...
    def unregister(self) -> None: ...

class TqdmCallback(_Callback):
    tqdm_class: type[Incomplete]
    def __init__(
        self, start: Incomplete | None = ..., pretask: Incomplete | None = ..., tqdm_class: type[Incomplete] = ..., **tqdm_kwargs
    ) -> None: ...
    def display(self) -> None: ...
