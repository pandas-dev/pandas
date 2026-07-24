from collections.abc import Callable
from typing import Any, overload
from typing_extensions import Self, TypeAlias

from tensorflow import Tensor

class Regularizer:
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...
    def __call__(self, x: Tensor, /) -> Tensor: ...

_Regularizer: TypeAlias = str | dict[str, Any] | Regularizer | None  # noqa: Y047

@overload
def get(identifier: None) -> None: ...
@overload
def get(identifier: str | dict[str, Any] | Regularizer) -> Regularizer: ...
@overload
def get(identifier: Callable[[Tensor], Tensor]) -> Callable[[Tensor], Tensor]: ...
def __getattr__(name: str): ...  # incomplete module
