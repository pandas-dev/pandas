from collections.abc import Callable
from typing import Any, overload
from typing_extensions import Self, TypeAlias

from tensorflow import Tensor
from tensorflow._aliases import DTypeLike, ShapeLike, TensorCompatible

class Initializer:
    def __call__(self, shape: ShapeLike, dtype: DTypeLike | None = None) -> Tensor: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...

class Constant(Initializer):
    def __init__(self, value: TensorCompatible = 0.0) -> None: ...

class GlorotNormal(Initializer):
    def __init__(self, seed: int | None = None) -> None: ...

class GlorotUniform(Initializer):
    def __init__(self, seed: int | None = None) -> None: ...

class TruncatedNormal(Initializer):
    def __init__(self, mean: TensorCompatible = 0.0, stddev: TensorCompatible = 0.05, seed: int | None = None) -> None: ...

class RandomNormal(Initializer):
    def __init__(self, mean: TensorCompatible = 0.0, stddev: TensorCompatible = 0.05, seed: int | None = None) -> None: ...

class RandomUniform(Initializer):
    def __init__(self, minval: TensorCompatible = -0.05, maxval: TensorCompatible = 0.05, seed: int | None = None) -> None: ...

class Zeros(Initializer): ...

constant = Constant
glorot_normal = GlorotNormal
glorot_uniform = GlorotUniform
truncated_normal = TruncatedNormal
zeros = Zeros

_Initializer: TypeAlias = (  # noqa: Y047
    str | Initializer | type[Initializer] | Callable[[ShapeLike], Tensor] | dict[str, Any] | None
)

@overload
def get(identifier: None) -> None: ...
@overload
def get(identifier: str | Initializer | dict[str, Any] | type[Initializer]) -> Initializer: ...
@overload
def get(identifier: Callable[[ShapeLike], Tensor]) -> Callable[[ShapeLike], Tensor]: ...
def __getattr__(name: str): ...  # incomplete module
