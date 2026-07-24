from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import Any, Final, Literal, TypeVar, overload
from typing_extensions import Self, TypeAlias, TypeGuard

from tensorflow import Tensor
from tensorflow._aliases import KerasSerializable, TensorCompatible
from tensorflow.keras.metrics import (
    binary_crossentropy as binary_crossentropy,
    categorical_crossentropy as categorical_crossentropy,
)

class Loss(ABC):
    reduction: _ReductionValues
    name: str | None
    def __init__(self, name: str | None = None, reduction: _ReductionValues = "sum_over_batch_size", dtype=None) -> None: ...
    @abstractmethod
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...
    def get_config(self) -> dict[str, Any]: ...
    def __call__(
        self, y_true: TensorCompatible, y_pred: TensorCompatible, sample_weight: TensorCompatible | None = None
    ) -> Tensor: ...

class BinaryCrossentropy(Loss):
    def __init__(
        self,
        from_logits: bool = False,
        label_smoothing: float = 0.0,
        axis: int = -1,
        reduction: _ReductionValues = "sum_over_batch_size",
        name: str | None = "binary_crossentropy",
        dtype=None,
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class BinaryFocalCrossentropy(Loss):
    def __init__(
        self,
        apply_class_balancing: bool = False,
        alpha: float = 0.25,
        gamma: float = 2.0,
        from_logits: bool = False,
        label_smoothing: float = 0.0,
        axis: int = -1,
        reduction: _ReductionValues = "sum_over_batch_size",
        name: str | None = "binary_focal_crossentropy",
        dtype=None,
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class CategoricalCrossentropy(Loss):
    def __init__(
        self,
        from_logits: bool = False,
        label_smoothing: float = 0.0,
        axis: int = -1,
        reduction: _ReductionValues = "sum_over_batch_size",
        name: str | None = "categorical_crossentropy",
        dtype=None,
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class CategoricalHinge(Loss):
    def __init__(
        self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "categorical_hinge", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class CosineSimilarity(Loss):
    def __init__(
        self,
        axis: int = -1,
        reduction: _ReductionValues = "sum_over_batch_size",
        name: str | None = "cosine_similarity",
        dtype=None,
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class Hinge(Loss):
    def __init__(self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "hinge", dtype=None) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class Huber(Loss):
    def __init__(
        self, delta: float = 1.0, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "huber_loss", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class KLDivergence(Loss):
    def __init__(
        self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "kl_divergence", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class LogCosh(Loss):
    def __init__(
        self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "log_cosh", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class MeanAbsoluteError(Loss):
    def __init__(
        self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "mean_absolute_error", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class MeanAbsolutePercentageError(Loss):
    def __init__(
        self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "mean_absolute_percentage_error", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class MeanSquaredError(Loss):
    def __init__(
        self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "mean_squared_error", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class MeanSquaredLogarithmicError(Loss):
    def __init__(
        self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "mean_squared_logarithmic_error", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class Poisson(Loss):
    def __init__(self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "poisson", dtype=None) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class SparseCategoricalCrossentropy(Loss):
    def __init__(
        self,
        from_logits: bool = False,
        ignore_class: int | None = None,
        reduction: _ReductionValues = "sum_over_batch_size",
        name: str = "sparse_categorical_crossentropy",
        dtype=None,
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class SquaredHinge(Loss):
    def __init__(
        self, reduction: _ReductionValues = "sum_over_batch_size", name: str | None = "squared_hinge", dtype=None
    ) -> None: ...
    def call(self, y_true: Tensor, y_pred: Tensor) -> Tensor: ...

class Reduction:
    AUTO: Final = "auto"
    NONE: Final = "none"
    SUM: Final = "sum"
    SUM_OVER_BATCH_SIZE: Final = "sum_over_batch_size"
    @classmethod
    def all(cls) -> tuple[_ReductionValues, ...]: ...
    @classmethod
    def validate(cls, key: object) -> TypeGuard[_ReductionValues]: ...

_ReductionValues: TypeAlias = Literal["auto", "none", "sum", "sum_over_batch_size"]

def categorical_hinge(y_true: TensorCompatible, y_pred: TensorCompatible) -> Tensor: ...
def huber(y_true: TensorCompatible, y_pred: TensorCompatible, delta: float = 1.0) -> Tensor: ...
def deserialize(name: str | dict[str, Any], custom_objects: dict[str, Any] | None = None) -> Loss: ...
def serialize(loss: KerasSerializable) -> dict[str, Any]: ...

_FuncT = TypeVar("_FuncT", bound=Callable[..., Any])

@overload
def get(identifier: None) -> None: ...
@overload
def get(identifier: str | dict[str, Any]) -> Loss: ...
@overload
def get(identifier: _FuncT) -> _FuncT: ...

# This is complete with respect to methods documented defined here,
# but many methods get re-exported here from tf.keras.metrics that aren't
# covered yet.
def __getattr__(name: str): ...  # incomplete module
