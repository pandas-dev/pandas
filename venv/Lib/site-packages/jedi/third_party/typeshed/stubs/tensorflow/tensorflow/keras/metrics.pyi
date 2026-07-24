from _typeshed import Incomplete
from abc import ABCMeta, abstractmethod
from collections.abc import Callable, Iterable, Sequence
from enum import Enum
from typing import Any, Literal, type_check_only
from typing_extensions import Self, TypeAlias

import tensorflow as tf
from tensorflow import Operation, Tensor
from tensorflow._aliases import DTypeLike, KerasSerializable, TensorCompatible
from tensorflow.keras.initializers import _Initializer

_Output: TypeAlias = Tensor | dict[str, Tensor]

class Metric(tf.keras.layers.Layer[tf.Tensor, tf.Tensor], metaclass=ABCMeta):
    def __init__(self, dtype: DTypeLike | None = None, name: str | None = None) -> None: ...
    def __new__(cls, *args: Any, **kwargs: Any) -> Self: ...
    def reset_state(self) -> None: ...
    @abstractmethod
    def update_state(
        self, y_true: TensorCompatible, y_pred: TensorCompatible, sample_weight: TensorCompatible | None = None
    ) -> Operation | None: ...
    @abstractmethod
    def result(self) -> _Output: ...
    # Metric inherits from keras.Layer, but its add_weight method is incompatible with the one from "Layer".
    def add_weight(  # type: ignore[override]
        self,
        name: str,
        shape: Iterable[int | None] | None = (),
        aggregation: tf.VariableAggregation = ...,
        synchronization: tf.VariableSynchronization = ...,
        initializer: _Initializer | None = None,
        dtype: DTypeLike | None = None,
    ) -> None: ...

class AUC(Metric):
    _from_logits: bool
    _num_labels: int
    num_labels: int | None
    def __init__(
        self,
        num_thresholds: int = 200,
        curve: Literal["ROC", "PR"] = "ROC",
        summation_method: Literal["interpolation", "minoring", "majoring"] = "interpolation",
        name: str | None = None,
        dtype: DTypeLike | None = None,
        thresholds: Sequence[float] | None = None,
        multi_label: bool = False,
        num_labels: int | None = None,
        label_weights: TensorCompatible | None = None,
        from_logits: bool = False,
    ) -> None: ...
    def update_state(
        self, y_true: TensorCompatible, y_pred: TensorCompatible, sample_weight: TensorCompatible | None = None
    ) -> Operation: ...
    def result(self) -> tf.Tensor: ...

class Precision(Metric):
    def __init__(
        self,
        thresholds: float | Sequence[float] | None = None,
        top_k: int | None = None,
        class_id: int | None = None,
        name: str | None = None,
        dtype: DTypeLike | None = None,
    ) -> None: ...
    def update_state(
        self, y_true: TensorCompatible, y_pred: TensorCompatible, sample_weight: TensorCompatible | None = None
    ) -> Operation: ...
    def result(self) -> tf.Tensor: ...

class Recall(Metric):
    def __init__(
        self,
        thresholds: float | Sequence[float] | None = None,
        top_k: int | None = None,
        class_id: int | None = None,
        name: str | None = None,
        dtype: DTypeLike | None = None,
    ) -> None: ...
    def update_state(
        self, y_true: TensorCompatible, y_pred: TensorCompatible, sample_weight: TensorCompatible | None = None
    ) -> Operation: ...
    def result(self) -> tf.Tensor: ...

class MeanMetricWrapper(Metric):
    def __init__(
        self, fn: Callable[[tf.Tensor, tf.Tensor], tf.Tensor], name: str | None = None, dtype: DTypeLike | None = None
    ) -> None: ...
    def update_state(
        self, y_true: TensorCompatible, y_pred: TensorCompatible, sample_weight: TensorCompatible | None = None
    ) -> Operation: ...
    def result(self) -> tf.Tensor: ...

class BinaryAccuracy(MeanMetricWrapper):
    def __init__(self, name: str | None = "binary_accuracy", dtype: DTypeLike | None = None, threshold: float = 0.5) -> None: ...

class Accuracy(MeanMetricWrapper):
    def __init__(self, name: str | None = "accuracy", dtype: DTypeLike | None = None) -> None: ...

class CategoricalAccuracy(MeanMetricWrapper):
    def __init__(self, name: str | None = "categorical_accuracy", dtype: DTypeLike | None = None) -> None: ...

class TopKCategoricalAccuracy(MeanMetricWrapper):
    def __init__(self, k: int = 5, name: str | None = "top_k_categorical_accuracy", dtype: DTypeLike | None = None) -> None: ...

class SparseTopKCategoricalAccuracy(MeanMetricWrapper):
    def __init__(
        self, k: int = 5, name: str | None = "sparse_top_k_categorical_accuracy", dtype: DTypeLike | None = None
    ) -> None: ...

class MeanSquaredError(MeanMetricWrapper):
    def __init__(self, name: str | None = "mean_squared_error", dtype: DTypeLike | None = None) -> None: ...

# TODO: Actually tensorflow.python.keras.utils.metrics_utils.Reduction, but that module
# is currently missing from the stub.
@type_check_only
class _Reduction(Enum):
    SUM = "sum"
    SUM_OVER_BATCH_SIZE = "sum_over_batch_size"
    WEIGHTED_MEAN = "weighted_mean"

class Reduce(Metric):
    reduction: _Reduction
    total: Incomplete
    count: Incomplete  # only defined for some reductions
    def __init__(self, reduction: _Reduction, name: str | None, dtype: DTypeLike | None = None) -> None: ...
    def update_state(self, values, sample_weight=None): ...  # type: ignore[override]
    def result(self) -> Tensor: ...

class Mean(Reduce):
    def __init__(self, name: str | None = "mean", dtype: DTypeLike | None = None) -> None: ...

def serialize(metric: KerasSerializable) -> dict[str, Any]: ...
def binary_crossentropy(
    y_true: TensorCompatible, y_pred: TensorCompatible, from_logits: bool = False, label_smoothing: float = 0.0, axis: int = -1
) -> Tensor: ...
def categorical_crossentropy(
    y_true: TensorCompatible, y_pred: TensorCompatible, from_logits: bool = False, label_smoothing: float = 0.0, axis: int = -1
) -> Tensor: ...
def __getattr__(name: str): ...  # incomplete module
