from abc import abstractmethod
from collections.abc import Sequence
from typing import Any
from typing_extensions import Self

import tensorflow as tf

class LearningRateSchedule:
    # At runtime these methods are abstract even though class is not ABC.
    @abstractmethod
    def __call__(self, step: int | tf.Tensor) -> float | tf.Tensor: ...
    @abstractmethod
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...

class PiecewiseConstantDecay(LearningRateSchedule):
    def __init__(
        self,
        boundaries: Sequence[tf.Tensor] | Sequence[float],
        values: Sequence[float] | Sequence[tf.Tensor],
        name: str = "PiecewiseConstant",
    ) -> None: ...
    def __call__(self, step: int | tf.Tensor) -> float | tf.Tensor: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...

class InverseTimeDecay(LearningRateSchedule):
    def __init__(
        self,
        initial_learning_rate: float | tf.Tensor,
        decay_steps: int,
        decay_rate: float,
        staircase: bool = False,
        name: str = "InverseTimeDecay",
    ) -> None: ...
    def __call__(self, step: int | tf.Tensor) -> float | tf.Tensor: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...

class PolynomialDecay(LearningRateSchedule):
    def __init__(
        self,
        initial_learning_rate: float | tf.Tensor,
        decay_steps: int,
        end_learning_rate: float | tf.Tensor = 0.0001,
        power: float = 1.0,
        cycle: bool = False,
        name: str = "PolynomialDecay",
    ) -> None: ...
    def __call__(self, step: int | tf.Tensor) -> float | tf.Tensor: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...

class CosineDecay(LearningRateSchedule):
    def __init__(
        self,
        initial_learning_rate: float | tf.Tensor,
        decay_steps: int,
        alpha: float | tf.Tensor = 0.0,
        name: str = "CosineDecay",
        warmup_target: int | tf.Tensor | None = None,  # float32 or float64 Tensor
        warmup_steps: int | tf.Tensor = 0,  # int32 or int64 Tensor
    ) -> None: ...
    def __call__(self, step: int | tf.Tensor) -> float | tf.Tensor: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...

class CosineDecayRestarts(LearningRateSchedule):
    def __init__(
        self,
        initial_learning_rate: float | tf.Tensor,
        first_decay_steps: int | tf.Tensor,
        t_mul: float | tf.Tensor = 2.0,
        m_mul: float | tf.Tensor = 1.0,
        alpha: float | tf.Tensor = 0.0,
        name: str = "SGDRDecay",
    ) -> None: ...
    def __call__(self, step: int | tf.Tensor) -> float | tf.Tensor: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...

class ExponentialDecay(LearningRateSchedule):
    def __init__(
        self,
        initial_learning_rate: float | tf.Tensor,
        decay_steps: int | tf.Tensor,
        decay_rate: float | tf.Tensor,
        staircase: bool = False,
        name: str = "ExponentialDecay",
    ) -> None: ...
    def __call__(self, step: int | tf.Tensor) -> float | tf.Tensor: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...

def deserialize(config: dict[str, Any], custom_objects: dict[str, type] | None = None) -> LearningRateSchedule: ...
def serialize(learning_rate_schedule: LearningRateSchedule) -> dict[str, Any]: ...
