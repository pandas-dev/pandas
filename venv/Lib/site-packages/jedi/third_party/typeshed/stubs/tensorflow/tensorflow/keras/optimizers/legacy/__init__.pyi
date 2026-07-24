from collections.abc import Callable, Iterable
from typing import Any
from typing_extensions import TypeAlias

import tensorflow as tf
from tensorflow._aliases import Gradients
from tensorflow.keras.optimizers import schedules as schedules
from tensorflow.python.trackable.base import Trackable

_LearningRate: TypeAlias = float | tf.Tensor | schedules.LearningRateSchedule | Callable[[], float | tf.Tensor]
_GradientAggregator: TypeAlias = Callable[[list[tuple[Gradients, tf.Variable]]], list[tuple[Gradients, tf.Variable]]] | None
_GradientTransformer: TypeAlias = (
    Iterable[Callable[[list[tuple[Gradients, tf.Variable]]], list[tuple[Gradients, tf.Variable]]]] | None
)

# kwargs here and in other optimizers can be given better type after Unpack[TypedDict], PEP 692, is supported.
class Optimizer(Trackable):
    _name: str
    _iterations: tf.Variable | None
    _weights: list[tf.Variable]
    gradient_aggregator: _GradientAggregator
    gradient_transformers: _GradientTransformer
    learning_rate: _LearningRate
    def __init__(
        self,
        name: str,
        gradient_aggregator: _GradientAggregator = None,
        gradient_transformers: _GradientTransformer = None,
        **kwargs: Any,
    ) -> None: ...

class Adam(Optimizer):
    def __init__(
        self,
        learning_rate: _LearningRate = 0.001,
        beta_1: float = 0.9,
        beta_2: float = 0.999,
        epsilon: float = 1e-07,
        amsgrad: bool = False,
        name: str = "Adam",
        **kwargs: Any,
    ) -> None: ...

class Adagrad(Optimizer):
    _initial_accumulator_value: float

    def __init__(
        self,
        learning_rate: _LearningRate = 0.001,
        initial_accumulator_value: float = 0.1,
        epsilon: float = 1e-7,
        name: str = "Adagrad",
        **kwargs: Any,
    ) -> None: ...

class SGD(Optimizer):
    def __init__(
        self, learning_rate: _LearningRate = 0.01, momentum: float = 0.0, nesterov: bool = False, name: str = "SGD", **kwargs: Any
    ) -> None: ...

def __getattr__(name: str): ...  # incomplete module
