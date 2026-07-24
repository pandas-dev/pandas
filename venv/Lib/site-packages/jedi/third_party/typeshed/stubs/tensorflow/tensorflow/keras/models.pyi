from _typeshed import Incomplete
from collections.abc import Callable, Container, Iterator
from pathlib import Path
from typing import Any, Literal
from typing_extensions import Self, TypeAlias, deprecated

import numpy as np
import numpy.typing as npt
import tensorflow as tf
from tensorflow import Variable
from tensorflow._aliases import ContainerGeneric, ShapeLike, TensorCompatible
from tensorflow.keras.layers import Layer, _InputT_contra, _OutputT_co
from tensorflow.keras.optimizers import Optimizer

_Loss: TypeAlias = str | tf.keras.losses.Loss | Callable[[TensorCompatible, TensorCompatible], tf.Tensor]
_Metric: TypeAlias = str | tf.keras.metrics.Metric | Callable[[TensorCompatible, TensorCompatible], tf.Tensor] | None

# Missing keras.src.backend.tensorflow.trainer.TensorFlowTrainer as a base class, which is not exposed by tensorflow
class Model(Layer[_InputT_contra, _OutputT_co]):
    _train_counter: tf.Variable
    _test_counter: tf.Variable
    optimizer: Optimizer | None
    # This is actually TensorFlowTrainer.loss
    @deprecated("Instead, use `model.compute_loss(x, y, y_pred, sample_weight)`.")
    def loss(self, y: TensorCompatible | None, y_pred: TensorCompatible | None, sample_weight=None) -> tf.Tensor | None: ...
    stop_training: bool

    def __new__(cls, *args: Any, **kwargs: Any) -> Model[_InputT_contra, _OutputT_co]: ...
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    def __setattr__(self, name: str, value: Any) -> None: ...
    def __reduce__(self): ...
    def build(self, input_shape: ShapeLike) -> None: ...
    def __call__(
        self, inputs: _InputT_contra, *, training: bool = False, mask: TensorCompatible | None = None
    ) -> _OutputT_co: ...
    def call(self, inputs: _InputT_contra, training: bool | None = None, mask: TensorCompatible | None = None) -> _OutputT_co: ...
    # Ideally loss/metrics/output would share the same structure but higher kinded types are not supported.
    def compile(
        self,
        optimizer: Optimizer | str = "rmsprop",
        loss: ContainerGeneric[_Loss] | None = None,
        loss_weights: ContainerGeneric[float] | None = None,
        metrics: ContainerGeneric[_Metric] | None = None,
        weighted_metrics: ContainerGeneric[_Metric] | None = None,
        run_eagerly: bool = False,
        steps_per_execution: int | Literal["auto"] = 1,
        jit_compile: bool | Literal["auto"] = "auto",
        auto_scale_loss: bool | None = True,
    ) -> None: ...
    @property
    def metrics(self) -> list[Incomplete]: ...
    @property
    def metrics_names(self) -> list[str]: ...
    @property
    def distribute_strategy(self) -> tf.distribute.Strategy: ...
    @property
    def run_eagerly(self) -> bool: ...
    @property
    def jit_compile(self) -> bool: ...
    @property
    def distribute_reduction_method(self) -> Incomplete | Literal["auto"]: ...
    def train_step(self, data: TensorCompatible): ...
    def compute_loss(
        self,
        x: TensorCompatible | None = None,
        y: TensorCompatible | None = None,
        y_pred: TensorCompatible | None = None,
        sample_weight=None,
        training: bool = True,
    ) -> tf.Tensor | None: ...
    def compute_metrics(
        self, x: TensorCompatible, y: TensorCompatible, y_pred: TensorCompatible, sample_weight=None
    ) -> dict[str, float]: ...
    def get_metrics_result(self) -> dict[str, float]: ...
    def make_train_function(self, force: bool = False) -> Callable[[tf.data.Iterator[Incomplete]], dict[str, float]]: ...
    def fit(
        self,
        x: TensorCompatible | dict[str, TensorCompatible] | tf.data.Dataset[Incomplete] | None = None,
        y: TensorCompatible | dict[str, TensorCompatible] | tf.data.Dataset[Incomplete] | None = None,
        batch_size: int | None = None,
        epochs: int = 1,
        verbose: Literal["auto", 0, 1, 2] = "auto",
        callbacks: list[tf.keras.callbacks.Callback] | None = None,
        validation_split: float = 0.0,
        validation_data: TensorCompatible | tf.data.Dataset[Any] | None = None,
        shuffle: bool = True,
        class_weight: dict[int, float] | None = None,
        sample_weight: npt.NDArray[np.float64] | None = None,
        initial_epoch: int = 0,
        steps_per_epoch: int | None = None,
        validation_steps: int | None = None,
        validation_batch_size: int | None = None,
        validation_freq: int | Container[int] = 1,
    ) -> tf.keras.callbacks.History: ...
    def test_step(self, data: TensorCompatible) -> dict[str, float]: ...
    def make_test_function(self, force: bool = False) -> Callable[[tf.data.Iterator[Incomplete]], dict[str, float]]: ...
    def evaluate(
        self,
        x: TensorCompatible | dict[str, TensorCompatible] | tf.data.Dataset[Incomplete] | None = None,
        y: TensorCompatible | dict[str, TensorCompatible] | tf.data.Dataset[Incomplete] | None = None,
        batch_size: int | None = None,
        verbose: Literal["auto", 0, 1, 2] = "auto",
        sample_weight: npt.NDArray[np.float64] | None = None,
        steps: int | None = None,
        callbacks: list[tf.keras.callbacks.Callback] | None = None,
        return_dict: bool = False,
        **kwargs: Any,
    ) -> float | list[float]: ...
    def predict_step(self, data: _InputT_contra) -> _OutputT_co: ...
    def make_predict_function(self, force: bool = False) -> Callable[[tf.data.Iterator[Incomplete]], _OutputT_co]: ...
    def predict(
        self,
        x: TensorCompatible | tf.data.Dataset[Incomplete],
        batch_size: int | None = None,
        verbose: Literal["auto", 0, 1, 2] = "auto",
        steps: int | None = None,
        callbacks: list[tf.keras.callbacks.Callback] | None = None,
    ) -> _OutputT_co: ...
    def reset_metrics(self) -> None: ...
    def train_on_batch(
        self,
        x: TensorCompatible | dict[str, TensorCompatible] | tf.data.Dataset[Incomplete],
        y: TensorCompatible | dict[str, TensorCompatible] | tf.data.Dataset[Incomplete] | None = None,
        sample_weight: npt.NDArray[np.float64] | None = None,
        class_weight: dict[int, float] | None = None,
        return_dict: bool = False,
    ) -> float | list[float]: ...
    def test_on_batch(
        self,
        x: TensorCompatible | dict[str, TensorCompatible] | tf.data.Dataset[Incomplete],
        y: TensorCompatible | dict[str, TensorCompatible] | tf.data.Dataset[Incomplete] | None = None,
        sample_weight: npt.NDArray[np.float64] | None = None,
        return_dict: bool = False,
    ) -> float | list[float]: ...
    def predict_on_batch(self, x: Iterator[_InputT_contra]) -> npt.NDArray[Incomplete]: ...
    @property
    def trainable_weights(self) -> list[Variable]: ...
    @property
    def non_trainable_weights(self) -> list[Variable]: ...
    def get_weights(self): ...
    def save(self, filepath: str | Path, overwrite: bool = True, zipped: bool | None = None) -> None: ...
    def save_weights(self, filepath: str | Path, overwrite: bool = True) -> None: ...
    # kwargs are from keras.saving.saving_api.load_weights
    def load_weights(self, filepath: str | Path, skip_mismatch: bool = False, *, by_name: bool = False) -> None: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any], custom_objects=None) -> Self: ...
    def to_json(self, **kwargs: Any) -> str: ...
    @property
    def weights(self) -> list[Variable]: ...
    def summary(
        self,
        line_length: None | int = None,
        positions: None | list[float] = None,
        print_fn: None | Callable[[str], None] = None,
        expand_nested: bool = False,
        show_trainable: bool = False,
        layer_range: None | list[str] | tuple[str, str] = None,
    ) -> None: ...
    @property
    def layers(self) -> list[Layer[Incomplete, Incomplete]]: ...
    def get_layer(self, name: str | None = None, index: int | None = None) -> Layer[Incomplete, Incomplete]: ...
    def get_compile_config(self) -> dict[str, Any]: ...
    def compile_from_config(self, config: dict[str, Any]) -> Self: ...
    def export(self, filepath: str | Path, format: str = "tf_saved_model", verbose: bool = True) -> None: ...

def __getattr__(name: str): ...  # incomplete module
