from _typeshed import Incomplete
from builtins import bool as _bool
from collections.abc import Generator, Mapping, Sequence
from contextlib import contextmanager
from types import TracebackType
from typing import overload
from typing_extensions import Self

import tensorflow as tf
from tensorflow import Tensor, UnconnectedGradients, Variable
from tensorflow._aliases import ContainerGradients, ContainerTensors, ContainerTensorsLike, Gradients, TensorLike

class ForwardAccumulator:
    def __init__(self, primals: Tensor, tangents: Tensor) -> None: ...
    def jvp(
        self, primals: Tensor, unconnected_gradients: tf.UnconnectedGradients = tf.UnconnectedGradients.NONE  # noqa: Y011
    ) -> Tensor | None: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, typ: type[BaseException] | None, value: BaseException | None, traceback: TracebackType | None) -> None: ...

class GradientTape:
    def __init__(self, persistent: _bool = False, watch_accessed_variables: _bool = True) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, typ: type[BaseException] | None, value: BaseException | None, traceback: TracebackType | None) -> None: ...
    # Higher kinded types would be nice here and these overloads are a way to simulate some of them.
    @overload
    def gradient(
        self,
        target: ContainerTensors,
        sources: TensorLike,
        output_gradients: list[Tensor] | None = None,
        unconnected_gradients: UnconnectedGradients = ...,
    ) -> Gradients: ...
    @overload
    def gradient(
        self,
        target: ContainerTensors,
        sources: Sequence[Tensor],
        output_gradients: list[Tensor] | None = None,
        unconnected_gradients: UnconnectedGradients = ...,
    ) -> list[Gradients]: ...
    @overload
    def gradient(
        self,
        target: ContainerTensors,
        sources: Mapping[str, Tensor],
        output_gradients: list[Tensor] | None = None,
        unconnected_gradients: UnconnectedGradients = ...,
    ) -> dict[str, Gradients]: ...
    @overload
    def gradient(
        self,
        target: ContainerTensors,
        sources: ContainerTensors,
        output_gradients: list[Tensor] | None = None,
        unconnected_gradients: UnconnectedGradients = ...,
    ) -> ContainerGradients: ...
    @contextmanager
    def stop_recording(self) -> Generator[None]: ...
    def reset(self) -> None: ...
    def watch(self, tensor: ContainerTensorsLike) -> None: ...
    def watched_variables(self) -> tuple[Variable, ...]: ...
    def __getattr__(self, name: str) -> Incomplete: ...
