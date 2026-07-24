# Commonly used type aliases.
# Everything in this module is private for stubs. There is no runtime equivalent.

from collections.abc import Iterable, Mapping, Sequence
from typing import Any, Protocol, TypeVar, type_check_only
from typing_extensions import TypeAlias

import numpy as np
import numpy.typing as npt
import tensorflow as tf
from tensorflow.dtypes import DType
from tensorflow.keras.layers import InputSpec

_T = TypeVar("_T")
ContainerGeneric: TypeAlias = Mapping[str, ContainerGeneric[_T]] | Sequence[ContainerGeneric[_T]] | _T

TensorLike: TypeAlias = tf.Tensor | tf.RaggedTensor | tf.SparseTensor
SparseTensorLike: TypeAlias = tf.Tensor | tf.SparseTensor
RaggedTensorLike: TypeAlias = tf.Tensor | tf.RaggedTensor
# _RaggedTensorLikeT = TypeVar("_RaggedTensorLikeT", tf.Tensor, tf.RaggedTensor)
Gradients: TypeAlias = tf.Tensor | tf.IndexedSlices

@type_check_only
class KerasSerializable1(Protocol):
    def get_config(self) -> dict[str, Any]: ...

@type_check_only
class KerasSerializable2(Protocol):
    __name__: str

KerasSerializable: TypeAlias = KerasSerializable1 | KerasSerializable2

TensorValue: TypeAlias = tf.Tensor  # Alias for a 0D Tensor
Integer: TypeAlias = TensorValue | int | IntArray | np.number[Any]  # Here IntArray are assumed to be 0D.
Float: TypeAlias = Integer | float | FloatArray
Slice: TypeAlias = tf.Tensor | tf.RaggedTensor | int | slice | None
FloatDataSequence: TypeAlias = Sequence[float] | Sequence[FloatDataSequence]
IntDataSequence: TypeAlias = Sequence[int] | Sequence[IntDataSequence]
StrDataSequence: TypeAlias = Sequence[str] | Sequence[StrDataSequence]
DataSequence: TypeAlias = FloatDataSequence | StrDataSequence | IntDataSequence
ScalarTensorCompatible: TypeAlias = tf.Tensor | str | float | np.ndarray[Any, Any] | np.number[Any]
UIntTensorCompatible: TypeAlias = tf.Tensor | int | UIntArray
FloatTensorCompatible: TypeAlias = tf.Tensor | int | IntArray | float | FloatArray | np.number[Any]
StringTensorCompatible: TypeAlias = tf.Tensor | str | npt.NDArray[np.str_] | Sequence[StringTensorCompatible]

TensorCompatible: TypeAlias = ScalarTensorCompatible | Sequence[TensorCompatible]
# _TensorCompatibleT = TypeVar("_TensorCompatibleT", bound=TensorCompatible)
# Sparse tensors are very annoying. Some operations work on them, but many do not.
# You will need to manually verify if an operation supports them. SparseTensorCompatible is intended to be a
# broader type than TensorCompatible and not all operations will support broader version. If unsure,
# use TensorCompatible instead.
SparseTensorCompatible: TypeAlias = TensorCompatible | tf.SparseTensor
# TensorFlow tries to convert anything passed as input. Meaning that even if, for example, only a Tensor of int32
# is allowed, a numpy array of strings that can be converted to int32 will work. Therefore having anything more specific
# then AnyArray might cause false positives, while AnyArray might cause false negatives.
TensorOrArray: TypeAlias = tf.Tensor | AnyArray

ShapeLike: TypeAlias = tf.TensorShape | Iterable[ScalarTensorCompatible | None] | int | tf.Tensor
DTypeLike: TypeAlias = DType | str | np.dtype[Any] | int

ContainerTensors: TypeAlias = ContainerGeneric[tf.Tensor]
ContainerTensorsLike: TypeAlias = ContainerGeneric[TensorLike]
ContainerTensorCompatible: TypeAlias = ContainerGeneric[TensorCompatible]
ContainerGradients: TypeAlias = ContainerGeneric[Gradients]
ContainerTensorShape: TypeAlias = ContainerGeneric[tf.TensorShape]
ContainerInputSpec: TypeAlias = ContainerGeneric[InputSpec]

AnyArray: TypeAlias = npt.NDArray[Any]
FloatArray: TypeAlias = npt.NDArray[np.float16 | np.float32 | np.float64]
UIntArray: TypeAlias = npt.NDArray[np.uint | np.uint8 | np.uint16 | np.uint32 | np.uint64]
SignedIntArray: TypeAlias = npt.NDArray[np.int_ | np.int8 | np.int16 | np.int32 | np.int64]
IntArray: TypeAlias = UIntArray | SignedIntArray
