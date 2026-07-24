from _typeshed import Incomplete
from collections.abc import Callable, Iterable, Sequence
from typing import Any, Generic, Literal, TypeVar, overload, type_check_only
from typing_extensions import Self, TypeAlias

import tensorflow as tf
from tensorflow import Tensor, Variable
from tensorflow._aliases import AnyArray, DataSequence, DTypeLike, Float, TensorCompatible, TensorLike
from tensorflow.keras.activations import _Activation
from tensorflow.keras.constraints import Constraint
from tensorflow.keras.initializers import _Initializer
from tensorflow.keras.regularizers import Regularizer, _Regularizer

_InputT_contra = TypeVar("_InputT_contra", contravariant=True)
_OutputT_co = TypeVar("_OutputT_co", covariant=True)

class InputSpec:
    dtype: str | None
    shape: tuple[int | None, ...]
    ndim: int | None
    max_ndim: int | None
    min_ndim: int | None
    axes: dict[int, int | None] | None
    def __init__(
        self,
        dtype: DTypeLike | None = None,
        shape: Iterable[int | None] | None = None,
        ndim: int | None = None,
        max_ndim: int | None = None,
        min_ndim: int | None = None,
        axes: dict[int, int | None] | None = None,
        allow_last_axis_squeeze: bool = False,
        name: str | None = None,
        optional: bool = False,
    ) -> None: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> type[Self]: ...

# Most layers have input and output type of just Tensor and when we support default type variables,
# maybe worth trying.
class Layer(tf.Module, Generic[_InputT_contra, _OutputT_co]):
    # The most general type is ContainerGeneric[InputSpec] as it really
    # depends on _InputT_contra. For most Layers it is just InputSpec
    # though. Maybe describable with HKT?
    input_spec: InputSpec | Any

    @property
    def trainable(self) -> bool: ...
    @trainable.setter
    def trainable(self, value: bool) -> None: ...
    def __init__(
        self,
        *,
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: DTypeLike | None = None,
        autocast: bool = True,
        name: str | None = None,
        # **kwargs
        input_dim: int | None = None,
        input_shape: Any = None,
    ) -> None: ...

    # *args/**kwargs are allowed, but have obscure footguns and tensorflow documentation discourages their usage.
    # First argument will automatically be cast to layer's compute dtype, but any other tensor arguments will not be.
    # Also various tensorflow tools/apis can misbehave if they encounter a layer with *args/**kwargs.
    def __call__(
        self, inputs: _InputT_contra, *, training: bool = False, mask: TensorCompatible | None = None
    ) -> _OutputT_co: ...
    def call(self, inputs: _InputT_contra, /) -> _OutputT_co: ...

    # input_shape's real type depends on _InputT_contra, but we can't express that without HKT.
    # For example _InputT_contra tf.Tensor -> tf.TensorShape, _InputT_contra dict[str, tf.Tensor] -> dict[str, tf.TensorShape].
    def build(self, input_shape: Any, /) -> None: ...
    @overload
    def compute_output_shape(self: Layer[tf.Tensor, tf.Tensor], input_shape: tf.TensorShape, /) -> tf.TensorShape: ...
    @overload
    def compute_output_shape(self, input_shape: Any, /) -> Any: ...
    def add_weight(
        self,
        shape: Iterable[int | None] | None = None,
        initializer: _Initializer | None = None,
        dtype: DTypeLike | None = None,
        trainable: bool = True,
        autocast: bool = True,
        regularizer: _Regularizer = None,
        constraint: _Constraint = None,
        aggregation: Literal["mean", "sum", "only_first_replica"] = "mean",
        name: str | None = None,
    ) -> tf.Variable: ...
    def add_loss(self, loss: tf.Tensor | Sequence[tf.Tensor] | Callable[[], tf.Tensor]) -> None: ...
    def count_params(self) -> int: ...
    @property
    def trainable_variables(self) -> list[Variable]: ...
    @property
    def non_trainable_variables(self) -> list[Variable]: ...
    @property
    def trainable_weights(self) -> list[Variable]: ...
    @property
    def non_trainable_weights(self) -> list[Variable]: ...
    @property
    def losses(self) -> list[Tensor]: ...
    def get_weights(self) -> list[AnyArray]: ...
    def set_weights(self, weights: Sequence[AnyArray]) -> None: ...
    def get_config(self) -> dict[str, Any]: ...
    @classmethod
    def from_config(cls, config: dict[str, Any]) -> Self: ...
    def __getattr__(self, name: str) -> Incomplete: ...

# Every layer has trainable, dtype, name, and dynamic. At runtime these
# are mainly handled with **kwargs, passed up and then validated.
# In actual implementation there's 12 allowed keyword arguments, but only
# 4 are documented and other 8 are mainly internal. The other 8 can be found
# https://github.com/keras-team/keras/blob/e6784e4302c7b8cd116b74a784f4b78d60e83c26/keras/engine/base_layer.py#L329
# PEP 692 support would be very helpful here and allow removing stubtest allowlist for
# all layer constructors.

# TODO: Replace last Any after adding tf.keras.mixed_precision.Policy.
_LayerDtype: TypeAlias = DTypeLike | dict[str, Any] | Any

_Constraint: TypeAlias = str | dict[str, Any] | Constraint | None

# IndexLookup is not exported by Keras
@type_check_only
class _IndexLookup(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        max_tokens: int | None,
        num_oov_indices: int,
        mask_token: str | None,
        oov_token: str,
        vocabulary_dtype: Literal["int64", "string"],
        vocabulary: str | None | TensorCompatible = None,
        idf_weights: TensorCompatible | None = None,
        invert: bool = False,
        output_mode: Literal["int", "count", "multi_hot", "one_hot", "tf_idf"] = "int",
        sparse: bool = False,
        pad_to_max_tokens: bool = False,
        name: str | None = None,
        *,
        # **kwargs
        vocabulary_size: int | None = None,
        has_input_vocabulary: bool = ...,
        trainable: bool | None = None,
        dtype: _LayerDtype | None = None,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        autocast: bool = True,
    ) -> None: ...
    def compute_output_signature(self, input_spec) -> tf.TensorSpec: ...
    def get_vocabulary(self, include_special_tokens: bool = True) -> list[Incomplete]: ...
    def vocabulary_size(self) -> int: ...

class StringLookup(_IndexLookup):
    def __init__(
        self,
        max_tokens: int | None = None,
        num_oov_indices: int = 1,
        mask_token: str | None = None,
        oov_token: str = "[UNK]",
        vocabulary: str | None | TensorCompatible = None,
        idf_weights: TensorCompatible | None = None,
        invert: bool = False,
        output_mode: Literal["int", "count", "multi_hot", "one_hot", "tf_idf"] = "int",
        pad_to_max_tokens: bool = False,
        sparse: bool = False,
        encoding: str = "utf-8",
        name: str | None = None,
        *,
        # **kwargs passed to IndexLookup
        vocabulary_size: int | None = None,
        has_input_vocabulary: bool = ...,
        trainable: bool | None = None,
        dtype: _LayerDtype | None = None,
        activity_regularizer: _Regularizer = None,
        autocast: bool = True,
    ) -> None: ...
    def adapt(self, data: tf.data.Dataset[TensorLike] | AnyArray | DataSequence, steps: Float | None = None) -> None: ...

class IntegerLookup(_IndexLookup):
    def __init__(
        self,
        max_tokens: int | None = None,
        num_oov_indices: int = 1,
        mask_token: int | None = None,
        oov_token: int = -1,
        vocabulary: str | None | TensorCompatible = None,
        vocabulary_dtype: Literal["int64"] = "int64",
        idf_weights: TensorCompatible | None = None,
        invert: bool = False,
        output_mode: Literal["int", "count", "multi_hot", "one_hot", "tf_idf"] = "int",
        sparse: bool = False,
        pad_to_max_tokens: bool = False,
        name: str | None = None,
        *,
        # **kwargs passed to IndexLookup
        vocabulary_size: int | None = None,
        has_input_vocabulary: bool = ...,
        trainable: bool | None = None,
        dtype: _LayerDtype | None = None,
        activity_regularizer: _Regularizer = None,
        autocast: bool = True,
    ) -> None: ...
    def adapt(self, data: tf.data.Dataset[TensorLike] | AnyArray | DataSequence, steps: Float | None = None) -> None: ...

# Layer's compute_output_shape commonly have instance as first argument name instead of self.
# This is an artifact of actual implementation commonly uses a decorator to define it.
# Layer.build has same weirdness sometimes. For both marked as positional only.
class Dense(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        units: int,
        activation: _Activation = None,
        use_bias: bool = True,
        kernel_initializer: _Initializer = "glorot_uniform",
        bias_initializer: _Initializer = "zeros",
        kernel_regularizer: _Regularizer = None,
        bias_regularizer: _Regularizer = None,
        activity_regularizer: _Regularizer = None,
        kernel_constraint: _Constraint = None,
        bias_constraint: _Constraint = None,
        lora_rank: int | None = None,
        *,
        # **kwargs passed to Layer
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class BatchNormalization(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        axis: int = -1,
        momentum: float = 0.99,
        epsilon: float = 0.001,
        center: bool = True,
        scale: bool = True,
        beta_initializer: _Initializer = "zeros",
        gamma_initializer: _Initializer = "ones",
        moving_mean_initializer: _Initializer = "zeros",
        moving_variance_initializer: _Initializer = "ones",
        beta_regularizer: _Regularizer = None,
        gamma_regularizer: _Regularizer = None,
        beta_constraint: _Constraint = None,
        gamma_constraint: _Constraint = None,
        synchronized: bool = False,
        *,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class ReLU(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        max_value: float | None = None,
        negative_slope: float | None = 0.0,
        threshold: float | None = 0.0,
        *,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class Dropout(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        rate: float,
        noise_shape: TensorCompatible | Sequence[int | None] | None = None,
        seed: int | None = None,
        *,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class Embedding(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        embeddings_initializer: _Initializer = "uniform",
        embeddings_regularizer: _Regularizer = None,
        embeddings_constraint: _Constraint = None,
        mask_zero: bool = False,
        weights=None,
        lora_rank: int | None = None,
        *,
        input_length: int | None = None,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class Conv2D(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        filters: int,
        kernel_size: int | Iterable[int],
        strides: int | Iterable[int] = (1, 1),
        padding: Literal["valid", "same"] = "valid",
        data_format: None | Literal["channels_last", "channels_first"] = None,
        dilation_rate: int | Iterable[int] = (1, 1),
        groups: int = 1,
        activation: _Activation = None,
        use_bias: bool = True,
        kernel_initializer: _Initializer = "glorot_uniform",
        bias_initializer: _Initializer = "zeros",
        kernel_regularizer: _Regularizer = None,
        bias_regularizer: _Regularizer = None,
        activity_regularizer: _Regularizer = None,
        kernel_constraint: _Constraint = None,
        bias_constraint: _Constraint = None,
        *,
        # **kwargs passed to Layer
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

Convolution2D = Conv2D

class Identity(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        *,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class LayerNormalization(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        axis: int = -1,
        epsilon: float = 0.001,
        center: bool = True,
        scale: bool = True,
        rms_scaling: bool = False,
        beta_initializer: _Initializer = "zeros",
        gamma_initializer: _Initializer = "ones",
        beta_regularizer: _Regularizer = None,
        gamma_regularizer: _Regularizer = None,
        beta_constraint: _Constraint = None,
        gamma_constraint: _Constraint = None,
        *,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class MultiHeadAttention(Layer[Any, tf.Tensor]):
    def __init__(
        self,
        num_heads: int,
        key_dim: int | None,
        value_dim: int | None = None,
        dropout: float = 0.0,
        use_bias: bool = True,
        output_shape: tuple[int, ...] | None = None,
        attention_axes: tuple[int, ...] | None = None,
        kernel_initializer: _Initializer = "glorot_uniform",
        bias_initializer: _Initializer = "zeros",
        kernel_regularizer: Regularizer | None = None,
        bias_regularizer: _Regularizer | None = None,
        activity_regularizer: _Regularizer | None = None,
        kernel_constraint: _Constraint | None = None,
        bias_constraint: _Constraint | None = None,
        seed: int | None = None,
        *,
        # **kwargs passed to Layer
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...
    @overload  # type: ignore[override]
    def __call__(
        self,
        query: tf.Tensor,
        value: tf.Tensor,
        key: tf.Tensor | None,
        attention_mask: tf.Tensor | None,
        return_attention_scores: Literal[False],
        training: bool,
        use_causal_mask: bool,
    ) -> tf.Tensor: ...
    @overload
    def __call__(
        self,
        query: tf.Tensor,
        value: tf.Tensor,
        key: tf.Tensor | None,
        attention_mask: tf.Tensor | None,
        return_attention_scores: Literal[True],
        training: bool,
        use_causal_mask: bool,
    ) -> tuple[tf.Tensor, tf.Tensor]: ...
    @overload
    def __call__(
        self,
        query: tf.Tensor,
        value: tf.Tensor,
        key: tf.Tensor | None = None,
        attention_mask: tf.Tensor | None = None,
        return_attention_scores: bool = False,
        training: bool = False,
        use_causal_mask: bool = False,
    ) -> tuple[tf.Tensor, tf.Tensor] | tf.Tensor: ...

class GaussianDropout(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        rate: float,
        seed: int | None = None,
        *,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class Activation(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        activation: _Activation = None,
        *,
        # **kwargs passed to Layer
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class GlobalAveragePooling2D(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        data_format: Literal["channels_last", "channels_first"] | None = None,
        keepdims: bool = False,
        *,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

class MaxPool2D(Layer[tf.Tensor, tf.Tensor]):
    def __init__(
        self,
        pool_size: int | tuple[int, int] = (2, 2),
        strides: int | tuple[int, int] | None = None,
        padding: Literal["valid", "same"] = "valid",
        data_format: Literal["channels_last", "channels_first"] | None = None,
        *,
        # **kwargs passed to Layer
        activity_regularizer: _Regularizer = None,
        trainable: bool = True,
        dtype: _LayerDtype | None = None,
        autocast: bool = True,
        name: str | None = None,
    ) -> None: ...

def __getattr__(name: str): ...  # incomplete module
