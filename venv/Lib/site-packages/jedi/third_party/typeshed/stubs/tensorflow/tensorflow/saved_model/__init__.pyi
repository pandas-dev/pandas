from _typeshed import Incomplete
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Generic, Literal, TypeVar
from typing_extensions import ParamSpec, TypeAlias

import tensorflow as tf
from tensorflow.python.training.tracking.autotrackable import AutoTrackable
from tensorflow.saved_model.experimental import VariablePolicy
from tensorflow.types.experimental import ConcreteFunction, PolymorphicFunction

_P = ParamSpec("_P")
_R_co = TypeVar("_R_co", covariant=True)

class Asset:
    @property
    def asset_path(self) -> tf.Tensor: ...
    def __init__(self, path: str | Path | tf.Tensor) -> None: ...

class LoadOptions:
    __slots__ = (
        "allow_partial_checkpoint",
        "experimental_io_device",
        "experimental_skip_checkpoint",
        "experimental_variable_policy",
        "experimental_load_function_aliases",
    )
    allow_partial_checkpoint: bool
    experimental_io_device: str | None
    experimental_skip_checkpoint: bool
    experimental_variable_policy: VariablePolicy | None
    experimental_load_function_aliases: bool

    def __init__(
        self,
        allow_partial_checkpoint: bool = False,
        experimental_io_device: str | None = None,
        experimental_skip_checkpoint: bool = False,
        experimental_variable_policy: (
            VariablePolicy | Literal["expand_distributed_variables", "save_variable_devices"] | None
        ) = None,
        experimental_load_function_aliases: bool = False,
    ) -> None: ...

class SaveOptions:
    __slots__ = (
        "namespace_whitelist",
        "save_debug_info",
        "function_aliases",
        "experimental_debug_stripper",
        "experimental_io_device",
        "experimental_variable_policy",
        "experimental_custom_gradients",
        "experimental_image_format",
        "experimental_skip_saver",
        "experimental_sharding_callback",
        "extra_tags",
    )
    namespace_whitelist: list[str]
    save_debug_info: bool
    function_aliases: dict[str, PolymorphicFunction[..., object]]
    experimental_debug_stripper: bool
    experimental_io_device: str
    experimental_variable_policy: VariablePolicy
    experimental_custom_gradients: bool
    experimental_image_format: bool
    experimental_skip_saver: bool
    experimental_sharding_callback: Incomplete | None
    extra_tags: Incomplete | None
    def __init__(
        self,
        namespace_whitelist: list[str] | None = None,
        save_debug_info: bool = False,
        function_aliases: Mapping[str, PolymorphicFunction[..., object]] | None = None,
        experimental_debug_stripper: bool = False,
        experimental_io_device: str | None = None,
        experimental_variable_policy: str | VariablePolicy | None = None,
        experimental_custom_gradients: bool = True,
        experimental_image_format: bool = False,
        experimental_skip_saver: bool = False,
        experimental_sharding_callback=None,
        extra_tags=None,
    ) -> None: ...

def contains_saved_model(export_dir: str | Path) -> bool: ...

class _LoadedAttributes(Generic[_P, _R_co]):
    signatures: Mapping[str, ConcreteFunction[_P, _R_co]]

class _LoadedModel(AutoTrackable, _LoadedAttributes[_P, _R_co]):
    variables: list[tf.Variable]
    trainable_variables: list[tf.Variable]
    # TF1 model artifact specific
    graph: tf.Graph

def load(
    export_dir: str, tags: str | Sequence[str] | None = None, options: LoadOptions | None = None
) -> _LoadedModel[..., Any]: ...

_TF_Function: TypeAlias = ConcreteFunction[..., object] | PolymorphicFunction[..., object]

def save(
    obj: tf.Module,
    export_dir: str,
    signatures: _TF_Function | Mapping[str, _TF_Function] | None = None,
    options: SaveOptions | None = None,
) -> None: ...

ASSETS_DIRECTORY: str = "assets"
ASSETS_KEY: str = "saved_model_assets"
CLASSIFY_INPUTS: str = "inputs"
CLASSIFY_METHOD_NAME: str = "tensorflow/serving/classify"
CLASSIFY_OUTPUT_CLASSES: str = "classes"
CLASSIFY_OUTPUT_SCORES: str = "scores"
DEBUG_DIRECTORY: str = "debug"
DEBUG_INFO_FILENAME_PB: str = "saved_model_debug_info.pb"
DEFAULT_SERVING_SIGNATURE_DEF_KEY: str = "serving_default"
GPU: str = "gpu"
PREDICT_INPUTS: str = "inputs"
PREDICT_METHOD_NAME: str = "tensorflow/serving/predict"
PREDICT_OUTPUTS: str = "outputs"
REGRESS_INPUTS: str = "inputs"
REGRESS_METHOD_NAME: str = "tensorflow/serving/regress"
REGRESS_OUTPUTS: str = "outputs"
SAVED_MODEL_FILENAME_PB: str = "saved_model.pb"
SAVED_MODEL_FILENAME_PBTXT: str = "saved_model.pbtxt"
SAVED_MODEL_SCHEMA_VERSION: int = 1
SERVING: str = "serve"
TPU: str = "tpu"
TRAINING: str = "train"
VARIABLES_DIRECTORY: str = "variables"
VARIABLES_FILENAME: str = "variables"
