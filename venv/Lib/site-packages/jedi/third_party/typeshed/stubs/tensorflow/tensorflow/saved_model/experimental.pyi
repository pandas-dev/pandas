from enum import Enum
from typing_extensions import Self

import tensorflow as tf
from tensorflow._aliases import Integer, TensorValue
from tensorflow.python.trackable.resource import CapturableResource

class Fingerprint:
    saved_model_checksum: TensorValue | None
    graph_def_program_hash: TensorValue | None = None
    signature_def_hash: TensorValue | None = None
    saved_object_graph_hash: TensorValue | None = None
    checkpoint_hash: TensorValue | None = None
    version: TensorValue | None = None
    # In practice it seems like any type is accepted, but that might cause issues later on.
    def __init__(
        self,
        saved_model_checksum: Integer | None = None,
        graph_def_program_hash: Integer | None = None,
        signature_def_hash: Integer | None = None,
        saved_object_graph_hash: Integer | None = None,
        checkpoint_hash: Integer | None = None,
        version: Integer | None = None,
    ) -> None: ...
    @classmethod
    def from_proto(cls, proto) -> Self: ...
    def singleprint(self) -> str: ...

class TrackableResource(CapturableResource):
    @property
    def resource_handle(self) -> tf.Tensor: ...
    def __init__(self, device: str = "") -> None: ...

class VariablePolicy(Enum):
    EXPAND_DISTRIBUTED_VARIABLES = "expand_distributed_variables"
    NONE = None
    SAVE_VARIABLE_DEVICES = "save_variable_devices"

def read_fingerprint(export_dir: str) -> Fingerprint: ...
