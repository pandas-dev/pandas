from _typeshed import Incomplete
from collections.abc import Callable
from typing import Any, TypeVar
from typing_extensions import Self

import numpy as np
import tensorflow as tf
from tensorflow.core.example.example_pb2 import Example as Example
from tensorflow.core.example.feature_pb2 import (
    BytesList as BytesList,
    Feature as Feature,
    Features as Features,
    FloatList as FloatList,
    Int64List as Int64List,
)
from tensorflow.core.protobuf.cluster_pb2 import ClusterDef as ClusterDef
from tensorflow.core.protobuf.tensorflow_server_pb2 import ServerDef as ServerDef
from tensorflow.python.trackable.base import Trackable
from tensorflow.python.training.tracking.autotrackable import AutoTrackable

class CheckpointOptions:
    __slots__ = (
        "experimental_io_device",
        "experimental_enable_async_checkpoint",
        "experimental_write_callbacks",
        "enable_async",
        "experimental_sharding_callback",
        "experimental_skip_slot_variables",
    )
    experimental_io_device: None | str
    experimental_enable_async_checkpoint: bool
    experimental_write_callbacks: None | list[Callable[[str], object] | Callable[[], object]]
    enable_async: bool
    experimental_sharding_callback: Incomplete  # should be ShardingCallback
    experimental_skip_slot_variables: bool

    def __init__(
        self,
        experimental_io_device: None | str = None,
        experimental_enable_async_checkpoint: bool = False,
        experimental_write_callbacks: None | list[Callable[[str], object] | Callable[[], object]] = None,
        enable_async: bool = False,
        experimental_skip_slot_variables: bool = False,
        experimental_sharding_callback=None,
    ) -> None: ...

_T = TypeVar("_T", bound=list[str] | tuple[str] | dict[int, str])

class ClusterSpec:
    def __init__(self, cluster: dict[str, _T] | ClusterDef | ClusterSpec) -> None: ...
    def as_dict(self) -> dict[str, list[str] | tuple[str] | dict[int, str]]: ...
    def num_tasks(self, job_name: str) -> int: ...

class _CheckpointLoadStatus:
    def assert_consumed(self) -> Self: ...
    def assert_existing_objects_matched(self) -> Self: ...
    def assert_nontrivial_match(self) -> Self: ...
    def expect_partial(self) -> Self: ...

class Checkpoint(AutoTrackable):
    def __init__(self, root: Trackable | None = None, **kwargs: Trackable) -> None: ...
    def read(self, save_path: str, options: CheckpointOptions | None = None) -> _CheckpointLoadStatus: ...
    def restore(self, save_path: str, options: CheckpointOptions | None = None) -> _CheckpointLoadStatus: ...
    def save(self, file_prefix: str, options: CheckpointOptions | None = None) -> str: ...
    def sync(self) -> None: ...
    def write(self, file_prefix: str, options: CheckpointOptions | None = None) -> str: ...

class CheckpointManager:
    def __init__(
        self,
        checkpoint: Checkpoint,
        directory: str,
        max_to_keep: int,
        keep_checkpoint_every_n_hours: int | None = None,
        checkpoint_name: str = "ckpt",
        step_counter: tf.Variable | None = None,
        checkpoint_interval: int | None = None,
        init_fn: Callable[[], object] | None = None,
    ) -> None: ...
    def _sweep(self) -> None: ...

def latest_checkpoint(checkpoint_dir: str, latest_filename: str | None = None) -> str: ...
def load_variable(ckpt_dir_or_file: str, name: str) -> np.ndarray[Any, Any]: ...
def list_variables(ckpt_dir_or_file: str) -> list[tuple[str, list[int]]]: ...
def __getattr__(name: str): ...  # incomplete module
