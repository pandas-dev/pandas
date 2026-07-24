from collections.abc import Callable, Sequence
from typing import Final, TypeVar

from tensorflow import Tensor
from tensorflow._aliases import TensorCompatible
from tensorflow.data import Dataset

AUTOTUNE: Final = -1
INFINITE_CARDINALITY: Final = -1
SHARD_HINT: Final = -1
UNKNOWN_CARDINALITY: Final = -2

_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")

def parallel_interleave(
    map_func: Callable[[_T1], Dataset[_T2]],
    cycle_length: int,
    block_length: int = 1,
    sloppy: bool | None = False,
    buffer_output_elements: int | None = None,
    prefetch_input_elements: int | None = None,
) -> Callable[[Dataset[_T1]], Dataset[_T2]]: ...
def enable_debug_mode() -> None: ...
def cardinality(dataset: Dataset[object]) -> Tensor: ...
def sample_from_datasets(
    datasets: Sequence[Dataset[_T1]],
    weights: TensorCompatible | None = None,
    seed: int | None = None,
    stop_on_empty_dataset: bool = False,
) -> Dataset[_T1]: ...
def __getattr__(name: str): ...  # incomplete module
