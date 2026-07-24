import typing_extensions
from typing import TypedDict, type_check_only

from tensorflow.config import PhysicalDevice

@type_check_only
class _MemoryInfo(TypedDict):
    current: int
    peak: int

def get_memory_info(device: str) -> _MemoryInfo: ...
def reset_memory_stats(device: str) -> None: ...
@typing_extensions.deprecated("This function is deprecated in favor of tf.config.experimental.get_memory_info")
def get_memory_usage(device: PhysicalDevice) -> int: ...
def get_memory_growth(device: PhysicalDevice) -> bool: ...
def set_memory_growth(device: PhysicalDevice, enable: bool) -> None: ...
def __getattr__(name: str): ...  # incomplete module
