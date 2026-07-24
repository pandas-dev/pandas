from typing import NamedTuple

from tensorflow.config import experimental as experimental

class PhysicalDevice(NamedTuple):
    name: str
    device_type: str

def list_physical_devices(device_type: None | str = None) -> list[PhysicalDevice]: ...
def get_visible_devices(device_type: None | str = None) -> list[PhysicalDevice]: ...
def set_visible_devices(devices: list[PhysicalDevice] | PhysicalDevice, device_type: None | str = None) -> None: ...
def __getattr__(name: str): ...  # incomplete module
