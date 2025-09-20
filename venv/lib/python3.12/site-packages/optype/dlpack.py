"""
Types and interfaces for DLPack, as used the Array API.
https://github.com/dmlc/dlpack
"""

import enum
import sys
from typing import Any, Protocol

if sys.version_info >= (3, 13):
    from types import CapsuleType
    from typing import TypeVar, runtime_checkable
else:
    from typing_extensions import CapsuleType, TypeVar, runtime_checkable

__all__ = "CanDLPack", "CanDLPackCompat", "CanDLPackDevice"


def __dir__() -> tuple[str, str, str]:
    return __all__


###

_TypeT_co = TypeVar("_TypeT_co", bound=enum.Enum | int, default=int, covariant=True)
_DeviceT_co = TypeVar("_DeviceT_co", bound=int, default=int, covariant=True)

###


class DLDeviceType(enum.IntEnum):
    # GPU device
    CPU = 1
    # Cuda GPU device
    CUDA = 2
    # Pinned CUDA GPU memory vy `cudaMallocHost`
    CPU_PINNED = 3
    # OpenCL devices
    OPENCL = 4
    # Vulkan buffer for next generation graphics
    VULKAN = 7
    # Metal for Apple GPU
    METAL = 8
    # Verilog simulation buffer
    VPI = 9
    # ROCm GPU's for AMD GPU's
    ROCM = 10
    # CUDA managed/unified memory allocated by `cudaMallocManaged`
    CUDA_MANAGED = 13
    # Unified shared memory allocated on a oneAPI non-partititioned
    # device. Call to oneAPI runtime is required to determine the device
    # type, the USM allocation type and the sycl context it is bound to.
    ONE_API = 14


class DLDataTypeCode(enum.IntEnum):
    # signed integer
    INT = 0
    # unsigned integer
    UINT = 1
    # IEEE floating point
    FLOAT = 2
    # Opaque handle type, reserved for testing purposes.
    OPAQUE_HANDLE = 3
    # bfloat16
    BFLOAT = 4
    # complex number (C/C++/Python layout: compact struct per complex number)
    COMPLEX = 5
    # boolean
    BOOL = 6


# NOTE: Because `__dlpack__` doesn't mutate the type, and the type parameters bind to
# the *co*variant `tuple`, they should be *co*variant; NOT *contra*variant!
# (This shows that PEP 695 claim -- variance can always be inferred -- is nonsense.)


# requires numpy>=2.1
@runtime_checkable
class CanDLPack(Protocol[_TypeT_co, _DeviceT_co]):  # type: ignore[misc] # pyright: ignore[reportInvalidTypeVarUse]
    def __dlpack__(
        self,
        /,
        *,
        stream: int | None = None,
        max_version: tuple[int, int] | None = None,
        dl_device: tuple[_TypeT_co, _DeviceT_co] | None = None,
        # NOTE: This should be `bool | None`, but because of an incorrect annotation in
        # `numpy.ndarray.__dlpack__` on `numpy < 2.2.0`, this is not possible.
        copy: bool | Any | None = None,
    ) -> CapsuleType: ...


@runtime_checkable
class CanDLPackCompat(Protocol):
    def __dlpack__(self, /, *, stream: None = None) -> CapsuleType: ...


@runtime_checkable
class CanDLPackDevice(Protocol[_TypeT_co, _DeviceT_co]):
    def __dlpack_device__(self, /) -> tuple[_TypeT_co, _DeviceT_co]: ...
