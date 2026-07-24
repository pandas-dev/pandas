import sys
from collections.abc import Callable, Mapping, Sequence
from typing import ClassVar, Final
from typing_extensions import TypeAlias

__docformat__: str

paFloat32: Final[int]
paInt32: Final[int]
paInt24: Final[int]
paInt16: Final[int]
paInt8: Final[int]
paUInt8: Final[int]
paCustomFormat: Final[int]

paInDevelopment: Final[int]
paDirectSound: Final[int]
paMME: Final[int]
paASIO: Final[int]
paSoundManager: Final[int]
paCoreAudio: Final[int]
paOSS: Final[int]
paALSA: Final[int]
paAL: Final[int]
paBeOS: Final[int]
paWDMKS: Final[int]
paJACK: Final[int]
paWASAPI: Final[int]
paNoDevice: Final[int]

paNoError: Final[int]
paNotInitialized: Final[int]
paUnanticipatedHostError: Final[int]
paInvalidChannelCount: Final[int]
paInvalidSampleRate: Final[int]
paInvalidDevice: Final[int]
paInvalidFlag: Final[int]
paSampleFormatNotSupported: Final[int]
paBadIODeviceCombination: Final[int]
paInsufficientMemory: Final[int]
paBufferTooBig: Final[int]
paBufferTooSmall: Final[int]
paNullCallback: Final[int]
paBadStreamPtr: Final[int]
paTimedOut: Final[int]
paInternalError: Final[int]
paDeviceUnavailable: Final[int]
paIncompatibleHostApiSpecificStreamInfo: Final[int]
paStreamIsStopped: Final[int]
paStreamIsNotStopped: Final[int]
paInputOverflowed: Final[int]
paOutputUnderflowed: Final[int]
paHostApiNotFound: Final[int]
paInvalidHostApi: Final[int]
paCanNotReadFromACallbackStream: Final[int]
paCanNotWriteToACallbackStream: Final[int]
paCanNotReadFromAnOutputOnlyStream: Final[int]
paCanNotWriteToAnInputOnlyStream: Final[int]
paIncompatibleStreamHostApi: Final[int]

paContinue: Final[int]
paComplete: Final[int]
paAbort: Final[int]

paInputUnderflow: Final[int]
paInputOverflow: Final[int]
paOutputUnderflow: Final[int]
paOutputOverflow: Final[int]
paPrimingOutput: Final[int]

paFramesPerBufferUnspecified: Final[int]

if sys.platform == "darwin":
    class PaMacCoreStreamInfo:
        paMacCoreChangeDeviceParameters: Final[int]
        paMacCoreFailIfConversionRequired: Final[int]
        paMacCoreConversionQualityMin: Final[int]
        paMacCoreConversionQualityMedium: Final[int]
        paMacCoreConversionQualityLow: Final[int]
        paMacCoreConversionQualityHigh: Final[int]
        paMacCoreConversionQualityMax: Final[int]
        paMacCorePlayNice: Final[int]
        paMacCorePro: Final[int]
        paMacCoreMinimizeCPUButPlayNice: Final[int]
        paMacCoreMinimizeCPU: Final[int]
        def __init__(self, flags: int | None = ..., channel_map: _ChannelMap | None = ...) -> None: ...
        def get_flags(self) -> int: ...
        def get_channel_map(self) -> _ChannelMap | None: ...

    _PaMacCoreStreamInfo: TypeAlias = PaMacCoreStreamInfo
else:
    _PaMacCoreStreamInfo: TypeAlias = None

# Auxiliary types
_ChannelMap: TypeAlias = Sequence[int]
_PaHostApiInfo: TypeAlias = Mapping[str, str | int]
_PaDeviceInfo: TypeAlias = Mapping[str, str | int | float]
_StreamCallback: TypeAlias = Callable[[bytes | None, int, Mapping[str, float], int], tuple[bytes | None, int]]

def get_format_from_width(width: int, unsigned: bool = ...) -> int: ...
def get_portaudio_version() -> int: ...
def get_portaudio_version_text() -> str: ...
def get_sample_size(format: int) -> int: ...

class Stream:
    def __init__(
        self,
        PA_manager: PyAudio,
        rate: int,
        channels: int,
        format: int,
        input: bool = ...,
        output: bool = ...,
        input_device_index: int | None = ...,
        output_device_index: int | None = ...,
        frames_per_buffer: int = ...,
        start: bool = ...,
        input_host_api_specific_stream_info: _PaMacCoreStreamInfo | None = ...,
        output_host_api_specific_stream_info: _PaMacCoreStreamInfo | None = ...,
        stream_callback: _StreamCallback | None = ...,
    ) -> None: ...
    def close(self) -> None: ...
    def get_cpu_load(self) -> float: ...
    def get_input_latency(self) -> float: ...
    def get_output_latency(self) -> float: ...
    def get_read_available(self) -> int: ...
    def get_time(self) -> float: ...
    def get_write_available(self) -> int: ...
    def is_active(self) -> bool: ...
    def is_stopped(self) -> bool: ...
    def read(self, num_frames: int, exception_on_overflow: bool = ...) -> bytes: ...
    def start_stream(self) -> None: ...
    def stop_stream(self) -> None: ...
    def write(self, frames: bytes, num_frames: int | None = ..., exception_on_underflow: bool = ...) -> None: ...

# Use an alias to workaround pyright complaints about recursive definitions in the PyAudio class
_Stream = Stream

class PyAudio:
    Stream: ClassVar[type[_Stream]]
    def __init__(self) -> None: ...
    def close(self, stream: _Stream) -> None: ...
    def get_default_host_api_info(self) -> _PaHostApiInfo: ...
    def get_default_input_device_info(self) -> _PaDeviceInfo: ...
    def get_default_output_device_info(self) -> _PaDeviceInfo: ...
    def get_device_count(self) -> int: ...
    def get_device_info_by_host_api_device_index(self, host_api_index: int, host_api_device_index: int) -> _PaDeviceInfo: ...
    def get_device_info_by_index(self, device_index: int) -> _PaDeviceInfo: ...
    def get_format_from_width(self, width: int, unsigned: bool = ...) -> int: ...
    def get_host_api_count(self) -> int: ...
    def get_host_api_info_by_index(self, host_api_index: int) -> _PaHostApiInfo: ...
    def get_host_api_info_by_type(self, host_api_type: int) -> _PaHostApiInfo: ...
    def get_sample_size(self, format: int) -> int: ...
    def is_format_supported(
        self,
        rate: int,
        input_device: int | None = ...,
        input_channels: int | None = ...,
        input_format: int | None = ...,
        output_device: int | None = ...,
        output_channels: int | None = ...,
        output_format: int | None = ...,
    ) -> bool: ...
    def open(
        self,
        rate: int,
        channels: int,
        format: int,
        input: bool = ...,
        output: bool = ...,
        input_device_index: int | None = ...,
        output_device_index: int | None = ...,
        frames_per_buffer: int = ...,
        start: bool = ...,
        input_host_api_specific_stream_info: _PaMacCoreStreamInfo | None = ...,
        output_host_api_specific_stream_info: _PaMacCoreStreamInfo | None = ...,
        stream_callback: _StreamCallback | None = ...,
    ) -> _Stream: ...
    def terminate(self) -> None: ...
