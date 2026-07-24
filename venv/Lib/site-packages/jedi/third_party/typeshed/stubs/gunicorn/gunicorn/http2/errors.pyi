from typing import Final

class HTTP2ErrorCode:
    NO_ERROR: Final = 0x0
    PROTOCOL_ERROR: Final = 0x1
    INTERNAL_ERROR: Final = 0x2
    FLOW_CONTROL_ERROR: Final = 0x3
    SETTINGS_TIMEOUT: Final = 0x4
    STREAM_CLOSED: Final = 0x5
    FRAME_SIZE_ERROR: Final = 0x6
    REFUSED_STREAM: Final = 0x7
    CANCEL: Final = 0x8
    COMPRESSION_ERROR: Final = 0x9
    CONNECT_ERROR: Final = 0xA
    ENHANCE_YOUR_CALM: Final = 0xB
    INADEQUATE_SECURITY: Final = 0xC
    HTTP_1_1_REQUIRED: Final = 0xD

class HTTP2Error(Exception):
    message: str
    error_code: int

    def __init__(self, message: str | None = None, error_code: int | None = None) -> None: ...

class HTTP2ProtocolError(HTTP2Error): ...
class HTTP2InternalError(HTTP2Error): ...
class HTTP2FlowControlError(HTTP2Error): ...
class HTTP2SettingsTimeout(HTTP2Error): ...
class HTTP2StreamClosed(HTTP2Error): ...
class HTTP2FrameSizeError(HTTP2Error): ...
class HTTP2RefusedStream(HTTP2Error): ...
class HTTP2Cancel(HTTP2Error): ...
class HTTP2CompressionError(HTTP2Error): ...
class HTTP2ConnectError(HTTP2Error): ...
class HTTP2EnhanceYourCalm(HTTP2Error): ...
class HTTP2InadequateSecurity(HTTP2Error): ...
class HTTP2RequiresHTTP11(HTTP2Error): ...

class HTTP2StreamError(HTTP2Error):
    stream_id: int

    def __init__(self, stream_id: int, message: str | None = None, error_code: int | None = None) -> None: ...

class HTTP2ConnectionError(HTTP2Error): ...
class HTTP2ConfigurationError(HTTP2Error): ...

class HTTP2NotAvailable(HTTP2Error):
    def __init__(self, message: str | None = None) -> None: ...

__all__ = [
    "HTTP2ErrorCode",
    "HTTP2Error",
    "HTTP2ProtocolError",
    "HTTP2InternalError",
    "HTTP2FlowControlError",
    "HTTP2SettingsTimeout",
    "HTTP2StreamClosed",
    "HTTP2FrameSizeError",
    "HTTP2RefusedStream",
    "HTTP2Cancel",
    "HTTP2CompressionError",
    "HTTP2ConnectError",
    "HTTP2EnhanceYourCalm",
    "HTTP2InadequateSecurity",
    "HTTP2RequiresHTTP11",
    "HTTP2StreamError",
    "HTTP2ConnectionError",
    "HTTP2ConfigurationError",
    "HTTP2NotAvailable",
]
