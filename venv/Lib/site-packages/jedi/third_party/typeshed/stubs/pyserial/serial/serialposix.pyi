import sys
from _typeshed import ReadableBuffer
from typing_extensions import Never

from serial.serialutil import SerialBase

class PlatformSpecificBase:
    BAUDRATE_CONSTANTS: dict[int, int]
    def set_low_latency_mode(self, low_latency_settings: bool) -> None: ...

CMSPAR: int
if sys.platform == "linux":
    TCGETS2: int
    TCSETS2: int
    BOTHER: int
    TIOCGRS485: int
    TIOCSRS485: int
    SER_RS485_ENABLED: int
    SER_RS485_RTS_ON_SEND: int
    SER_RS485_RTS_AFTER_SEND: int
    SER_RS485_RX_DURING_TX: int

    class PlatformSpecific(PlatformSpecificBase): ...

elif sys.platform == "cygwin":
    class PlatformSpecific(PlatformSpecificBase): ...

elif sys.platform == "darwin":
    IOSSIOSPEED: int

    class PlatformSpecific(PlatformSpecificBase):
        osx_version: list[str]
        TIOCSBRK: int
        TIOCCBRK: int

else:
    class PlatformSpecific(PlatformSpecificBase): ...

TIOCMGET: int
TIOCMBIS: int
TIOCMBIC: int
TIOCMSET: int
TIOCM_DTR: int
TIOCM_RTS: int
TIOCM_CTS: int
TIOCM_CAR: int
TIOCM_RNG: int
TIOCM_DSR: int
TIOCM_CD: int
TIOCM_RI: int
TIOCINQ: int
TIOCOUTQ: int
TIOCM_zero_str: bytes
TIOCM_RTS_str: bytes
TIOCM_DTR_str: bytes
TIOCSBRK: int
TIOCCBRK: int

class Serial(SerialBase, PlatformSpecific):
    fd: int | None
    pipe_abort_read_w: int | None
    pipe_abort_read_r: int | None
    pipe_abort_write_w: int | None
    pipe_abort_write_r: int | None
    def open(self) -> None: ...
    @property
    def in_waiting(self) -> int: ...
    def read(self, size: int = 1) -> bytes: ...
    def cancel_read(self) -> None: ...
    def cancel_write(self) -> None: ...
    def write(self, b: ReadableBuffer, /) -> int | None: ...
    def reset_input_buffer(self) -> None: ...
    def reset_output_buffer(self) -> None: ...
    def send_break(self, duration: float = 0.25) -> None: ...
    @property
    def cts(self) -> bool: ...
    @property
    def dsr(self) -> bool: ...
    @property
    def ri(self) -> bool: ...
    @property
    def cd(self) -> bool: ...
    @property
    def out_waiting(self) -> int: ...
    def set_input_flow_control(self, enable: bool = True) -> None: ...
    def set_output_flow_control(self, enable: bool = True) -> None: ...
    def nonblocking(self) -> None: ...

class PosixPollSerial(Serial): ...

class VTIMESerial(Serial):
    @property
    def cancel_read(self) -> Never: ...
