import sys

from serial.serialutil import *

if sys.platform == "win32":
    from serial.serialwin32 import Serial as Serial
else:
    from serial.serialposix import PosixPollSerial as PosixPollSerial, Serial as Serial, VTIMESerial as VTIMESerial
# TODO: java? cli? These platforms raise flake8-pyi Y008. Should they be included with a noqa?

__version__: str
VERSION: str
protocol_handler_packages: list[str]

def serial_for_url(
    url: str | None,
    baudrate: int = ...,
    bytesize: int = ...,
    parity: str = ...,
    stopbits: float = ...,
    timeout: float | None = ...,
    xonxoff: bool = ...,
    rtscts: bool = ...,
    write_timeout: float | None = ...,
    dsrdtr: bool = ...,
    inter_byte_timeout: float | None = ...,
    exclusive: float | None = ...,
    *,
    do_not_open: bool = ...,
) -> Serial: ...
