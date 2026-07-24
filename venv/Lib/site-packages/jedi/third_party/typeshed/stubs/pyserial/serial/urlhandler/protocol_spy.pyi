from collections.abc import Generator
from typing import TextIO, type_check_only

import serial

def sixteen(data: bytes) -> Generator[tuple[str, str] | tuple[None, None]]: ...
def hexdump(data: bytes) -> Generator[tuple[int, str]]: ...
@type_check_only
class _Formatter:
    def rx(self, data: bytes) -> None: ...
    def tx(self, data: bytes) -> None: ...
    def control(self, name: str, value: str) -> None: ...

class FormatRaw(_Formatter):
    output: TextIO
    color: bool
    rx_color: str
    tx_color: str
    def __init__(self, output: TextIO, color: bool) -> None: ...

class FormatHexdump(_Formatter):
    start_time: float
    output: TextIO
    color: bool
    rx_color: str
    tx_color: str
    control_color: str
    def __init__(self, output: TextIO, color: bool) -> None: ...
    def write_line(self, timestamp: float, label: str, value: str, value2: str = "") -> None: ...

class Serial(serial.Serial):
    formatter: FormatRaw | FormatHexdump | None
    show_all: bool
    def from_url(self, url: str) -> str: ...
