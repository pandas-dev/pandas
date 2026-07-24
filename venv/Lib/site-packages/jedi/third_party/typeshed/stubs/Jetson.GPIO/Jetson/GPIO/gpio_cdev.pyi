import ctypes
from dataclasses import dataclass
from typing import Final, Literal

from .gpio_pin_data import ChannelInfo

GPIO_HIGH: Final = 1

GPIOHANDLE_REQUEST_INPUT: Final = 0x1
GPIOHANDLE_REQUEST_OUTPUT: Final = 0x2

GPIOEVENT_REQUEST_RISING_EDGE: Final = 0x1
GPIOEVENT_REQUEST_FALLING_EDGE: Final = 0x2
GPIOEVENT_REQUEST_BOTH_EDGES: Final = 0x3

GPIO_GET_CHIPINFO_IOCTL: Final = 0x8044B401
GPIO_GET_LINEINFO_IOCTL: Final = 0xC048B402
GPIO_GET_LINEHANDLE_IOCTL: Final = 0xC16CB403
GPIOHANDLE_GET_LINE_VALUES_IOCTL: Final = 0xC040B408
GPIOHANDLE_SET_LINE_VALUES_IOCTL: Final = 0xC040B409
GPIO_GET_LINEEVENT_IOCTL: Final = 0xC030B404

class gpiochip_info(ctypes.Structure):
    name: str
    label: str
    lines: int

class gpiohandle_request(ctypes.Structure):
    lineoffsets: list[int]
    flags: int
    default_values: list[int]
    consumer_label: str
    lines: int
    fd: int

class gpiohandle_data(ctypes.Structure):
    values: list[int]

class gpioline_info(ctypes.Structure):
    line_offset: int
    flags: int
    name: str
    consumer: str

class gpioline_info_changed(ctypes.Structure):
    line_info: gpioline_info
    timestamp: int
    event_type: int
    padding: list[int]

class gpioevent_request(ctypes.Structure):
    lineoffset: int
    handleflags: int
    eventflags: int
    consumer_label: str
    fd: int

class gpioevent_data(ctypes.Structure):
    timestamp: int
    id: int

class GPIOError(IOError): ...

def chip_open(gpio_chip: str) -> int: ...
def chip_check_info(label: str, gpio_device: str) -> int | None: ...
def chip_open_by_label(label: str) -> int: ...
def close_chip(chip_fd: int) -> None: ...
def open_line(ch_info: ChannelInfo, request: int) -> None: ...
def close_line(line_handle: int) -> None: ...
def request_handle(line_offset: int, direction: Literal[0, 1], initial: Literal[0, 1], consumer: str) -> gpiohandle_request: ...
def request_event(line_offset: int, edge: int, consumer: str) -> gpioevent_request: ...
def get_value(line_handle: int) -> int: ...
def set_value(line_handle: int, value: int) -> None: ...
@dataclass
class PadCtlRegister:
    is_gpio: bool
    is_input: bool
    is_tristate: bool
    def __init__(self, value: int) -> None: ...
    @property
    def is_bidi(self) -> bool: ...

def check_pinmux(ch_info: ChannelInfo, direction: int) -> None: ...
