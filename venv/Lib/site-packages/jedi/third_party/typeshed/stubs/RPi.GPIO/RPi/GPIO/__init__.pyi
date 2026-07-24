from collections.abc import Callable
from typing import Final, Literal, TypedDict, type_check_only
from typing_extensions import TypeAlias

@type_check_only
class _RPi_Info(TypedDict):
    P1_REVISION: int
    REVISION: str
    TYPE: str
    MANUFACTURER: str
    PROCESSOR: str
    RAM: str

VERSION: str
RPI_INFO: _RPi_Info
RPI_REVISION: int

HIGH: Literal[1]
LOW: Literal[0]

OUT: Final = 0
IN: Final = 1
HARD_PWM: Final = 43
SERIAL: Final = 40
I2C: Final = 42
SPI: Final = 41
UNKNOWN: Final = -1

BOARD: Final = 10
BCM: Final = 11

PUD_OFF: Final = 20
PUD_UP: Final = 22
PUD_DOWN: Final = 21

RISING: Final = 31
FALLING: Final = 32
BOTH: Final = 33

_EventCallback: TypeAlias = Callable[[int], object]

def setup(
    channel: int | list[int] | tuple[int, ...], direction: Literal[0, 1], pull_up_down: int = 20, initial: int = -1
) -> None: ...
def cleanup(channel: int | list[int] | tuple[int, ...] = -666) -> None: ...
def output(
    channel: int | list[int] | tuple[int, ...],
    value: Literal[0, 1] | bool | list[Literal[0, 1] | bool] | tuple[Literal[0, 1] | bool, ...],
    /,
) -> None: ...
def input(channel: int, /) -> bool: ...
def setmode(mode: Literal[10, 11], /) -> None: ...
def getmode() -> Literal[10, 11] | None: ...
def add_event_detect(channel: int, edge: int, callback: _EventCallback | None = None, bouncetime: int = -666) -> None: ...
def remove_event_detect(channel: int, /) -> None: ...
def event_detected(channel: int, /) -> bool: ...
def add_event_callback(channel: int, callback: _EventCallback) -> None: ...
def wait_for_edge(channel: int, edge: int, bouncetime: int = -666, timeout: int = -1) -> int | None: ...
def gpio_function(channel: int, /) -> int: ...
def setwarnings(gpio_warnings: bool, /) -> None: ...

class PWM:
    def __init__(self, channel: int, frequency: float, /) -> None: ...
    def start(self, dutycycle: float, /) -> None: ...
    def ChangeDutyCycle(self, dutycycle: float, /) -> None: ...
    def ChangeFrequency(self, frequency: float, /) -> None: ...
    def stop(self) -> None: ...
