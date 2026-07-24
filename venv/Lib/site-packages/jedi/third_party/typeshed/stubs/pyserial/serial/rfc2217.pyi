import logging
from _typeshed import ReadableBuffer
from collections.abc import Callable, Generator
from typing import Any

from serial.serialutil import SerialBase

LOGGER_LEVELS: dict[str, int]
SE: bytes
NOP: bytes
DM: bytes
BRK: bytes
IP: bytes
AO: bytes
AYT: bytes
EC: bytes
EL: bytes
GA: bytes
SB: bytes
WILL: bytes
WONT: bytes
DO: bytes
DONT: bytes
IAC: bytes
IAC_DOUBLED: bytes
BINARY: bytes
ECHO: bytes
SGA: bytes
COM_PORT_OPTION: bytes
SET_BAUDRATE: bytes
SET_DATASIZE: bytes
SET_PARITY: bytes
SET_STOPSIZE: bytes
SET_CONTROL: bytes
NOTIFY_LINESTATE: bytes
NOTIFY_MODEMSTATE: bytes
FLOWCONTROL_SUSPEND: bytes
FLOWCONTROL_RESUME: bytes
SET_LINESTATE_MASK: bytes
SET_MODEMSTATE_MASK: bytes
PURGE_DATA: bytes
SERVER_SET_BAUDRATE: bytes
SERVER_SET_DATASIZE: bytes
SERVER_SET_PARITY: bytes
SERVER_SET_STOPSIZE: bytes
SERVER_SET_CONTROL: bytes
SERVER_NOTIFY_LINESTATE: bytes
SERVER_NOTIFY_MODEMSTATE: bytes
SERVER_FLOWCONTROL_SUSPEND: bytes
SERVER_FLOWCONTROL_RESUME: bytes
SERVER_SET_LINESTATE_MASK: bytes
SERVER_SET_MODEMSTATE_MASK: bytes
SERVER_PURGE_DATA: bytes
RFC2217_ANSWER_MAP: dict[bytes, bytes]
SET_CONTROL_REQ_FLOW_SETTING: bytes
SET_CONTROL_USE_NO_FLOW_CONTROL: bytes
SET_CONTROL_USE_SW_FLOW_CONTROL: bytes
SET_CONTROL_USE_HW_FLOW_CONTROL: bytes
SET_CONTROL_REQ_BREAK_STATE: bytes
SET_CONTROL_BREAK_ON: bytes
SET_CONTROL_BREAK_OFF: bytes
SET_CONTROL_REQ_DTR: bytes
SET_CONTROL_DTR_ON: bytes
SET_CONTROL_DTR_OFF: bytes
SET_CONTROL_REQ_RTS: bytes
SET_CONTROL_RTS_ON: bytes
SET_CONTROL_RTS_OFF: bytes
SET_CONTROL_REQ_FLOW_SETTING_IN: bytes
SET_CONTROL_USE_NO_FLOW_CONTROL_IN: bytes
SET_CONTROL_USE_SW_FLOW_CONTOL_IN: bytes
SET_CONTROL_USE_HW_FLOW_CONTOL_IN: bytes
SET_CONTROL_USE_DCD_FLOW_CONTROL: bytes
SET_CONTROL_USE_DTR_FLOW_CONTROL: bytes
SET_CONTROL_USE_DSR_FLOW_CONTROL: bytes
LINESTATE_MASK_TIMEOUT: int
LINESTATE_MASK_SHIFTREG_EMPTY: int
LINESTATE_MASK_TRANSREG_EMPTY: int
LINESTATE_MASK_BREAK_DETECT: int
LINESTATE_MASK_FRAMING_ERROR: int
LINESTATE_MASK_PARTIY_ERROR: int
LINESTATE_MASK_OVERRUN_ERROR: int
LINESTATE_MASK_DATA_READY: int
MODEMSTATE_MASK_CD: int
MODEMSTATE_MASK_RI: int
MODEMSTATE_MASK_DSR: int
MODEMSTATE_MASK_CTS: int
MODEMSTATE_MASK_CD_CHANGE: int
MODEMSTATE_MASK_RI_CHANGE: int
MODEMSTATE_MASK_DSR_CHANGE: int
MODEMSTATE_MASK_CTS_CHANGE: int
PURGE_RECEIVE_BUFFER: bytes
PURGE_TRANSMIT_BUFFER: bytes
PURGE_BOTH_BUFFERS: bytes
RFC2217_PARITY_MAP: dict[str, int]
RFC2217_REVERSE_PARITY_MAP: dict[int, str]
RFC2217_STOPBIT_MAP: dict[int | float, int]
RFC2217_REVERSE_STOPBIT_MAP: dict[int, int | float]
M_NORMAL: int
M_IAC_SEEN: int
M_NEGOTIATE: int
REQUESTED: str
ACTIVE: str
INACTIVE: str
REALLY_INACTIVE: str

class TelnetOption:
    connection: Serial
    name: str
    option: bytes
    send_yes: bytes
    send_no: bytes
    ack_yes: bytes
    ack_no: bytes
    state: str
    active: bool
    activation_callback: Callable[[], Any]

    def __init__(
        self,
        connection: Serial,
        name: str,
        option: bytes,
        send_yes: bytes,
        send_no: bytes,
        ack_yes: bytes,
        ack_no: bytes,
        initial_state: str,
        activation_callback: Callable[[], Any] | None = None,
    ) -> None: ...
    def process_incoming(self, command: bytes) -> None: ...

class TelnetSubnegotiation:
    connection: Serial
    name: str
    option: bytes
    value: bytes | None
    ack_option: bytes
    state: str
    def __init__(self, connection: Serial, name: str, option: bytes, ack_option: bytes | None = None) -> None: ...
    def set(self, value: bytes) -> None: ...
    def is_ready(self) -> bool: ...
    @property
    def active(self) -> bool: ...
    def wait(self, timeout: float = 3) -> None: ...
    def check_answer(self, suboption: bytes) -> None: ...

class Serial(SerialBase):
    logger: logging.Logger | None
    def open(self) -> None: ...
    def from_url(self, url: str) -> tuple[str, int]: ...
    @property
    def in_waiting(self) -> int: ...
    def read(self, size: int = 1) -> bytes: ...
    def write(self, b: ReadableBuffer, /) -> int | None: ...
    def reset_input_buffer(self) -> None: ...
    def reset_output_buffer(self) -> None: ...
    @property
    def cts(self) -> bool: ...
    @property
    def dsr(self) -> bool: ...
    @property
    def ri(self) -> bool: ...
    @property
    def cd(self) -> bool: ...
    def telnet_send_option(self, action: bytes, option: bytes) -> None: ...
    def rfc2217_send_subnegotiation(self, option: bytes, value: bytes = b"") -> None: ...
    def rfc2217_send_purge(self, value: bytes) -> None: ...
    def rfc2217_set_control(self, value: bytes) -> None: ...
    def rfc2217_flow_server_ready(self) -> None: ...
    def get_modem_state(self) -> int: ...

class PortManager:
    serial: Serial
    connection: Serial
    logger: logging.Logger | None
    mode: int
    suboption: bytes | None
    telnet_command: bytes | None
    modemstate_mask: int
    last_modemstate: int | None
    linstate_mask: int
    def __init__(self, serial_port: Serial, connection: Serial, logger: logging.Logger | None = None) -> None: ...
    def telnet_send_option(self, action: bytes, option: bytes) -> None: ...
    def rfc2217_send_subnegotiation(self, option: bytes, value: bytes = b"") -> None: ...
    def check_modem_lines(self, force_notification: bool = False) -> None: ...
    def escape(self, data: bytes) -> Generator[bytes]: ...
    def filter(self, data: bytes) -> Generator[bytes]: ...
