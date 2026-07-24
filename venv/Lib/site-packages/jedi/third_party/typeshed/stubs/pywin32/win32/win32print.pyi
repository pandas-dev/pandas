from _typeshed import Incomplete
from typing import Literal

import _win32typing

def OpenPrinter(printer: str, Defaults: Incomplete | None = ..., /) -> _win32typing.PyPrinterHANDLE: ...
def GetPrinter(hPrinter: _win32typing.PyPrinterHANDLE, Level: int = ..., /): ...
def SetPrinter(hPrinter: _win32typing.PyPrinterHANDLE, Level, pPrinter, Command, /) -> None: ...
def ClosePrinter(hPrinter: _win32typing.PyPrinterHANDLE, /) -> None: ...
def AddPrinterConnection(printer: str, /): ...
def DeletePrinterConnection(printer: str, /): ...
def EnumPrinters(flags, name: str | None = ..., level: int = ..., /): ...
def GetDefaultPrinter() -> str: ...
def GetDefaultPrinterW() -> str: ...
def SetDefaultPrinter(printer: str, /): ...
def SetDefaultPrinterW(Printer: str, /): ...
def StartDocPrinter(
    hprinter: _win32typing.PyPrinterHANDLE | int, level: Literal[1], tuple: tuple[str, str, str | None], /
) -> int: ...
def EndDocPrinter(hPrinter: _win32typing.PyPrinterHANDLE, /): ...
def AbortPrinter(hPrinter: _win32typing.PyPrinterHANDLE, /) -> None: ...
def StartPagePrinter(hprinter: _win32typing.PyPrinterHANDLE, /) -> None: ...
def EndPagePrinter(hprinter: _win32typing.PyPrinterHANDLE, /) -> None: ...
def StartDoc(hdc: int, docinfo, /): ...
def EndDoc(hdc: int, /) -> None: ...
def AbortDoc(hdc: int, /) -> None: ...
def StartPage(hdc: int, /) -> None: ...
def EndPage(hdc: int, /) -> None: ...
def WritePrinter(hprinter: int | _win32typing.PyPrinterHANDLE, buf: bytes | bytearray | memoryview, /) -> int: ...
def EnumJobs(hPrinter: _win32typing.PyPrinterHANDLE, FirstJob, NoJobs, Level=..., /): ...
def GetJob(hPrinter: _win32typing.PyPrinterHANDLE, JobID, Level: int = ..., /): ...
def SetJob(hPrinter: _win32typing.PyPrinterHANDLE, JobID, Level, JobInfo, Command, /): ...
def DocumentProperties(
    HWnd: int,
    hPrinter: _win32typing.PyPrinterHANDLE,
    DeviceName: str,
    DevModeOutput: _win32typing.PyDEVMODE,
    DevModeInput: _win32typing.PyDEVMODE,
    Mode,
    /,
): ...
def EnumPrintProcessors(Server: str | None = ..., Environment: str | None = ..., /) -> tuple[str, ...]: ...
def EnumPrintProcessorDatatypes(ServerName: str, PrintProcessorName: str, /) -> tuple[str, ...]: ...
def EnumPrinterDrivers(Server: str | None = ..., Environment: str | None = ..., Level=..., /) -> tuple[Incomplete, ...]: ...
def EnumForms(hprinter: _win32typing.PyPrinterHANDLE, /) -> tuple[_win32typing.FORM_INFO_1, ...]: ...
def AddForm(hprinter: _win32typing.PyPrinterHANDLE, Form, /) -> None: ...
def DeleteForm(hprinter: _win32typing.PyPrinterHANDLE, FormName: str, /) -> None: ...
def GetForm(hprinter: _win32typing.PyPrinterHANDLE, FormName: str, /) -> None: ...
def SetForm(hprinter: _win32typing.PyPrinterHANDLE, FormName: str, Form, /) -> None: ...
def AddJob(hprinter: _win32typing.PyPrinterHANDLE, /) -> None: ...
def ScheduleJob(hprinter: _win32typing.PyPrinterHANDLE, JobId, /) -> None: ...
def DeviceCapabilities(Device: str, Port: str, Capability, DEVMODE: _win32typing.PyDEVMODE | None = ..., /) -> None: ...
def GetDeviceCaps(hdc: int, Index, /): ...
def EnumMonitors(Name: str, Level, /) -> tuple[Incomplete, ...]: ...
def EnumPorts(Name: str, Level, /) -> tuple[Incomplete, ...]: ...
def GetPrintProcessorDirectory(Name: str, Environment: str, /) -> str: ...
def GetPrinterDriverDirectory(Name: str, Environment: str, /) -> str: ...
def AddPrinter(Name: str, Level, pPrinter, /) -> _win32typing.PyPrinterHANDLE: ...
def DeletePrinter(hPrinter: _win32typing.PyPrinterHANDLE, /) -> None: ...
def DeletePrinterDriver(Server: str, Environment: str, DriverName: str, /) -> None: ...
def DeletePrinterDriverEx(Server: str, Environment: str, DriverName: str, DeleteFlag, VersionFlag, /) -> None: ...
def FlushPrinter(Printer: _win32typing.PyPrinterHANDLE, Buf, Sleep, /): ...

DEF_PRIORITY: int
DI_APPBANDING: int
DI_ROPS_READ_DESTINATION: int
DPD_DELETE_ALL_FILES: int
DPD_DELETE_SPECIFIC_VERSION: int
DPD_DELETE_UNUSED_FILES: int
DSPRINT_PENDING: int
DSPRINT_PUBLISH: int
DSPRINT_REPUBLISH: int
DSPRINT_UNPUBLISH: int
DSPRINT_UPDATE: int
FORM_BUILTIN: int
FORM_PRINTER: int
FORM_USER: int
JOB_ACCESS_ADMINISTER: int
JOB_ACCESS_READ: int
JOB_ALL_ACCESS: int
JOB_CONTROL_CANCEL: int
JOB_CONTROL_DELETE: int
JOB_CONTROL_LAST_PAGE_EJECTED: int
JOB_CONTROL_PAUSE: int
JOB_CONTROL_RESTART: int
JOB_CONTROL_RESUME: int
JOB_CONTROL_SENT_TO_PRINTER: int
JOB_EXECUTE: int
JOB_INFO_1: int
JOB_POSITION_UNSPECIFIED: int
JOB_READ: int
JOB_STATUS_BLOCKED_DEVQ: int
JOB_STATUS_COMPLETE: int
JOB_STATUS_DELETED: int
JOB_STATUS_DELETING: int
JOB_STATUS_ERROR: int
JOB_STATUS_OFFLINE: int
JOB_STATUS_PAPEROUT: int
JOB_STATUS_PAUSED: int
JOB_STATUS_PRINTED: int
JOB_STATUS_PRINTING: int
JOB_STATUS_RESTART: int
JOB_STATUS_SPOOLING: int
JOB_STATUS_USER_INTERVENTION: int
JOB_WRITE: int
MAX_PRIORITY: int
MIN_PRIORITY: int
PORT_STATUS_DOOR_OPEN: int
PORT_STATUS_NO_TONER: int
PORT_STATUS_OFFLINE: int
PORT_STATUS_OUTPUT_BIN_FULL: int
PORT_STATUS_OUT_OF_MEMORY: int
PORT_STATUS_PAPER_JAM: int
PORT_STATUS_PAPER_OUT: int
PORT_STATUS_PAPER_PROBLEM: int
PORT_STATUS_POWER_SAVE: int
PORT_STATUS_TONER_LOW: int
PORT_STATUS_TYPE_ERROR: int
PORT_STATUS_TYPE_INFO: int
PORT_STATUS_TYPE_WARNING: int
PORT_STATUS_USER_INTERVENTION: int
PORT_STATUS_WARMING_UP: int
PORT_TYPE_NET_ATTACHED: int
PORT_TYPE_READ: int
PORT_TYPE_REDIRECTED: int
PORT_TYPE_WRITE: int
PRINTER_ACCESS_ADMINISTER: int
PRINTER_ACCESS_USE: int
PRINTER_ALL_ACCESS: int
PRINTER_ATTRIBUTE_DEFAULT: int
PRINTER_ATTRIBUTE_DIRECT: int
PRINTER_ATTRIBUTE_DO_COMPLETE_FIRST: int
PRINTER_ATTRIBUTE_ENABLE_BIDI: int
PRINTER_ATTRIBUTE_ENABLE_DEVQ: int
PRINTER_ATTRIBUTE_FAX: int
PRINTER_ATTRIBUTE_HIDDEN: int
PRINTER_ATTRIBUTE_KEEPPRINTEDJOBS: int
PRINTER_ATTRIBUTE_LOCAL: int
PRINTER_ATTRIBUTE_NETWORK: int
PRINTER_ATTRIBUTE_PUBLISHED: int
PRINTER_ATTRIBUTE_QUEUED: int
PRINTER_ATTRIBUTE_RAW_ONLY: int
PRINTER_ATTRIBUTE_SHARED: int
PRINTER_ATTRIBUTE_TS: int
PRINTER_ATTRIBUTE_WORK_OFFLINE: int
PRINTER_CONTROL_PAUSE: int
PRINTER_CONTROL_PURGE: int
PRINTER_CONTROL_RESUME: int
PRINTER_CONTROL_SET_STATUS: int
PRINTER_ENUM_CONNECTIONS: int
PRINTER_ENUM_CONTAINER: int
PRINTER_ENUM_DEFAULT: int
PRINTER_ENUM_EXPAND: int
PRINTER_ENUM_ICON1: int
PRINTER_ENUM_ICON2: int
PRINTER_ENUM_ICON3: int
PRINTER_ENUM_ICON4: int
PRINTER_ENUM_ICON5: int
PRINTER_ENUM_ICON6: int
PRINTER_ENUM_ICON7: int
PRINTER_ENUM_ICON8: int
PRINTER_ENUM_LOCAL: int
PRINTER_ENUM_NAME: int
PRINTER_ENUM_NETWORK: int
PRINTER_ENUM_REMOTE: int
PRINTER_ENUM_SHARED: int
PRINTER_EXECUTE: int
PRINTER_INFO_1: int
PRINTER_READ: int
PRINTER_STATUS_BUSY: int
PRINTER_STATUS_DOOR_OPEN: int
PRINTER_STATUS_ERROR: int
PRINTER_STATUS_INITIALIZING: int
PRINTER_STATUS_IO_ACTIVE: int
PRINTER_STATUS_MANUAL_FEED: int
PRINTER_STATUS_NOT_AVAILABLE: int
PRINTER_STATUS_NO_TONER: int
PRINTER_STATUS_OFFLINE: int
PRINTER_STATUS_OUTPUT_BIN_FULL: int
PRINTER_STATUS_OUT_OF_MEMORY: int
PRINTER_STATUS_PAGE_PUNT: int
PRINTER_STATUS_PAPER_JAM: int
PRINTER_STATUS_PAPER_OUT: int
PRINTER_STATUS_PAPER_PROBLEM: int
PRINTER_STATUS_PAUSED: int
PRINTER_STATUS_PENDING_DELETION: int
PRINTER_STATUS_POWER_SAVE: int
PRINTER_STATUS_PRINTING: int
PRINTER_STATUS_PROCESSING: int
PRINTER_STATUS_SERVER_UNKNOWN: int
PRINTER_STATUS_TONER_LOW: int
PRINTER_STATUS_USER_INTERVENTION: int
PRINTER_STATUS_WAITING: int
PRINTER_STATUS_WARMING_UP: int
PRINTER_WRITE: int
SERVER_ACCESS_ADMINISTER: int
SERVER_ACCESS_ENUMERATE: int
SERVER_ALL_ACCESS: int
SERVER_EXECUTE: int
SERVER_READ: int
SERVER_WRITE: int
