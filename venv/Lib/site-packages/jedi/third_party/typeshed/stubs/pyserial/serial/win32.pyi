import sys
from _typeshed import Incomplete
from ctypes import Structure, Union, _CField, _NamedFuncPointer, _Pointer, c_int64, c_ulong, c_void_p
from ctypes.wintypes import DWORD
from typing_extensions import TypeAlias

if sys.platform == "win32":
    def is_64bit() -> bool: ...

    ULONG_PTR: type[c_int64 | c_ulong]

    class _SECURITY_ATTRIBUTES(Structure):
        nLength: _CField[Incomplete, Incomplete, Incomplete]
        lpSecurityDescriptor: _CField[Incomplete, Incomplete, Incomplete]
        bInheritHandle: _CField[Incomplete, Incomplete, Incomplete]

    LPSECURITY_ATTRIBUTES: type[_Pointer[_SECURITY_ATTRIBUTES]]
    CreateEvent: _NamedFuncPointer
    CreateFile: _NamedFuncPointer
    # The following are included in __all__ but their existence is not guaranteed as
    # they are defined in a try/except block. Their aliases above are always defined.
    CreateEventW: _NamedFuncPointer
    CreateFileW: _NamedFuncPointer

    class _OVERLAPPED(Structure):
        Internal: _CField[Incomplete, Incomplete, Incomplete]
        InternalHigh: _CField[Incomplete, Incomplete, Incomplete]
        Offset: _CField[Incomplete, Incomplete, Incomplete]
        OffsetHigh: _CField[Incomplete, Incomplete, Incomplete]
        Pointer: _CField[Incomplete, Incomplete, Incomplete]
        hEvent: _CField[Incomplete, Incomplete, Incomplete]

    OVERLAPPED: TypeAlias = _OVERLAPPED

    class _COMSTAT(Structure):
        fCtsHold: _CField[Incomplete, Incomplete, Incomplete]
        fDsrHold: _CField[Incomplete, Incomplete, Incomplete]
        fRlsdHold: _CField[Incomplete, Incomplete, Incomplete]
        fXoffHold: _CField[Incomplete, Incomplete, Incomplete]
        fXoffSent: _CField[Incomplete, Incomplete, Incomplete]
        fEof: _CField[Incomplete, Incomplete, Incomplete]
        fTxim: _CField[Incomplete, Incomplete, Incomplete]
        fReserved: _CField[Incomplete, Incomplete, Incomplete]
        cbInQue: _CField[Incomplete, Incomplete, Incomplete]
        cbOutQue: _CField[Incomplete, Incomplete, Incomplete]

    COMSTAT: TypeAlias = _COMSTAT

    class _DCB(Structure):
        DCBlength: _CField[Incomplete, Incomplete, Incomplete]
        BaudRate: _CField[Incomplete, Incomplete, Incomplete]
        fBinary: _CField[Incomplete, Incomplete, Incomplete]
        fParity: _CField[Incomplete, Incomplete, Incomplete]
        fOutxCtsFlow: _CField[Incomplete, Incomplete, Incomplete]
        fOutxDsrFlow: _CField[Incomplete, Incomplete, Incomplete]
        fDtrControl: _CField[Incomplete, Incomplete, Incomplete]
        fDsrSensitivity: _CField[Incomplete, Incomplete, Incomplete]
        fTXContinueOnXoff: _CField[Incomplete, Incomplete, Incomplete]
        fOutX: _CField[Incomplete, Incomplete, Incomplete]
        fInX: _CField[Incomplete, Incomplete, Incomplete]
        fErrorChar: _CField[Incomplete, Incomplete, Incomplete]
        fNull: _CField[Incomplete, Incomplete, Incomplete]
        fRtsControl: _CField[Incomplete, Incomplete, Incomplete]
        fAbortOnError: _CField[Incomplete, Incomplete, Incomplete]
        fDummy2: _CField[Incomplete, Incomplete, Incomplete]
        wReserved: _CField[Incomplete, Incomplete, Incomplete]
        XonLim: _CField[Incomplete, Incomplete, Incomplete]
        XoffLim: _CField[Incomplete, Incomplete, Incomplete]
        ByteSize: _CField[Incomplete, Incomplete, Incomplete]
        Parity: _CField[Incomplete, Incomplete, Incomplete]
        StopBits: _CField[Incomplete, Incomplete, Incomplete]
        XonChar: _CField[Incomplete, Incomplete, Incomplete]
        XoffChar: _CField[Incomplete, Incomplete, Incomplete]
        ErrorChar: _CField[Incomplete, Incomplete, Incomplete]
        EofChar: _CField[Incomplete, Incomplete, Incomplete]
        EvtChar: _CField[Incomplete, Incomplete, Incomplete]
        wReserved1: _CField[Incomplete, Incomplete, Incomplete]

    DCB: TypeAlias = _DCB

    class _COMMTIMEOUTS(Structure):
        ReadIntervalTimeout: _CField[Incomplete, Incomplete, Incomplete]
        ReadTotalTimeoutMultiplier: _CField[Incomplete, Incomplete, Incomplete]
        ReadTotalTimeoutConstant: _CField[Incomplete, Incomplete, Incomplete]
        WriteTotalTimeoutMultiplier: _CField[Incomplete, Incomplete, Incomplete]
        WriteTotalTimeoutConstant: _CField[Incomplete, Incomplete, Incomplete]

    COMMTIMEOUTS: TypeAlias = _COMMTIMEOUTS

    GetLastError: _NamedFuncPointer
    LPOVERLAPPED: type[_Pointer[_OVERLAPPED]]
    LPDWORD: type[_Pointer[DWORD]]
    GetOverlappedResult: _NamedFuncPointer
    ResetEvent: _NamedFuncPointer
    LPCVOID = c_void_p
    WriteFile: _NamedFuncPointer
    LPVOID = c_void_p
    ReadFile: _NamedFuncPointer
    CloseHandle: _NamedFuncPointer
    ClearCommBreak: _NamedFuncPointer
    LPCOMSTAT: type[_Pointer[_COMSTAT]]
    ClearCommError: _NamedFuncPointer
    SetupComm: _NamedFuncPointer
    EscapeCommFunction: _NamedFuncPointer
    GetCommModemStatus: _NamedFuncPointer
    LPDCB: type[_Pointer[_DCB]]
    GetCommState: _NamedFuncPointer
    LPCOMMTIMEOUTS: type[_Pointer[_COMMTIMEOUTS]]
    GetCommTimeouts: _NamedFuncPointer
    PurgeComm: _NamedFuncPointer
    SetCommBreak: _NamedFuncPointer
    SetCommMask: _NamedFuncPointer
    SetCommState: _NamedFuncPointer
    SetCommTimeouts: _NamedFuncPointer
    WaitForSingleObject: _NamedFuncPointer
    WaitCommEvent: _NamedFuncPointer
    CancelIoEx: _NamedFuncPointer

    ONESTOPBIT: int
    TWOSTOPBITS: int
    NOPARITY: int
    ODDPARITY: int
    EVENPARITY: int
    RTS_CONTROL_HANDSHAKE: int
    RTS_CONTROL_ENABLE: int
    DTR_CONTROL_HANDSHAKE: int
    DTR_CONTROL_ENABLE: int
    MS_DSR_ON: int
    EV_RING: int
    EV_PERR: int
    EV_ERR: int
    SETXOFF: int
    EV_RXCHAR: int
    GENERIC_WRITE: int
    PURGE_TXCLEAR: int
    FILE_FLAG_OVERLAPPED: int
    EV_DSR: int
    MAXDWORD: int
    EV_RLSD: int
    ERROR_IO_PENDING: int
    MS_CTS_ON: int
    EV_EVENT1: int
    EV_RX80FULL: int
    PURGE_RXABORT: int
    FILE_ATTRIBUTE_NORMAL: int
    PURGE_TXABORT: int
    SETXON: int
    OPEN_EXISTING: int
    MS_RING_ON: int
    EV_TXEMPTY: int
    EV_RXFLAG: int
    MS_RLSD_ON: int
    GENERIC_READ: int
    EV_EVENT2: int
    EV_CTS: int
    EV_BREAK: int
    PURGE_RXCLEAR: int

    class N11_OVERLAPPED4DOLLAR_48E(Union):
        Offset: _CField[Incomplete, Incomplete, Incomplete]
        OffsetHigh: _CField[Incomplete, Incomplete, Incomplete]
        Pointer: _CField[Incomplete, Incomplete, Incomplete]

    class N11_OVERLAPPED4DOLLAR_484DOLLAR_49E(Structure):
        Offset: _CField[Incomplete, Incomplete, Incomplete]
        OffsetHigh: _CField[Incomplete, Incomplete, Incomplete]

    PVOID: TypeAlias = c_void_p

    __all__ = [
        "GetLastError",
        "MS_CTS_ON",
        "FILE_ATTRIBUTE_NORMAL",
        "DTR_CONTROL_ENABLE",
        "_COMSTAT",
        "MS_RLSD_ON",
        "GetOverlappedResult",
        "SETXON",
        "PURGE_TXABORT",
        "PurgeComm",
        "N11_OVERLAPPED4DOLLAR_48E",
        "EV_RING",
        "ONESTOPBIT",
        "SETXOFF",
        "PURGE_RXABORT",
        "GetCommState",
        "RTS_CONTROL_ENABLE",
        "_DCB",
        "CreateEvent",
        "_COMMTIMEOUTS",
        "_SECURITY_ATTRIBUTES",
        "EV_DSR",
        "EV_PERR",
        "EV_RXFLAG",
        "OPEN_EXISTING",
        "DCB",
        "FILE_FLAG_OVERLAPPED",
        "EV_CTS",
        "SetupComm",
        "LPOVERLAPPED",
        "EV_TXEMPTY",
        "ClearCommBreak",
        "LPSECURITY_ATTRIBUTES",
        "SetCommBreak",
        "SetCommTimeouts",
        "COMMTIMEOUTS",
        "ODDPARITY",
        "EV_RLSD",
        "GetCommModemStatus",
        "EV_EVENT2",
        "PURGE_TXCLEAR",
        "EV_BREAK",
        "EVENPARITY",
        "LPCVOID",
        "COMSTAT",
        "ReadFile",
        "PVOID",
        "_OVERLAPPED",
        "WriteFile",
        "GetCommTimeouts",
        "ResetEvent",
        "EV_RXCHAR",
        "LPCOMSTAT",
        "ClearCommError",
        "ERROR_IO_PENDING",
        "EscapeCommFunction",
        "GENERIC_READ",
        "RTS_CONTROL_HANDSHAKE",
        "OVERLAPPED",
        "DTR_CONTROL_HANDSHAKE",
        "PURGE_RXCLEAR",
        "GENERIC_WRITE",
        "LPDCB",
        "CreateEventW",
        "SetCommMask",
        "EV_EVENT1",
        "SetCommState",
        "LPVOID",
        "CreateFileW",
        "LPDWORD",
        "EV_RX80FULL",
        "TWOSTOPBITS",
        "LPCOMMTIMEOUTS",
        "MAXDWORD",
        "MS_DSR_ON",
        "MS_RING_ON",
        "N11_OVERLAPPED4DOLLAR_484DOLLAR_49E",
        "EV_ERR",
        "ULONG_PTR",
        "CreateFile",
        "NOPARITY",
        "CloseHandle",
    ]
