from _typeshed import Incomplete

import _win32typing
from win32.lib.pywintypes import error as error

def GetNamedPipeHandleState(
    hPipe: int, bGetCollectionData=..., /
) -> tuple[Incomplete, Incomplete, Incomplete, Incomplete, str]: ...
def SetNamedPipeHandleState(
    hPipe: int, Mode: int, MaxCollectionCount: None | Incomplete, CollectDataTimeout: None | Incomplete, /
) -> None: ...
def ConnectNamedPipe(hPipe: int, overlapped: _win32typing.PyOVERLAPPED | None = ..., /): ...
def TransactNamedPipe(
    pipeName,
    writeData: str,
    buffer_bufSize: _win32typing.PyOVERLAPPEDReadBuffer,
    overlapped: _win32typing.PyOVERLAPPED | None = ...,
    /,
) -> str: ...
def CallNamedPipe(pipeName, data: str, bufSize, timeOut, /) -> str: ...
def CreatePipe(sa: _win32typing.PySECURITY_ATTRIBUTES, nSize: int, /) -> tuple[int, int]: ...
def FdCreatePipe(sa: _win32typing.PySECURITY_ATTRIBUTES, nSize, mode, /) -> tuple[Incomplete, Incomplete]: ...
def CreateNamedPipe(
    pipeName: str,
    openMode,
    pipeMode,
    nMaxInstances,
    nOutBufferSize,
    nInBufferSize,
    nDefaultTimeOut,
    sa: _win32typing.PySECURITY_ATTRIBUTES,
    /,
) -> int: ...
def DisconnectNamedPipe(hFile: int, /) -> None: ...
def GetOverlappedResult(hFile: int, overlapped: _win32typing.PyOVERLAPPED, bWait: int | bool, /) -> int: ...
def WaitNamedPipe(pipeName: str, timeout, /) -> None: ...
def GetNamedPipeInfo(hNamedPipe: int, /) -> tuple[int, int, int, int]: ...
def PeekNamedPipe(hPipe: int, size: int, /) -> tuple[str, int, Incomplete]: ...
def GetNamedPipeClientProcessId(hPipe: int, /): ...
def GetNamedPipeServerProcessId(hPipe: int, /): ...
def GetNamedPipeClientSessionId(hPipe: int, /): ...
def GetNamedPipeServerSessionId(hPipe: int, /): ...
def popen(cmdstring: str, mode: str, /): ...
def popen2(*args): ...  # incomplete
def popen3(*args): ...  # incomplete
def popen4(*args): ...  # incomplete

FILE_FLAG_FIRST_PIPE_INSTANCE: int
PIPE_ACCEPT_REMOTE_CLIENTS: int
PIPE_REJECT_REMOTE_CLIENTS: int
NMPWAIT_NOWAIT: int
NMPWAIT_USE_DEFAULT_WAIT: int
NMPWAIT_WAIT_FOREVER: int
PIPE_ACCESS_DUPLEX: int
PIPE_ACCESS_INBOUND: int
PIPE_ACCESS_OUTBOUND: int
PIPE_NOWAIT: int
PIPE_READMODE_BYTE: int
PIPE_READMODE_MESSAGE: int
PIPE_TYPE_BYTE: int
PIPE_TYPE_MESSAGE: int
PIPE_UNLIMITED_INSTANCES: int
PIPE_WAIT: int
UNICODE: int
