from _typeshed import Incomplete

import _win32typing
from win32.lib.pywintypes import error as error

def CreatePhonebookEntry(hWnd: int, fileName: str | None = ..., /) -> None: ...
def Dial(
    dialExtensions, fileName: str, RasDialParams: _win32typing.RASDIALPARAMS, callback, /
) -> tuple[Incomplete, Incomplete]: ...
def EditPhonebookEntry(hWnd: int, fileName: str, entryName: str | None = ..., /) -> None: ...
def EnumConnections(): ...
def EnumEntries(reserved: str | None = ..., fileName: str | None = ..., /) -> None: ...
def GetConnectStatus(hrasconn, /) -> tuple[Incomplete, Incomplete, str, str]: ...
def GetEntryDialParams(
    fileName: str, entryName: str, /
) -> tuple[tuple[Incomplete, Incomplete, Incomplete, Incomplete, Incomplete, Incomplete], bool]: ...
def GetErrorString(error, /) -> str: ...
def HangUp(hras, /) -> None: ...
def IsHandleValid(hras: int | None, /) -> bool: ...
def SetEntryDialParams(fileName: str, RasDialParams, bSavePassword, /) -> None: ...
def RASDIALEXTENSIONS() -> _win32typing.RASDIALEXTENSIONS: ...

RASCS_AllDevicesConnected: int
RASCS_AuthAck: int
RASCS_AuthCallback: int
RASCS_AuthChangePassword: int
RASCS_Authenticate: int
RASCS_Authenticated: int
RASCS_AuthLinkSpeed: int
RASCS_AuthNotify: int
RASCS_AuthProject: int
RASCS_AuthRetry: int
RASCS_CallbackComplete: int
RASCS_CallbackSetByCaller: int
RASCS_ConnectDevice: int
RASCS_Connected: int
RASCS_DeviceConnected: int
RASCS_Disconnected: int
RASCS_Interactive: int
RASCS_LogonNetwork: int
RASCS_OpenPort: int
RASCS_PasswordExpired: int
RASCS_PortOpened: int
RASCS_PrepareForCallback: int
RASCS_Projected: int
RASCS_ReAuthenticate: int
RASCS_RetryAuthentication: int
RASCS_StartAuthentication: int
RASCS_WaitForCallback: int
RASCS_WaitForModemReset: int

def GetEapUserIdentity(phoneBook: str | None, entry: str, flags: int, hwnd: _win32typing.PyHANDLE | int | None = None, /): ...

RASEAPF_Logon: int
RASEAPF_NonInteractive: int
RASEAPF_Preview: int
