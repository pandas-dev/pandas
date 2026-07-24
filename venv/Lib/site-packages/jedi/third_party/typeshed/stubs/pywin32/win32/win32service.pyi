from _typeshed import Incomplete
from collections.abc import Iterable

import _win32typing
from win32.lib.pywintypes import error as error

def GetThreadDesktop(ThreadId: int, /) -> _win32typing.PyHDESK: ...
def EnumWindowStations() -> tuple[tuple[str, Incomplete], ...]: ...
def GetUserObjectInformation(Handle: int, _type, /) -> None: ...
def SetUserObjectInformation(Handle: int, info, _type, /) -> None: ...
def OpenWindowStation(szWinSta, Inherit, DesiredAccess, /) -> _win32typing.PyHWINSTA: ...
def OpenDesktop(szDesktop, Flags, Inherit, DesiredAccess, /) -> _win32typing.PyHDESK: ...
def CreateDesktop(
    Desktop, Flags, DesiredAccess, SecurityAttributes: _win32typing.PySECURITY_ATTRIBUTES, /
) -> _win32typing.PyHDESK: ...
def OpenInputDesktop(Flags, Inherit, DesiredAccess, /) -> _win32typing.PyHDESK: ...
def GetProcessWindowStation() -> _win32typing.PyHWINSTA: ...
def CreateWindowStation(
    WindowStation, Flags, DesiredAccess, SecurityAttributes: _win32typing.PySECURITY_ATTRIBUTES, /
) -> _win32typing.PyHWINSTA: ...
def EnumServicesStatus(
    hSCManager: _win32typing.PySC_HANDLE | int, ServiceType: int = ..., ServiceState: int = ..., /
) -> tuple[Incomplete, ...]: ...
def EnumServicesStatusEx(
    SCManager: _win32typing.PySC_HANDLE, ServiceType, ServiceState, InfoLevel, GroupName: Incomplete | None = ..., /
) -> tuple[Incomplete, ...]: ...
def EnumDependentServices(hService: _win32typing.PySC_HANDLE, ServiceState, /) -> tuple[Incomplete, ...]: ...
def QueryServiceConfig(hService: _win32typing.PySC_HANDLE, /): ...
def StartService(hService: _win32typing.PySC_HANDLE, args: Iterable[str] | None, /) -> None: ...
def OpenService(scHandle: _win32typing.PySC_HANDLE, name: str, desiredAccess, /) -> _win32typing.PySC_HANDLE: ...
def OpenSCManager(machineName: str | None, dbName: str | None, desiredAccess: int, /) -> _win32typing.PySC_HANDLE: ...
def CloseServiceHandle(scHandle: _win32typing.PySC_HANDLE, /) -> None: ...
def QueryServiceStatus(hService: _win32typing.PySC_HANDLE, /) -> _win32typing.SERVICE_STATUS: ...
def QueryServiceStatusEx(hService: _win32typing.PySC_HANDLE, /) -> _win32typing.SERVICE_STATUS: ...
def SetServiceObjectSecurity(
    Handle: _win32typing.PySC_HANDLE, SecurityInformation, SecurityDescriptor: _win32typing.PySECURITY_DESCRIPTOR, /
) -> None: ...
def QueryServiceObjectSecurity(
    Handle: _win32typing.PySC_HANDLE, SecurityInformation, /
) -> _win32typing.PySECURITY_DESCRIPTOR: ...
def GetServiceKeyName(hSCManager: _win32typing.PySC_HANDLE, DisplayName, /): ...
def GetServiceDisplayName(hSCManager: _win32typing.PySC_HANDLE, ServiceName, /): ...
def SetServiceStatus(
    scHandle, serviceStatus: _win32typing.SERVICE_STATUS | tuple[int, int, int, int, int, int, int], /
) -> None: ...
def ControlService(scHandle: _win32typing.PySC_HANDLE, code, /) -> _win32typing.SERVICE_STATUS: ...
def DeleteService(scHandle: _win32typing.PySC_HANDLE, /) -> None: ...
def CreateService(
    scHandle: _win32typing.PySC_HANDLE,
    name: str,
    displayName: str,
    desiredAccess: int,
    serviceType: int,
    startType: int,
    errorControl: int,
    binaryFile: str,
    loadOrderGroup: str | None,
    bFetchTag: bool,
    serviceDeps: Iterable[Incomplete] | None,
    acctName: str | None,
    password: str | None,
    /,
) -> _win32typing.PySC_HANDLE: ...
def ChangeServiceConfig(
    hService: _win32typing.PySC_HANDLE,
    serviceType: int,
    startType: int,
    errorControl: int,
    binaryFile: str | None,
    loadOrderGroup: str | None,
    bFetchTag: bool,
    serviceDeps: Iterable[Incomplete] | None,
    acctName: str | None,
    password: str | None,
    displayName: str | None,
    /,
): ...
def LockServiceDatabase(sc_handle: _win32typing.PySC_HANDLE, /): ...
def UnlockServiceDatabase(lock, /): ...
def QueryServiceLockStatus(hSCManager: _win32typing.PySC_HANDLE, /) -> tuple[Incomplete, str, Incomplete]: ...
def ChangeServiceConfig2(hService: _win32typing.PySC_HANDLE, InfoLevel, info, /) -> None: ...
def QueryServiceConfig2(hService: _win32typing.PySC_HANDLE, InfoLevel, /): ...

DBT_CONFIGCHANGECANCELED: int
DBT_CONFIGCHANGED: int
DBT_CUSTOMEVENT: int
DBT_DEVICEARRIVAL: int
DBT_DEVICEQUERYREMOVE: int
DBT_DEVICEQUERYREMOVEFAILED: int
DBT_DEVICEREMOVECOMPLETE: int
DBT_DEVICEREMOVEPENDING: int
DBT_DEVICETYPESPECIFIC: int
DBT_QUERYCHANGECONFIG: int
DF_ALLOWOTHERACCOUNTHOOK: int
SC_ACTION_NONE: int
SC_ACTION_REBOOT: int
SC_ACTION_RESTART: int
SC_ACTION_RUN_COMMAND: int
SC_ENUM_PROCESS_INFO: int
SC_GROUP_IDENTIFIER: int
SC_MANAGER_ALL_ACCESS: int
SC_MANAGER_CONNECT: int
SC_MANAGER_CREATE_SERVICE: int
SC_MANAGER_ENUMERATE_SERVICE: int
SC_MANAGER_LOCK: int
SC_MANAGER_MODIFY_BOOT_CONFIG: int
SC_MANAGER_QUERY_LOCK_STATUS: int
SERVICE_ACCEPT_HARDWAREPROFILECHANGE: int
SERVICE_ACCEPT_NETBINDCHANGE: int
SERVICE_ACCEPT_PARAMCHANGE: int
SERVICE_ACCEPT_PAUSE_CONTINUE: int
SERVICE_ACCEPT_POWEREVENT: int
SERVICE_ACCEPT_PRESHUTDOWN: int
SERVICE_ACCEPT_SESSIONCHANGE: int
SERVICE_ACCEPT_SHUTDOWN: int
SERVICE_ACCEPT_STOP: int
SERVICE_ACTIVE: int
SERVICE_ALL_ACCESS: int
SERVICE_AUTO_START: int
SERVICE_BOOT_START: int
SERVICE_CHANGE_CONFIG: int
SERVICE_CONFIG_DELAYED_AUTO_START_INFO: int
SERVICE_CONFIG_DESCRIPTION: int
SERVICE_CONFIG_FAILURE_ACTIONS: int
SERVICE_CONFIG_FAILURE_ACTIONS_FLAG: int
SERVICE_CONFIG_PRESHUTDOWN_INFO: int
SERVICE_CONFIG_REQUIRED_PRIVILEGES_INFO: int
SERVICE_CONFIG_SERVICE_SID_INFO: int
SERVICE_CONTINUE_PENDING: int
SERVICE_CONTROL_CONTINUE: int
SERVICE_CONTROL_DEVICEEVENT: int
SERVICE_CONTROL_HARDWAREPROFILECHANGE: int
SERVICE_CONTROL_INTERROGATE: int
SERVICE_CONTROL_NETBINDADD: int
SERVICE_CONTROL_NETBINDDISABLE: int
SERVICE_CONTROL_NETBINDENABLE: int
SERVICE_CONTROL_NETBINDREMOVE: int
SERVICE_CONTROL_PARAMCHANGE: int
SERVICE_CONTROL_PAUSE: int
SERVICE_CONTROL_POWEREVENT: int
SERVICE_CONTROL_PRESHUTDOWN: int
SERVICE_CONTROL_SESSIONCHANGE: int
SERVICE_CONTROL_SHUTDOWN: int
SERVICE_CONTROL_STOP: int
SERVICE_DEMAND_START: int
SERVICE_DISABLED: int
SERVICE_DRIVER: int
SERVICE_ENUMERATE_DEPENDENTS: int
SERVICE_ERROR_CRITICAL: int
SERVICE_ERROR_IGNORE: int
SERVICE_ERROR_NORMAL: int
SERVICE_ERROR_SEVERE: int
SERVICE_FILE_SYSTEM_DRIVER: int
SERVICE_INACTIVE: int
SERVICE_INTERACTIVE_PROCESS: int
SERVICE_INTERROGATE: int
SERVICE_KERNEL_DRIVER: int
SERVICE_NO_CHANGE: int
SERVICE_PAUSE_CONTINUE: int
SERVICE_PAUSE_PENDING: int
SERVICE_PAUSED: int
SERVICE_QUERY_CONFIG: int
SERVICE_QUERY_STATUS: int
SERVICE_RUNNING: int
SERVICE_SID_TYPE_NONE: int
SERVICE_SID_TYPE_RESTRICTED: int
SERVICE_SID_TYPE_UNRESTRICTED: int
SERVICE_SPECIFIC_ERROR: int
SERVICE_START: int
SERVICE_START_PENDING: int
SERVICE_STATE_ALL: int
SERVICE_STOP: int
SERVICE_STOP_PENDING: int
SERVICE_STOPPED: int
SERVICE_SYSTEM_START: int
SERVICE_USER_DEFINED_CONTROL: int
SERVICE_WIN32: int
SERVICE_WIN32_OWN_PROCESS: int
SERVICE_WIN32_SHARE_PROCESS: int
UOI_FLAGS: int
UOI_NAME: int
UOI_TYPE: int
UOI_USER_SID: int
WSF_VISIBLE: int
HDESKType = _win32typing.PyHDESK
HWINSTAType = _win32typing.PyHWINSTA
UNICODE: int
