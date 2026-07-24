import ctypes
import sys
from _typeshed import Incomplete
from collections.abc import Generator
from ctypes.wintypes import DWORD

from serial.tools.list_ports_common import ListPortInfo

if sys.platform == "win32":

    def ValidHandle(
        value: type[ctypes._CData] | None, func: ctypes._FuncPointer, arguments: tuple[ctypes._CData, ...]
    ) -> ctypes._CData: ...

    NULL: int
    HDEVINFO = ctypes.c_void_p
    LPCTSTR = ctypes.c_wchar_p
    PCTSTR = ctypes.c_wchar_p
    PTSTR = ctypes.c_wchar_p
    LPDWORD: ctypes._Pointer[DWORD]
    PDWORD: ctypes._Pointer[DWORD]
    LPBYTE = ctypes.c_void_p
    PBYTE = ctypes.c_void_p
    ACCESS_MASK = DWORD
    REGSAM = ACCESS_MASK

    class GUID(ctypes.Structure):
        Data1: ctypes._CField[Incomplete, Incomplete, Incomplete]
        Data2: ctypes._CField[Incomplete, Incomplete, Incomplete]
        Data3: ctypes._CField[Incomplete, Incomplete, Incomplete]
        Data4: ctypes._CField[Incomplete, Incomplete, Incomplete]

    class SP_DEVINFO_DATA(ctypes.Structure):
        cbSize: ctypes._CField[Incomplete, Incomplete, Incomplete]
        ClassGuid: ctypes._CField[Incomplete, Incomplete, Incomplete]
        DevInst: ctypes._CField[Incomplete, Incomplete, Incomplete]
        Reserved: ctypes._CField[Incomplete, Incomplete, Incomplete]

    PSP_DEVINFO_DATA: type[ctypes._Pointer[SP_DEVINFO_DATA]]
    PSP_DEVICE_INTERFACE_DETAIL_DATA = ctypes.c_void_p
    setupapi: ctypes.WinDLL
    SetupDiDestroyDeviceInfoList: ctypes._NamedFuncPointer
    SetupDiClassGuidsFromName: ctypes._NamedFuncPointer
    SetupDiEnumDeviceInfo: ctypes._NamedFuncPointer
    SetupDiGetClassDevs: ctypes._NamedFuncPointer
    SetupDiGetDeviceRegistryProperty: ctypes._NamedFuncPointer
    SetupDiGetDeviceInstanceId: ctypes._NamedFuncPointer
    SetupDiOpenDevRegKey: ctypes._NamedFuncPointer
    advapi32: ctypes.WinDLL
    RegCloseKey: ctypes._NamedFuncPointer
    RegQueryValueEx: ctypes._NamedFuncPointer
    cfgmgr32: ctypes.WinDLL
    CM_Get_Parent: ctypes._NamedFuncPointer
    CM_Get_Device_IDW: ctypes._NamedFuncPointer
    CM_MapCrToWin32Err: ctypes._NamedFuncPointer
    DIGCF_PRESENT: int
    DIGCF_DEVICEINTERFACE: int
    INVALID_HANDLE_VALUE: int
    ERROR_INSUFFICIENT_BUFFER: int
    ERROR_NOT_FOUND: int
    SPDRP_HARDWAREID: int
    SPDRP_FRIENDLYNAME: int
    SPDRP_LOCATION_PATHS: int
    SPDRP_MFG: int
    DICS_FLAG_GLOBAL: int
    DIREG_DEV: int
    KEY_READ: int
    MAX_USB_DEVICE_TREE_TRAVERSAL_DEPTH: int

    def get_parent_serial_number(
        child_devinst: ctypes._CData,
        child_vid: int | None,
        child_pid: int | None,
        depth: int = 0,
        last_serial_number: str | None = None,
    ) -> str: ...
    def iterate_comports() -> Generator[ListPortInfo]: ...
    def comports(include_links: bool = False) -> list[ListPortInfo]: ...
