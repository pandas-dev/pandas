import builtins
import ctypes
import sys
from _typeshed import Incomplete
from types import TracebackType
from typing_extensions import Self

if sys.platform == "win32":
    def format_system_message(errno: int) -> str | None: ...

    class WindowsError(builtins.WindowsError):
        def __init__(self, value: int | None = None) -> None: ...
        @property
        def message(self) -> str: ...
        @property
        def code(self) -> int: ...

    def handle_nonzero_success(result: int) -> None: ...
    GMEM_MOVEABLE: int
    GlobalAlloc: Incomplete
    GlobalLock: Incomplete
    GlobalUnlock: Incomplete
    GlobalSize: Incomplete
    CreateFileMapping: Incomplete
    MapViewOfFile: Incomplete
    UnmapViewOfFile: Incomplete
    RtlMoveMemory: Incomplete

    class MemoryMap:
        name: str
        length: int
        security_attributes: Incomplete | None
        pos: int
        filemap: Incomplete
        view: Incomplete
        def __init__(self, name: str, length: int, security_attributes=None) -> None: ...
        def __enter__(self) -> Self: ...
        def seek(self, pos: int) -> None: ...
        def write(self, msg: bytes) -> None: ...
        def read(self, n: int) -> bytes: ...
        def __exit__(
            self, exc_type: type[BaseException] | None, exc_val: BaseException | None, tb: TracebackType | None
        ) -> None: ...

    READ_CONTROL: int
    STANDARD_RIGHTS_REQUIRED: int
    STANDARD_RIGHTS_READ: int
    STANDARD_RIGHTS_WRITE: int
    STANDARD_RIGHTS_EXECUTE: int
    STANDARD_RIGHTS_ALL: int
    POLICY_VIEW_LOCAL_INFORMATION: int
    POLICY_VIEW_AUDIT_INFORMATION: int
    POLICY_GET_PRIVATE_INFORMATION: int
    POLICY_TRUST_ADMIN: int
    POLICY_CREATE_ACCOUNT: int
    POLICY_CREATE_SECRET: int
    POLICY_CREATE_PRIVILEGE: int
    POLICY_SET_DEFAULT_QUOTA_LIMITS: int
    POLICY_SET_AUDIT_REQUIREMENTS: int
    POLICY_AUDIT_LOG_ADMIN: int
    POLICY_SERVER_ADMIN: int
    POLICY_LOOKUP_NAMES: int
    POLICY_NOTIFICATION: int
    POLICY_ALL_ACCESS: int
    POLICY_READ: int
    POLICY_WRITE: int
    POLICY_EXECUTE: int

    class TokenAccess:
        TOKEN_QUERY: int

    class TokenInformationClass:
        TokenUser: int

    class TOKEN_USER(ctypes.Structure):
        num: int
        SID: Incomplete
        ATTRIBUTES: Incomplete

    class SECURITY_DESCRIPTOR(ctypes.Structure):
        SECURITY_DESCRIPTOR_CONTROL: Incomplete
        REVISION: int
        Revision: int
        Sbz1: Incomplete
        Control: Incomplete
        Owner: Incomplete
        Group: Incomplete
        Sacl: Incomplete
        Dacl: Incomplete

    class SECURITY_ATTRIBUTES(ctypes.Structure):
        nLength: int
        lpSecurityDescriptor: int
        bInheritHandle: bool
        def __init__(self, *args, **kwargs) -> None: ...
        @property
        def descriptor(self): ...
        @descriptor.setter
        def descriptor(self, value) -> None: ...

    def GetTokenInformation(token, information_class): ...
    def OpenProcessToken(proc_handle, access): ...
    def get_current_user() -> TOKEN_USER: ...
    def get_security_attributes_for_user(user: TOKEN_USER | None = None) -> SECURITY_ATTRIBUTES: ...
