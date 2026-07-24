from _typeshed import Incomplete
from typing import TypedDict, type_check_only

@type_check_only
class _DirtyErrorDict(TypedDict):
    error_type: str
    message: str
    details: dict[str, Incomplete]

class DirtyError(Exception):
    message: str
    details: dict[str, Incomplete]

    def __init__(self, message: str, details: dict[str, Incomplete] | None = None) -> None: ...
    def to_dict(self) -> _DirtyErrorDict: ...
    @classmethod
    def from_dict(cls, data: dict[str, Incomplete]) -> DirtyError: ...

class DirtyTimeoutError(DirtyError):
    timeout: float | None

    def __init__(self, message: str = "Operation timed out", timeout: float | None = None) -> None: ...

class DirtyConnectionError(DirtyError):
    socket_path: str | None

    def __init__(self, message: str = "Connection failed", socket_path: str | None = None) -> None: ...

class DirtyWorkerError(DirtyError):
    worker_id: int | None
    traceback: str | None

    def __init__(self, message: str, worker_id: int | None = None, traceback: str | None = None) -> None: ...

class DirtyAppError(DirtyError):
    app_path: str | None
    action: str | None
    traceback: str | None

    def __init__(
        self, message: str, app_path: str | None = None, action: str | None = None, traceback: str | None = None
    ) -> None: ...

class DirtyAppNotFoundError(DirtyAppError):
    app_path: str

    def __init__(self, app_path: str) -> None: ...

class DirtyNoWorkersAvailableError(DirtyError):
    app_path: str

    def __init__(self, app_path: str, message: str | None = None) -> None: ...

class DirtyProtocolError(DirtyError):
    def __init__(self, message: str = "Protocol error", raw_data: str | bytes | None = None) -> None: ...
