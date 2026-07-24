from _typeshed import Incomplete
from typing import Literal, TypedDict, type_check_only
from typing_extensions import NotRequired

from gunicorn.arbiter import Arbiter

@type_check_only
class _Worker(TypedDict):
    pid: int
    age: int
    booted: bool
    last_heartbeat: float
    aborted: NotRequired[bool]
    apps: NotRequired[list[str]]

@type_check_only
class _ShowWorkersReturnType(TypedDict):
    workers: list[_Worker]
    count: int

@type_check_only
class _App(TypedDict):
    import_path: str
    worker_count: int | None
    current_workers: int
    worker_pids: list[int]

@type_check_only
class _ShowDirtyReturnType(TypedDict):
    enabled: bool
    pid: int | None
    workers: list[_Worker]
    apps: list[_App]

@type_check_only
class _ShowStatsReturnType(TypedDict):
    uptime: float | None
    pid: int
    workers_current: int
    workers_target: int
    workers_spawned: int
    workers_killed: int
    reloads: int
    dirty_arbiter_pid: int | None

@type_check_only
class _ListenerInfo(TypedDict):
    address: str
    fd: int
    type: Literal["unix", "tcp", "tcp6", "unknown"]

@type_check_only
class _ShowListenersReturnType(TypedDict):
    listeners: list[_ListenerInfo]
    count: int

@type_check_only
class _WorkerAddReturnType(TypedDict):
    added: int
    previous: int
    total: int

@type_check_only
class _WorkerRemovedReturnType(TypedDict):
    removed: int
    previous: int
    total: int

@type_check_only
class _WorkerKillSucessReturnType(TypedDict):
    success: Literal[True]
    killed: int

@type_check_only
class _WorkerKillFailedReturnType(TypedDict):
    success: Literal[False]
    error: str

@type_check_only
class _ReloadReturnType(TypedDict):
    status: Literal["reloading"]

@type_check_only
class _ReopenReturnType(TypedDict):
    status: Literal["reopening"]

@type_check_only
class _ShutdownReturnType(TypedDict):
    status: Literal["shutting_down"]
    mode: str

@type_check_only
class _HelpReturnType(TypedDict):
    commands: dict[str, str]

class CommandHandlers:
    arbiter: Arbiter
    def __init__(self, arbiter: Arbiter) -> None: ...
    def show_workers(self) -> _ShowWorkersReturnType: ...
    def show_dirty(self) -> _ShowDirtyReturnType: ...
    def show_config(self) -> dict[str, Incomplete]: ...
    def show_stats(self) -> _ShowStatsReturnType: ...
    def show_listeners(self) -> _ShowListenersReturnType: ...
    def worker_add(self, count: int = 1) -> _WorkerAddReturnType: ...
    def worker_remove(self, count: int = 1) -> _WorkerRemovedReturnType: ...
    def worker_kill(self, pid: int) -> _WorkerKillSucessReturnType | _WorkerKillFailedReturnType: ...
    def dirty_add(self, count: int = 1) -> dict[str, Incomplete]: ...
    def dirty_remove(self, count: int = 1) -> dict[str, Incomplete]: ...
    def reload(self) -> _ReloadReturnType: ...
    def reopen(self) -> _ReopenReturnType: ...
    def shutdown(self, mode: str = "graceful") -> _ShutdownReturnType: ...
    def show_all(self) -> dict[str, Incomplete]: ...
    def help(self) -> _HelpReturnType: ...
