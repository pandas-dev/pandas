from _typeshed import Incomplete
from collections.abc import Callable
from typing import Any, Protocol, TypeVar, type_check_only

_T = TypeVar("_T")

# at runtime, socketio.namespace.BaseNamespace, but socketio isn't py.typed
@type_check_only
class _BaseNamespace(Protocol):
    def is_asyncio_based(self) -> bool: ...
    def trigger_event(self, event: str, *args): ...

# at runtime, socketio.namespace.BaseNamespace, but socketio isn't py.typed
class _Namespace(_BaseNamespace, Protocol):
    def emit(
        self,
        event: str,
        data=None,
        to=None,
        room: str | None = None,
        skip_sid=None,
        namespace: str | None = None,
        callback: Callable[..., Incomplete] | None = None,
        ignore_queue: bool = False,
    ): ...
    def send(
        self,
        data,
        to=None,
        room: str | None = None,
        skip_sid=None,
        namespace: str | None = None,
        callback: Callable[..., Incomplete] | None = None,
        ignore_queue: bool = False,
    ) -> None: ...
    def call(
        self, event: str, data=None, to=None, sid=None, namespace: str | None = None, timeout=None, ignore_queue: bool = False
    ): ...
    def enter_room(self, sid, room: str, namespace: str | None = None): ...
    def leave_room(self, sid, room: str, namespace: str | None = None): ...
    def close_room(self, room: str, namespace: str | None = None): ...
    def rooms(self, sid, namespace: str | None = None): ...
    def get_session(self, sid, namespace: str | None = None): ...
    def save_session(self, sid, session, namespace: str | None = None): ...
    def session(self, sid, namespace: str | None = None): ...
    def disconnect(self, sid, namespace: str | None = None): ...

class Namespace(_Namespace):
    def __init__(self, namespace: str | None = None) -> None: ...
    def trigger_event(self, event: str, *args): ...
    def emit(  # type: ignore[override]
        self,
        event: str,
        data=None,
        room: str | None = None,
        include_self: bool = True,
        namespace: str | None = None,
        callback: Callable[..., _T] | None = None,
    ) -> _T | tuple[str, int]: ...
    def send(  # type: ignore[override]
        self,
        data,
        room: str | None = None,
        include_self: bool = True,
        namespace: str | None = None,
        callback: Callable[..., Any] | None = None,
    ) -> None: ...
    def close_room(self, room: str, namespace: str | None = None) -> None: ...
