import datetime
from collections.abc import Awaitable
from typing import Any

from asgiref.typing import ASGIReceiveCallable, ASGISendCallable
from channels.consumer import _ChannelScope
from channels.utils import _ChannelApplication
from django.contrib.sessions.backends.base import SessionBase

class CookieMiddleware:
    inner: _ChannelApplication

    def __init__(self, inner: _ChannelApplication) -> None: ...

    # Returns the same type as the provided _ChannelApplication.
    async def __call__(self, scope: _ChannelScope, receive: ASGIReceiveCallable, send: ASGISendCallable) -> Any: ...
    @classmethod
    def set_cookie(
        cls,
        message: dict[str, Any],
        key: str,
        value: str = "",
        max_age: int | None = None,
        expires: str | datetime.datetime | None = None,
        path: str = "/",
        domain: str | None = None,
        secure: bool = False,
        httponly: bool = False,
        samesite: str = "lax",
    ) -> None: ...
    @classmethod
    def delete_cookie(cls, message: dict[str, Any], key: str, path: str = "/", domain: str | None = None) -> None: ...

class InstanceSessionWrapper:
    save_message_types: list[str]
    cookie_response_message_types: list[str]
    cookie_name: str
    session_store: SessionBase
    scope: _ChannelScope
    activated: bool
    real_send: ASGISendCallable

    def __init__(self, scope: _ChannelScope, send: ASGISendCallable) -> None: ...
    async def resolve_session(self) -> None: ...
    async def send(self, message: dict[str, Any]) -> Awaitable[None]: ...
    async def save_session(self) -> None: ...

class SessionMiddleware:
    inner: _ChannelApplication

    def __init__(self, inner: _ChannelApplication) -> None: ...
    async def __call__(
        self, scope: _ChannelScope, receive: ASGIReceiveCallable, send: ASGISendCallable
    ) -> _ChannelApplication: ...

def SessionMiddlewareStack(inner: _ChannelApplication) -> _ChannelApplication: ...
