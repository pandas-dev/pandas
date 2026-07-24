from _typeshed import Incomplete
from collections.abc import Awaitable
from typing import Any, ClassVar, Protocol, TypedDict, type_check_only

from asgiref.typing import ASGIReceiveCallable, ASGISendCallable, Scope, WebSocketScope
from channels.auth import UserLazyObject
from channels.layers import BaseChannelLayer
from django.contrib.sessions.backends.base import SessionBase
from django.utils.functional import LazyObject

# _LazySession is a LazyObject that wraps a SessionBase instance.
# We subclass both for type checking purposes to expose SessionBase attributes,
# and suppress mypy's "misc" error with `# type: ignore[misc]`.
@type_check_only
class _LazySession(SessionBase, LazyObject[Incomplete]):  # type: ignore[misc]
    _wrapped: SessionBase

@type_check_only
class _URLRoute(TypedDict):
    # Values extracted from Django's URLPattern matching,
    # passed through ASGI scope routing.
    # `args` and `kwargs` are the result of pattern matching against the URL path.
    args: tuple[Any, ...]
    kwargs: dict[str, Any]

# Channel Scope definition
@type_check_only
class _ChannelScope(WebSocketScope, total=False):
    # Channels specific
    channel: str
    url_route: _URLRoute
    path_remaining: str

    # Auth specific
    cookies: dict[str, str]
    session: _LazySession
    user: UserLazyObject | None

# Accepts any ASGI message dict with a required "type" key (str),
# but allows additional arbitrary keys for flexibility.
def get_handler_name(message: dict[str, Any]) -> str: ...
@type_check_only
class _ASGIApplicationProtocol(Protocol):
    consumer_class: AsyncConsumer

    # Accepts any initialization kwargs passed to the consumer class.
    # Typed as `Any` to allow flexibility in subclass-specific arguments.
    consumer_initkwargs: Any

    def __call__(self, scope: Scope, receive: ASGIReceiveCallable, send: ASGISendCallable) -> Awaitable[None]: ...

class AsyncConsumer:
    channel_layer_alias: ClassVar[str]

    scope: _ChannelScope
    channel_layer: BaseChannelLayer
    channel_name: str
    channel_receive: ASGIReceiveCallable
    base_send: ASGISendCallable

    async def __call__(self, scope: _ChannelScope, receive: ASGIReceiveCallable, send: ASGISendCallable) -> None: ...
    async def dispatch(self, message: dict[str, Any]) -> None: ...
    async def send(self, message: dict[str, Any]) -> None: ...

    # initkwargs will be used to instantiate the consumer instance.
    @classmethod
    def as_asgi(cls, **initkwargs: Any) -> _ASGIApplicationProtocol: ...

class SyncConsumer(AsyncConsumer):

    # Since we're overriding asynchronous methods with synchronous ones,
    # we need to use `# type: ignore[override]` to suppress mypy errors.
    async def dispatch(self, message: dict[str, Any]) -> None: ...  # type: ignore[override]
    def send(self, message: dict[str, Any]) -> None: ...  # type: ignore[override]
