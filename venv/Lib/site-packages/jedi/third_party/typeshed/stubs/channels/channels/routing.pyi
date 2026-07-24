from typing import Any, type_check_only

from asgiref.typing import ASGIReceiveCallable, ASGISendCallable, Scope
from django.urls.resolvers import URLPattern

from .consumer import _ASGIApplicationProtocol, _ChannelScope
from .utils import _ChannelApplication

def get_default_application() -> ProtocolTypeRouter: ...

class ProtocolTypeRouter:
    application_mapping: dict[str, _ChannelApplication]

    def __init__(self, application_mapping: dict[str, Any]) -> None: ...
    async def __call__(self, scope: Scope, receive: ASGIReceiveCallable, send: ASGISendCallable) -> None: ...

@type_check_only
class _ExtendedURLPattern(URLPattern):
    callback: _ASGIApplicationProtocol | URLRouter

class URLRouter:
    routes: list[_ExtendedURLPattern | URLRouter]

    def __init__(self, routes: list[_ExtendedURLPattern | URLRouter]) -> None: ...
    async def __call__(self, scope: Scope, receive: ASGIReceiveCallable, send: ASGISendCallable) -> None: ...

class ChannelNameRouter:
    application_mapping: dict[str, _ChannelApplication]

    def __init__(self, application_mapping: dict[str, _ChannelApplication]) -> None: ...
    async def __call__(self, scope: _ChannelScope, receive: ASGIReceiveCallable, send: ASGISendCallable) -> None: ...
