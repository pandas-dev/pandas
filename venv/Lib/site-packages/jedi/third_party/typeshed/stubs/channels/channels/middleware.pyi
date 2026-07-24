from asgiref.typing import ASGIReceiveCallable, ASGISendCallable

from .consumer import _ChannelScope
from .utils import _ChannelApplication

class BaseMiddleware:
    inner: _ChannelApplication

    def __init__(self, inner: _ChannelApplication) -> None: ...
    async def __call__(
        self, scope: _ChannelScope, receive: ASGIReceiveCallable, send: ASGISendCallable
    ) -> _ChannelApplication: ...
