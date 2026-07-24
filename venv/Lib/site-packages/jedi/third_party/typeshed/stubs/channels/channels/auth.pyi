from _typeshed import Incomplete

from asgiref.typing import ASGIReceiveCallable, ASGISendCallable
from channels.middleware import BaseMiddleware
from django.contrib.auth.backends import BaseBackend
from django.contrib.auth.base_user import AbstractBaseUser
from django.contrib.auth.models import AnonymousUser
from django.utils.functional import LazyObject

from .consumer import _ChannelScope
from .utils import _ChannelApplication

async def get_user(scope: _ChannelScope) -> AbstractBaseUser | AnonymousUser: ...
async def login(scope: _ChannelScope, user: AbstractBaseUser, backend: BaseBackend | None = None) -> None: ...
async def logout(scope: _ChannelScope) -> None: ...

# Inherits AbstractBaseUser to improve autocomplete and show this is a lazy proxy for a user.
# At runtime, it's just a LazyObject that wraps the actual user instance.
class UserLazyObject(AbstractBaseUser, LazyObject[Incomplete]): ...

class AuthMiddleware(BaseMiddleware):
    def populate_scope(self, scope: _ChannelScope) -> None: ...
    async def resolve_scope(self, scope: _ChannelScope) -> None: ...
    async def __call__(
        self, scope: _ChannelScope, receive: ASGIReceiveCallable, send: ASGISendCallable
    ) -> _ChannelApplication: ...

def AuthMiddlewareStack(inner: _ChannelApplication) -> _ChannelApplication: ...
