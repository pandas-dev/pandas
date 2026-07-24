from collections.abc import Callable
from typing import Any, ClassVar
from typing_extensions import TypeAlias

from channels.routing import ProtocolTypeRouter
from channels.utils import _ChannelApplication
from django.contrib.staticfiles.handlers import ASGIStaticFilesHandler
from django.test.testcases import TransactionTestCase

DaphneProcess: TypeAlias = Any  # TODO: temporary hack for daphne.testing.DaphneProcess; remove once daphne provides types

_StaticWrapper: TypeAlias = Callable[[ProtocolTypeRouter], _ChannelApplication]

def make_application(*, static_wrapper: _StaticWrapper | None) -> Any: ...
def set_database_connection() -> None: ...

class ChannelsLiveServerTestCase(TransactionTestCase):
    host: ClassVar[str] = "localhost"
    ProtocolServerProcess: ClassVar[type[DaphneProcess]] = ...
    static_wrapper: ClassVar[type[ASGIStaticFilesHandler]] = ...
    serve_static: ClassVar[bool] = True

    @property
    def live_server_url(self) -> str: ...
    @property
    def live_server_ws_url(self) -> str: ...
