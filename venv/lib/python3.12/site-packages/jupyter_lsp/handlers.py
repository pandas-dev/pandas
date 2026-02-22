""" tornado handler for managing and communicating with language servers
"""

from typing import Optional, Text

from jupyter_core.utils import ensure_async
from jupyter_server.base.handlers import APIHandler, JupyterHandler
from jupyter_server.utils import url_path_join as ujoin
from tornado import web
from tornado.websocket import WebSocketHandler

try:
    from jupyter_server.auth.decorator import authorized
except ImportError:

    def authorized(method):  # type: ignore
        """A no-op fallback for `jupyter_server 1.x`"""
        return method


try:
    from jupyter_server.base.websocket import WebSocketMixin
except ImportError:
    from jupyter_server.base.zmqhandlers import WebSocketMixin

from .manager import LanguageServerManager
from .schema import SERVERS_RESPONSE
from .specs.utils import censored_spec

AUTH_RESOURCE = "lsp"


class BaseHandler(APIHandler):
    manager = None  # type: LanguageServerManager

    def initialize(self, manager: LanguageServerManager):
        self.manager = manager


class BaseJupyterHandler(JupyterHandler):
    manager = None  # type: LanguageServerManager

    def initialize(self, manager: LanguageServerManager):
        self.manager = manager


class LanguageServerWebSocketHandler(  # type: ignore
    WebSocketMixin, WebSocketHandler, BaseJupyterHandler
):
    """Setup tornado websocket to route to language server sessions.

    The logic of `get` and `pre_get` methods is derived from jupyter-server ws handlers,
    and should be kept in sync to follow best practice established by upstream; see:
    https://github.com/jupyter-server/jupyter_server/blob/v2.12.5/jupyter_server/services/kernels/websocket.py#L36
    """

    auth_resource = AUTH_RESOURCE

    language_server: Optional[Text] = None

    async def pre_get(self):
        """Handle a pre_get."""
        # authenticate first
        # authenticate the request before opening the websocket
        user = self.current_user
        if user is None:
            self.log.warning("Couldn't authenticate WebSocket connection")
            raise web.HTTPError(403)

        if not hasattr(self, "authorizer"):
            return

        # authorize the user.
        is_authorized = await ensure_async(
            self.authorizer.is_authorized(self, user, "execute", AUTH_RESOURCE)
        )
        if not is_authorized:
            raise web.HTTPError(403)

    async def get(self, *args, **kwargs):
        """Get an event socket."""
        await self.pre_get()
        res = super().get(*args, **kwargs)
        if res is not None:
            await res

    async def open(self, language_server):
        await self.manager.ready()
        self.language_server = language_server
        self.manager.subscribe(self)
        self.log.debug("[{}] Opened a handler".format(self.language_server))
        super().open()

    async def on_message(self, message):
        self.log.debug("[{}] Handling a message".format(self.language_server))
        await self.manager.on_client_message(message, self)

    def on_close(self):
        self.manager.unsubscribe(self)
        self.log.debug("[{}] Closed a handler".format(self.language_server))


class LanguageServersHandler(BaseHandler):
    """Reports the status of all current servers

    Response should conform to schema in schema/servers.schema.json
    """

    auth_resource = AUTH_RESOURCE
    validator = SERVERS_RESPONSE

    @web.authenticated
    @authorized
    async def get(self):
        """finish with the JSON representations of the sessions"""
        await self.manager.ready()

        response = {
            "version": 2,
            "sessions": {
                language_server: session.to_json()
                for language_server, session in self.manager.sessions.items()
            },
            "specs": {
                key: censored_spec(spec)
                for key, spec in self.manager.all_language_servers.items()
            },
        }

        errors = list(self.validator.iter_errors(response))

        if errors:  # pragma: no cover
            self.log.warning("{} validation errors: {}".format(len(errors), errors))

        self.finish(response)


def add_handlers(nbapp):
    """Add Language Server routes to the notebook server web application"""
    lsp_url = ujoin(nbapp.base_url, "lsp")
    re_langservers = "(?P<language_server>.*)"

    opts = {"manager": nbapp.language_server_manager}

    nbapp.web_app.add_handlers(
        ".*",
        [
            (ujoin(lsp_url, "status"), LanguageServersHandler, opts),
            (
                ujoin(lsp_url, "ws", re_langservers),
                LanguageServerWebSocketHandler,
                opts,
            ),
        ],
    )
