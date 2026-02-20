"""HTTP handler to shut down the Jupyter server."""

from tornado import ioloop, web

from jupyter_server.auth.decorator import authorized
from jupyter_server.base.handlers import JupyterHandler

AUTH_RESOURCE = "server"


class ShutdownHandler(JupyterHandler):
    """A shutdown API handler."""

    auth_resource = AUTH_RESOURCE

    @web.authenticated
    @authorized
    async def post(self):
        """Shut down the server."""
        self.log.info("Shutting down on /api/shutdown request.")

        if self.serverapp:
            await self.serverapp._cleanup()

        ioloop.IOLoop.current().stop()


default_handlers = [
    (r"/api/shutdown", ShutdownHandler),
]
