"""Tornado handlers for kernels.

Preliminary documentation at https://github.com/ipython/ipython/wiki/IPEP-16%3A-Notebook-multi-directory-dashboard-and-URL-mapping#kernels-api
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import json

try:
    from jupyter_client.jsonutil import json_default
except ImportError:
    from jupyter_client.jsonutil import date_default as json_default

from jupyter_core.utils import ensure_async
from tornado import web

from jupyter_server.auth.decorator import authorized
from jupyter_server.utils import url_escape, url_path_join

from ...base.handlers import APIHandler
from .websocket import KernelWebsocketHandler

AUTH_RESOURCE = "kernels"


class KernelsAPIHandler(APIHandler):
    """A kernels API handler."""

    auth_resource = AUTH_RESOURCE


class MainKernelHandler(KernelsAPIHandler):
    """The root kernel handler."""

    @web.authenticated
    @authorized
    async def get(self):
        """Get the list of running kernels."""
        km = self.kernel_manager
        kernels = await ensure_async(km.list_kernels())
        self.finish(json.dumps(kernels, default=json_default))

    @web.authenticated
    @authorized
    async def post(self):
        """Start a kernel."""
        km = self.kernel_manager
        model = self.get_json_body()
        if model is None:
            model = {"name": km.default_kernel_name}
        else:
            model.setdefault("name", km.default_kernel_name)

        kernel_id = await ensure_async(
            km.start_kernel(  # type:ignore[has-type]
                kernel_name=model["name"], path=model.get("path")
            )
        )
        model = await ensure_async(km.kernel_model(kernel_id))
        location = url_path_join(self.base_url, "api", "kernels", url_escape(kernel_id))
        self.set_header("Location", location)
        self.set_status(201)
        self.finish(json.dumps(model, default=json_default))


class KernelHandler(KernelsAPIHandler):
    """A kernel API handler."""

    @web.authenticated
    @authorized
    async def get(self, kernel_id):
        """Get a kernel model."""
        km = self.kernel_manager
        model = await ensure_async(km.kernel_model(kernel_id))
        self.finish(json.dumps(model, default=json_default))

    @web.authenticated
    @authorized
    async def delete(self, kernel_id):
        """Remove a kernel."""
        km = self.kernel_manager
        await ensure_async(km.shutdown_kernel(kernel_id))
        self.set_status(204)
        self.finish()


class KernelActionHandler(KernelsAPIHandler):
    """A kernel action API handler."""

    @web.authenticated
    @authorized
    async def post(self, kernel_id, action):
        """Interrupt or restart a kernel."""
        km = self.kernel_manager
        if action == "interrupt":
            await ensure_async(km.interrupt_kernel(kernel_id))  # type:ignore[func-returns-value]
            self.set_status(204)
        if action == "restart":
            try:
                await km.restart_kernel(kernel_id)
            except Exception:
                message = "Exception restarting kernel"
                self.log.error(message, exc_info=True)
                self.write(json.dumps({"message": message, "traceback": ""}))
                self.set_status(500)
            else:
                model = await ensure_async(km.kernel_model(kernel_id))
                self.write(json.dumps(model, default=json_default))
        self.finish()


# -----------------------------------------------------------------------------
# URL to handler mappings
# -----------------------------------------------------------------------------
_kernel_id_regex = r"(?P<kernel_id>\w+-\w+-\w+-\w+-\w+)"
_kernel_action_regex = r"(?P<action>restart|interrupt)"

default_handlers = [
    (r"/api/kernels", MainKernelHandler),
    (r"/api/kernels/%s" % _kernel_id_regex, KernelHandler),
    (
        rf"/api/kernels/{_kernel_id_regex}/{_kernel_action_regex}",
        KernelActionHandler,
    ),
    (r"/api/kernels/%s/channels" % _kernel_id_regex, KernelWebsocketHandler),
]
