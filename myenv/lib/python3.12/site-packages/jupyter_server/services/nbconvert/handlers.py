"""API Handlers for nbconvert."""

import asyncio
import json

from anyio.to_thread import run_sync
from tornado import web

from jupyter_server.auth.decorator import authorized

from ...base.handlers import APIHandler

AUTH_RESOURCE = "nbconvert"


class NbconvertRootHandler(APIHandler):
    """The nbconvert root API handler."""

    auth_resource = AUTH_RESOURCE
    _exporter_lock: asyncio.Lock

    def initialize(self, **kwargs):
        """Initialize an nbconvert root handler."""
        super().initialize(**kwargs)
        # share lock across instances of this handler class
        if not hasattr(self.__class__, "_exporter_lock"):
            self.__class__._exporter_lock = asyncio.Lock()
        self._exporter_lock = self.__class__._exporter_lock

    @web.authenticated
    @authorized
    async def get(self):
        """Get the list of nbconvert exporters."""
        try:
            from nbconvert.exporters import base
        except ImportError as e:
            raise web.HTTPError(500, "Could not import nbconvert: %s" % e) from e
        res = {}
        # Some exporters use the filesystem when instantiating, delegate that
        # to a thread so we don't block the event loop for it.
        exporters = await run_sync(base.get_export_names)
        async with self._exporter_lock:
            for exporter_name in exporters:
                try:
                    exporter_class = await run_sync(base.get_exporter, exporter_name)
                except ValueError:
                    # I think the only way this will happen is if the entrypoint
                    # is uninstalled while this method is running
                    continue
                # XXX: According to the docs, it looks like this should be set to None
                # if the exporter shouldn't be exposed to the front-end and a friendly
                # name if it should. However, none of the built-in exports have it defined.
                # if not exporter_class.export_from_notebook:
                #    continue
                res[exporter_name] = {
                    "output_mimetype": exporter_class.output_mimetype,
                }

        self.finish(json.dumps(res))


default_handlers = [
    (r"/api/nbconvert", NbconvertRootHandler),
]
