"""Tornado handlers for plugin management."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import dataclasses
import json

from jupyter_server.base.handlers import APIHandler
from tornado import web

from jupyterlab.extensions.manager import PluginManager


class PluginHandler(APIHandler):
    def initialize(self, manager: PluginManager):
        super().initialize()
        self.manager = manager

    @web.authenticated
    async def get(self):
        """GET query returns info on plugins locks"""
        # note: this is informative only - validation is server-side
        locks = await self.manager.plugin_locks()
        self.set_status(200)
        self.finish(json.dumps(locks))

    @web.authenticated
    async def post(self):
        """POST query performs an action on a specific plugin

        Body arguments:
            {
                "cmd": Action to perform - ["enable", "disable"]
                "plugin_name": Plugin name
            }
        """
        data = self.get_json_body()
        cmd = data["cmd"]
        name = data["plugin_name"]
        if cmd not in ("enable", "disable") or not name:
            raise web.HTTPError(
                422,
                f"Could not process instruction {cmd!r} with plugin name {name!r}",
            )

        ret_value = None
        try:
            if cmd == "enable":
                ret_value = await self.manager.enable(name)
            elif cmd == "disable":
                ret_value = await self.manager.disable(name)
        except Exception as e:
            raise web.HTTPError(500, str(e)) from e

        if ret_value.status == "error":
            self.set_status(500)
        else:
            self.set_status(201)
        self.finish(json.dumps(dataclasses.asdict(ret_value)))


# The path for lab plugins handler.
plugins_handler_path = r"/lab/api/plugins"
