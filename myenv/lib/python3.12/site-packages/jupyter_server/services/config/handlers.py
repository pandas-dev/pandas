"""Tornado handlers for frontend config storage."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
import json

from tornado import web

from jupyter_server.auth.decorator import authorized

from ...base.handlers import APIHandler

AUTH_RESOURCE = "config"


class ConfigHandler(APIHandler):
    """A config API handler."""

    auth_resource = AUTH_RESOURCE

    @web.authenticated
    @authorized
    def get(self, section_name):
        """Get config by section name."""
        self.set_header("Content-Type", "application/json")
        self.finish(json.dumps(self.config_manager.get(section_name)))

    @web.authenticated
    @authorized
    def put(self, section_name):
        """Set a config section by name."""
        data = self.get_json_body()  # Will raise 400 if content is not valid JSON
        self.config_manager.set(section_name, data)
        self.set_status(204)

    @web.authenticated
    @authorized
    def patch(self, section_name):
        """Update a config section by name."""
        new_data = self.get_json_body()
        section = self.config_manager.update(section_name, new_data)
        self.finish(json.dumps(section))


# URL to handler mappings

section_name_regex = r"(?P<section_name>\w+)"

default_handlers = [
    (r"/api/config/%s" % section_name_regex, ConfigHandler),
]
