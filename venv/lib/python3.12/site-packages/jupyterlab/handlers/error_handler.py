"""An error handler for JupyterLab."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from jupyter_server.base.handlers import JupyterHandler
from jupyter_server.extension.handler import ExtensionHandlerMixin
from tornado import web

TEMPLATE = """
<!DOCTYPE HTML>
<html>
<head>
    <meta charset="utf-8">
    <title>JupyterLab Error</title>
</head>
<body>
<h1>JupyterLab Error<h1>
{messages}
</body>
"""


class ErrorHandler(ExtensionHandlerMixin, JupyterHandler):
    def initialize(self, messages=None, name=None):
        super().initialize(name=name)
        self.messages = messages

    @web.authenticated
    @web.removeslash
    def get(self):
        msgs = [f"<h2>{msg}</h2>" for msg in self.messages]
        self.write(TEMPLATE.format(messages="\n".join(msgs)))
