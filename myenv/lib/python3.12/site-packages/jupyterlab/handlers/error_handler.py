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
%s
</body>
"""


class ErrorHandler(ExtensionHandlerMixin, JupyterHandler):
    def initialize(self, messages=None, name=None):
        super().initialize(name=name)
        self.messages = messages

    @web.authenticated
    @web.removeslash
    def get(self):
        msgs = ["<h2>%s</h2>" % msg for msg in self.messages]
        self.write(TEMPLATE % "\n".join(msgs))
