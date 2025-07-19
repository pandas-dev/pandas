"""Tornado handlers for logging out of the Jupyter Server."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from ..base.handlers import JupyterHandler
from .decorator import allow_unauthenticated


class LogoutHandler(JupyterHandler):
    """An auth logout handler."""

    @allow_unauthenticated
    def get(self):
        """Handle a logout."""
        self.identity_provider.clear_login_cookie(self)
        if self.login_available:
            message = {"info": "Successfully logged out."}
        else:
            message = {"warning": "Cannot log out. Jupyter Server authentication is disabled."}
        self.write(self.render_template("logout.html", message=message))


default_handlers = [(r"/logout", LogoutHandler)]
