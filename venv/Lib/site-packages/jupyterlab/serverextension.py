# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from jupyter_server.utils import url_path_join
from tornado.web import RedirectHandler


def load_jupyter_server_extension(serverapp):
    from .labapp import LabApp

    """Temporary server extension shim when using
    old notebook server.
    """
    extension = LabApp()
    extension.serverapp = serverapp
    extension.load_config_file()
    extension.update_config(serverapp.config)
    extension.parse_command_line(serverapp.extra_args)
    extension.handlers.extend(
        [
            (
                r"/static/favicons/favicon.ico",
                RedirectHandler,
                {"url": url_path_join(serverapp.base_url, "static/base/images/favicon.ico")},
            ),
            (
                r"/static/favicons/favicon-busy-1.ico",
                RedirectHandler,
                {"url": url_path_join(serverapp.base_url, "static/base/images/favicon-busy-1.ico")},
            ),
            (
                r"/static/favicons/favicon-busy-2.ico",
                RedirectHandler,
                {"url": url_path_join(serverapp.base_url, "static/base/images/favicon-busy-2.ico")},
            ),
            (
                r"/static/favicons/favicon-busy-3.ico",
                RedirectHandler,
                {"url": url_path_join(serverapp.base_url, "static/base/images/favicon-busy-3.ico")},
            ),
            (
                r"/static/favicons/favicon-file.ico",
                RedirectHandler,
                {"url": url_path_join(serverapp.base_url, "static/base/images/favicon-file.ico")},
            ),
            (
                r"/static/favicons/favicon-notebook.ico",
                RedirectHandler,
                {
                    "url": url_path_join(
                        serverapp.base_url, "static/base/images/favicon-notebook.ico"
                    )
                },
            ),
            (
                r"/static/favicons/favicon-terminal.ico",
                RedirectHandler,
                {
                    "url": url_path_join(
                        serverapp.base_url, "static/base/images/favicon-terminal.ico"
                    )
                },
            ),
            (
                r"/static/logo/logo.png",
                RedirectHandler,
                {"url": url_path_join(serverapp.base_url, "static/base/images/logo.png")},
            ),
        ]
    )
    extension.initialize()
