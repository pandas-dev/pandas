"""A JupyterHub EntryPoint that defaults to use JupyterLab"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import os

from jupyter_server.serverapp import ServerApp
from traitlets import default

from .labapp import LabApp

if not os.environ.get("JUPYTERHUB_SINGLEUSER_APP"):
    # setting this env prior to import of jupyterhub.singleuser avoids unnecessary import of notebook
    os.environ["JUPYTERHUB_SINGLEUSER_APP"] = "jupyter_server.serverapp.ServerApp"

try:
    from jupyterhub.singleuser.mixins import make_singleuser_app
except ImportError:
    # backward-compat with jupyterhub < 1.3
    from jupyterhub.singleuser import SingleUserNotebookApp as SingleUserServerApp
else:
    SingleUserServerApp = make_singleuser_app(ServerApp)


class SingleUserLabApp(SingleUserServerApp):
    @default("default_url")
    def _default_url(self):
        return "/lab"

    def find_server_extensions(self):
        """unconditionally enable jupyterlab server extension

        never called if using legacy SingleUserNotebookApp
        """
        super().find_server_extensions()
        self.jpserver_extensions[LabApp.get_extension_package()] = True


def main(argv=None):
    return SingleUserLabApp.launch_instance(argv)


if __name__ == "__main__":
    main()
