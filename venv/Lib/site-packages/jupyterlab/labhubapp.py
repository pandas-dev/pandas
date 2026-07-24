"""A JupyterHub EntryPoint that defaults to use JupyterLab"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import os
import sys

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
    try:
        from jupyterhub.singleuser import SingleUserNotebookApp as SingleUserServerApp
    except ImportError as e:
        # jupyterhub is not installed at all
        venv_info = sys.prefix
        is_venv = sys.base_prefix != sys.prefix
        venv_type = "virtual environment" if is_venv else "Python environment"

        error_msg = (
            f"JupyterHub is not installed and is required to run this application.\n\n"
            f"Current {venv_type}: {venv_info}\n\n"
            f"Python sys.path entries searched:\n"
        )
        for path in sys.path:
            error_msg += f"  - {path}\n"
        error_msg += (
            f"\nTo fix this issue, install jupyterhub:\n"
            f"  pip install jupyterhub\n\n"
            f"Original error: {e}"
        )
        raise ImportError(error_msg) from e
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
