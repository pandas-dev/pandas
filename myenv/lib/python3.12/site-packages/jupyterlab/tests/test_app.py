# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""A lab app that runs a sub process for a demo or a test."""

import atexit
import json
import os
import shutil
import sys
import tempfile

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

from os import path as osp
from os.path import join as pjoin
from stat import S_IRGRP, S_IROTH, S_IRUSR
from tempfile import TemporaryDirectory
from unittest.mock import patch

import jupyter_core
import jupyterlab_server
from ipykernel.kernelspec import write_kernel_spec
from jupyter_server.serverapp import ServerApp
from jupyterlab_server.process_app import ProcessApp
from traitlets import default

HERE = osp.realpath(osp.dirname(__file__))


def _create_template_dir():
    template_dir = tempfile.mkdtemp(prefix="mock_static")
    index_filepath = osp.join(template_dir, "index.html")
    with open(index_filepath, "w") as fid:
        fid.write(
            """
<!DOCTYPE HTML>
<html>
<head>
    <meta charset="utf-8">
    <title>{% block title %}Jupyter Lab Test{% endblock %}</title>
    <meta http-equiv="X-UA-Compatible" content="IE=edge" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    {% block meta %}
    {% endblock %}
</head>
<body>
  <h1>JupyterLab Test Application</h1>
  <div id="site">
    {% block site %}
    {% endblock site %}
  </div>
  {% block after_site %}
  {% endblock after_site %}
</body>
</html>"""
        )
    return template_dir


def _create_static_dir():
    static_dir = tempfile.mkdtemp(prefix="mock_static")
    return static_dir


def _create_schemas_dir():
    """Create a temporary directory for schemas."""
    root_dir = tempfile.mkdtemp(prefix="mock_schemas")
    extension_dir = osp.join(root_dir, "@jupyterlab", "apputils-extension")
    os.makedirs(extension_dir)

    # Get schema content.
    schema_package = jupyterlab_server.__name__
    schema_path = "tests/schemas/@jupyterlab/apputils-extension/themes.json"
    themes = files(schema_package).joinpath(schema_path).read_bytes()
    with open(osp.join(extension_dir, "themes.json"), "w") as fid:
        fid.write(themes.decode("utf-8"))
    atexit.register(lambda: shutil.rmtree(root_dir, True))
    return root_dir


def _create_user_settings_dir():
    """Create a temporary directory for workspaces."""
    root_dir = tempfile.mkdtemp(prefix="mock_user_settings")
    atexit.register(lambda: shutil.rmtree(root_dir, True))
    return root_dir


def _create_workspaces_dir():
    """Create a temporary directory for workspaces."""
    root_dir = tempfile.mkdtemp(prefix="mock_workspaces")
    atexit.register(lambda: shutil.rmtree(root_dir, True))
    return root_dir


class TestEnv:
    """Set Jupyter path variables to a temporary directory

    Useful as a context manager or with explicit start/stop
    """

    def start(self):
        self.test_dir = td = TemporaryDirectory()
        self.env_patch = patch.dict(
            os.environ,
            {
                "JUPYTER_CONFIG_DIR": pjoin(td.name, "jupyter"),
                "JUPYTER_DATA_DIR": pjoin(td.name, "jupyter_data"),
                "JUPYTER_RUNTIME_DIR": pjoin(td.name, "jupyter_runtime"),
                "IPYTHONDIR": pjoin(td.name, "ipython"),
            },
        )
        self.env_patch.start()
        self.path_patch = patch.multiple(
            jupyter_core.paths,
            SYSTEM_JUPYTER_PATH=[pjoin(td.name, "share", "jupyter")],
            ENV_JUPYTER_PATH=[pjoin(td.name, "env", "share", "jupyter")],
            SYSTEM_CONFIG_PATH=[pjoin(td.name, "etc", "jupyter")],
            ENV_CONFIG_PATH=[pjoin(td.name, "env", "etc", "jupyter")],
        )
        self.path_patch.start()

    def stop(self):
        self.env_patch.stop()
        self.path_patch.stop()
        try:
            self.test_dir.cleanup()
        except OSError:
            pass

    def __enter__(self):
        self.start()
        return self.test_dir.name

    def __exit__(self, *exc_info):
        self.stop()


class ProcessTestApp(ProcessApp):
    """A process app for running tests, includes a mock contents directory."""

    allow_origin = "*"

    def initialize_templates(self):
        self.static_paths = [_create_static_dir()]
        self.template_paths = [_create_template_dir()]

    def initialize_settings(self):
        self.env_patch = TestEnv()
        self.env_patch.start()
        ProcessApp.__init__(self)

        self.settings["allow_origin"] = ProcessTestApp.allow_origin

        self.static_dir = self.static_paths[0]
        self.template_dir = self.template_paths[0]
        self.schemas_dir = _create_schemas_dir()
        self.user_settings_dir = _create_user_settings_dir()
        self.workspaces_dir = _create_workspaces_dir()

        self._install_default_kernels()
        self.settings["kernel_manager"].default_kernel_name = "echo"

        super().initialize_settings()

    def _install_kernel(self, kernel_name, kernel_spec):
        """Install a kernel spec to the data directory.

        Parameters
        ----------
        kernel_name: str
            Name of the kernel.
        kernel_spec: dict
            The kernel spec for the kernel
        """
        paths = jupyter_core.paths
        kernel_dir = pjoin(paths.jupyter_data_dir(), "kernels", kernel_name)
        os.makedirs(kernel_dir)
        with open(pjoin(kernel_dir, "kernel.json"), "w") as f:
            f.write(json.dumps(kernel_spec))

    def _install_default_kernels(self):
        # Install echo and ipython kernels - should be done after env patch
        self._install_kernel(
            kernel_name="echo",
            kernel_spec={
                "argv": [
                    sys.executable,
                    "-m",
                    "jupyterlab.tests.echo_kernel",
                    "-f",
                    "{connection_file}",
                ],
                "display_name": "Echo Kernel",
                "language": "echo",
            },
        )

        paths = jupyter_core.paths
        ipykernel_dir = pjoin(paths.jupyter_data_dir(), "kernels", "ipython")
        write_kernel_spec(ipykernel_dir)

    def _process_finished(self, future):
        self.serverapp.http_server.stop()
        self.serverapp.io_loop.stop()
        self.env_patch.stop()
        try:
            os._exit(future.result())
        except Exception as e:
            self.log.error(str(e))
            os._exit(1)


class RootedServerApp(ServerApp):
    @default("root_dir")
    def _default_root_dir(self):
        """Create a temporary directory with some file structure."""
        root_dir = tempfile.mkdtemp(prefix="mock_root")
        os.mkdir(osp.join(root_dir, "src"))
        with open(osp.join(root_dir, "src", "temp.txt"), "w") as fid:
            fid.write("hello")

        readonly_filepath = osp.join(root_dir, "src", "readonly-temp.txt")
        with open(readonly_filepath, "w") as fid:
            fid.write("hello from a readonly file")

        os.chmod(readonly_filepath, S_IRUSR | S_IRGRP | S_IROTH)
        atexit.register(lambda: shutil.rmtree(root_dir, True))
        return root_dir
