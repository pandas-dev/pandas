"""pytest configuration for jupyterlite-core"""

import os
import shutil
import socket
import subprocess
import sys
import time
import warnings
from pathlib import Path

import pytest

from jupyterlite_core.constants import ALL_APP_ARCHIVES, NPM_SOURCE_DATE_EPOCH

HERE = Path(__file__).parent
FIXTURES = HERE / "fixtures"
WHEELS = [*FIXTURES.glob("*.whl")]
CONDA_PKGS = [*FIXTURES.glob("*.tar.bz2"), *FIXTURES.glob("*.conda")]


CI = os.environ.get("CI", None)
DARWIN = sys.platform.startswith("darwin")
LINUX = sys.platform.startswith("linux")
PYPY = "__pypy__" in sys.builtin_module_names


@pytest.fixture
def an_empty_lite_dir(tmp_path):
    lite_dir = tmp_path / "a_lite_dir"
    lite_dir.mkdir()

    yield lite_dir

    # for windows flake
    for retry in range(5):
        if lite_dir.exists():
            try:
                shutil.rmtree(lite_dir)
            except Exception as error:  # pragma: no covers
                warnings.warn(
                    f"Attempt {retry}: failed to clean up {lite_dir}: {error}",
                    stacklevel=2,
                )
                time.sleep(5)


@pytest.fixture
def source_date_epoch():
    """get a SOURCE_DATE_EPOCH

    loosely derived from https://reproducible-builds.org/docs/source-date-epoch/#python
    """
    now = int(time.time())
    print("SOURCE_DATE_EPOCH is", now)
    return f"{now}"


@pytest.fixture(params=sorted(ALL_APP_ARCHIVES))
def a_lite_app_archive(request):
    return request.param


@pytest.fixture
def the_npm_source_date_epoch():
    return NPM_SOURCE_DATE_EPOCH


@pytest.fixture
def a_simple_lite_ipynb():
    from nbformat.v4 import new_notebook, writes

    nb = new_notebook(
        metadata={
            "jupyter-lite": {
                "jupyter-config-data": {
                    "federated_extensions": [
                        {
                            "extension": "./extension",
                            "load": "static/remoteEntry.abc123.js",
                            "name": "@org/pkg",
                        }
                    ],
                    "disabledExtensions": ["@org/pkg"],
                    "settingsOverrides": {},
                }
            }
        }
    )
    return writes(nb)


@pytest.fixture
def an_unused_port():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("127.0.0.1", 0))
    sock.listen(1)
    port = sock.getsockname()[1]
    sock.close()
    return port


@pytest.fixture
def a_fixture_server(an_unused_port):
    p = subprocess.Popen(
        ["python", "-m", "http.server", "-b", "127.0.0.1", f"{an_unused_port}"],
        cwd=str(FIXTURES),
    )
    url = f"http://localhost:{an_unused_port}"
    yield url
    p.terminate()
