"""test configuration for jupyterlite-pyodide-kernel"""

from pathlib import Path
import sys
import pytest
import jupyterlite_core.tests.conftest
from jupyterlite_core.tests.conftest import (
    a_fixture_server,
    an_empty_lite_dir,
    an_unused_port,
)


from jupyterlite_pyodide_kernel.constants import (
    PYODIDE_VERSION,
    PYODIDE_KERNEL_NPM_NAME,
)

__all__ = [
    "a_fixture_server",
    "a_pyodide_server",
    "a_pyodide_tarball",
    "an_empty_lite_dir",
    "an_unused_port",
    "index_cmd",
]

HERE = Path(__file__).parent
FIXTURES = HERE / "fixtures"

# patch the upstream fixtures with local ones
jupyterlite_core.tests.conftest.FIXTURES = FIXTURES


IN_TREE_EXTENSION = (HERE / "../labextension").resolve()
SHARE = Path(sys.prefix) / "share/jupyter/labextensions"
IN_SHARE_EXTENSION = (SHARE / PYODIDE_KERNEL_NPM_NAME).resolve()
# prefer testing an in-tree extension, if available
PYODIDE_KERNEL_EXTENSION = (
    IN_TREE_EXTENSION if IN_TREE_EXTENSION.exists() else IN_SHARE_EXTENSION
)

WHEELS = [*FIXTURES.glob("*.whl")]

PYODIDE_GH = "https://github.com/pyodide/pyodide/releases/download"
PYODIDE_TARBALL = f"pyodide-core-{PYODIDE_VERSION}.tar.bz2"
PYODIDE_URL = f"{PYODIDE_GH}/{PYODIDE_VERSION}/{PYODIDE_TARBALL}"
PYODIDE_FIXTURE = FIXTURES / ".pyodide" / PYODIDE_VERSION / PYODIDE_TARBALL


@pytest.fixture
def index_cmd():
    """get the command line arguments for indexing a folder."""
    return ("jupyter", "piplite", "index")


@pytest.fixture
def a_pyodide_tarball():
    """maybe fetch the pyodide archive"""
    if not PYODIDE_FIXTURE.exists():  # pragma: no cover
        import shutil
        import urllib.request

        PYODIDE_FIXTURE.parent.mkdir(exist_ok=True, parents=True)
        with urllib.request.urlopen(PYODIDE_URL) as response:
            with PYODIDE_FIXTURE.open("wb") as fd:
                shutil.copyfileobj(response, fd)

    unpacked = PYODIDE_FIXTURE.parent / "pyodide/pyodide"

    if not unpacked.is_dir():  # pragma: no cover
        from jupyterlite_core.addons.base import BaseAddon
        from jupyterlite_core.manager import LiteManager

        manager = LiteManager()
        BaseAddon(manager=manager).extract_one(PYODIDE_FIXTURE, unpacked.parent)

    return PYODIDE_FIXTURE


@pytest.fixture
def a_pyodide_server(an_unused_port, a_pyodide_tarball):  # pragma: no cover
    """serve up the pyodide archive"""
    import subprocess

    root = a_pyodide_tarball.parent

    p = subprocess.Popen(
        ["python", "-m", "http.server", "-b", "127.0.0.1", f"{an_unused_port}"],
        cwd=str(root),
    )
    url = f"http://localhost:{an_unused_port}"
    yield url
    p.terminate()
