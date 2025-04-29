"""tests of various mechanisms of using pyodide"""

import os
import shutil
from pathlib import Path

import pytest

from jupyterlite_pyodide_kernel.constants import PYODIDE_URL_ENV_VAR


@pytest.mark.parametrize(
    "approach,path,kind",
    [
        ["cli", "local", "archive"],
        ["cli", "local", "folder"],
        ["cli", "remote", "archive"],
        ["env", "local", "archive"],
        ["env", "local", "folder"],
        ["env", "remote", "archive"],
        ["wellknown", None, None],
    ],
)
def test_pyodide(
    an_empty_lite_dir,
    script_runner,
    a_pyodide_tarball,
    a_pyodide_server,
    approach,
    path,
    kind,
):  # pragma: no cover
    """can we fetch a pyodide archive, or use a local copy?"""
    env = dict(os.environ)
    pargs = []

    if approach == "wellknown":
        static = an_empty_lite_dir / "static"
        static.mkdir(parents=True, exist_ok=True)
        shutil.copytree(
            a_pyodide_tarball.parent / "pyodide/pyodide",
            static / "pyodide",
        )
    else:
        url = a_pyodide_tarball

        if path == "remote":
            url = f"{a_pyodide_server}/{Path(url).name}"
        elif kind == "folder":
            url = str(Path(url).parent / "pyodide")

        if approach == "env":
            env[PYODIDE_URL_ENV_VAR] = str(url)
        else:
            pargs += ["--pyodide", url]

    kwargs = dict(cwd=str(an_empty_lite_dir), env=env)

    status = script_runner.run(["jupyter", "lite", "status", *pargs], **kwargs)
    assert status.success, "status did NOT succeed"

    build = script_runner.run(["jupyter", "lite", "build", *pargs], **kwargs)
    assert build.success, "the build did NOT succeed"

    pyodide_path = an_empty_lite_dir / "_output/static/pyodide/pyodide.js"
    assert pyodide_path.exists(), "pyodide.js does not exist"

    check = script_runner.run(["jupyter", "lite", "check", *pargs], **kwargs)
    assert check.success, "the check did NOT succeed"
