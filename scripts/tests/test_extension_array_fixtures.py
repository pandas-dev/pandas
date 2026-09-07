"""
Guard for the pandas base ExtensionArray test suite.

Third-party ExtensionArray authors run the tests in ``pandas/tests/extension/base``
against their own array by subclassing them, and copy
``pandas/tests/extension/conftest.py`` to get the fixtures those tests request.
That only works while the base tests stay clear of fixtures defined in
``pandas/conftest.py``, which pytest does not make available outside the pandas
test tree.

The check below copies the extension conftest and one of our own extension test
modules into a temporary directory -- so ``pandas/conftest.py`` is out of scope,
exactly as it is for a downstream user -- and runs them there.

See GH#56735.
"""

import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

import pandas as pd

extension_tests_path = Path(pd.__file__).parent / "tests" / "extension"


@pytest.mark.skipif(
    not (extension_tests_path / "conftest.py").exists(),
    reason="pandas installed without its test suite",
)
def test_base_extension_tests_run_without_pandas_conftest(tmp_path) -> None:
    for name in ["conftest.py", "test_datetime.py"]:
        shutil.copy(extension_tests_path / name, tmp_path / name)

    # PYTEST_ADDOPTS may name options that only pandas/conftest.py defines
    env = {key: value for key, value in os.environ.items() if key != "PYTEST_ADDOPTS"}
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", str(tmp_path)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        env=env,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(
            "The base ExtensionArray tests depend on fixtures that are not "
            "defined in pandas/tests/extension/conftest.py. Either define the "
            "fixture there as well, or use pytest.mark.parametrize instead.\n\n"
            f"{result.stdout}\n{result.stderr}"
        )
