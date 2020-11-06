from __future__ import annotations

from pathlib import Path
import subprocess

import pytest

FILE_0 = "cat_0 = Categorical()\ncat_1 = pd.Categorical()"
FILE_1 = "cat_0 = pd.Categorical()\ncat_1 = Categorical()"
FILE_2 = "cat_0 = Categorical()\ncat_1 = Categorical()"
FILE_3 = "cat_0 = pd.Categorical()\ncat_1 = pd.Categorical()"


@pytest.mark.parametrize("content", [FILE_0, FILE_1])
def test_inconsistent_usage(tmpdir, content: str) -> None:
    tmpfile = Path(tmpdir / "tmpfile.py")
    tmpfile.touch()
    tmpfile.write_text(content)
    output = subprocess.run(
        ["python", "scripts/check_for_inconsistent_pandas_namespace.py", str(tmpfile)],
        stderr=subprocess.PIPE,
    )

    # check stderr
    result = output.stderr.decode()
    expected = "Found both `pd.Categorical` and `Categorical` in"
    assert expected in result

    # check return code
    result = output.returncode
    expected = 1
    assert result == expected


@pytest.mark.parametrize("content", [FILE_2, FILE_3])
def test_consistent_usage(tmpdir, content: str) -> None:
    tmpfile = Path(tmpdir / "tmpfile.py")
    tmpfile.touch()
    tmpfile.write_text(content)
    output = subprocess.run(
        ["python", "scripts/check_for_inconsistent_pandas_namespace.py", str(tmpfile)],
        stderr=subprocess.PIPE,
    )

    # check stderr
    result = output.stderr.decode()
    expected = ""
    assert expected == result

    # check return code
    result = output.returncode
    expected = 0
    assert result == expected
