"""Test utility functions for hooks."""

from pathlib import Path

import pytest

from numpydoc.hooks import utils


@pytest.mark.parametrize(
    ["reason_file", "files", "expected_reason"],
    [
        (None, None, "current directory"),
        (None, ["x.py"], "file system root"),
        (".git", ["x.py"], "version control"),
        ("pyproject.toml", ["x.py"], "pyproject.toml"),
        ("setup.cfg", ["x.py"], "setup.cfg"),
    ],
)
def test_find_project_root(tmp_path, request, reason_file, files, expected_reason):
    """Test the process of finding the project root."""
    if reason_file:
        (tmp_path / reason_file).touch()

    if files:
        expected_dir = (
            Path(tmp_path.anchor) if expected_reason == "file system root" else tmp_path
        )

        for file in files:
            (tmp_path / file).touch()
    else:
        expected_dir = request.config.rootdir

    root_dir, reason = utils.find_project_root(files if not files else [tmp_path])
    assert reason == expected_reason
    assert root_dir == expected_dir
