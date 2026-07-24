"""tests for more kinds of contents"""

import json

import pytest


@pytest.mark.parametrize(
    "allow_hidden,expect_success,expect_content,extra_ignore",
    [
        [True, True, True, []],
        [False, False, False, []],
        [False, True, False, [r"/\.binder/"]],
    ],
)
def test_contents_with_dot(  # noqa: PLR0913
    allow_hidden,
    expect_success,
    expect_content,
    extra_ignore,
    an_empty_lite_dir,
    script_runner,
):
    """Can hidden files be exposed with contents?"""
    config = {
        "LiteBuildConfig": {
            "ignore_sys_prefix": True,
            "extra_ignore_contents": extra_ignore,
        },
        "ContentsManager": {"allow_hidden": allow_hidden},
    }
    print("config", config)
    (an_empty_lite_dir / "jupyter_lite_config.json").write_text(json.dumps(config))
    dot_binder = an_empty_lite_dir / ".binder"
    dot_binder.mkdir()
    postbuild = dot_binder / "postBuild"
    postbuild.write_text("#!/usr/bin/env bash\necho ok")

    result = script_runner.run(
        ["jupyter", "lite", "build", "--contents", "."],
        cwd=str(an_empty_lite_dir),
    )
    if expect_success:
        assert result.success
    else:
        assert not result.success
        assert "jupyter_lite_config" in result.stdout

    out = an_empty_lite_dir / "_output"
    root_contents_json = out / "api/contents/all.json"
    out_postbuild = out / "files/.binder/postBuild"
    hidden_contents_json = out / "api/contents/.binder/all.json"

    if expect_content:
        root_contents = json.loads(root_contents_json.read_text(encoding="utf-8"))
        assert len(root_contents["content"]) == 1, root_contents
        assert out_postbuild.exists()

        hidden_contents = json.loads(hidden_contents_json.read_text(encoding="utf-8"))

        postbuild_content = hidden_contents["content"][0]
        assert postbuild_content["name"] == "postBuild", postbuild_content
        assert postbuild_content["path"] == ".binder/postBuild", postbuild_content
    else:
        assert not root_contents_json.exists()


def test_contents_with_space(
    an_empty_lite_dir,
    script_runner,
):
    dir_name = "dir with spaces"
    contents_dir = an_empty_lite_dir / "contents" / dir_name
    contents_dir.mkdir(parents=True)
    file_name = "file name with spaces"
    contents_file = contents_dir / file_name
    contents_file.touch()

    result = script_runner.run(
        ["jupyter", "lite", "build", "--contents", "contents"],
        cwd=str(an_empty_lite_dir),
    )
    assert result.success

    out = an_empty_lite_dir / "_output"
    root_contents_json = out / "api/contents/all.json"
    contents_json = out / f"api/contents/{dir_name}/all.json"

    root_contents = json.loads(root_contents_json.read_text(encoding="utf-8"))
    assert len(root_contents["content"]) == 1, root_contents
    contents = json.loads(contents_json.read_text(encoding="utf-8"))
    content = contents["content"][0]
    assert content["name"] == file_name
    assert content["path"] == f"{dir_name}/{file_name}"


def test_contents_missing_jupyter_server(
    an_empty_lite_dir,
    script_runner,
    monkeypatch,
):
    """
    Test that an error is raised when contents are provided but jupyter_server is not installed
    """
    # Create a test file to be used as contents
    test_contents = an_empty_lite_dir / "test_contents"
    test_contents.mkdir()
    (test_contents / "test_file.txt").write_text("Test content")

    # Set environment variable to simulate jupyter_server not being installed
    monkeypatch.setenv("JUPYTERLITE_NO_JUPYTER_SERVER", "true")

    # Run the build command with contents
    result = script_runner.run(
        ["jupyter", "lite", "build", "--contents", "test_contents"],
        cwd=str(an_empty_lite_dir),
    )

    # The build should fail
    assert not result.success

    # Check if the expected error message is in the output
    expected_error = (
        "jupyter-server is not installed. You cannot add custom content to jupyterlite."
    )
    assert expected_error in result.stdout or expected_error in result.stderr


def test_contents_resolved_relative_to_lite_dir(
    an_empty_lite_dir,
    script_runner,
):
    """Contents paths from config should resolve relative to lite_dir, not CWD.

    When running ``jupyter lite build --lite-dir <dir>`` from a different working directory,
    ``"contents": ["."]`` in the config file should resolve to ``lite_dir``, not to CWD.
    """
    # Create a content file in the lite dir
    (an_empty_lite_dir / "notebook.ipynb").write_text("{}")
    config = {
        "LiteBuildConfig": {
            "ignore_sys_prefix": True,
            "contents": ["."],
        },
    }
    (an_empty_lite_dir / "jupyter_lite_config.json").write_text(json.dumps(config))

    # Run from a *different* directory, pointing --lite-dir at an_empty_lite_dir
    other_dir = an_empty_lite_dir.parent
    result = script_runner.run(
        ["jupyter", "lite", "build", "--lite-dir", str(an_empty_lite_dir)],
        cwd=str(other_dir),
    )
    assert result.success

    out = an_empty_lite_dir / "_output"
    root_contents_json = out / "api/contents/all.json"
    root_contents = json.loads(root_contents_json.read_text(encoding="utf-8"))
    assert len(root_contents["content"]) == 1, root_contents
    assert root_contents["content"][0]["name"] == "notebook.ipynb"


def test_workspaces_resolved_relative_to_lite_dir(
    an_empty_lite_dir,
    script_runner,
):
    """Workspace paths from config should resolve relative to lite_dir, not CWD.

    When running ``jupyter lite build --lite-dir <dir>`` from a different working directory,
    the ``"workspaces"`` config option should resolve to ``lite_dir``, not to CWD.
    """
    workspace = an_empty_lite_dir / "default.jupyterlab-workspace"
    workspace.write_text(
        json.dumps({"data": {}, "metadata": {"id": "default"}}),
        encoding="utf-8",
    )
    config = {
        "LiteBuildConfig": {
            "ignore_sys_prefix": True,
            "workspaces": ["default.jupyterlab-workspace"],
        },
    }
    (an_empty_lite_dir / "jupyter_lite_config.json").write_text(json.dumps(config))

    other_dir = an_empty_lite_dir.parent
    result = script_runner.run(
        ["jupyter", "lite", "build", "--lite-dir", str(an_empty_lite_dir)],
        cwd=str(other_dir),
    )
    assert result.success

    out = an_empty_lite_dir / "_output"
    workspaces_json = out / "api/workspaces/all.json"
    workspaces = json.loads(workspaces_json.read_text(encoding="utf-8"))
    assert "default" in workspaces, workspaces


def test_workspaces_path_in_py_config_resolved_relative_to_lite_dir(
    an_empty_lite_dir,
    script_runner,
):
    """Path-based workspaces in py config should resolve relative to lite_dir."""
    workspace = an_empty_lite_dir / "default.jupyterlab-workspace"
    workspace.write_text(
        json.dumps({"data": {}, "metadata": {"id": "default"}}),
        encoding="utf-8",
    )
    (an_empty_lite_dir / "jupyter_lite_config.py").write_text(
        "\n".join(
            [
                "from pathlib import Path",
                "c = get_config()",
                "c.LiteBuildConfig.ignore_sys_prefix = True",
                "c.LiteBuildConfig.workspaces = [Path('default.jupyterlab-workspace')]",
            ]
        ),
        encoding="utf-8",
    )

    other_dir = an_empty_lite_dir.parent
    result = script_runner.run(
        ["jupyter", "lite", "build", "--lite-dir", str(an_empty_lite_dir)],
        cwd=str(other_dir),
    )
    assert result.success

    out = an_empty_lite_dir / "_output"
    workspaces_json = out / "api/workspaces/all.json"
    workspaces = json.loads(workspaces_json.read_text(encoding="utf-8"))
    assert "default" in workspaces, workspaces


def test_contents_path_in_py_config_resolved_relative_to_lite_dir(
    an_empty_lite_dir,
    script_runner,
):
    """Path-based contents in py config should resolve relative to lite_dir."""
    # Use a subdirectory so the .py config file itself is not picked up as content
    content_dir = an_empty_lite_dir / "my_contents"
    content_dir.mkdir()
    (content_dir / "notebook.ipynb").write_text("{}")
    (an_empty_lite_dir / "jupyter_lite_config.py").write_text(
        "\n".join(
            [
                "from pathlib import Path",
                "c = get_config()",
                "c.LiteBuildConfig.ignore_sys_prefix = True",
                "c.LiteBuildConfig.contents = [Path('my_contents')]",
            ]
        ),
        encoding="utf-8",
    )

    other_dir = an_empty_lite_dir.parent
    result = script_runner.run(
        ["jupyter", "lite", "build", "--lite-dir", str(an_empty_lite_dir)],
        cwd=str(other_dir),
    )
    assert result.success

    out = an_empty_lite_dir / "_output"
    root_contents_json = out / "api/contents/all.json"
    root_contents = json.loads(root_contents_json.read_text(encoding="utf-8"))
    assert len(root_contents["content"]) == 1, root_contents
    assert root_contents["content"][0]["name"] == "notebook.ipynb"
