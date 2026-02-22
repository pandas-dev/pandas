from __future__ import annotations

import importlib
import json
import pathlib
import platform
import sys

import click
import pytest
import yaml
from click.testing import CliRunner

import dask
import dask.cli
from dask._compatibility import importlib_metadata


def test_config_get_no_key():
    runner = CliRunner()
    result = runner.invoke(dask.cli.config_get)
    assert result.exit_code == 2
    expected = (
        "Usage: get [OPTIONS] KEY\n"
        "Try 'get --help' for help.\n\n"
        "Error: Missing argument 'KEY'.\n"
    )
    assert result.output == expected


def test_config_get_value():
    runner = CliRunner()
    result = runner.invoke(dask.cli.config_get, ["array"])
    assert result.exit_code == 0
    assert result.output.startswith("backend:")
    assert len(result.output.splitlines()) > 2


def test_config_get_bad_value():
    runner = CliRunner()
    result = runner.invoke(dask.cli.config_get, ["bad_key"])
    assert result.exit_code != 0
    assert result.output.startswith("Section not found")


def test_config_get_none():
    with dask.config.set({"foo.bar": None}):
        runner = CliRunner()
        result = runner.invoke(dask.cli.config_get, ["foo.bar"])
        assert result.exit_code == 0
        assert result.output == "None\n"


@pytest.fixture
def tmp_conf_dir(tmpdir, monkeypatch):
    # Set temporary DASK_CONFIG in dask.config and handle reloading of module
    # before and after test so module initialization takes place setting attrs
    # like PATH, config, paths, and module level constants
    monkeypatch.setenv("DASK_CONFIG", str(tmpdir))
    originals = dask.config.__dict__.copy()
    dask.config = importlib.reload(dask.config)
    dask.config.paths = [str(tmpdir)]
    try:
        yield pathlib.Path(tmpdir)
    finally:
        dask.config = importlib.reload(dask.config)
        dask.config.__dict__.update(originals)


@pytest.mark.parametrize("value", ("333MiB", 2, [1, 2], {"foo": "bar"}, None))
@pytest.mark.parametrize("empty_config", (True, False))
@pytest.mark.parametrize("file", (None, "bar.yaml", "foo/bar.yaml"))
@pytest.mark.parametrize("existing_key_config_file", (True, False))
def test_config_set_value(
    tmp_conf_dir,
    value,
    empty_config,
    file,
    existing_key_config_file,
):
    # Potentially custom file location, CLI defaulting to the default dir / 'dask.yaml'
    config_file = tmp_conf_dir / (file or "dask.yaml")

    if not empty_config:
        expected_conf = {"dataframe": {"foo": "bar"}}
        config_file.parent.mkdir(parents=True, exist_ok=True)
        config_file.write_text(yaml.dump(expected_conf))
    else:
        expected_conf = dict()
        assert not config_file.exists()

    cmd = ["fizz.buzz", str(value)]

    # default dask config location when file not specified or first file containing key
    if file:
        cmd.extend(["--file", str(config_file)])

    runner = CliRunner()

    if existing_key_config_file and not file:
        # Now we want to ensure the config file being used is this one b/c it already
        # has this key and hasn't explicitly provided a different file to use.
        config_file = pathlib.Path(tmp_conf_dir) / "existing_conf.yaml"
        cmd_ = ["fizz.buzz", "foobar", "--file", str(config_file)]
        runner.invoke(dask.cli.config_set, cmd_, catch_exceptions=False)

        # Expected config is now just what's in this file, not the default file
        expected_conf = {"fizz": {"buzz": "foobar"}}

    result = runner.invoke(dask.cli.config_set, cmd, catch_exceptions=False)

    expected = f"Updated [fizz.buzz] to [{value}], config saved to {config_file}\n"
    assert expected == result.output

    actual_conf = yaml.safe_load(config_file.read_text())

    expected_conf.update({"fizz": {"buzz": value}})
    assert expected_conf == actual_conf


def test_config_find(tmp_conf_dir):
    runner = CliRunner()

    # no config files
    result = runner.invoke(dask.cli.config_find, ["fizz.buzz"], catch_exceptions=False)
    expected = (
        f"Unable to find [fizz.buzz] in any of the following paths:\n{tmp_conf_dir}\n"
    )
    assert result.output == expected

    conf1 = tmp_conf_dir / "conf1.yaml"
    conf2 = tmp_conf_dir / "conf2.yaml"
    conf3 = tmp_conf_dir / "conf3.yaml"

    conf1.write_text(yaml.dump({"fizz": {"buzz": 1}}))
    conf2.write_text(yaml.dump({"fizz": {"buzz": 2}}))
    conf3.write_text(yaml.dump({"foo": {"bar": 1}}))  # shouldn't show up

    result = runner.invoke(dask.cli.config_find, ["fizz.buzz"], catch_exceptions=False)
    expected = (
        "Found [fizz.buzz] in the following files:\n"
        f"{conf1}  [fizz.buzz=1]\n"
        f"{conf2}  [fizz.buzz=2]\n"
    )
    assert result.output == expected


def test_config_list():
    runner = CliRunner()
    result = runner.invoke(dask.cli.config_list)
    assert result.exit_code == 0
    assert "array:" in result.output


def test_version():
    runner = CliRunner()
    result = runner.invoke(dask.cli.cli, ["--version"])
    assert result.exit_code == 0
    assert result.output == f"cli, version {dask.__version__}\n"


def test_info_versions():
    runner = CliRunner()
    result = runner.invoke(dask.cli.versions)
    assert result.exit_code == 0

    # $ dask info versions
    # will print to stdout a json like struct, so result.output can be
    # loaded with json.
    table = json.loads(result.output)

    assert table["Python"] == ".".join(str(x) for x in sys.version_info[:3])
    assert table["dask"] == dask.__version__
    assert table["Platform"] == platform.uname().system

    try:
        from distributed import __version__ as distributed_version
    except ImportError:
        distributed_version = None

    assert table["distributed"] == distributed_version


@click.group()
def dummy_cli():
    pass


def bad_command():
    pass


@click.command(name="good")
def good_command():
    pass


@click.command(name="good")
def good_command_2():
    pass


def test_register_command_ep():
    from dask.cli import _register_command_ep

    bad_ep = importlib_metadata.EntryPoint(
        name="bad",
        value="dask.tests.test_cli:bad_command",
        group="dask_cli",
    )

    good_ep = importlib_metadata.EntryPoint(
        name="good",
        value="dask.tests.test_cli:good_command",
        group="dask_cli",
    )

    class ErrorEP:
        @property
        def name(self):
            return "foo"

        def load(self):
            raise ImportError("Entrypoint could not be imported")

    with pytest.warns(UserWarning, match="must be instances of"):
        _register_command_ep(dummy_cli, bad_ep)

    with pytest.warns(UserWarning, match="exception occurred"):
        _register_command_ep(dummy_cli, ErrorEP())

    _register_command_ep(dummy_cli, good_ep)
    assert "good" in dummy_cli.commands
    assert dummy_cli.commands["good"] is good_command


@click.group()
def dummy_cli_2():
    pass


def test_repeated_name_registration_warn():
    from dask.cli import _register_command_ep

    one = importlib_metadata.EntryPoint(
        name="one",
        value="dask.tests.test_cli:good_command",
        group="dask_cli",
    )

    two = importlib_metadata.EntryPoint(
        name="two",
        value="dask.tests.test_cli:good_command_2",
        group="dask_cli",
    )

    _register_command_ep(dummy_cli_2, one)
    with pytest.warns(UserWarning, match="While registering the command with name"):
        _register_command_ep(dummy_cli_2, two)
