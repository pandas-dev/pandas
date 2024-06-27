from __future__ import annotations

import pathlib
import warnings
from functools import reduce
from typing import Any

import click
import yaml

import dask
from dask import __version__
from dask._compatibility import importlib_metadata

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
    "max_content_width": 88,
}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
def cli():
    """Dask command line interface."""
    pass


@cli.command()
def docs():
    """Open Dask documentation (https://docs.dask.org/) in a web browser."""
    import webbrowser

    webbrowser.open("https://docs.dask.org")


@cli.group()
def info():
    """Information about your dask installation."""
    pass


@info.command()
def versions():
    """Print versions of Dask related projects."""
    from dask.utils import show_versions

    show_versions()


@cli.group()
def config():
    """Dask config settings"""
    pass


@config.command(name="get")
@click.argument("key", required=True)
def config_get(key=None):
    """Print config key, or the whole config."""
    try:
        data = reduce(lambda d, k: d[k], key.split("."), dask.config.config)
    except (KeyError, TypeError):
        click.echo(click.style(f"Section not found: {key}", fg="red"), err=True)
        exit(1)

    if isinstance(data, (list, dict)):
        click.echo_via_pager(yaml.dump(data))
    elif data is None:
        click.echo("None")
    else:
        click.echo(data)


@config.command(name="find")
@click.argument("key", required=True)
def config_find(key):
    """Find a Dask config key by searching config location paths"""
    paths = list(dask.config.paths_containing_key(key))
    if paths:
        click.echo(f"Found [{key}] in the following files:")
        max_len = len(str(max(paths, key=lambda v: len(str(v)))))
        for path in paths:
            config = next(dask.config.collect_yaml([path]))
            value = dask.config.get(key, config=config)
            spacing = " " * (max_len - len(str(path)))
            click.echo(f"{path} {spacing} [{key}={value}]")
    else:
        click.echo(f"Unable to find [{key}] in any of the following paths:")
        click.echo("\n".join(map(str, dask.config.paths)))


@config.command(name="set")
@click.argument("key", required=True)
@click.argument("value", required=True)
@click.option(
    "--file",
    default=None,
    required=False,
    type=pathlib.Path,
    help=(
        "File to use, defaulting to the first config file containing the key, "
        "'DASK_CONFIG' environment variable location or finally, "
        "'~/.config/dask/dask.yaml'"
    ),
)
def config_set(key, value, file):
    """Set a Dask config key to a new value"""
    value = dask.config.interpret_value(value)
    _, path = save_config(key, value, file)
    click.echo(f"Updated [{key}] to [{value}], config saved to {path}")


def save_config(
    key: str,
    value: Any,
    config_file: pathlib.Path | None = None,
) -> tuple[dict, pathlib.Path]:
    """
    Save new config values to dask config file,
    return new config and path to file.
    """
    # If explicit config file wasn't supplied, try finding the first existing
    # config file w/ this key, falling back to the default config path
    if config_file is None:
        config_file = next(
            dask.config.paths_containing_key(key),
            pathlib.Path(dask.config.PATH) / "dask.yaml",
        )
    config_file.parent.mkdir(exist_ok=True, parents=True)

    if config_file.exists():
        config = yaml.safe_load(config_file.read_text()) or {}
    else:
        config = {}

    dask.config.set({key: value}, config=config)
    try:
        config_file.write_text(yaml.dump(config))
    except OSError as e:
        raise RuntimeError(
            f"""

For some reason we couldn't write config to {config_file}.
Perhaps you don't have permissions here?

You can change the directory where Dask writes config files
to somewhere else using the DASK_CONFIG environment variable.

For example, you could set the following:

    export DASK_CONFIG=~/.dask
"""
        ) from e

    return config, config_file


@config.command(name="list")
def config_list():
    """Print the whole config."""
    click.echo_via_pager(yaml.dump(dask.config.config))


def _register_command_ep(interface, entry_point):
    """Add `entry_point` command to `interface`.

    Parameters
    ----------
    interface : click.Command or click.Group
        The click interface to augment with `entry_point`.
    entry_point : importlib.metadata.EntryPoint
        The entry point which loads to a ``click.Command`` or
        ``click.Group`` instance to be added as a sub-command or
        sub-group in `interface`.

    """
    try:
        command = entry_point.load()
    except Exception as e:
        warnings.warn(
            f"While registering the command with name '{entry_point.name}', an "
            f"exception occurred; {e}."
        )
        return
    if not isinstance(command, (click.Command, click.Group)):
        warnings.warn(
            "entry points in 'dask_cli' must be instances of "
            f"click.Command or click.Group, not {type(command)}."
        )
        return

    if command.name in interface.commands:
        warnings.warn(
            f"While registering the command with name '{command.name}', an "
            "existing command or group; the original has been overwritten."
        )

    interface.add_command(command)


def run_cli():
    """Run the dask command line interface."""

    # discover "dask_cli" entry points and try to register them to the
    # top level `cli`.
    for ep in importlib_metadata.entry_points(group="dask_cli"):
        _register_command_ep(cli, ep)

    cli()
