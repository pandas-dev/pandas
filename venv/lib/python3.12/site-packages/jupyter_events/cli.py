"""The cli for jupyter events."""
from __future__ import annotations

import json
import pathlib
import platform

import click
from jsonschema import ValidationError
from rich.console import Console
from rich.json import JSON
from rich.markup import escape
from rich.padding import Padding
from rich.style import Style

from jupyter_events.schema import EventSchema, EventSchemaFileAbsent, EventSchemaLoadingError

WIN = platform.system() == "Windows"


class RC:
    """Return code enum."""

    OK = 0
    INVALID = 1
    UNPARSABLE = 2
    NOT_FOUND = 3


class EMOJI:
    """Terminal emoji enum"""

    X = "XX" if WIN else "\u274c"
    OK = "OK" if WIN else "\u2714"


console = Console()
error_console = Console(stderr=True)


@click.group()
@click.version_option()
def main() -> None:
    """A simple CLI tool to quickly validate JSON schemas against
    Jupyter Event's custom validator.

    You can see Jupyter Event's meta-schema here:

        https://raw.githubusercontent.com/jupyter/jupyter_events/main/jupyter_events/schemas/event-metaschema.yml
    """


@click.command()
@click.argument("schema")
@click.pass_context
def validate(ctx: click.Context, schema: str) -> int:
    """Validate a SCHEMA against Jupyter Event's meta schema.

    SCHEMA can be a JSON/YAML string or filepath to a schema.
    """
    console.rule("Validating the following schema", style=Style(color="blue"))

    _schema = None
    try:
        # attempt to read schema as a serialized string
        _schema = EventSchema._load_schema(schema)
    except EventSchemaLoadingError:
        # pass here to avoid printing traceback of this exception if next block
        # excepts
        pass

    # if not a serialized schema string, try to interpret it as a path to schema file
    if _schema is None:
        schema_path = pathlib.Path(schema)
        try:
            _schema = EventSchema._load_schema(schema_path)
        except (EventSchemaLoadingError, EventSchemaFileAbsent) as e:
            # no need for full tracestack for user error exceptions. just print
            # the error message and return
            error_console.print(f"[bold red]ERROR[/]: {e}")
            return ctx.exit(RC.UNPARSABLE)

    # Print what was found.
    schema_json = JSON(json.dumps(_schema))
    console.print(Padding(schema_json, (1, 0, 1, 4)))
    # Now validate this schema against the meta-schema.
    try:
        EventSchema(_schema)
        console.rule("Results", style=Style(color="green"))
        out = Padding(f"[green]{EMOJI.OK}[white] Nice work! This schema is valid.", (1, 0, 1, 0))
        console.print(out)
        return ctx.exit(RC.OK)
    except ValidationError as err:
        error_console.rule("Results", style=Style(color="red"))
        error_console.print(f"[red]{EMOJI.X} [white]The schema failed to validate.")
        error_console.print("\nWe found the following error with your schema:")
        out = escape(str(err))  # type:ignore[assignment]
        error_console.print(Padding(out, (1, 0, 1, 4)))
        return ctx.exit(RC.INVALID)


main.add_command(validate)
