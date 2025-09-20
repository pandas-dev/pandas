# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
# pylint: disable=broad-except,redefined-builtin,redefined-outer-name
from functools import wraps
from typing import (
    Any,
    Callable,
    Dict,
    Literal,
    Optional,
    Tuple,
)

import click
from click import Context

from pyiceberg import __version__
from pyiceberg.catalog import URI, Catalog, load_catalog
from pyiceberg.cli.output import ConsoleOutput, JsonOutput, Output
from pyiceberg.exceptions import NoSuchNamespaceError, NoSuchPropertyException, NoSuchTableError
from pyiceberg.table import TableProperties
from pyiceberg.table.refs import SnapshotRef, SnapshotRefType
from pyiceberg.utils.properties import property_as_int


def catch_exception() -> Callable:  # type: ignore
    def decorator(func: Callable) -> Callable:  # type: ignore
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any):  # type: ignore
            try:
                return func(*args, **kwargs)
            except Exception as e:
                ctx: Context = click.get_current_context(silent=True)
                _, output = _catalog_and_output(ctx)
                output.exception(e)
                ctx.exit(1)

        return wrapper

    return decorator


@click.group()
@click.option("--catalog")
@click.option("--verbose", type=click.BOOL)
@click.option("--output", type=click.Choice(["text", "json"]), default="text")
@click.option("--ugi")
@click.option("--uri")
@click.option("--credential")
@click.pass_context
def run(
    ctx: Context,
    catalog: Optional[str],
    verbose: bool,
    output: str,
    ugi: Optional[str],
    uri: Optional[str],
    credential: Optional[str],
) -> None:
    properties = {}
    if ugi:
        properties["ugi"] = ugi
    if uri:
        properties[URI] = uri
    if credential:
        properties["credential"] = credential

    ctx.ensure_object(dict)
    if output == "text":
        ctx.obj["output"] = ConsoleOutput(verbose=verbose)
    else:
        ctx.obj["output"] = JsonOutput(verbose=verbose)

    try:
        ctx.obj["catalog"] = load_catalog(catalog, **properties)
    except Exception as e:
        ctx.obj["output"].exception(e)
        ctx.exit(1)

    if not isinstance(ctx.obj["catalog"], Catalog):
        ctx.obj["output"].exception(
            ValueError("Could not determine catalog type from uri. REST (http/https) and Hive (thrift) is supported")
        )
        ctx.exit(1)


def _catalog_and_output(ctx: Context) -> Tuple[Catalog, Output]:
    """Small helper to set the types."""
    return ctx.obj["catalog"], ctx.obj["output"]


@run.command()
@click.pass_context
@click.argument("parent", required=False)
@catch_exception()
def list(ctx: Context, parent: Optional[str]) -> None:  # pylint: disable=redefined-builtin
    """List tables or namespaces."""
    catalog, output = _catalog_and_output(ctx)

    identifiers = []
    if parent:
        # Do we have tables under parent namespace?
        identifiers = catalog.list_tables(parent)
    if not identifiers:
        # List hierarchical namespaces if parent, root namespaces otherwise.
        identifiers = catalog.list_namespaces(parent or ())
    output.identifiers(identifiers)


@run.command()
@click.option("--entity", type=click.Choice(["any", "namespace", "table"]), default="any")
@click.argument("identifier")
@click.pass_context
@catch_exception()
def describe(ctx: Context, entity: Literal["name", "namespace", "table"], identifier: str) -> None:
    """Describe a namespace or a table."""
    catalog, output = _catalog_and_output(ctx)
    identifier_tuple = Catalog.identifier_to_tuple(identifier)

    is_namespace = False
    if entity in {"namespace", "any"} and len(identifier_tuple) > 0:
        try:
            namespace_properties = catalog.load_namespace_properties(identifier_tuple)
            output.describe_properties(namespace_properties)
            is_namespace = True
        except NoSuchNamespaceError as exc:
            if entity != "any" or len(identifier_tuple) == 1:  # type: ignore
                raise exc

    is_table = False
    if entity in {"table", "any"} and len(identifier_tuple) > 1:
        try:
            catalog_table = catalog.load_table(identifier)
            output.describe_table(catalog_table)
            is_table = True
        except NoSuchTableError as exc:
            if entity != "any":
                raise exc

    if is_namespace is False and is_table is False:
        raise NoSuchTableError(f"Table or namespace does not exist: {identifier}")


@run.command()
@click.argument("identifier")
@click.option("--history", is_flag=True)
@click.pass_context
@catch_exception()
def files(ctx: Context, identifier: str, history: bool) -> None:
    """List all the files of the table."""
    catalog, output = _catalog_and_output(ctx)

    catalog_table = catalog.load_table(identifier)
    output.files(catalog_table, history)


@run.command()
@click.argument("identifier")
@click.pass_context
@catch_exception()
def schema(ctx: Context, identifier: str) -> None:
    """Get the schema of the table."""
    catalog, output = _catalog_and_output(ctx)
    table = catalog.load_table(identifier)
    output.schema(table.schema())


@run.command()
@click.argument("identifier")
@click.pass_context
@catch_exception()
def spec(ctx: Context, identifier: str) -> None:
    """Return the partition spec of the table."""
    catalog, output = _catalog_and_output(ctx)
    table = catalog.load_table(identifier)
    output.spec(table.spec())


@run.command()
@click.argument("identifier")
@click.pass_context
@catch_exception()
def uuid(ctx: Context, identifier: str) -> None:
    """Return the UUID of the table."""
    catalog, output = _catalog_and_output(ctx)
    metadata = catalog.load_table(identifier).metadata
    output.uuid(metadata.table_uuid)


@run.command()
@click.argument("identifier")
@click.pass_context
@catch_exception()
def location(ctx: Context, identifier: str) -> None:
    """Return the location of the table."""
    catalog, output = _catalog_and_output(ctx)
    table = catalog.load_table(identifier)
    output.text(table.location())


@run.command()
@click.pass_context
@catch_exception()
def version(ctx: Context) -> None:
    """Print pyiceberg version."""
    ctx.obj["output"].version(__version__)


@run.group()
def create() -> None:
    """Operation to create a namespace."""


@create.command()
@click.argument("identifier")
@click.pass_context
@catch_exception()
def namespace(ctx: Context, identifier: str) -> None:
    """Create a namespace."""
    catalog, output = _catalog_and_output(ctx)

    catalog.create_namespace(identifier)
    output.text(f"Created namespace: {identifier}")


@run.group()
def drop() -> None:
    """Operations to drop a namespace or table."""


@drop.command()
@click.argument("identifier")
@click.pass_context
@catch_exception()
def table(ctx: Context, identifier: str) -> None:  # noqa: F811
    """Drop a table."""
    catalog, output = _catalog_and_output(ctx)

    catalog.drop_table(identifier)
    output.text(f"Dropped table: {identifier}")


@drop.command()  # type: ignore
@click.argument("identifier")
@click.pass_context
@catch_exception()
def namespace(ctx: Context, identifier: str) -> None:  # noqa: F811
    """Drop a namespace."""
    catalog, output = _catalog_and_output(ctx)

    catalog.drop_namespace(identifier)
    output.text(f"Dropped namespace: {identifier}")


@run.command()
@click.argument("from_identifier")
@click.argument("to_identifier")
@click.pass_context
@catch_exception()
def rename(ctx: Context, from_identifier: str, to_identifier: str) -> None:
    """Rename a table."""
    catalog, output = _catalog_and_output(ctx)

    catalog.rename_table(from_identifier, to_identifier)
    output.text(f"Renamed table from {from_identifier} to {to_identifier}")


@run.group()
def properties() -> None:
    """Properties on tables/namespaces."""


@properties.group()
def get() -> None:
    """Fetch properties on tables/namespaces."""


@get.command("namespace")
@click.argument("identifier")
@click.argument("property_name", required=False)
@click.pass_context
@catch_exception()
def get_namespace(ctx: Context, identifier: str, property_name: str) -> None:
    """Fetch properties on a namespace."""
    catalog, output = _catalog_and_output(ctx)
    identifier_tuple = Catalog.identifier_to_tuple(identifier)

    namespace_properties = catalog.load_namespace_properties(identifier_tuple)

    if property_name:
        if property_value := namespace_properties.get(property_name):
            output.text(property_value)
        else:
            raise NoSuchPropertyException(f"Could not find property {property_name} on namespace {identifier}")
    else:
        output.describe_properties(namespace_properties)


@get.command("table")
@click.argument("identifier")
@click.argument("property_name", required=False)
@click.pass_context
@catch_exception()
def get_table(ctx: Context, identifier: str, property_name: str) -> None:
    """Fetch properties on a table."""
    catalog, output = _catalog_and_output(ctx)
    identifier_tuple = Catalog.identifier_to_tuple(identifier)

    metadata = catalog.load_table(identifier_tuple).metadata

    if property_name:
        if property_value := metadata.properties.get(property_name):
            output.text(property_value)
        else:
            raise NoSuchPropertyException(f"Could not find property {property_name} on table {identifier}")
    else:
        output.describe_properties(metadata.properties)


@properties.group()
def set() -> None:
    """Set a property on tables/namespaces."""


@set.command()  # type: ignore
@click.argument("identifier")
@click.argument("property_name")
@click.argument("property_value")
@click.pass_context
@catch_exception()
def namespace(ctx: Context, identifier: str, property_name: str, property_value: str) -> None:  # noqa: F811
    """Set a property on a namespace."""
    catalog, output = _catalog_and_output(ctx)

    catalog.update_namespace_properties(identifier, updates={property_name: property_value})
    output.text(f"Updated {property_name} on {identifier}")


@set.command()  # type: ignore
@click.argument("identifier")
@click.argument("property_name")
@click.argument("property_value")
@click.pass_context
@catch_exception()
def table(ctx: Context, identifier: str, property_name: str, property_value: str) -> None:  # noqa: F811
    """Set a property on a table."""
    catalog, output = _catalog_and_output(ctx)
    identifier_tuple = Catalog.identifier_to_tuple(identifier)

    table = catalog.load_table(identifier_tuple)
    with table.transaction() as tx:
        tx.set_properties({property_name: property_value})
    output.text(f"Set {property_name}={property_value} on {identifier}")


@properties.group()
def remove() -> None:
    """Remove a property from tables/namespaces."""


@remove.command()  # type: ignore
@click.argument("identifier")
@click.argument("property_name")
@click.pass_context
@catch_exception()
def namespace(ctx: Context, identifier: str, property_name: str) -> None:  # noqa: F811
    """Remove a property from a namespace."""
    catalog, output = _catalog_and_output(ctx)

    result = catalog.update_namespace_properties(identifier, removals={property_name})

    if result.removed == [property_name]:
        output.text(f"Property {property_name} removed from {identifier}")
    else:
        raise NoSuchPropertyException(f"Property {property_name} does not exist on {identifier}")


@remove.command()  # type: ignore
@click.argument("identifier")
@click.argument("property_name")
@click.pass_context
@catch_exception()
def table(ctx: Context, identifier: str, property_name: str) -> None:  # noqa: F811
    """Remove a property from a table."""
    catalog, output = _catalog_and_output(ctx)
    table = catalog.load_table(identifier)
    if property_name in table.metadata.properties:
        with table.transaction() as tx:
            tx.remove_properties(property_name)
        output.text(f"Property {property_name} removed from {identifier}")
    else:
        raise NoSuchPropertyException(f"Property {property_name} does not exist on {identifier}")


@run.command()
@click.argument("identifier")
@click.option("--type", required=False)
@click.option("--verbose", type=click.BOOL)
@click.pass_context
@catch_exception()
def list_refs(ctx: Context, identifier: str, type: str, verbose: bool) -> None:
    """List all the refs in the provided table."""
    catalog, output = _catalog_and_output(ctx)
    table = catalog.load_table(identifier)
    refs = table.refs()
    if type:
        type = type.lower()
        if type not in {SnapshotRefType.BRANCH, SnapshotRefType.TAG}:
            raise ValueError(f"Type must be either branch or tag, got: {type}")

    relevant_refs = [
        (ref_name, ref.snapshot_ref_type, _retention_properties(ref, table.properties))
        for (ref_name, ref) in refs.items()
        if not type or ref.snapshot_ref_type == type
    ]

    output.describe_refs(relevant_refs)


def _retention_properties(ref: SnapshotRef, table_properties: Dict[str, str]) -> Dict[str, str]:
    retention_properties = {}
    if ref.snapshot_ref_type == SnapshotRefType.BRANCH:
        default_min_snapshots_to_keep = property_as_int(
            table_properties,
            TableProperties.MIN_SNAPSHOTS_TO_KEEP,
            TableProperties.MIN_SNAPSHOTS_TO_KEEP_DEFAULT,
        )

        default_max_snapshot_age_ms = property_as_int(
            table_properties,
            TableProperties.MAX_SNAPSHOT_AGE_MS,
            TableProperties.MAX_SNAPSHOT_AGE_MS_DEFAULT,
        )

        retention_properties["min_snapshots_to_keep"] = (
            str(ref.min_snapshots_to_keep) if ref.min_snapshots_to_keep else str(default_min_snapshots_to_keep)
        )
        retention_properties["max_snapshot_age_ms"] = (
            str(ref.max_snapshot_age_ms) if ref.max_snapshot_age_ms else str(default_max_snapshot_age_ms)
        )
    else:
        retention_properties["min_snapshots_to_keep"] = "N/A"
        retention_properties["max_snapshot_age_ms"] = "N/A"

    retention_properties["max_ref_age_ms"] = str(ref.max_ref_age_ms) if ref.max_ref_age_ms else "forever"

    return retention_properties
