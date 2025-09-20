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
import json
from abc import ABC, abstractmethod
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
)
from uuid import UUID

from rich.console import Console
from rich.table import Table as RichTable
from rich.tree import Tree

from pyiceberg.partitioning import PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.table import Table
from pyiceberg.table.metadata import TableMetadata
from pyiceberg.table.refs import SnapshotRefType
from pyiceberg.typedef import IcebergBaseModel, Identifier, Properties


class Output(ABC):
    """Output interface for exporting."""

    @abstractmethod
    def exception(self, ex: Exception) -> None: ...

    @abstractmethod
    def identifiers(self, identifiers: List[Identifier]) -> None: ...

    @abstractmethod
    def describe_table(self, table: Table) -> None: ...

    @abstractmethod
    def files(self, table: Table, history: bool) -> None: ...

    @abstractmethod
    def describe_properties(self, properties: Properties) -> None: ...

    @abstractmethod
    def text(self, response: str) -> None: ...

    @abstractmethod
    def schema(self, schema: Schema) -> None: ...

    @abstractmethod
    def spec(self, spec: PartitionSpec) -> None: ...

    @abstractmethod
    def uuid(self, uuid: Optional[UUID]) -> None: ...

    @abstractmethod
    def version(self, version: str) -> None: ...

    @abstractmethod
    def describe_refs(self, refs: List[Tuple[str, SnapshotRefType, Dict[str, str]]]) -> None: ...


class ConsoleOutput(Output):
    """Writes to the console."""

    verbose: bool

    def __init__(self, **properties: Any) -> None:
        self.verbose = properties.get("verbose", False)

    @property
    def _table(self) -> RichTable:
        return RichTable.grid(padding=(0, 2))

    def exception(self, ex: Exception) -> None:
        if self.verbose:
            Console(stderr=True).print_exception()
        else:
            Console(stderr=True).print(ex)

    def identifiers(self, identifiers: List[Identifier]) -> None:
        table = self._table
        for identifier in identifiers:
            table.add_row(".".join(identifier))

        Console().print(table)

    def describe_table(self, table: Table) -> None:
        metadata = table.metadata
        table_properties = self._table

        for key, value in metadata.properties.items():
            table_properties.add_row(key, value)

        schema_tree = Tree(f"Schema, id={table.metadata.current_schema_id}")
        for field in table.schema().fields:
            schema_tree.add(str(field))

        snapshot_tree = Tree("Snapshots")
        for snapshot in metadata.snapshots:
            snapshot_tree.add(f"Snapshot {snapshot.snapshot_id}, schema {snapshot.schema_id}: {snapshot.manifest_list}")

        output_table = self._table
        output_table.add_row("Table format version", str(metadata.format_version))
        output_table.add_row("Metadata location", table.metadata_location)
        output_table.add_row("Table UUID", str(table.metadata.table_uuid))
        output_table.add_row("Last Updated", str(metadata.last_updated_ms))
        output_table.add_row("Partition spec", str(table.spec()))
        output_table.add_row("Sort order", str(table.sort_order()))
        output_table.add_row("Current schema", schema_tree)
        output_table.add_row("Current snapshot", str(table.current_snapshot()))
        output_table.add_row("Snapshots", snapshot_tree)
        output_table.add_row("Properties", table_properties)
        Console().print(output_table)

    def files(self, table: Table, history: bool) -> None:
        if history:
            snapshots = table.metadata.snapshots
        else:
            if snapshot := table.current_snapshot():
                snapshots = [snapshot]
            else:
                snapshots = []

        snapshot_tree = Tree(f"Snapshots: {'.'.join(table.name())}")
        io = table.io

        for snapshot in snapshots:
            list_tree = snapshot_tree.add(
                f"Snapshot {snapshot.snapshot_id}, schema {snapshot.schema_id}: {snapshot.manifest_list}"
            )

            manifest_list = snapshot.manifests(io)
            for manifest in manifest_list:
                manifest_tree = list_tree.add(f"Manifest: {manifest.manifest_path}")
                for manifest_entry in manifest.fetch_manifest_entry(io, discard_deleted=False):
                    manifest_tree.add(f"Datafile: {manifest_entry.data_file.file_path}")
        Console().print(snapshot_tree)

    def describe_properties(self, properties: Properties) -> None:
        output_table = self._table
        for k, v in properties.items():
            output_table.add_row(k, v)
        Console().print(output_table)

    def text(self, response: str) -> None:
        Console(soft_wrap=True).print(response)

    def schema(self, schema: Schema) -> None:
        output_table = self._table
        for field in schema.fields:
            output_table.add_row(field.name, str(field.field_type), field.doc or "")
        Console().print(output_table)

    def spec(self, spec: PartitionSpec) -> None:
        Console().print(str(spec))

    def uuid(self, uuid: Optional[UUID]) -> None:
        Console().print(str(uuid) if uuid else "missing")

    def version(self, version: str) -> None:
        Console().print(version)

    def describe_refs(self, ref_details: List[Tuple[str, SnapshotRefType, Dict[str, str]]]) -> None:
        refs_table = RichTable(title="Snapshot Refs")
        refs_table.add_column("Ref")
        refs_table.add_column("Type")
        refs_table.add_column("Max ref age ms")
        refs_table.add_column("Min snapshots to keep")
        refs_table.add_column("Max snapshot age ms")
        for name, type, ref_detail in ref_details:
            refs_table.add_row(
                name, type, ref_detail["max_ref_age_ms"], ref_detail["min_snapshots_to_keep"], ref_detail["max_snapshot_age_ms"]
            )
        Console().print(refs_table)


class JsonOutput(Output):
    """Writes json to stdout."""

    verbose: bool

    def __init__(self, **properties: Any) -> None:
        self.verbose = properties.get("verbose", False)

    def _out(self, d: Any) -> None:
        print(json.dumps(d))

    def exception(self, ex: Exception) -> None:
        self._out({"type": ex.__class__.__name__, "message": str(ex)})

    def identifiers(self, identifiers: List[Identifier]) -> None:
        self._out([".".join(identifier) for identifier in identifiers])

    def describe_table(self, table: Table) -> None:
        class FauxTable(IcebergBaseModel):
            """Just to encode it using Pydantic."""

            identifier: Identifier
            metadata_location: str
            metadata: TableMetadata

        print(
            FauxTable(
                identifier=table.name(), metadata=table.metadata, metadata_location=table.metadata_location
            ).model_dump_json()
        )

    def describe_properties(self, properties: Properties) -> None:
        self._out(properties)

    def text(self, response: str) -> None:
        print(json.dumps(response))

    def schema(self, schema: Schema) -> None:
        print(schema.model_dump_json())

    def files(self, table: Table, history: bool) -> None:
        pass

    def spec(self, spec: PartitionSpec) -> None:
        print(spec.model_dump_json())

    def uuid(self, uuid: Optional[UUID]) -> None:
        self._out({"uuid": str(uuid) if uuid else "missing"})

    def version(self, version: str) -> None:
        self._out({"version": version})

    def describe_refs(self, refs: List[Tuple[str, SnapshotRefType, Dict[str, str]]]) -> None:
        self._out(
            [
                {"name": name, "type": type, detail_key: detail_val}
                for name, type, detail in refs
                for detail_key, detail_val in detail.items()
            ]
        )
