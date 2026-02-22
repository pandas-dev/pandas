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
from __future__ import annotations

import itertools
import uuid
from abc import ABC, abstractmethod
from datetime import datetime
from functools import singledispatch
from typing import TYPE_CHECKING, Annotated, Any, Generic, Literal, TypeVar, cast

from pydantic import Field, field_validator, model_serializer, model_validator

from pyiceberg.exceptions import CommitFailedException
from pyiceberg.partitioning import PARTITION_FIELD_ID_START, PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.table.metadata import SUPPORTED_TABLE_FORMAT_VERSION, TableMetadata, TableMetadataUtil
from pyiceberg.table.refs import MAIN_BRANCH, SnapshotRef, SnapshotRefType
from pyiceberg.table.snapshots import (
    MetadataLogEntry,
    Snapshot,
    SnapshotLogEntry,
)
from pyiceberg.table.sorting import SortOrder
from pyiceberg.table.statistics import (
    PartitionStatisticsFile,
    StatisticsFile,
    filter_statistics_by_snapshot_id,
)
from pyiceberg.typedef import (
    IcebergBaseModel,
    Properties,
)
from pyiceberg.types import (
    transform_dict_value_to_str,
)
from pyiceberg.utils.datetime import datetime_to_millis
from pyiceberg.utils.properties import property_as_int

if TYPE_CHECKING:
    from pydantic.functional_serializers import ModelWrapSerializerWithoutInfo

    from pyiceberg.table import Transaction

U = TypeVar("U")


class UpdateTableMetadata(ABC, Generic[U]):
    _transaction: Transaction

    def __init__(self, transaction: Transaction) -> None:
        self._transaction = transaction

    @abstractmethod
    def _commit(self) -> UpdatesAndRequirements: ...

    def commit(self) -> None:
        self._transaction._apply(*self._commit())

    def __exit__(self, _: Any, value: Any, traceback: Any) -> None:
        """Close and commit the change."""
        self.commit()

    def __enter__(self) -> U:
        """Update the table."""
        return self  # type: ignore


class AssignUUIDUpdate(IcebergBaseModel):
    action: Literal["assign-uuid"] = Field(default="assign-uuid")
    uuid: uuid.UUID


class UpgradeFormatVersionUpdate(IcebergBaseModel):
    action: Literal["upgrade-format-version"] = Field(default="upgrade-format-version")
    format_version: int = Field(alias="format-version")


class AddSchemaUpdate(IcebergBaseModel):
    action: Literal["add-schema"] = Field(default="add-schema")
    schema_: Schema = Field(alias="schema")


class SetCurrentSchemaUpdate(IcebergBaseModel):
    action: Literal["set-current-schema"] = Field(default="set-current-schema")
    schema_id: int = Field(
        alias="schema-id", description="Schema ID to set as current, or -1 to set last added schema", default=-1
    )


class AddPartitionSpecUpdate(IcebergBaseModel):
    action: Literal["add-spec"] = Field(default="add-spec")
    spec: PartitionSpec


class SetDefaultSpecUpdate(IcebergBaseModel):
    action: Literal["set-default-spec"] = Field(default="set-default-spec")
    spec_id: int = Field(
        alias="spec-id", description="Partition spec ID to set as the default, or -1 to set last added spec", default=-1
    )


class AddSortOrderUpdate(IcebergBaseModel):
    action: Literal["add-sort-order"] = Field(default="add-sort-order")
    sort_order: SortOrder = Field(alias="sort-order")


class SetDefaultSortOrderUpdate(IcebergBaseModel):
    action: Literal["set-default-sort-order"] = Field(default="set-default-sort-order")
    sort_order_id: int = Field(
        alias="sort-order-id", description="Sort order ID to set as the default, or -1 to set last added sort order", default=-1
    )


class AddSnapshotUpdate(IcebergBaseModel):
    action: Literal["add-snapshot"] = Field(default="add-snapshot")
    snapshot: Snapshot


class SetSnapshotRefUpdate(IcebergBaseModel):
    action: Literal["set-snapshot-ref"] = Field(default="set-snapshot-ref")
    ref_name: str = Field(alias="ref-name")
    type: Literal[SnapshotRefType.TAG, SnapshotRefType.BRANCH]
    snapshot_id: int = Field(alias="snapshot-id")
    max_ref_age_ms: Annotated[int | None, Field(alias="max-ref-age-ms", default=None)]
    max_snapshot_age_ms: Annotated[int | None, Field(alias="max-snapshot-age-ms", default=None)]
    min_snapshots_to_keep: Annotated[int | None, Field(alias="min-snapshots-to-keep", default=None)]


class RemoveSnapshotsUpdate(IcebergBaseModel):
    action: Literal["remove-snapshots"] = Field(default="remove-snapshots")
    snapshot_ids: list[int] = Field(alias="snapshot-ids")


class RemoveSnapshotRefUpdate(IcebergBaseModel):
    action: Literal["remove-snapshot-ref"] = Field(default="remove-snapshot-ref")
    ref_name: str = Field(alias="ref-name")


class SetLocationUpdate(IcebergBaseModel):
    action: Literal["set-location"] = Field(default="set-location")
    location: str


class SetPropertiesUpdate(IcebergBaseModel):
    action: Literal["set-properties"] = Field(default="set-properties")
    updates: dict[str, str]

    @field_validator("updates", mode="before")
    def transform_properties_dict_value_to_str(cls, properties: Properties) -> dict[str, str]:
        return transform_dict_value_to_str(properties)


class RemovePropertiesUpdate(IcebergBaseModel):
    action: Literal["remove-properties"] = Field(default="remove-properties")
    removals: list[str]


class SetStatisticsUpdate(IcebergBaseModel):
    action: Literal["set-statistics"] = Field(default="set-statistics")
    statistics: StatisticsFile
    snapshot_id: int | None = Field(
        None,
        alias="snapshot-id",
        description="snapshot-id is **DEPRECATED for REMOVAL** since it contains redundant information. "
        "Use `statistics.snapshot-id` field instead.",
    )

    @model_validator(mode="before")
    def validate_snapshot_id(cls, data: dict[str, Any]) -> dict[str, Any]:
        stats = data["statistics"]
        if isinstance(stats, StatisticsFile):
            snapshot_id = stats.snapshot_id
        elif isinstance(stats, dict):
            snapshot_id = cast(int, stats.get("snapshot-id"))
        else:
            snapshot_id = None

        data["snapshot_id"] = snapshot_id

        return data


class RemoveStatisticsUpdate(IcebergBaseModel):
    action: Literal["remove-statistics"] = Field(default="remove-statistics")
    snapshot_id: int = Field(alias="snapshot-id")


class RemovePartitionSpecsUpdate(IcebergBaseModel):
    action: Literal["remove-partition-specs"] = Field(default="remove-partition-specs")
    spec_ids: list[int] = Field(alias="spec-ids")


class RemoveSchemasUpdate(IcebergBaseModel):
    action: Literal["remove-schemas"] = Field(default="remove-schemas")
    schema_ids: list[int] = Field(alias="schema-ids")


class SetPartitionStatisticsUpdate(IcebergBaseModel):
    action: Literal["set-partition-statistics"] = Field(default="set-partition-statistics")
    partition_statistics: PartitionStatisticsFile


class RemovePartitionStatisticsUpdate(IcebergBaseModel):
    action: Literal["remove-partition-statistics"] = Field(default="remove-partition-statistics")
    snapshot_id: int = Field(alias="snapshot-id")


TableUpdate = Annotated[
    AssignUUIDUpdate
    | UpgradeFormatVersionUpdate
    | AddSchemaUpdate
    | SetCurrentSchemaUpdate
    | AddPartitionSpecUpdate
    | SetDefaultSpecUpdate
    | AddSortOrderUpdate
    | SetDefaultSortOrderUpdate
    | AddSnapshotUpdate
    | SetSnapshotRefUpdate
    | RemoveSnapshotsUpdate
    | RemoveSnapshotRefUpdate
    | SetLocationUpdate
    | SetPropertiesUpdate
    | RemovePropertiesUpdate
    | SetStatisticsUpdate
    | RemoveStatisticsUpdate
    | RemovePartitionSpecsUpdate
    | RemoveSchemasUpdate
    | SetPartitionStatisticsUpdate
    | RemovePartitionStatisticsUpdate,
    Field(discriminator="action"),
]


class _TableMetadataUpdateContext:
    _updates: list[TableUpdate]

    def __init__(self) -> None:
        self._updates = []

    def add_update(self, update: TableUpdate) -> None:
        self._updates.append(update)

    def is_added_snapshot(self, snapshot_id: int) -> bool:
        return any(
            update.snapshot.snapshot_id == snapshot_id for update in self._updates if isinstance(update, AddSnapshotUpdate)
        )

    def is_added_schema(self, schema_id: int) -> bool:
        return any(update.schema_.schema_id == schema_id for update in self._updates if isinstance(update, AddSchemaUpdate))

    def is_added_partition_spec(self, spec_id: int) -> bool:
        return any(update.spec.spec_id == spec_id for update in self._updates if isinstance(update, AddPartitionSpecUpdate))

    def is_added_sort_order(self, sort_order_id: int) -> bool:
        return any(
            update.sort_order.order_id == sort_order_id for update in self._updates if isinstance(update, AddSortOrderUpdate)
        )

    def has_changes(self) -> bool:
        return len(self._updates) > 0


@singledispatch
def _apply_table_update(update: TableUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    """Apply a table update to the table metadata.

    Args:
        update: The update to be applied.
        base_metadata: The base metadata to be updated.
        context: Contains previous updates and other change tracking information in the current transaction.

    Returns:
        The updated metadata.

    """
    raise NotImplementedError(f"Unsupported table update: {update}")


@_apply_table_update.register(AssignUUIDUpdate)
def _(update: AssignUUIDUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    if update.uuid == base_metadata.table_uuid:
        return base_metadata

    context.add_update(update)
    return base_metadata.model_copy(update={"table_uuid": update.uuid})


@_apply_table_update.register(SetLocationUpdate)
def _(update: SetLocationUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    context.add_update(update)
    return base_metadata.model_copy(update={"location": update.location})


@_apply_table_update.register(UpgradeFormatVersionUpdate)
def _(
    update: UpgradeFormatVersionUpdate,
    base_metadata: TableMetadata,
    context: _TableMetadataUpdateContext,
) -> TableMetadata:
    if update.format_version > SUPPORTED_TABLE_FORMAT_VERSION:
        raise ValueError(f"Unsupported table format version: {update.format_version}")
    elif update.format_version < base_metadata.format_version:
        raise ValueError(f"Cannot downgrade v{base_metadata.format_version} table to v{update.format_version}")
    elif update.format_version == base_metadata.format_version:
        return base_metadata

    updated_metadata = base_metadata.model_copy(update={"format_version": update.format_version})

    context.add_update(update)
    return TableMetadataUtil._construct_without_validation(updated_metadata)


@_apply_table_update.register(SetPropertiesUpdate)
def _(update: SetPropertiesUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    if len(update.updates) == 0:
        return base_metadata

    properties = dict(base_metadata.properties)
    properties.update(update.updates)

    context.add_update(update)
    return base_metadata.model_copy(update={"properties": properties})


@_apply_table_update.register(RemovePropertiesUpdate)
def _(update: RemovePropertiesUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    if len(update.removals) == 0:
        return base_metadata

    properties = dict(base_metadata.properties)
    for key in update.removals:
        properties.pop(key)

    context.add_update(update)
    return base_metadata.model_copy(update={"properties": properties})


@_apply_table_update.register(AddSchemaUpdate)
def _(update: AddSchemaUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    metadata_updates: dict[str, Any] = {
        "last_column_id": max(base_metadata.last_column_id, update.schema_.highest_field_id),
        "schemas": base_metadata.schemas + [update.schema_],
    }

    context.add_update(update)
    return base_metadata.model_copy(update=metadata_updates)


@_apply_table_update.register(SetCurrentSchemaUpdate)
def _(update: SetCurrentSchemaUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    new_schema_id = update.schema_id
    if new_schema_id == -1:
        # The last added schema should be in base_metadata.schemas at this point
        new_schema_id = max(schema.schema_id for schema in base_metadata.schemas)
        if not context.is_added_schema(new_schema_id):
            raise ValueError("Cannot set current schema to last added schema when no schema has been added")

    if new_schema_id == base_metadata.current_schema_id:
        return base_metadata

    schema = base_metadata.schema_by_id(new_schema_id)
    if schema is None:
        raise ValueError(f"Schema with id {new_schema_id} does not exist")

    context.add_update(update)
    return base_metadata.model_copy(update={"current_schema_id": new_schema_id})


@_apply_table_update.register(AddPartitionSpecUpdate)
def _(update: AddPartitionSpecUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    for spec in base_metadata.partition_specs:
        # Only raise in case of a discrepancy
        if spec.spec_id == update.spec.spec_id and spec != update.spec:
            raise ValueError(f"Partition spec with id {spec.spec_id} already exists: {spec}")

    metadata_updates: dict[str, Any] = {
        "partition_specs": base_metadata.partition_specs + [update.spec],
        "last_partition_id": max(
            max([field.field_id for field in update.spec.fields], default=0),
            base_metadata.last_partition_id or PARTITION_FIELD_ID_START - 1,
        ),
    }

    context.add_update(update)
    return base_metadata.model_copy(update=metadata_updates)


@_apply_table_update.register(SetDefaultSpecUpdate)
def _(update: SetDefaultSpecUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    new_spec_id = update.spec_id
    if new_spec_id == -1:
        new_spec_id = max(spec.spec_id for spec in base_metadata.partition_specs)
        if not context.is_added_partition_spec(new_spec_id):
            raise ValueError("Cannot set current partition spec to last added one when no partition spec has been added")
    if new_spec_id == base_metadata.default_spec_id:
        return base_metadata
    found_spec_id = False
    for spec in base_metadata.partition_specs:
        found_spec_id = spec.spec_id == new_spec_id
        if found_spec_id:
            break

    if not found_spec_id:
        raise ValueError(f"Failed to find spec with id {new_spec_id}")

    context.add_update(update)
    return base_metadata.model_copy(update={"default_spec_id": new_spec_id})


@_apply_table_update.register(AddSnapshotUpdate)
def _(update: AddSnapshotUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    if len(base_metadata.schemas) == 0:
        raise ValueError("Attempting to add a snapshot before a schema is added")
    elif len(base_metadata.partition_specs) == 0:
        raise ValueError("Attempting to add a snapshot before a partition spec is added")
    elif len(base_metadata.sort_orders) == 0:
        raise ValueError("Attempting to add a snapshot before a sort order is added")
    elif base_metadata.snapshot_by_id(update.snapshot.snapshot_id) is not None:
        raise ValueError(f"Snapshot with id {update.snapshot.snapshot_id} already exists")
    elif (
        base_metadata.format_version == 2
        and update.snapshot.sequence_number is not None
        and update.snapshot.sequence_number <= base_metadata.last_sequence_number
        and update.snapshot.parent_snapshot_id is not None
    ):
        raise ValueError(
            f"Cannot add snapshot with sequence number {update.snapshot.sequence_number} "
            f"older than last sequence number {base_metadata.last_sequence_number}"
        )
    elif base_metadata.format_version >= 3 and update.snapshot.first_row_id is None:
        raise ValueError("Cannot add snapshot without first row id")
    elif (
        base_metadata.format_version >= 3
        and update.snapshot.first_row_id is not None
        and base_metadata.next_row_id is not None
        and update.snapshot.first_row_id < base_metadata.next_row_id
    ):
        raise ValueError(
            f"Cannot add a snapshot with first row id smaller than the table's next-row-id "
            f"{update.snapshot.first_row_id} < {base_metadata.next_row_id}"
        )

    context.add_update(update)
    return base_metadata.model_copy(
        update={
            "last_updated_ms": update.snapshot.timestamp_ms,
            "last_sequence_number": update.snapshot.sequence_number,
            "snapshots": base_metadata.snapshots + [update.snapshot],
            "next_row_id": base_metadata.next_row_id + update.snapshot.added_rows
            if base_metadata.format_version >= 3
            and base_metadata.next_row_id is not None
            and update.snapshot.added_rows is not None
            else None,
        }
    )


@_apply_table_update.register(SetSnapshotRefUpdate)
def _(update: SetSnapshotRefUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    snapshot_ref = SnapshotRef(
        snapshot_id=update.snapshot_id,
        snapshot_ref_type=update.type,
        min_snapshots_to_keep=update.min_snapshots_to_keep,
        max_snapshot_age_ms=update.max_snapshot_age_ms,
        max_ref_age_ms=update.max_ref_age_ms,
    )

    existing_ref = base_metadata.refs.get(update.ref_name)
    if existing_ref is not None and existing_ref == snapshot_ref:
        return base_metadata

    snapshot = base_metadata.snapshot_by_id(snapshot_ref.snapshot_id)
    if snapshot is None:
        raise ValueError(f"Cannot set {update.ref_name} to unknown snapshot {snapshot_ref.snapshot_id}")

    metadata_updates: dict[str, Any] = {}
    if context.is_added_snapshot(snapshot_ref.snapshot_id):
        metadata_updates["last_updated_ms"] = snapshot.timestamp_ms

    if update.ref_name == MAIN_BRANCH:
        metadata_updates["current_snapshot_id"] = snapshot_ref.snapshot_id
        if "last_updated_ms" not in metadata_updates:
            metadata_updates["last_updated_ms"] = datetime_to_millis(datetime.now().astimezone())

        metadata_updates["snapshot_log"] = base_metadata.snapshot_log + [
            SnapshotLogEntry(
                snapshot_id=snapshot_ref.snapshot_id,
                timestamp_ms=metadata_updates["last_updated_ms"],
            )
        ]

    metadata_updates["refs"] = {**base_metadata.refs, update.ref_name: snapshot_ref}
    context.add_update(update)
    return base_metadata.model_copy(update=metadata_updates)


@_apply_table_update.register(RemoveSnapshotsUpdate)
def _(update: RemoveSnapshotsUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    for remove_snapshot_id in update.snapshot_ids:
        if not any(snapshot.snapshot_id == remove_snapshot_id for snapshot in base_metadata.snapshots):
            raise ValueError(f"Snapshot with snapshot id {remove_snapshot_id} does not exist: {base_metadata.snapshots}")

    snapshots = [
        (
            snapshot.model_copy(update={"parent_snapshot_id": None})
            if snapshot.parent_snapshot_id in update.snapshot_ids
            else snapshot
        )
        for snapshot in base_metadata.snapshots
        if snapshot.snapshot_id not in update.snapshot_ids
    ]
    snapshot_log = [
        snapshot_log_entry
        for snapshot_log_entry in base_metadata.snapshot_log
        if snapshot_log_entry.snapshot_id not in update.snapshot_ids
    ]

    remove_ref_updates = (
        RemoveSnapshotRefUpdate(ref_name=ref_name)
        for ref_name, ref in base_metadata.refs.items()
        if ref.snapshot_id in update.snapshot_ids
    )
    remove_statistics_updates = (
        RemoveStatisticsUpdate(snapshot_id=statistics_file.snapshot_id)
        for statistics_file in base_metadata.statistics
        if statistics_file.snapshot_id in update.snapshot_ids
    )
    updates = itertools.chain(remove_ref_updates, remove_statistics_updates)
    new_metadata = base_metadata
    for upd in updates:
        new_metadata = _apply_table_update(upd, new_metadata, context)

    context.add_update(update)
    return new_metadata.model_copy(update={"snapshots": snapshots, "snapshot_log": snapshot_log})


@_apply_table_update.register(RemoveSnapshotRefUpdate)
def _(update: RemoveSnapshotRefUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    if update.ref_name not in base_metadata.refs:
        return base_metadata

    existing_ref = base_metadata.refs[update.ref_name]
    if base_metadata.snapshot_by_id(existing_ref.snapshot_id) is None:
        raise ValueError(f"Cannot remove {update.ref_name} ref with unknown snapshot {existing_ref.snapshot_id}")

    current_snapshot_id = None if update.ref_name == MAIN_BRANCH else base_metadata.current_snapshot_id

    metadata_refs = {ref_name: ref for ref_name, ref in base_metadata.refs.items() if ref_name != update.ref_name}
    context.add_update(update)
    return base_metadata.model_copy(update={"refs": metadata_refs, "current_snapshot_id": current_snapshot_id})


@_apply_table_update.register(AddSortOrderUpdate)
def _(update: AddSortOrderUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    for sort in base_metadata.sort_orders:
        # Only raise in case of a discrepancy
        if sort.order_id == update.sort_order.order_id and sort != update.sort_order:
            raise ValueError(f"Sort-order with id {sort.order_id} already exists: {sort}")

    context.add_update(update)
    return base_metadata.model_copy(
        update={
            "sort_orders": base_metadata.sort_orders + [update.sort_order],
        }
    )


@_apply_table_update.register(SetDefaultSortOrderUpdate)
def _(
    update: SetDefaultSortOrderUpdate,
    base_metadata: TableMetadata,
    context: _TableMetadataUpdateContext,
) -> TableMetadata:
    new_sort_order_id = update.sort_order_id
    if new_sort_order_id == -1:
        # The last added sort order should be in base_metadata.sort_orders at this point
        new_sort_order_id = max(sort_order.order_id for sort_order in base_metadata.sort_orders)
        if not context.is_added_sort_order(new_sort_order_id):
            raise ValueError("Cannot set current sort order to the last added one when no sort order has been added")

    if new_sort_order_id == base_metadata.default_sort_order_id:
        return base_metadata

    sort_order = base_metadata.sort_order_by_id(new_sort_order_id)
    if sort_order is None:
        raise ValueError(f"Sort order with id {new_sort_order_id} does not exist")

    context.add_update(update)
    return base_metadata.model_copy(update={"default_sort_order_id": new_sort_order_id})


@_apply_table_update.register(SetStatisticsUpdate)
def _(update: SetStatisticsUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    statistics = filter_statistics_by_snapshot_id(base_metadata.statistics, update.statistics.snapshot_id)
    context.add_update(update)

    return base_metadata.model_copy(update={"statistics": statistics + [update.statistics]})


@_apply_table_update.register(RemoveStatisticsUpdate)
def _(update: RemoveStatisticsUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    if not any(stat.snapshot_id == update.snapshot_id for stat in base_metadata.statistics):
        raise ValueError(f"Statistics with snapshot id {update.snapshot_id} does not exist")

    statistics = filter_statistics_by_snapshot_id(base_metadata.statistics, update.snapshot_id)
    context.add_update(update)

    return base_metadata.model_copy(update={"statistics": statistics})


@_apply_table_update.register(RemovePartitionSpecsUpdate)
def _(update: RemovePartitionSpecsUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    for remove_spec_id in update.spec_ids:
        if not any(spec.spec_id == remove_spec_id for spec in base_metadata.partition_specs):
            raise ValueError(f"Partition spec with id {remove_spec_id} does not exist")

    if base_metadata.default_spec_id in update.spec_ids:
        raise ValueError(f"Cannot remove default partition spec: {base_metadata.default_spec_id}")

    partition_specs = [spec for spec in base_metadata.partition_specs if spec.spec_id not in update.spec_ids]

    context.add_update(update)
    return base_metadata.model_copy(update={"partition_specs": partition_specs})


@_apply_table_update.register(RemoveSchemasUpdate)
def _(update: RemoveSchemasUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    # This method should error if any schemas do not exist.
    # It should error if the default schema is being removed.
    # Otherwise, remove the schemas listed in update.schema_ids.
    for remove_schema_id in update.schema_ids:
        if not any(schema.schema_id == remove_schema_id for schema in base_metadata.schemas):
            raise ValueError(f"Schema with schema id {remove_schema_id} does not exist")
        if base_metadata.current_schema_id == remove_schema_id:
            raise ValueError(f"Cannot remove current schema with id {remove_schema_id}")

    schemas = [schema for schema in base_metadata.schemas if schema.schema_id not in update.schema_ids]
    context.add_update(update)

    return base_metadata.model_copy(update={"schemas": schemas})


@_apply_table_update.register(SetPartitionStatisticsUpdate)
def _(update: SetPartitionStatisticsUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext) -> TableMetadata:
    partition_statistics = filter_statistics_by_snapshot_id(
        base_metadata.partition_statistics, update.partition_statistics.snapshot_id
    )
    context.add_update(update)

    return base_metadata.model_copy(update={"partition_statistics": partition_statistics + [update.partition_statistics]})


@_apply_table_update.register(RemovePartitionStatisticsUpdate)
def _(
    update: RemovePartitionStatisticsUpdate, base_metadata: TableMetadata, context: _TableMetadataUpdateContext
) -> TableMetadata:
    if not any(part_stat.snapshot_id == update.snapshot_id for part_stat in base_metadata.partition_statistics):
        raise ValueError(f"Partition Statistics with snapshot id {update.snapshot_id} does not exist")

    statistics = filter_statistics_by_snapshot_id(base_metadata.partition_statistics, update.snapshot_id)
    context.add_update(update)

    return base_metadata.model_copy(update={"partition_statistics": statistics})


def update_table_metadata(
    base_metadata: TableMetadata,
    updates: tuple[TableUpdate, ...],
    enforce_validation: bool = False,
    metadata_location: str | None = None,
) -> TableMetadata:
    """Update the table metadata with the given updates in one transaction.

    Args:
        base_metadata: The base metadata to be updated.
        updates: The updates in one transaction.
        enforce_validation: Whether to trigger validation after applying the updates.
        metadata_location: Current metadata location of the table

    Returns:
        The metadata with the updates applied.
    """
    context = _TableMetadataUpdateContext()
    new_metadata = base_metadata

    for update in updates:
        new_metadata = _apply_table_update(update, new_metadata, context)

    # Update last_updated_ms if it was not updated by update operations
    if context.has_changes():
        if metadata_location:
            new_metadata = _update_table_metadata_log(new_metadata, metadata_location, base_metadata.last_updated_ms)
        if base_metadata.last_updated_ms == new_metadata.last_updated_ms:
            new_metadata = new_metadata.model_copy(update={"last_updated_ms": datetime_to_millis(datetime.now().astimezone())})

    # Check correctness of partition spec, and sort order
    new_metadata.spec().check_compatible(new_metadata.schema())
    new_metadata.sort_order().check_compatible(new_metadata.schema())

    if enforce_validation:
        return TableMetadataUtil.parse_obj(new_metadata.model_dump())
    else:
        return new_metadata.model_copy(deep=True)


def _update_table_metadata_log(base_metadata: TableMetadata, metadata_location: str, last_updated_ms: int) -> TableMetadata:
    from pyiceberg.table import TableProperties

    """
    Update the metadata log of the table.

    Args:
        base_metadata: The base metadata to be updated.
        metadata_location: Current metadata location of the table
        last_updated_ms: The timestamp of the last update of table metadata

    Returns:
        The metadata with the updates applied to metadata-log.
    """
    max_metadata_log_entries = max(
        1,
        property_as_int(
            base_metadata.properties,
            TableProperties.METADATA_PREVIOUS_VERSIONS_MAX,
            TableProperties.METADATA_PREVIOUS_VERSIONS_MAX_DEFAULT,
        ),  # type: ignore
    )
    previous_metadata_log = base_metadata.metadata_log
    if len(base_metadata.metadata_log) >= max_metadata_log_entries:  # type: ignore
        remove_index = len(base_metadata.metadata_log) - max_metadata_log_entries + 1  # type: ignore
        previous_metadata_log = base_metadata.metadata_log[remove_index:]
    metadata_updates: dict[str, Any] = {
        "metadata_log": previous_metadata_log + [MetadataLogEntry(metadata_file=metadata_location, timestamp_ms=last_updated_ms)]
    }
    return base_metadata.model_copy(update=metadata_updates)


class ValidatableTableRequirement(IcebergBaseModel):
    type: str

    @abstractmethod
    def validate(self, base_metadata: TableMetadata | None) -> None:
        """Validate the requirement against the base metadata.

        Args:
            base_metadata: The base metadata to be validated against.

        Raises:
            CommitFailedException: When the requirement is not met.
        """
        ...


class AssertCreate(ValidatableTableRequirement):
    """The table must not already exist; used for create transactions."""

    type: Literal["assert-create"] = Field(default="assert-create")

    def validate(self, base_metadata: TableMetadata | None) -> None:
        if base_metadata is not None:
            raise CommitFailedException("Table already exists")


class AssertTableUUID(ValidatableTableRequirement):
    """The table UUID must match the requirement's `uuid`."""

    type: Literal["assert-table-uuid"] = Field(default="assert-table-uuid")
    uuid: uuid.UUID

    def validate(self, base_metadata: TableMetadata | None) -> None:
        if base_metadata is None:
            raise CommitFailedException("Requirement failed: current table metadata is missing")
        elif self.uuid != base_metadata.table_uuid:
            raise CommitFailedException(f"Table UUID does not match: {self.uuid} != {base_metadata.table_uuid}")


class AssertRefSnapshotId(ValidatableTableRequirement):
    """The table branch or tag identified by the requirement's `ref` must reference the requirement's `snapshot-id`.

    if `snapshot-id` is `null` or missing, the ref must not already exist.
    """

    type: Literal["assert-ref-snapshot-id"] = Field(default="assert-ref-snapshot-id")
    ref: str = Field(...)
    snapshot_id: int | None = Field(default=None, alias="snapshot-id")

    @model_serializer(mode="wrap")
    def serialize_model(self, handler: ModelWrapSerializerWithoutInfo) -> dict[str, Any]:
        partial_result = handler(self)
        # Ensure "snapshot-id" is always present, even if value is None
        return {**partial_result, "snapshot-id": self.snapshot_id}

    def validate(self, base_metadata: TableMetadata | None) -> None:
        if base_metadata is None:
            raise CommitFailedException("Requirement failed: current table metadata is missing")
        elif len(base_metadata.snapshots) == 0 and self.ref != MAIN_BRANCH:
            raise CommitFailedException(
                f"Requirement failed: Table has no snapshots and can only be written to the {MAIN_BRANCH} BRANCH."
            )
        elif snapshot_ref := base_metadata.refs.get(self.ref):
            ref_type = snapshot_ref.snapshot_ref_type
            if self.snapshot_id is None:
                raise CommitFailedException(f"Requirement failed: {ref_type} {self.ref} was created concurrently")
            elif self.snapshot_id != snapshot_ref.snapshot_id:
                raise CommitFailedException(
                    f"Requirement failed: {ref_type} {self.ref} has changed: "
                    f"expected id {self.snapshot_id}, found {snapshot_ref.snapshot_id}"
                )
        elif self.snapshot_id is not None:
            raise CommitFailedException(f"Requirement failed: branch or tag {self.ref} is missing, expected {self.snapshot_id}")


class AssertLastAssignedFieldId(ValidatableTableRequirement):
    """The table's last assigned column id must match the requirement's `last-assigned-field-id`."""

    type: Literal["assert-last-assigned-field-id"] = Field(default="assert-last-assigned-field-id")
    last_assigned_field_id: int = Field(..., alias="last-assigned-field-id")

    def validate(self, base_metadata: TableMetadata | None) -> None:
        if base_metadata is None:
            raise CommitFailedException("Requirement failed: current table metadata is missing")
        elif base_metadata.last_column_id != self.last_assigned_field_id:
            raise CommitFailedException(
                f"Requirement failed: last assigned field id has changed: "
                f"expected {self.last_assigned_field_id}, found {base_metadata.last_column_id}"
            )


class AssertCurrentSchemaId(ValidatableTableRequirement):
    """The table's current schema id must match the requirement's `current-schema-id`."""

    type: Literal["assert-current-schema-id"] = Field(default="assert-current-schema-id")
    current_schema_id: int = Field(..., alias="current-schema-id")

    def validate(self, base_metadata: TableMetadata | None) -> None:
        if base_metadata is None:
            raise CommitFailedException("Requirement failed: current table metadata is missing")
        elif self.current_schema_id != base_metadata.current_schema_id:
            raise CommitFailedException(
                f"Requirement failed: current schema id has changed: "
                f"expected {self.current_schema_id}, found {base_metadata.current_schema_id}"
            )


class AssertLastAssignedPartitionId(ValidatableTableRequirement):
    """The table's last assigned partition id must match the requirement's `last-assigned-partition-id`."""

    type: Literal["assert-last-assigned-partition-id"] = Field(default="assert-last-assigned-partition-id")
    last_assigned_partition_id: int | None = Field(..., alias="last-assigned-partition-id")

    def validate(self, base_metadata: TableMetadata | None) -> None:
        if base_metadata is None:
            raise CommitFailedException("Requirement failed: current table metadata is missing")
        elif base_metadata.last_partition_id != self.last_assigned_partition_id:
            raise CommitFailedException(
                f"Requirement failed: last assigned partition id has changed: "
                f"expected {self.last_assigned_partition_id}, found {base_metadata.last_partition_id}"
            )


class AssertDefaultSpecId(ValidatableTableRequirement):
    """The table's default spec id must match the requirement's `default-spec-id`."""

    type: Literal["assert-default-spec-id"] = Field(default="assert-default-spec-id")
    default_spec_id: int = Field(..., alias="default-spec-id")

    def validate(self, base_metadata: TableMetadata | None) -> None:
        if base_metadata is None:
            raise CommitFailedException("Requirement failed: current table metadata is missing")
        elif self.default_spec_id != base_metadata.default_spec_id:
            raise CommitFailedException(
                f"Requirement failed: default spec id has changed: "
                f"expected {self.default_spec_id}, found {base_metadata.default_spec_id}"
            )


class AssertDefaultSortOrderId(ValidatableTableRequirement):
    """The table's default sort order id must match the requirement's `default-sort-order-id`."""

    type: Literal["assert-default-sort-order-id"] = Field(default="assert-default-sort-order-id")
    default_sort_order_id: int = Field(..., alias="default-sort-order-id")

    def validate(self, base_metadata: TableMetadata | None) -> None:
        if base_metadata is None:
            raise CommitFailedException("Requirement failed: current table metadata is missing")
        elif self.default_sort_order_id != base_metadata.default_sort_order_id:
            raise CommitFailedException(
                f"Requirement failed: default sort order id has changed: "
                f"expected {self.default_sort_order_id}, found {base_metadata.default_sort_order_id}"
            )


TableRequirement = Annotated[
    AssertCreate
    | AssertTableUUID
    | AssertRefSnapshotId
    | AssertLastAssignedFieldId
    | AssertCurrentSchemaId
    | AssertLastAssignedPartitionId
    | AssertDefaultSpecId
    | AssertDefaultSortOrderId,
    Field(discriminator="type"),
]

UpdatesAndRequirements = tuple[tuple[TableUpdate, ...], tuple[TableRequirement, ...]]
