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

import time
import warnings
from collections import defaultdict
from collections.abc import Iterable, Mapping
from enum import Enum
from typing import TYPE_CHECKING, Any

from pydantic import Field, PrivateAttr, model_serializer

from pyiceberg.io import FileIO
from pyiceberg.manifest import DataFile, DataFileContent, ManifestFile, _manifests
from pyiceberg.partitioning import UNPARTITIONED_PARTITION_SPEC, PartitionSpec
from pyiceberg.schema import Schema

if TYPE_CHECKING:
    from pyiceberg.table.metadata import TableMetadata
from pyiceberg.typedef import IcebergBaseModel

ADDED_DATA_FILES = "added-data-files"
ADDED_DELETE_FILES = "added-delete-files"
ADDED_EQUALITY_DELETES = "added-equality-deletes"
ADDED_FILE_SIZE = "added-files-size"
ADDED_POSITION_DELETES = "added-position-deletes"
ADDED_POSITION_DELETE_FILES = "added-position-delete-files"
ADDED_RECORDS = "added-records"
DELETED_DATA_FILES = "deleted-data-files"
DELETED_RECORDS = "deleted-records"
ADDED_EQUALITY_DELETE_FILES = "added-equality-delete-files"
REMOVED_DELETE_FILES = "removed-delete-files"
REMOVED_EQUALITY_DELETES = "removed-equality-deletes"
REMOVED_EQUALITY_DELETE_FILES = "removed-equality-delete-files"
REMOVED_FILE_SIZE = "removed-files-size"
REMOVED_POSITION_DELETES = "removed-position-deletes"
REMOVED_POSITION_DELETE_FILES = "removed-position-delete-files"
TOTAL_EQUALITY_DELETES = "total-equality-deletes"
TOTAL_POSITION_DELETES = "total-position-deletes"
TOTAL_DATA_FILES = "total-data-files"
TOTAL_DELETE_FILES = "total-delete-files"
TOTAL_RECORDS = "total-records"
TOTAL_FILE_SIZE = "total-files-size"
CHANGED_PARTITION_COUNT_PROP = "changed-partition-count"
CHANGED_PARTITION_PREFIX = "partitions."
PARTITION_SUMMARY_PROP = "partition-summaries-included"
OPERATION = "operation"

INITIAL_SEQUENCE_NUMBER = 0


class Operation(Enum):
    """Describes the operation.

    Possible operation values are:
        - append: Only data files were added and no files were removed.
        - replace: Data and delete files were added and removed without changing table data;
            i.e., compaction, changing the data file format, or relocating data files.
        - overwrite: Data and delete files were added and removed in a logical overwrite operation.
        - delete: Data files were removed and their contents logically deleted and/or delete files
            were added to delete rows.
    """

    APPEND = "append"
    REPLACE = "replace"
    OVERWRITE = "overwrite"
    DELETE = "delete"

    def __repr__(self) -> str:
        """Return the string representation of the Operation class."""
        return f"Operation.{self.name}"


class UpdateMetrics:
    added_file_size: int
    removed_file_size: int
    added_data_files: int
    removed_data_files: int
    added_eq_delete_files: int
    removed_eq_delete_files: int
    added_pos_delete_files: int
    removed_pos_delete_files: int
    added_delete_files: int
    removed_delete_files: int
    added_records: int
    deleted_records: int
    added_pos_deletes: int
    removed_pos_deletes: int
    added_eq_deletes: int
    removed_eq_deletes: int

    def __init__(self) -> None:
        self.added_file_size = 0
        self.removed_file_size = 0
        self.added_data_files = 0
        self.removed_data_files = 0
        self.added_eq_delete_files = 0
        self.removed_eq_delete_files = 0
        self.added_pos_delete_files = 0
        self.removed_pos_delete_files = 0
        self.added_delete_files = 0
        self.removed_delete_files = 0
        self.added_records = 0
        self.deleted_records = 0
        self.added_pos_deletes = 0
        self.removed_pos_deletes = 0
        self.added_eq_deletes = 0
        self.removed_eq_deletes = 0

    def add_file(self, data_file: DataFile) -> None:
        self.added_file_size += data_file.file_size_in_bytes

        if data_file.content == DataFileContent.DATA:
            self.added_data_files += 1
            self.added_records += data_file.record_count
        elif data_file.content == DataFileContent.POSITION_DELETES:
            self.added_delete_files += 1
            self.added_pos_delete_files += 1
            self.added_pos_deletes += data_file.record_count
        elif data_file.content == DataFileContent.EQUALITY_DELETES:
            self.added_delete_files += 1
            self.added_eq_delete_files += 1
            self.added_eq_deletes += data_file.record_count
        else:
            raise ValueError(f"Unknown data file content: {data_file.content}")

    def remove_file(self, data_file: DataFile) -> None:
        self.removed_file_size += data_file.file_size_in_bytes

        if data_file.content == DataFileContent.DATA:
            self.removed_data_files += 1
            self.deleted_records += data_file.record_count
        elif data_file.content == DataFileContent.POSITION_DELETES:
            self.removed_delete_files += 1
            self.removed_pos_delete_files += 1
            self.removed_pos_deletes += data_file.record_count
        elif data_file.content == DataFileContent.EQUALITY_DELETES:
            self.removed_delete_files += 1
            self.removed_eq_delete_files += 1
            self.removed_eq_deletes += data_file.record_count
        else:
            raise ValueError(f"Unknown data file content: {data_file.content}")

    def to_dict(self) -> dict[str, str]:
        properties: dict[str, str] = {}
        set_when_positive(properties, self.added_file_size, ADDED_FILE_SIZE)
        set_when_positive(properties, self.removed_file_size, REMOVED_FILE_SIZE)
        set_when_positive(properties, self.added_data_files, ADDED_DATA_FILES)
        set_when_positive(properties, self.removed_data_files, DELETED_DATA_FILES)
        set_when_positive(properties, self.added_eq_delete_files, ADDED_EQUALITY_DELETE_FILES)
        set_when_positive(properties, self.removed_eq_delete_files, REMOVED_EQUALITY_DELETE_FILES)
        set_when_positive(properties, self.added_pos_delete_files, ADDED_POSITION_DELETE_FILES)
        set_when_positive(properties, self.removed_pos_delete_files, REMOVED_POSITION_DELETE_FILES)
        set_when_positive(properties, self.added_delete_files, ADDED_DELETE_FILES)
        set_when_positive(properties, self.removed_delete_files, REMOVED_DELETE_FILES)
        set_when_positive(properties, self.added_records, ADDED_RECORDS)
        set_when_positive(properties, self.deleted_records, DELETED_RECORDS)
        set_when_positive(properties, self.added_pos_deletes, ADDED_POSITION_DELETES)
        set_when_positive(properties, self.removed_pos_deletes, REMOVED_POSITION_DELETES)
        set_when_positive(properties, self.added_eq_deletes, ADDED_EQUALITY_DELETES)
        set_when_positive(properties, self.removed_eq_deletes, REMOVED_EQUALITY_DELETES)
        return properties


class Summary(IcebergBaseModel, Mapping[str, str]):
    """A class that stores the summary information for a Snapshot.

    The snapshot summary’s operation field is used by some operations,
    like snapshot expiration, to skip processing certain snapshots.
    """

    operation: Operation = Field()
    _additional_properties: dict[str, str] = PrivateAttr()

    def __init__(self, operation: Operation | None = None, **data: Any) -> None:
        if operation is None:
            warnings.warn("Encountered invalid snapshot summary: operation is missing, defaulting to overwrite", stacklevel=2)
            operation = Operation.OVERWRITE
        super().__init__(operation=operation, **data)
        self._additional_properties = data

    def __getitem__(self, __key: str) -> Any | None:  # type: ignore
        """Return a key as it is a map."""
        if __key.lower() == "operation":
            return self.operation
        else:
            return self._additional_properties.get(__key)

    def __setitem__(self, key: str, value: Any) -> None:
        """Set a key as it is a map."""
        if key.lower() == "operation":
            self.operation = value
        else:
            self._additional_properties[key] = value

    def __len__(self) -> int:
        """Return the number of keys in the summary."""
        # Operation is required
        return 1 + len(self._additional_properties)

    @model_serializer
    def ser_model(self) -> dict[str, str]:
        return {
            "operation": str(self.operation.value),
            **self._additional_properties,
        }

    @property
    def additional_properties(self) -> dict[str, str]:
        return self._additional_properties

    def __repr__(self) -> str:
        """Return the string representation of the Summary class."""
        repr_properties = f", **{repr(self._additional_properties)}" if self._additional_properties else ""
        return f"Summary({repr(self.operation)}{repr_properties})"

    def __eq__(self, other: Any) -> bool:
        """Compare if the summary is equal to another summary."""
        return (
            self.operation == other.operation and self.additional_properties == other.additional_properties
            if isinstance(other, Summary)
            else False
        )


class Snapshot(IcebergBaseModel):
    snapshot_id: int = Field(alias="snapshot-id")
    parent_snapshot_id: int | None = Field(alias="parent-snapshot-id", default=None)
    sequence_number: int | None = Field(alias="sequence-number", default=INITIAL_SEQUENCE_NUMBER)
    timestamp_ms: int = Field(alias="timestamp-ms", default_factory=lambda: int(time.time() * 1000))
    manifest_list: str = Field(alias="manifest-list", description="Location of the snapshot's manifest list file")
    summary: Summary | None = Field(default=None)
    schema_id: int | None = Field(alias="schema-id", default=None)
    first_row_id: int | None = Field(
        alias="first-row-id", default=None, description="assigned to the first row in the first data file in the first manifest"
    )
    added_rows: int | None = Field(
        alias="added-rows", default=None, description="The upper bound of the number of rows with assigned row IDs"
    )

    def __str__(self) -> str:
        """Return the string representation of the Snapshot class."""
        operation = f"{self.summary.operation}: " if self.summary else ""
        parent_id = f", parent_id={self.parent_snapshot_id}" if self.parent_snapshot_id else ""
        schema_id = f", schema_id={self.schema_id}" if self.schema_id is not None else ""
        result_str = f"{operation}id={self.snapshot_id}{parent_id}{schema_id}"
        return result_str

    def __repr__(self) -> str:
        """Return the string representation of the Snapshot class."""
        fields = [
            f"snapshot_id={self.snapshot_id}",
            f"parent_snapshot_id={self.parent_snapshot_id}",
            f"sequence_number={self.sequence_number}",
            f"timestamp_ms={self.timestamp_ms}",
            f"manifest_list='{self.manifest_list}'",
            f"summary={repr(self.summary)}" if self.summary else None,
            f"schema_id={self.schema_id}" if self.schema_id is not None else None,
            f"first_row_id={self.first_row_id}" if self.first_row_id is not None else None,
            f"added_rows={self.added_rows}" if self.added_rows is not None else None,
        ]
        filtered_fields = [field for field in fields if field is not None]
        return f"Snapshot({', '.join(filtered_fields)})"

    def manifests(self, io: FileIO) -> list[ManifestFile]:
        """Return the manifests for the given snapshot."""
        return list(_manifests(io, self.manifest_list))


class MetadataLogEntry(IcebergBaseModel):
    metadata_file: str = Field(alias="metadata-file")
    timestamp_ms: int = Field(alias="timestamp-ms")


class SnapshotLogEntry(IcebergBaseModel):
    snapshot_id: int = Field(alias="snapshot-id")
    timestamp_ms: int = Field(alias="timestamp-ms")


class SnapshotSummaryCollector:
    metrics: UpdateMetrics
    partition_metrics: defaultdict[str, UpdateMetrics]
    max_changed_partitions_for_summaries: int

    def __init__(self, partition_summary_limit: int = 0) -> None:
        self.metrics = UpdateMetrics()
        self.partition_metrics = defaultdict(UpdateMetrics)
        self.max_changed_partitions_for_summaries = partition_summary_limit

    def set_partition_summary_limit(self, limit: int) -> None:
        self.max_changed_partitions_for_summaries = limit

    def add_file(self, data_file: DataFile, schema: Schema, partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC) -> None:
        self.metrics.add_file(data_file)
        if len(data_file.partition) > 0:
            self.update_partition_metrics(partition_spec=partition_spec, file=data_file, is_add_file=True, schema=schema)

    def remove_file(
        self, data_file: DataFile, schema: Schema, partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC
    ) -> None:
        self.metrics.remove_file(data_file)
        if len(data_file.partition) > 0:
            self.update_partition_metrics(partition_spec=partition_spec, file=data_file, is_add_file=False, schema=schema)

    def update_partition_metrics(self, partition_spec: PartitionSpec, file: DataFile, is_add_file: bool, schema: Schema) -> None:
        partition_path = partition_spec.partition_to_path(file.partition, schema)
        partition_metrics: UpdateMetrics = self.partition_metrics[partition_path]

        if is_add_file:
            partition_metrics.add_file(file)
        else:
            partition_metrics.remove_file(file)

    def build(self) -> dict[str, str]:
        properties = self.metrics.to_dict()
        changed_partitions_size = len(self.partition_metrics)
        set_when_positive(properties, changed_partitions_size, CHANGED_PARTITION_COUNT_PROP)
        if changed_partitions_size <= self.max_changed_partitions_for_summaries:
            if changed_partitions_size > 0:
                properties[PARTITION_SUMMARY_PROP] = "true"
            for partition_path, update_metrics_partition in self.partition_metrics.items():
                if (summary := self._partition_summary(update_metrics_partition)) and len(summary) != 0:
                    properties[CHANGED_PARTITION_PREFIX + partition_path] = summary

        return properties

    def _partition_summary(self, update_metrics: UpdateMetrics) -> str:
        return ",".join([f"{prop}={val}" for prop, val in update_metrics.to_dict().items()])


def update_snapshot_summaries(summary: Summary, previous_summary: Mapping[str, str] | None = None) -> Summary:
    if summary.operation not in {Operation.APPEND, Operation.OVERWRITE, Operation.DELETE}:
        raise ValueError(f"Operation not implemented: {summary.operation}")

    if not previous_summary:
        previous_summary = {
            TOTAL_DATA_FILES: "0",
            TOTAL_DELETE_FILES: "0",
            TOTAL_RECORDS: "0",
            TOTAL_FILE_SIZE: "0",
            TOTAL_POSITION_DELETES: "0",
            TOTAL_EQUALITY_DELETES: "0",
        }

    def _update_totals(total_property: str, added_property: str, removed_property: str) -> None:
        if previous_total_str := previous_summary.get(total_property):
            try:
                new_total = int(previous_total_str)
                if new_total >= 0 and (added := summary.get(added_property)):
                    new_total += int(added)
                if new_total >= 0 and (removed := summary.get(removed_property)):
                    new_total -= int(removed)
            except ValueError as e:
                raise ValueError(f"Could not parse summary property {total_property} to an int: {previous_total_str}") from e

            if new_total >= 0:
                summary[total_property] = str(new_total)

    _update_totals(
        total_property=TOTAL_DATA_FILES,
        added_property=ADDED_DATA_FILES,
        removed_property=DELETED_DATA_FILES,
    )
    _update_totals(
        total_property=TOTAL_DELETE_FILES,
        added_property=ADDED_DELETE_FILES,
        removed_property=REMOVED_DELETE_FILES,
    )
    _update_totals(
        total_property=TOTAL_RECORDS,
        added_property=ADDED_RECORDS,
        removed_property=DELETED_RECORDS,
    )
    _update_totals(
        total_property=TOTAL_FILE_SIZE,
        added_property=ADDED_FILE_SIZE,
        removed_property=REMOVED_FILE_SIZE,
    )
    _update_totals(
        total_property=TOTAL_POSITION_DELETES,
        added_property=ADDED_POSITION_DELETES,
        removed_property=REMOVED_POSITION_DELETES,
    )
    _update_totals(
        total_property=TOTAL_EQUALITY_DELETES,
        added_property=ADDED_EQUALITY_DELETES,
        removed_property=REMOVED_EQUALITY_DELETES,
    )

    return summary


def set_when_positive(properties: dict[str, str], num: int, property_name: str) -> None:
    if num > 0:
        properties[property_name] = str(num)


def ancestors_of(current_snapshot: Snapshot | None, table_metadata: TableMetadata) -> Iterable[Snapshot]:
    """Get the ancestors of and including the given snapshot."""
    snapshot = current_snapshot
    while snapshot is not None:
        yield snapshot
        if snapshot.parent_snapshot_id is None:
            break
        snapshot = table_metadata.snapshot_by_id(snapshot.parent_snapshot_id)


def ancestors_between(from_snapshot: Snapshot | None, to_snapshot: Snapshot, table_metadata: TableMetadata) -> Iterable[Snapshot]:
    """Get the ancestors of and including the given snapshot between the to and from snapshots."""
    if from_snapshot is not None:
        for snapshot in ancestors_of(to_snapshot, table_metadata):
            yield snapshot
            if snapshot == from_snapshot:
                break
    else:
        yield from ancestors_of(to_snapshot, table_metadata)


def latest_ancestor_before_timestamp(table_metadata: TableMetadata, timestamp_ms: int) -> Snapshot | None:
    """Find the latest ancestor snapshot whose timestamp is before the provided timestamp.

    Args:
        table_metadata: The table metadata for a table
        timestamp_ms: lookup snapshots strictly before this timestamp

    Returns:
        The latest ancestor snapshot older than the timestamp, or None if not found.
    """
    result: Snapshot | None = None
    result_timestamp: int = 0

    for ancestor in ancestors_of(table_metadata.current_snapshot(), table_metadata):
        if timestamp_ms > ancestor.timestamp_ms > result_timestamp:
            result = ancestor
            result_timestamp = ancestor.timestamp_ms

    return result
