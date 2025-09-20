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

from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Set, Tuple

from pyiceberg.conversions import from_bytes
from pyiceberg.manifest import DataFileContent, ManifestContent, ManifestFile, PartitionFieldSummary
from pyiceberg.partitioning import PartitionSpec
from pyiceberg.table.snapshots import Snapshot, ancestors_of
from pyiceberg.types import PrimitiveType
from pyiceberg.utils.concurrent import ExecutorFactory
from pyiceberg.utils.singleton import _convert_to_hashable_type

if TYPE_CHECKING:
    import pyarrow as pa

    from pyiceberg.table import Table


class InspectTable:
    tbl: Table

    def __init__(self, tbl: Table) -> None:
        self.tbl = tbl

        try:
            import pyarrow as pa  # noqa
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError("For metadata operations PyArrow needs to be installed") from e

    def _get_snapshot(self, snapshot_id: Optional[int] = None) -> Snapshot:
        if snapshot_id is not None:
            if snapshot := self.tbl.metadata.snapshot_by_id(snapshot_id):
                return snapshot
            else:
                raise ValueError(f"Cannot find snapshot with ID {snapshot_id}")

        if snapshot := self.tbl.metadata.current_snapshot():
            return snapshot
        else:
            raise ValueError("Cannot get a snapshot as the table does not have any.")

    def snapshots(self) -> "pa.Table":
        import pyarrow as pa

        snapshots_schema = pa.schema(
            [
                pa.field("committed_at", pa.timestamp(unit="ms"), nullable=False),
                pa.field("snapshot_id", pa.int64(), nullable=False),
                pa.field("parent_id", pa.int64(), nullable=True),
                pa.field("operation", pa.string(), nullable=True),
                pa.field("manifest_list", pa.string(), nullable=False),
                pa.field("summary", pa.map_(pa.string(), pa.string()), nullable=True),
            ]
        )
        snapshots = []
        for snapshot in self.tbl.metadata.snapshots:
            if summary := snapshot.summary:
                operation = summary.operation.value
                additional_properties = snapshot.summary.additional_properties
            else:
                operation = None
                additional_properties = None

            snapshots.append(
                {
                    "committed_at": datetime.fromtimestamp(snapshot.timestamp_ms / 1000.0, tz=timezone.utc),
                    "snapshot_id": snapshot.snapshot_id,
                    "parent_id": snapshot.parent_snapshot_id,
                    "operation": str(operation),
                    "manifest_list": snapshot.manifest_list,
                    "summary": additional_properties,
                }
            )

        return pa.Table.from_pylist(
            snapshots,
            schema=snapshots_schema,
        )

    def entries(self, snapshot_id: Optional[int] = None) -> "pa.Table":
        import pyarrow as pa

        from pyiceberg.io.pyarrow import schema_to_pyarrow

        schema = self.tbl.metadata.schema()

        readable_metrics_struct = []

        def _readable_metrics_struct(bound_type: PrimitiveType) -> pa.StructType:
            pa_bound_type = schema_to_pyarrow(bound_type)
            return pa.struct(
                [
                    pa.field("column_size", pa.int64(), nullable=True),
                    pa.field("value_count", pa.int64(), nullable=True),
                    pa.field("null_value_count", pa.int64(), nullable=True),
                    pa.field("nan_value_count", pa.int64(), nullable=True),
                    pa.field("lower_bound", pa_bound_type, nullable=True),
                    pa.field("upper_bound", pa_bound_type, nullable=True),
                ]
            )

        for field in self.tbl.metadata.schema().fields:
            readable_metrics_struct.append(
                pa.field(schema.find_column_name(field.field_id), _readable_metrics_struct(field.field_type), nullable=False)
            )

        partition_record = self.tbl.metadata.specs_struct()
        pa_record_struct = schema_to_pyarrow(partition_record)

        entries_schema = pa.schema(
            [
                pa.field("status", pa.int8(), nullable=False),
                pa.field("snapshot_id", pa.int64(), nullable=False),
                pa.field("sequence_number", pa.int64(), nullable=False),
                pa.field("file_sequence_number", pa.int64(), nullable=False),
                pa.field(
                    "data_file",
                    pa.struct(
                        [
                            pa.field("content", pa.int8(), nullable=False),
                            pa.field("file_path", pa.string(), nullable=False),
                            pa.field("file_format", pa.string(), nullable=False),
                            pa.field("partition", pa_record_struct, nullable=False),
                            pa.field("record_count", pa.int64(), nullable=False),
                            pa.field("file_size_in_bytes", pa.int64(), nullable=False),
                            pa.field("column_sizes", pa.map_(pa.int32(), pa.int64()), nullable=True),
                            pa.field("value_counts", pa.map_(pa.int32(), pa.int64()), nullable=True),
                            pa.field("null_value_counts", pa.map_(pa.int32(), pa.int64()), nullable=True),
                            pa.field("nan_value_counts", pa.map_(pa.int32(), pa.int64()), nullable=True),
                            pa.field("lower_bounds", pa.map_(pa.int32(), pa.binary()), nullable=True),
                            pa.field("upper_bounds", pa.map_(pa.int32(), pa.binary()), nullable=True),
                            pa.field("key_metadata", pa.binary(), nullable=True),
                            pa.field("split_offsets", pa.list_(pa.int64()), nullable=True),
                            pa.field("equality_ids", pa.list_(pa.int32()), nullable=True),
                            pa.field("sort_order_id", pa.int32(), nullable=True),
                        ]
                    ),
                    nullable=False,
                ),
                pa.field("readable_metrics", pa.struct(readable_metrics_struct), nullable=True),
            ]
        )

        entries = []
        snapshot = self._get_snapshot(snapshot_id)
        for manifest in snapshot.manifests(self.tbl.io):
            for entry in manifest.fetch_manifest_entry(io=self.tbl.io, discard_deleted=False):
                column_sizes = entry.data_file.column_sizes or {}
                value_counts = entry.data_file.value_counts or {}
                null_value_counts = entry.data_file.null_value_counts or {}
                nan_value_counts = entry.data_file.nan_value_counts or {}
                lower_bounds = entry.data_file.lower_bounds or {}
                upper_bounds = entry.data_file.upper_bounds or {}
                readable_metrics = {
                    schema.find_column_name(field.field_id): {
                        "column_size": column_sizes.get(field.field_id),
                        "value_count": value_counts.get(field.field_id),
                        "null_value_count": null_value_counts.get(field.field_id),
                        "nan_value_count": nan_value_counts.get(field.field_id),
                        # Makes them readable
                        "lower_bound": from_bytes(field.field_type, lower_bound)
                        if (lower_bound := lower_bounds.get(field.field_id))
                        else None,
                        "upper_bound": from_bytes(field.field_type, upper_bound)
                        if (upper_bound := upper_bounds.get(field.field_id))
                        else None,
                    }
                    for field in self.tbl.metadata.schema().fields
                }

                partition = entry.data_file.partition
                partition_record_dict = {
                    field.name: partition[pos]
                    for pos, field in enumerate(self.tbl.metadata.specs()[manifest.partition_spec_id].fields)
                }

                entries.append(
                    {
                        "status": entry.status.value,
                        "snapshot_id": entry.snapshot_id,
                        "sequence_number": entry.sequence_number,
                        "file_sequence_number": entry.file_sequence_number,
                        "data_file": {
                            "content": entry.data_file.content,
                            "file_path": entry.data_file.file_path,
                            "file_format": entry.data_file.file_format,
                            "partition": partition_record_dict,
                            "record_count": entry.data_file.record_count,
                            "file_size_in_bytes": entry.data_file.file_size_in_bytes,
                            "column_sizes": dict(entry.data_file.column_sizes),
                            "value_counts": dict(entry.data_file.value_counts or {}),
                            "null_value_counts": dict(entry.data_file.null_value_counts or {}),
                            "nan_value_counts": dict(entry.data_file.nan_value_counts or {}),
                            "lower_bounds": entry.data_file.lower_bounds,
                            "upper_bounds": entry.data_file.upper_bounds,
                            "key_metadata": entry.data_file.key_metadata,
                            "split_offsets": entry.data_file.split_offsets,
                            "equality_ids": entry.data_file.equality_ids,
                            "sort_order_id": entry.data_file.sort_order_id,
                            "spec_id": entry.data_file.spec_id,
                        },
                        "readable_metrics": readable_metrics,
                    }
                )

        return pa.Table.from_pylist(
            entries,
            schema=entries_schema,
        )

    def refs(self) -> "pa.Table":
        import pyarrow as pa

        ref_schema = pa.schema(
            [
                pa.field("name", pa.string(), nullable=False),
                pa.field("type", pa.dictionary(pa.int32(), pa.string()), nullable=False),
                pa.field("snapshot_id", pa.int64(), nullable=False),
                pa.field("max_reference_age_in_ms", pa.int64(), nullable=True),
                pa.field("min_snapshots_to_keep", pa.int32(), nullable=True),
                pa.field("max_snapshot_age_in_ms", pa.int64(), nullable=True),
            ]
        )

        ref_results = []
        for ref in self.tbl.metadata.refs:
            if snapshot_ref := self.tbl.metadata.refs.get(ref):
                ref_results.append(
                    {
                        "name": ref,
                        "type": snapshot_ref.snapshot_ref_type.upper(),
                        "snapshot_id": snapshot_ref.snapshot_id,
                        "max_reference_age_in_ms": snapshot_ref.max_ref_age_ms,
                        "min_snapshots_to_keep": snapshot_ref.min_snapshots_to_keep,
                        "max_snapshot_age_in_ms": snapshot_ref.max_snapshot_age_ms,
                    }
                )

        return pa.Table.from_pylist(ref_results, schema=ref_schema)

    def partitions(self, snapshot_id: Optional[int] = None) -> "pa.Table":
        import pyarrow as pa

        from pyiceberg.io.pyarrow import schema_to_pyarrow

        table_schema = pa.schema(
            [
                pa.field("record_count", pa.int64(), nullable=False),
                pa.field("file_count", pa.int32(), nullable=False),
                pa.field("total_data_file_size_in_bytes", pa.int64(), nullable=False),
                pa.field("position_delete_record_count", pa.int64(), nullable=False),
                pa.field("position_delete_file_count", pa.int32(), nullable=False),
                pa.field("equality_delete_record_count", pa.int64(), nullable=False),
                pa.field("equality_delete_file_count", pa.int32(), nullable=False),
                pa.field("last_updated_at", pa.timestamp(unit="ms"), nullable=True),
                pa.field("last_updated_snapshot_id", pa.int64(), nullable=True),
            ]
        )

        partition_record = self.tbl.metadata.specs_struct()
        has_partitions = len(partition_record.fields) > 0

        if has_partitions:
            pa_record_struct = schema_to_pyarrow(partition_record)
            partitions_schema = pa.schema(
                [
                    pa.field("partition", pa_record_struct, nullable=False),
                    pa.field("spec_id", pa.int32(), nullable=False),
                ]
            )

            table_schema = pa.unify_schemas([partitions_schema, table_schema])

        snapshot = self._get_snapshot(snapshot_id)
        executor = ExecutorFactory.get_or_create()
        local_partitions_maps = executor.map(self._process_manifest, snapshot.manifests(self.tbl.io))

        partitions_map: Dict[Tuple[str, Any], Any] = {}
        for local_map in local_partitions_maps:
            for partition_record_key, partition_row in local_map.items():
                if partition_record_key not in partitions_map:
                    partitions_map[partition_record_key] = partition_row
                else:
                    existing = partitions_map[partition_record_key]
                    existing["record_count"] += partition_row["record_count"]
                    existing["file_count"] += partition_row["file_count"]
                    existing["total_data_file_size_in_bytes"] += partition_row["total_data_file_size_in_bytes"]
                    existing["position_delete_record_count"] += partition_row["position_delete_record_count"]
                    existing["position_delete_file_count"] += partition_row["position_delete_file_count"]
                    existing["equality_delete_record_count"] += partition_row["equality_delete_record_count"]
                    existing["equality_delete_file_count"] += partition_row["equality_delete_file_count"]

                    if partition_row["last_updated_at"] and (
                        not existing["last_updated_at"] or partition_row["last_updated_at"] > existing["last_updated_at"]
                    ):
                        existing["last_updated_at"] = partition_row["last_updated_at"]
                        existing["last_updated_snapshot_id"] = partition_row["last_updated_snapshot_id"]

        return pa.Table.from_pylist(
            partitions_map.values(),
            schema=table_schema,
        )

    def _process_manifest(self, manifest: ManifestFile) -> Dict[Tuple[str, Any], Any]:
        partitions_map: Dict[Tuple[str, Any], Any] = {}
        for entry in manifest.fetch_manifest_entry(io=self.tbl.io):
            partition = entry.data_file.partition
            partition_record_dict = {
                field.name: partition[pos]
                for pos, field in enumerate(self.tbl.metadata.specs()[manifest.partition_spec_id].fields)
            }
            entry_snapshot = self.tbl.snapshot_by_id(entry.snapshot_id) if entry.snapshot_id is not None else None

            partition_record_key = _convert_to_hashable_type(partition_record_dict)
            if partition_record_key not in partitions_map:
                partitions_map[partition_record_key] = {
                    "partition": partition_record_dict,
                    "spec_id": entry.data_file.spec_id,
                    "record_count": 0,
                    "file_count": 0,
                    "total_data_file_size_in_bytes": 0,
                    "position_delete_record_count": 0,
                    "position_delete_file_count": 0,
                    "equality_delete_record_count": 0,
                    "equality_delete_file_count": 0,
                    "last_updated_at": entry_snapshot.timestamp_ms if entry_snapshot else None,
                    "last_updated_snapshot_id": entry_snapshot.snapshot_id if entry_snapshot else None,
                }

            partition_row = partitions_map[partition_record_key]

            if entry_snapshot is not None:
                if (
                    partition_row["last_updated_at"] is None
                    or partition_row["last_updated_snapshot_id"] < entry_snapshot.timestamp_ms
                ):
                    partition_row["last_updated_at"] = entry_snapshot.timestamp_ms
                    partition_row["last_updated_snapshot_id"] = entry_snapshot.snapshot_id

            if entry.data_file.content == DataFileContent.DATA:
                partition_row["record_count"] += entry.data_file.record_count
                partition_row["file_count"] += 1
                partition_row["total_data_file_size_in_bytes"] += entry.data_file.file_size_in_bytes
            elif entry.data_file.content == DataFileContent.POSITION_DELETES:
                partition_row["position_delete_record_count"] += entry.data_file.record_count
                partition_row["position_delete_file_count"] += 1
            elif entry.data_file.content == DataFileContent.EQUALITY_DELETES:
                partition_row["equality_delete_record_count"] += entry.data_file.record_count
                partition_row["equality_delete_file_count"] += 1
            else:
                raise ValueError(f"Unknown DataFileContent ({entry.data_file.content})")

        return partitions_map

    def _get_manifests_schema(self) -> "pa.Schema":
        import pyarrow as pa

        partition_summary_schema = pa.struct(
            [
                pa.field("contains_null", pa.bool_(), nullable=False),
                pa.field("contains_nan", pa.bool_(), nullable=True),
                pa.field("lower_bound", pa.string(), nullable=True),
                pa.field("upper_bound", pa.string(), nullable=True),
            ]
        )

        manifest_schema = pa.schema(
            [
                pa.field("content", pa.int8(), nullable=False),
                pa.field("path", pa.string(), nullable=False),
                pa.field("length", pa.int64(), nullable=False),
                pa.field("partition_spec_id", pa.int32(), nullable=False),
                pa.field("added_snapshot_id", pa.int64(), nullable=False),
                pa.field("added_data_files_count", pa.int32(), nullable=False),
                pa.field("existing_data_files_count", pa.int32(), nullable=False),
                pa.field("deleted_data_files_count", pa.int32(), nullable=False),
                pa.field("added_delete_files_count", pa.int32(), nullable=False),
                pa.field("existing_delete_files_count", pa.int32(), nullable=False),
                pa.field("deleted_delete_files_count", pa.int32(), nullable=False),
                pa.field("partition_summaries", pa.list_(partition_summary_schema), nullable=False),
            ]
        )
        return manifest_schema

    def _get_all_manifests_schema(self) -> "pa.Schema":
        import pyarrow as pa

        all_manifests_schema = self._get_manifests_schema()
        all_manifests_schema = all_manifests_schema.append(pa.field("reference_snapshot_id", pa.int64(), nullable=False))
        return all_manifests_schema

    def _generate_manifests_table(self, snapshot: Optional[Snapshot], is_all_manifests_table: bool = False) -> "pa.Table":
        import pyarrow as pa

        def _partition_summaries_to_rows(
            spec: PartitionSpec, partition_summaries: List[PartitionFieldSummary]
        ) -> List[Dict[str, Any]]:
            rows = []
            for i, field_summary in enumerate(partition_summaries):
                field = spec.fields[i]
                partition_field_type = spec.partition_type(self.tbl.schema()).fields[i].field_type
                lower_bound = (
                    (
                        field.transform.to_human_string(
                            partition_field_type, from_bytes(partition_field_type, field_summary.lower_bound)
                        )
                    )
                    if field_summary.lower_bound
                    else None
                )
                upper_bound = (
                    (
                        field.transform.to_human_string(
                            partition_field_type, from_bytes(partition_field_type, field_summary.upper_bound)
                        )
                    )
                    if field_summary.upper_bound
                    else None
                )
                rows.append(
                    {
                        "contains_null": field_summary.contains_null,
                        "contains_nan": field_summary.contains_nan,
                        "lower_bound": lower_bound,
                        "upper_bound": upper_bound,
                    }
                )
            return rows

        specs = self.tbl.metadata.specs()
        manifests = []
        if snapshot:
            for manifest in snapshot.manifests(self.tbl.io):
                is_data_file = manifest.content == ManifestContent.DATA
                is_delete_file = manifest.content == ManifestContent.DELETES
                manifest_row = {
                    "content": manifest.content,
                    "path": manifest.manifest_path,
                    "length": manifest.manifest_length,
                    "partition_spec_id": manifest.partition_spec_id,
                    "added_snapshot_id": manifest.added_snapshot_id,
                    "added_data_files_count": manifest.added_files_count if is_data_file else 0,
                    "existing_data_files_count": manifest.existing_files_count if is_data_file else 0,
                    "deleted_data_files_count": manifest.deleted_files_count if is_data_file else 0,
                    "added_delete_files_count": manifest.added_files_count if is_delete_file else 0,
                    "existing_delete_files_count": manifest.existing_files_count if is_delete_file else 0,
                    "deleted_delete_files_count": manifest.deleted_files_count if is_delete_file else 0,
                    "partition_summaries": _partition_summaries_to_rows(specs[manifest.partition_spec_id], manifest.partitions)
                    if manifest.partitions
                    else [],
                }
                if is_all_manifests_table:
                    manifest_row["reference_snapshot_id"] = snapshot.snapshot_id
                manifests.append(manifest_row)

        return pa.Table.from_pylist(
            manifests,
            schema=self._get_all_manifests_schema() if is_all_manifests_table else self._get_manifests_schema(),
        )

    def manifests(self) -> "pa.Table":
        return self._generate_manifests_table(self.tbl.current_snapshot())

    def metadata_log_entries(self) -> "pa.Table":
        import pyarrow as pa

        from pyiceberg.table.snapshots import MetadataLogEntry

        table_schema = pa.schema(
            [
                pa.field("timestamp", pa.timestamp(unit="ms"), nullable=False),
                pa.field("file", pa.string(), nullable=False),
                pa.field("latest_snapshot_id", pa.int64(), nullable=True),
                pa.field("latest_schema_id", pa.int32(), nullable=True),
                pa.field("latest_sequence_number", pa.int64(), nullable=True),
            ]
        )

        def metadata_log_entry_to_row(metadata_entry: MetadataLogEntry) -> Dict[str, Any]:
            latest_snapshot = self.tbl.snapshot_as_of_timestamp(metadata_entry.timestamp_ms)
            return {
                "timestamp": metadata_entry.timestamp_ms,
                "file": metadata_entry.metadata_file,
                "latest_snapshot_id": latest_snapshot.snapshot_id if latest_snapshot else None,
                "latest_schema_id": latest_snapshot.schema_id if latest_snapshot else None,
                "latest_sequence_number": latest_snapshot.sequence_number if latest_snapshot else None,
            }

        # similar to MetadataLogEntriesTable in Java
        # https://github.com/apache/iceberg/blob/8a70fe0ff5f241aec8856f8091c77fdce35ad256/core/src/main/java/org/apache/iceberg/MetadataLogEntriesTable.java#L62-L66
        metadata_log_entries = self.tbl.metadata.metadata_log + [
            MetadataLogEntry(metadata_file=self.tbl.metadata_location, timestamp_ms=self.tbl.metadata.last_updated_ms)
        ]

        return pa.Table.from_pylist(
            [metadata_log_entry_to_row(entry) for entry in metadata_log_entries],
            schema=table_schema,
        )

    def history(self) -> "pa.Table":
        import pyarrow as pa

        history_schema = pa.schema(
            [
                pa.field("made_current_at", pa.timestamp(unit="ms"), nullable=False),
                pa.field("snapshot_id", pa.int64(), nullable=False),
                pa.field("parent_id", pa.int64(), nullable=True),
                pa.field("is_current_ancestor", pa.bool_(), nullable=False),
            ]
        )

        ancestors_ids = {snapshot.snapshot_id for snapshot in ancestors_of(self.tbl.current_snapshot(), self.tbl.metadata)}

        history = []
        metadata = self.tbl.metadata

        for snapshot_entry in metadata.snapshot_log:
            snapshot = metadata.snapshot_by_id(snapshot_entry.snapshot_id)

            history.append(
                {
                    "made_current_at": datetime.fromtimestamp(snapshot_entry.timestamp_ms / 1000.0, tz=timezone.utc),
                    "snapshot_id": snapshot_entry.snapshot_id,
                    "parent_id": snapshot.parent_snapshot_id if snapshot else None,
                    "is_current_ancestor": snapshot_entry.snapshot_id in ancestors_ids,
                }
            )

        return pa.Table.from_pylist(history, schema=history_schema)

    def _get_files_from_manifest(
        self, manifest_list: ManifestFile, data_file_filter: Optional[Set[DataFileContent]] = None
    ) -> "pa.Table":
        import pyarrow as pa

        files: list[dict[str, Any]] = []
        schema = self.tbl.metadata.schema()
        io = self.tbl.io

        for manifest_entry in manifest_list.fetch_manifest_entry(io):
            data_file = manifest_entry.data_file
            if data_file_filter and data_file.content not in data_file_filter:
                continue
            column_sizes = data_file.column_sizes or {}
            value_counts = data_file.value_counts or {}
            null_value_counts = data_file.null_value_counts or {}
            nan_value_counts = data_file.nan_value_counts or {}
            lower_bounds = data_file.lower_bounds or {}
            upper_bounds = data_file.upper_bounds or {}
            readable_metrics = {
                schema.find_column_name(field.field_id): {
                    "column_size": column_sizes.get(field.field_id),
                    "value_count": value_counts.get(field.field_id),
                    "null_value_count": null_value_counts.get(field.field_id),
                    "nan_value_count": nan_value_counts.get(field.field_id),
                    "lower_bound": from_bytes(field.field_type, lower_bound)
                    if (lower_bound := lower_bounds.get(field.field_id))
                    else None,
                    "upper_bound": from_bytes(field.field_type, upper_bound)
                    if (upper_bound := upper_bounds.get(field.field_id))
                    else None,
                }
                for field in self.tbl.metadata.schema().fields
            }
            partition = data_file.partition
            partition_record_dict = {
                field.name: partition[pos]
                for pos, field in enumerate(self.tbl.metadata.specs()[manifest_list.partition_spec_id].fields)
            }
            files.append(
                {
                    "content": data_file.content,
                    "file_path": data_file.file_path,
                    "file_format": data_file.file_format,
                    "spec_id": data_file.spec_id,
                    "partition": partition_record_dict,
                    "record_count": data_file.record_count,
                    "file_size_in_bytes": data_file.file_size_in_bytes,
                    "column_sizes": dict(data_file.column_sizes) if data_file.column_sizes is not None else None,
                    "value_counts": dict(data_file.value_counts) if data_file.value_counts is not None else None,
                    "null_value_counts": dict(data_file.null_value_counts) if data_file.null_value_counts is not None else None,
                    "nan_value_counts": dict(data_file.nan_value_counts) if data_file.nan_value_counts is not None else None,
                    "lower_bounds": dict(data_file.lower_bounds) if data_file.lower_bounds is not None else None,
                    "upper_bounds": dict(data_file.upper_bounds) if data_file.upper_bounds is not None else None,
                    "key_metadata": data_file.key_metadata,
                    "split_offsets": data_file.split_offsets,
                    "equality_ids": data_file.equality_ids,
                    "sort_order_id": data_file.sort_order_id,
                    "readable_metrics": readable_metrics,
                }
            )
        return pa.Table.from_pylist(
            files,
            schema=self._get_files_schema(),
        )

    def _get_files_schema(self) -> "pa.Schema":
        import pyarrow as pa

        from pyiceberg.io.pyarrow import schema_to_pyarrow

        schema = self.tbl.metadata.schema()
        readable_metrics_struct = []

        def _readable_metrics_struct(bound_type: PrimitiveType) -> pa.StructType:
            pa_bound_type = schema_to_pyarrow(bound_type)
            return pa.struct(
                [
                    pa.field("column_size", pa.int64(), nullable=True),
                    pa.field("value_count", pa.int64(), nullable=True),
                    pa.field("null_value_count", pa.int64(), nullable=True),
                    pa.field("nan_value_count", pa.int64(), nullable=True),
                    pa.field("lower_bound", pa_bound_type, nullable=True),
                    pa.field("upper_bound", pa_bound_type, nullable=True),
                ]
            )

        partition_record = self.tbl.metadata.specs_struct()
        pa_record_struct = schema_to_pyarrow(partition_record)

        for field in self.tbl.metadata.schema().fields:
            readable_metrics_struct.append(
                pa.field(schema.find_column_name(field.field_id), _readable_metrics_struct(field.field_type), nullable=False)
            )

        files_schema = pa.schema(
            [
                pa.field("content", pa.int8(), nullable=False),
                pa.field("file_path", pa.string(), nullable=False),
                pa.field("file_format", pa.dictionary(pa.int32(), pa.string()), nullable=False),
                pa.field("spec_id", pa.int32(), nullable=False),
                pa.field("partition", pa_record_struct, nullable=False),
                pa.field("record_count", pa.int64(), nullable=False),
                pa.field("file_size_in_bytes", pa.int64(), nullable=False),
                pa.field("column_sizes", pa.map_(pa.int32(), pa.int64()), nullable=True),
                pa.field("value_counts", pa.map_(pa.int32(), pa.int64()), nullable=True),
                pa.field("null_value_counts", pa.map_(pa.int32(), pa.int64()), nullable=True),
                pa.field("nan_value_counts", pa.map_(pa.int32(), pa.int64()), nullable=True),
                pa.field("lower_bounds", pa.map_(pa.int32(), pa.binary()), nullable=True),
                pa.field("upper_bounds", pa.map_(pa.int32(), pa.binary()), nullable=True),
                pa.field("key_metadata", pa.binary(), nullable=True),
                pa.field("split_offsets", pa.list_(pa.int64()), nullable=True),
                pa.field("equality_ids", pa.list_(pa.int32()), nullable=True),
                pa.field("sort_order_id", pa.int32(), nullable=True),
                pa.field("readable_metrics", pa.struct(readable_metrics_struct), nullable=True),
            ]
        )
        return files_schema

    def _files(self, snapshot_id: Optional[int] = None, data_file_filter: Optional[Set[DataFileContent]] = None) -> "pa.Table":
        import pyarrow as pa

        if not snapshot_id and not self.tbl.metadata.current_snapshot():
            return self._get_files_schema().empty_table()

        snapshot = self._get_snapshot(snapshot_id)
        io = self.tbl.io

        executor = ExecutorFactory.get_or_create()
        results = list(
            executor.map(
                lambda manifest_list: self._get_files_from_manifest(manifest_list, data_file_filter), snapshot.manifests(io)
            )
        )
        return pa.concat_tables(results)

    def files(self, snapshot_id: Optional[int] = None) -> "pa.Table":
        return self._files(snapshot_id)

    def data_files(self, snapshot_id: Optional[int] = None) -> "pa.Table":
        return self._files(snapshot_id, {DataFileContent.DATA})

    def delete_files(self, snapshot_id: Optional[int] = None) -> "pa.Table":
        return self._files(snapshot_id, {DataFileContent.POSITION_DELETES, DataFileContent.EQUALITY_DELETES})

    def all_manifests(self) -> "pa.Table":
        import pyarrow as pa

        snapshots = self.tbl.snapshots()
        if not snapshots:
            return pa.Table.from_pylist([], schema=self._get_all_manifests_schema())

        executor = ExecutorFactory.get_or_create()
        manifests_by_snapshots: Iterator["pa.Table"] = executor.map(
            lambda args: self._generate_manifests_table(*args), [(snapshot, True) for snapshot in snapshots]
        )
        return pa.concat_tables(manifests_by_snapshots)

    def _all_files(self, data_file_filter: Optional[Set[DataFileContent]] = None) -> "pa.Table":
        import pyarrow as pa

        snapshots = self.tbl.snapshots()
        if not snapshots:
            return pa.Table.from_pylist([], schema=self._get_files_schema())

        executor = ExecutorFactory.get_or_create()
        manifest_lists = executor.map(lambda snapshot: snapshot.manifests(self.tbl.io), snapshots)

        unique_manifests = {(manifest.manifest_path, manifest) for manifest_list in manifest_lists for manifest in manifest_list}

        file_lists = executor.map(
            lambda args: self._get_files_from_manifest(*args), [(manifest, data_file_filter) for _, manifest in unique_manifests]
        )

        return pa.concat_tables(file_lists)

    def all_files(self) -> "pa.Table":
        return self._all_files()

    def all_data_files(self) -> "pa.Table":
        return self._all_files({DataFileContent.DATA})

    def all_delete_files(self) -> "pa.Table":
        return self._all_files({DataFileContent.POSITION_DELETES, DataFileContent.EQUALITY_DELETES})
