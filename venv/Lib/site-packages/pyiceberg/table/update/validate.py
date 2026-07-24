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
from collections.abc import Iterator

from pyiceberg.exceptions import ValidationException
from pyiceberg.expressions import BooleanExpression
from pyiceberg.expressions.visitors import ROWS_CANNOT_MATCH, _InclusiveMetricsEvaluator
from pyiceberg.manifest import ManifestContent, ManifestEntry, ManifestEntryStatus, ManifestFile
from pyiceberg.schema import Schema
from pyiceberg.table import Table
from pyiceberg.table.snapshots import Operation, Snapshot, ancestors_between
from pyiceberg.typedef import Record

VALIDATE_DATA_FILES_EXIST_OPERATIONS: set[Operation] = {Operation.OVERWRITE, Operation.REPLACE, Operation.DELETE}
VALIDATE_ADDED_DATA_FILES_OPERATIONS: set[Operation] = {Operation.APPEND, Operation.OVERWRITE}


def _validation_history(
    table: Table,
    from_snapshot: Snapshot,
    to_snapshot: Snapshot,
    matching_operations: set[Operation],
    manifest_content_filter: ManifestContent,
) -> tuple[list[ManifestFile], set[int]]:
    """Return newly added manifests and snapshot IDs between the starting snapshot and parent snapshot.

    Args:
        table: Table to get the history from
        from_snapshot: Parent snapshot to get the history from
        to_snapshot: Starting snapshot
        matching_operations: Operations to match on
        manifest_content_filter: Manifest content type to filter

    Raises:
        ValidationException: If no matching snapshot is found or only one snapshot is found

    Returns:
        List of manifest files and set of snapshots ID's matching conditions
    """
    manifests_files: list[ManifestFile] = []
    snapshots: set[int] = set()

    last_snapshot = None
    for snapshot in ancestors_between(from_snapshot, to_snapshot, table.metadata):
        last_snapshot = snapshot
        summary = snapshot.summary
        if summary is None:
            raise ValidationException(f"No summary found for snapshot {snapshot}!")
        if summary.operation not in matching_operations:
            continue

        snapshots.add(snapshot.snapshot_id)
        # TODO: Maybe do the IO in a separate thread at some point, and collect at the bottom (we can easily merge the sets
        manifests_files.extend(
            [
                manifest
                for manifest in snapshot.manifests(table.io)
                if manifest.added_snapshot_id == snapshot.snapshot_id and manifest.content == manifest_content_filter
            ]
        )

    if last_snapshot is not None and last_snapshot.snapshot_id != from_snapshot.snapshot_id:
        raise ValidationException("No matching snapshot found.")

    return manifests_files, snapshots


def _filter_manifest_entries(
    entry: ManifestEntry,
    snapshot_ids: set[int],
    data_filter: BooleanExpression | None,
    partition_set: dict[int, set[Record]] | None,
    entry_status: ManifestEntryStatus | None,
    schema: Schema,
) -> bool:
    """Filter manifest entries based on data filter and partition set.

    Args:
        entry: Manifest entry to filter
        snapshot_ids: set of snapshot ids to match data files
        data_filter: Optional filter to match data files
        partition_set: Optional set of partitions to match data files
        entry_status: Optional status to match data files
        schema: schema for filtering

    Returns:
        True if the entry should be included, False otherwise
    """
    if entry.snapshot_id not in snapshot_ids:
        return False

    if entry_status is not None and entry.status != entry_status:
        return False

    if data_filter is not None:
        evaluator = _InclusiveMetricsEvaluator(schema, data_filter)
        if evaluator.eval(entry.data_file) is ROWS_CANNOT_MATCH:
            return False

    if partition_set is not None:
        partition = entry.data_file.partition
        spec_id = entry.data_file.spec_id
        if spec_id not in partition_set or partition not in partition_set[spec_id]:
            return False

    return True


def _deleted_data_files(
    table: Table,
    starting_snapshot: Snapshot,
    data_filter: BooleanExpression | None,
    partition_set: dict[int, set[Record]] | None,
    parent_snapshot: Snapshot | None,
) -> Iterator[ManifestEntry]:
    """Find deleted data files matching a filter since a starting snapshot.

    Args:
        table: Table to validate
        starting_snapshot: Snapshot current at the start of the operation
        data_filter: Expression used to find deleted data files
        partition_set: dict of {spec_id: set[partition]} to filter on
        parent_snapshot: Ending snapshot on the branch being validated

    Returns:
        List of conflicting manifest-entries
    """
    # if there is no current table state, no files have been deleted
    if parent_snapshot is None:
        return

    manifests, snapshot_ids = _validation_history(
        table,
        parent_snapshot,
        starting_snapshot,
        VALIDATE_DATA_FILES_EXIST_OPERATIONS,
        ManifestContent.DATA,
    )

    for manifest in manifests:
        for entry in manifest.fetch_manifest_entry(table.io, discard_deleted=False):
            if _filter_manifest_entries(
                entry, snapshot_ids, data_filter, partition_set, ManifestEntryStatus.DELETED, table.schema()
            ):
                yield entry


def _validate_deleted_data_files(
    table: Table,
    starting_snapshot: Snapshot,
    data_filter: BooleanExpression | None,
    parent_snapshot: Snapshot,
) -> None:
    """Validate that no files matching a filter have been deleted from the table since a starting snapshot.

    Args:
        table: Table to validate
        starting_snapshot: Snapshot current at the start of the operation
        data_filter: Expression used to find deleted data files
        parent_snapshot: Ending snapshot on the branch being validated

    """
    conflicting_entries = _deleted_data_files(table, starting_snapshot, data_filter, None, parent_snapshot)
    if any(conflicting_entries):
        conflicting_snapshots = {entry.snapshot_id for entry in conflicting_entries}
        raise ValidationException(f"Deleted data files were found matching the filter for snapshots {conflicting_snapshots}!")


def _added_data_files(
    table: Table,
    starting_snapshot: Snapshot,
    data_filter: BooleanExpression | None,
    partition_set: dict[int, set[Record]] | None,
    parent_snapshot: Snapshot | None,
) -> Iterator[ManifestEntry]:
    """Return manifest entries for data files added between the starting snapshot and parent snapshot.

    Args:
        table: Table to get the history from
        starting_snapshot: Starting snapshot to get the history from
        data_filter: Optional filter to match data files
        partition_set: Optional set of partitions to match data files
        parent_snapshot: Parent snapshot to get the history from

    Returns:
        Iterator of manifest entries for added data files matching the conditions
    """
    if parent_snapshot is None:
        return

    manifests, snapshot_ids = _validation_history(
        table,
        parent_snapshot,
        starting_snapshot,
        VALIDATE_ADDED_DATA_FILES_OPERATIONS,
        ManifestContent.DATA,
    )

    for manifest in manifests:
        for entry in manifest.fetch_manifest_entry(table.io):
            if _filter_manifest_entries(entry, snapshot_ids, data_filter, partition_set, None, table.schema()):
                yield entry


def _validate_added_data_files(
    table: Table,
    starting_snapshot: Snapshot,
    data_filter: BooleanExpression | None,
    parent_snapshot: Snapshot | None,
) -> None:
    """Validate that no files matching a filter have been added to the table since a starting snapshot.

    Args:
        table: Table to validate
        starting_snapshot: Snapshot current at the start of the operation
        data_filter: Expression used to find added data files
        parent_snapshot: Ending snapshot on the branch being validated

    """
    conflicting_entries = _added_data_files(table, starting_snapshot, data_filter, None, parent_snapshot)
    if any(conflicting_entries):
        conflicting_snapshots = {entry.snapshot_id for entry in conflicting_entries if entry.snapshot_id is not None}
        raise ValidationException(f"Added data files were found matching the filter for snapshots {conflicting_snapshots}!")
