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
from abc import abstractmethod
from collections import defaultdict
from collections.abc import Callable
from datetime import datetime
from functools import cached_property
from typing import TYPE_CHECKING, Generic

from pyiceberg.avro.codecs import AvroCompressionCodec
from pyiceberg.expressions import (
    AlwaysFalse,
    BooleanExpression,
    Or,
)
from pyiceberg.expressions.visitors import (
    ROWS_MIGHT_NOT_MATCH,
    ROWS_MUST_MATCH,
    _InclusiveMetricsEvaluator,
    _StrictMetricsEvaluator,
    inclusive_projection,
    manifest_evaluator,
)
from pyiceberg.io import FileIO, OutputFile
from pyiceberg.manifest import (
    DataFile,
    DataFileContent,
    ManifestContent,
    ManifestEntry,
    ManifestEntryStatus,
    ManifestFile,
    ManifestWriter,
    write_manifest,
    write_manifest_list,
)
from pyiceberg.partitioning import (
    PartitionSpec,
)
from pyiceberg.table.refs import MAIN_BRANCH, SnapshotRefType
from pyiceberg.table.snapshots import (
    Operation,
    Snapshot,
    SnapshotSummaryCollector,
    Summary,
    ancestors_of,
    latest_ancestor_before_timestamp,
    update_snapshot_summaries,
)
from pyiceberg.table.update import (
    AddSnapshotUpdate,
    AssertRefSnapshotId,
    RemoveSnapshotRefUpdate,
    RemoveSnapshotsUpdate,
    SetSnapshotRefUpdate,
    TableRequirement,
    TableUpdate,
    U,
    UpdatesAndRequirements,
    UpdateTableMetadata,
)
from pyiceberg.typedef import (
    EMPTY_DICT,
    KeyDefaultDict,
)
from pyiceberg.utils.bin_packing import ListPacker
from pyiceberg.utils.concurrent import ExecutorFactory
from pyiceberg.utils.datetime import datetime_to_millis
from pyiceberg.utils.properties import property_as_bool, property_as_int

if TYPE_CHECKING:
    from pyiceberg.table import Transaction


def _new_manifest_file_name(num: int, commit_uuid: uuid.UUID) -> str:
    return f"{commit_uuid}-m{num}.avro"


def _new_manifest_list_file_name(snapshot_id: int, attempt: int, commit_uuid: uuid.UUID) -> str:
    # Mimics the behavior in Java:
    # https://github.com/apache/iceberg/blob/c862b9177af8e2d83122220764a056f3b96fd00c/core/src/main/java/org/apache/iceberg/SnapshotProducer.java#L491
    return f"snap-{snapshot_id}-{attempt}-{commit_uuid}.avro"


class _SnapshotProducer(UpdateTableMetadata[U], Generic[U]):
    commit_uuid: uuid.UUID
    _io: FileIO
    _operation: Operation
    _snapshot_id: int
    _parent_snapshot_id: int | None
    _added_data_files: list[DataFile]
    _manifest_num_counter: itertools.count[int]
    _deleted_data_files: set[DataFile]
    _compression: AvroCompressionCodec
    _target_branch: str | None

    def __init__(
        self,
        operation: Operation,
        transaction: Transaction,
        io: FileIO,
        commit_uuid: uuid.UUID | None = None,
        snapshot_properties: dict[str, str] = EMPTY_DICT,
        branch: str | None = MAIN_BRANCH,
    ) -> None:
        super().__init__(transaction)
        self.commit_uuid = commit_uuid or uuid.uuid4()
        self._io = io
        self._operation = operation
        self._snapshot_id = self._transaction.table_metadata.new_snapshot_id()
        self._added_data_files = []
        self._deleted_data_files = set()
        self.snapshot_properties = snapshot_properties
        self._manifest_num_counter = itertools.count(0)
        from pyiceberg.table import TableProperties

        self._compression = self._transaction.table_metadata.properties.get(  # type: ignore
            TableProperties.WRITE_AVRO_COMPRESSION, TableProperties.WRITE_AVRO_COMPRESSION_DEFAULT
        )
        self._target_branch = self._validate_target_branch(branch=branch)
        self._parent_snapshot_id = (
            snapshot.snapshot_id if (snapshot := self._transaction.table_metadata.snapshot_by_name(self._target_branch)) else None
        )

    def _validate_target_branch(self, branch: str | None) -> str | None:
        # if branch is none, write will be written into a staging snapshot
        if branch is not None:
            if branch in self._transaction.table_metadata.refs:
                ref = self._transaction.table_metadata.refs[branch]
                if ref.snapshot_ref_type != SnapshotRefType.BRANCH:
                    raise ValueError(f"{branch} is a tag, not a branch. Tags cannot be targets for producing snapshots")
        return branch

    def append_data_file(self, data_file: DataFile) -> _SnapshotProducer[U]:
        self._added_data_files.append(data_file)
        return self

    def delete_data_file(self, data_file: DataFile) -> _SnapshotProducer[U]:
        self._deleted_data_files.add(data_file)
        return self

    def _calculate_added_rows(self, manifests: list[ManifestFile]) -> int:
        """Calculate the number of added rows from a list of manifest files."""
        added_rows = 0
        for manifest in manifests:
            if manifest.added_snapshot_id is None or manifest.added_snapshot_id == self._snapshot_id:
                if manifest.added_rows_count is None:
                    raise ValueError(
                        "Cannot determine number of added rows in snapshot because "
                        f"the entry for manifest {manifest.manifest_path} is missing the field `added-rows-count`"
                    )
                added_rows += manifest.added_rows_count
        return added_rows

    @abstractmethod
    def _deleted_entries(self) -> list[ManifestEntry]: ...

    @abstractmethod
    def _existing_manifests(self) -> list[ManifestFile]: ...

    def _process_manifests(self, manifests: list[ManifestFile]) -> list[ManifestFile]:
        """To perform any post-processing on the manifests before writing them to the new snapshot."""
        return manifests

    def _manifests(self) -> list[ManifestFile]:
        def _write_added_manifest() -> list[ManifestFile]:
            if self._added_data_files:
                with write_manifest(
                    format_version=self._transaction.table_metadata.format_version,
                    spec=self._transaction.table_metadata.spec(),
                    schema=self._transaction.table_metadata.schema(),
                    output_file=self.new_manifest_output(),
                    snapshot_id=self._snapshot_id,
                    avro_compression=self._compression,
                ) as writer:
                    for data_file in self._added_data_files:
                        writer.add(
                            ManifestEntry.from_args(
                                status=ManifestEntryStatus.ADDED,
                                snapshot_id=self._snapshot_id,
                                sequence_number=None,
                                file_sequence_number=None,
                                data_file=data_file,
                            )
                        )
                return [writer.to_manifest_file()]
            else:
                return []

        def _write_delete_manifest() -> list[ManifestFile]:
            # Check if we need to mark the files as deleted
            deleted_entries = self._deleted_entries()
            if len(deleted_entries) > 0:
                deleted_manifests = []
                partition_groups: dict[int, list[ManifestEntry]] = defaultdict(list)
                for deleted_entry in deleted_entries:
                    partition_groups[deleted_entry.data_file.spec_id].append(deleted_entry)
                for spec_id, entries in partition_groups.items():
                    with write_manifest(
                        format_version=self._transaction.table_metadata.format_version,
                        spec=self._transaction.table_metadata.specs()[spec_id],
                        schema=self._transaction.table_metadata.schema(),
                        output_file=self.new_manifest_output(),
                        snapshot_id=self._snapshot_id,
                        avro_compression=self._compression,
                    ) as writer:
                        for entry in entries:
                            writer.add_entry(entry)
                    deleted_manifests.append(writer.to_manifest_file())
                return deleted_manifests
            else:
                return []

        executor = ExecutorFactory.get_or_create()

        added_manifests = executor.submit(_write_added_manifest)
        delete_manifests = executor.submit(_write_delete_manifest)
        existing_manifests = executor.submit(self._existing_manifests)

        return self._process_manifests(added_manifests.result() + delete_manifests.result() + existing_manifests.result())

    def _summary(self, snapshot_properties: dict[str, str] = EMPTY_DICT) -> Summary:
        from pyiceberg.table import TableProperties

        # avoid copying metadata for each data file
        table_metadata = self._transaction.table_metadata

        partition_summary_limit = int(
            table_metadata.properties.get(
                TableProperties.WRITE_PARTITION_SUMMARY_LIMIT, TableProperties.WRITE_PARTITION_SUMMARY_LIMIT_DEFAULT
            )
        )
        ssc = SnapshotSummaryCollector(partition_summary_limit=partition_summary_limit)

        for data_file in self._added_data_files:
            ssc.add_file(
                data_file=data_file,
                partition_spec=table_metadata.spec(),
                schema=table_metadata.schema(),
            )

        if len(self._deleted_data_files) > 0:
            specs = table_metadata.specs()
            for data_file in self._deleted_data_files:
                ssc.remove_file(
                    data_file=data_file,
                    partition_spec=specs[data_file.spec_id],
                    schema=table_metadata.schema(),
                )

        previous_snapshot = (
            table_metadata.snapshot_by_id(self._parent_snapshot_id) if self._parent_snapshot_id is not None else None
        )

        return update_snapshot_summaries(
            summary=Summary(operation=self._operation, **ssc.build(), **snapshot_properties),
            previous_summary=previous_snapshot.summary if previous_snapshot is not None else None,
        )

    def _commit(self) -> UpdatesAndRequirements:
        new_manifests = self._manifests()
        next_sequence_number = self._transaction.table_metadata.next_sequence_number()

        summary = self._summary(self.snapshot_properties)
        file_name = _new_manifest_list_file_name(
            snapshot_id=self._snapshot_id,
            attempt=0,
            commit_uuid=self.commit_uuid,
        )
        location_provider = self._transaction._table.location_provider()
        manifest_list_file_path = location_provider.new_metadata_location(file_name)

        with write_manifest_list(
            format_version=self._transaction.table_metadata.format_version,
            output_file=self._io.new_output(manifest_list_file_path),
            snapshot_id=self._snapshot_id,
            parent_snapshot_id=self._parent_snapshot_id,
            sequence_number=next_sequence_number,
            avro_compression=self._compression,
        ) as writer:
            writer.add_manifests(new_manifests)

        first_row_id: int | None = None

        if self._transaction.table_metadata.format_version >= 3:
            first_row_id = self._transaction.table_metadata.next_row_id

        snapshot = Snapshot(
            snapshot_id=self._snapshot_id,
            parent_snapshot_id=self._parent_snapshot_id,
            manifest_list=manifest_list_file_path,
            sequence_number=next_sequence_number,
            summary=summary,
            schema_id=self._transaction.table_metadata.current_schema_id,
            first_row_id=first_row_id,
        )

        add_snapshot_update = AddSnapshotUpdate(snapshot=snapshot)

        if self._target_branch is None:
            return (
                (add_snapshot_update,),
                (),
            )
        else:
            return (
                (
                    add_snapshot_update,
                    SetSnapshotRefUpdate(
                        snapshot_id=self._snapshot_id,
                        parent_snapshot_id=self._parent_snapshot_id,
                        ref_name=self._target_branch,
                        type=SnapshotRefType.BRANCH,
                    ),
                ),
                (
                    AssertRefSnapshotId(
                        snapshot_id=self._transaction.table_metadata.refs[self._target_branch].snapshot_id
                        if self._target_branch in self._transaction.table_metadata.refs
                        else None,
                        ref=self._target_branch,
                    ),
                ),
            )

    @property
    def snapshot_id(self) -> int:
        return self._snapshot_id

    def spec(self, spec_id: int) -> PartitionSpec:
        return self._transaction.table_metadata.specs()[spec_id]

    def new_manifest_writer(self, spec: PartitionSpec) -> ManifestWriter:
        return write_manifest(
            format_version=self._transaction.table_metadata.format_version,
            spec=spec,
            schema=self._transaction.table_metadata.schema(),
            output_file=self.new_manifest_output(),
            snapshot_id=self._snapshot_id,
            avro_compression=self._compression,
        )

    def new_manifest_output(self) -> OutputFile:
        location_provider = self._transaction._table.location_provider()
        file_name = _new_manifest_file_name(num=next(self._manifest_num_counter), commit_uuid=self.commit_uuid)
        file_path = location_provider.new_metadata_location(file_name)
        return self._io.new_output(file_path)

    def fetch_manifest_entry(self, manifest: ManifestFile, discard_deleted: bool = True) -> list[ManifestEntry]:
        return manifest.fetch_manifest_entry(io=self._io, discard_deleted=discard_deleted)


class _DeleteFiles(_SnapshotProducer["_DeleteFiles"]):
    """Will delete manifest entries from the current snapshot based on the predicate.

    This will produce a DELETE snapshot:
        Data files were removed and their contents logically deleted and/or delete
        files were added to delete rows.

    From the specification
    """

    _predicate: BooleanExpression
    _case_sensitive: bool

    def __init__(
        self,
        operation: Operation,
        transaction: Transaction,
        io: FileIO,
        branch: str | None = MAIN_BRANCH,
        commit_uuid: uuid.UUID | None = None,
        snapshot_properties: dict[str, str] = EMPTY_DICT,
    ):
        super().__init__(operation, transaction, io, commit_uuid, snapshot_properties, branch)
        self._predicate = AlwaysFalse()
        self._case_sensitive = True

    def _commit(self) -> UpdatesAndRequirements:
        # Only produce a commit when there is something to delete
        if self.files_affected:
            return super()._commit()
        else:
            return (), ()

    def _build_partition_projection(self, spec_id: int) -> BooleanExpression:
        schema = self._transaction.table_metadata.schema()
        spec = self._transaction.table_metadata.specs()[spec_id]
        project = inclusive_projection(schema, spec, self._case_sensitive)
        return project(self._predicate)

    @cached_property
    def partition_filters(self) -> KeyDefaultDict[int, BooleanExpression]:
        return KeyDefaultDict(self._build_partition_projection)

    def _build_manifest_evaluator(self, spec_id: int) -> Callable[[ManifestFile], bool]:
        schema = self._transaction.table_metadata.schema()
        spec = self._transaction.table_metadata.specs()[spec_id]
        return manifest_evaluator(spec, schema, self.partition_filters[spec_id], self._case_sensitive)

    def delete_by_predicate(self, predicate: BooleanExpression, case_sensitive: bool = True) -> None:
        self._predicate = Or(self._predicate, predicate)
        self._case_sensitive = case_sensitive

    @cached_property
    def _compute_deletes(self) -> tuple[list[ManifestFile], list[ManifestEntry], bool]:
        """Computes all the delete operation and cache it when nothing changes.

        Returns:
            - List of existing manifests that are not affected by the delete operation.
            - The manifest-entries that are deleted based on the metadata.
            - Flag indicating that rewrites of data-files are needed.
        """
        schema = self._transaction.table_metadata.schema()

        def _copy_with_new_status(entry: ManifestEntry, status: ManifestEntryStatus) -> ManifestEntry:
            return ManifestEntry.from_args(
                status=status,
                # When a file is replaced or deleted from the dataset, its manifest entry fields store the
                # snapshot ID in which the file was deleted and status 2 (deleted).
                snapshot_id=self.snapshot_id if status == ManifestEntryStatus.DELETED else entry.snapshot_id,
                sequence_number=entry.sequence_number,
                file_sequence_number=entry.file_sequence_number,
                data_file=entry.data_file,
            )

        manifest_evaluators: dict[int, Callable[[ManifestFile], bool]] = KeyDefaultDict(self._build_manifest_evaluator)
        strict_metrics_evaluator = _StrictMetricsEvaluator(schema, self._predicate, case_sensitive=self._case_sensitive).eval
        inclusive_metrics_evaluator = _InclusiveMetricsEvaluator(
            schema, self._predicate, case_sensitive=self._case_sensitive
        ).eval

        existing_manifests = []
        total_deleted_entries = []
        partial_rewrites_needed = False
        self._deleted_data_files = set()

        # Determine the snapshot to read manifests from for deletion
        # Should be the current tip of the _target_branch
        parent_snapshot_id_for_delete_source = self._parent_snapshot_id
        if parent_snapshot_id_for_delete_source is not None:
            snapshot = self._transaction.table_metadata.snapshot_by_id(parent_snapshot_id_for_delete_source)
            if snapshot:  # Ensure snapshot is found
                for manifest_file in snapshot.manifests(io=self._io):
                    if manifest_file.content == ManifestContent.DATA:
                        if not manifest_evaluators[manifest_file.partition_spec_id](manifest_file):
                            # If the manifest isn't relevant, we can just keep it in the manifest-list
                            existing_manifests.append(manifest_file)
                        else:
                            # It is relevant, let's check out the content
                            deleted_entries = []
                            existing_entries = []
                            for entry in manifest_file.fetch_manifest_entry(io=self._io, discard_deleted=True):
                                if strict_metrics_evaluator(entry.data_file) == ROWS_MUST_MATCH:
                                    # Based on the metadata, it can be dropped right away
                                    deleted_entries.append(_copy_with_new_status(entry, ManifestEntryStatus.DELETED))
                                    self._deleted_data_files.add(entry.data_file)
                                else:
                                    # Based on the metadata, we cannot determine if it can be deleted
                                    existing_entries.append(_copy_with_new_status(entry, ManifestEntryStatus.EXISTING))
                                    if inclusive_metrics_evaluator(entry.data_file) != ROWS_MIGHT_NOT_MATCH:
                                        partial_rewrites_needed = True

                            if len(deleted_entries) > 0:
                                total_deleted_entries += deleted_entries

                                # Rewrite the manifest
                                if len(existing_entries) > 0:
                                    with write_manifest(
                                        format_version=self._transaction.table_metadata.format_version,
                                        spec=self._transaction.table_metadata.specs()[manifest_file.partition_spec_id],
                                        schema=self._transaction.table_metadata.schema(),
                                        output_file=self.new_manifest_output(),
                                        snapshot_id=self._snapshot_id,
                                        avro_compression=self._compression,
                                    ) as writer:
                                        for existing_entry in existing_entries:
                                            writer.add_entry(existing_entry)
                                    existing_manifests.append(writer.to_manifest_file())
                            else:
                                existing_manifests.append(manifest_file)
                    else:
                        existing_manifests.append(manifest_file)

        return existing_manifests, total_deleted_entries, partial_rewrites_needed

    def _existing_manifests(self) -> list[ManifestFile]:
        return self._compute_deletes[0]

    def _deleted_entries(self) -> list[ManifestEntry]:
        return self._compute_deletes[1]

    @property
    def rewrites_needed(self) -> bool:
        """Indicate if data files need to be rewritten."""
        return self._compute_deletes[2]

    @property
    def files_affected(self) -> bool:
        """Indicate if any manifest-entries can be dropped."""
        return len(self._deleted_entries()) > 0


class _FastAppendFiles(_SnapshotProducer["_FastAppendFiles"]):
    def _existing_manifests(self) -> list[ManifestFile]:
        """To determine if there are any existing manifest files.

        A fast append will add another ManifestFile to the ManifestList.
        All the existing manifest files are considered existing.
        """
        existing_manifests = []

        if self._parent_snapshot_id is not None:
            previous_snapshot = self._transaction.table_metadata.snapshot_by_id(self._parent_snapshot_id)

            if previous_snapshot is None:
                raise ValueError(f"Snapshot could not be found: {self._parent_snapshot_id}")

            for manifest in previous_snapshot.manifests(io=self._io):
                if manifest.has_added_files() or manifest.has_existing_files() or manifest.added_snapshot_id == self._snapshot_id:
                    existing_manifests.append(manifest)

        return existing_manifests

    def _deleted_entries(self) -> list[ManifestEntry]:
        """To determine if we need to record any deleted manifest entries.

        In case of an append, nothing is deleted.
        """
        return []


class _MergeAppendFiles(_FastAppendFiles):
    _target_size_bytes: int
    _min_count_to_merge: int
    _merge_enabled: bool

    def __init__(
        self,
        operation: Operation,
        transaction: Transaction,
        io: FileIO,
        branch: str | None = MAIN_BRANCH,
        commit_uuid: uuid.UUID | None = None,
        snapshot_properties: dict[str, str] = EMPTY_DICT,
    ) -> None:
        from pyiceberg.table import TableProperties

        super().__init__(operation, transaction, io, commit_uuid, snapshot_properties, branch)
        self._target_size_bytes = property_as_int(
            self._transaction.table_metadata.properties,
            TableProperties.MANIFEST_TARGET_SIZE_BYTES,
            TableProperties.MANIFEST_TARGET_SIZE_BYTES_DEFAULT,
        )  # type: ignore
        self._min_count_to_merge = property_as_int(
            self._transaction.table_metadata.properties,
            TableProperties.MANIFEST_MIN_MERGE_COUNT,
            TableProperties.MANIFEST_MIN_MERGE_COUNT_DEFAULT,
        )  # type: ignore
        self._merge_enabled = property_as_bool(
            self._transaction.table_metadata.properties,
            TableProperties.MANIFEST_MERGE_ENABLED,
            TableProperties.MANIFEST_MERGE_ENABLED_DEFAULT,
        )

    def _process_manifests(self, manifests: list[ManifestFile]) -> list[ManifestFile]:
        """To perform any post-processing on the manifests before writing them to the new snapshot.

        In _MergeAppendFiles, we merge manifests based on the target size and the minimum count to merge
        if automatic merge is enabled.
        """
        unmerged_data_manifests = [manifest for manifest in manifests if manifest.content == ManifestContent.DATA]
        unmerged_deletes_manifests = [manifest for manifest in manifests if manifest.content == ManifestContent.DELETES]

        data_manifest_merge_manager = _ManifestMergeManager(
            target_size_bytes=self._target_size_bytes,
            min_count_to_merge=self._min_count_to_merge,
            merge_enabled=self._merge_enabled,
            snapshot_producer=self,
        )

        return data_manifest_merge_manager.merge_manifests(unmerged_data_manifests) + unmerged_deletes_manifests


class _OverwriteFiles(_SnapshotProducer["_OverwriteFiles"]):
    """Overwrites data from the table. This will produce an OVERWRITE snapshot.

    Data and delete files were added and removed in a logical overwrite operation.
    """

    def _existing_manifests(self) -> list[ManifestFile]:
        """Determine if there are any existing manifest files."""
        existing_files = []

        if snapshot := self._transaction.table_metadata.snapshot_by_name(name=self._target_branch):
            for manifest_file in snapshot.manifests(io=self._io):
                entries = manifest_file.fetch_manifest_entry(io=self._io, discard_deleted=True)
                found_deleted_data_files = [entry.data_file for entry in entries if entry.data_file in self._deleted_data_files]

                if len(found_deleted_data_files) == 0:
                    existing_files.append(manifest_file)
                else:
                    # We have to rewrite the manifest file without the deleted data files
                    if any(entry.data_file not in found_deleted_data_files for entry in entries):
                        with write_manifest(
                            format_version=self._transaction.table_metadata.format_version,
                            spec=self._transaction.table_metadata.specs()[manifest_file.partition_spec_id],
                            schema=self._transaction.table_metadata.schema(),
                            output_file=self.new_manifest_output(),
                            snapshot_id=self._snapshot_id,
                            avro_compression=self._compression,
                        ) as writer:
                            for entry in entries:
                                if entry.data_file not in found_deleted_data_files:
                                    writer.add_entry(
                                        ManifestEntry.from_args(
                                            status=ManifestEntryStatus.EXISTING,
                                            snapshot_id=entry.snapshot_id,
                                            sequence_number=entry.sequence_number,
                                            file_sequence_number=entry.file_sequence_number,
                                            data_file=entry.data_file,
                                        )
                                    )
                        existing_files.append(writer.to_manifest_file())
        return existing_files

    def _deleted_entries(self) -> list[ManifestEntry]:
        """To determine if we need to record any deleted entries.

        With a full overwrite all the entries are considered deleted.
        With partial overwrites we have to use the predicate to evaluate
        which entries are affected.
        """
        if self._parent_snapshot_id is not None:
            previous_snapshot = self._transaction.table_metadata.snapshot_by_id(self._parent_snapshot_id)
            if previous_snapshot is None:
                # This should never happen since you cannot overwrite an empty table
                raise ValueError(f"Could not find the previous snapshot: {self._parent_snapshot_id}")

            executor = ExecutorFactory.get_or_create()

            def _get_entries(manifest: ManifestFile) -> list[ManifestEntry]:
                return [
                    ManifestEntry.from_args(
                        status=ManifestEntryStatus.DELETED,
                        snapshot_id=entry.snapshot_id,
                        sequence_number=entry.sequence_number,
                        file_sequence_number=entry.file_sequence_number,
                        data_file=entry.data_file,
                    )
                    for entry in manifest.fetch_manifest_entry(self._io, discard_deleted=True)
                    if entry.data_file.content == DataFileContent.DATA and entry.data_file in self._deleted_data_files
                ]

            list_of_entries = executor.map(_get_entries, previous_snapshot.manifests(self._io))
            return list(itertools.chain(*list_of_entries))
        else:
            return []


class UpdateSnapshot:
    _transaction: Transaction
    _io: FileIO
    _branch: str | None
    _snapshot_properties: dict[str, str]

    def __init__(
        self,
        transaction: Transaction,
        io: FileIO,
        branch: str | None = MAIN_BRANCH,
        snapshot_properties: dict[str, str] = EMPTY_DICT,
    ) -> None:
        self._transaction = transaction
        self._io = io
        self._snapshot_properties = snapshot_properties
        self._branch = branch

    def fast_append(self) -> _FastAppendFiles:
        return _FastAppendFiles(
            operation=Operation.APPEND,
            transaction=self._transaction,
            io=self._io,
            branch=self._branch,
            snapshot_properties=self._snapshot_properties,
        )

    def merge_append(self) -> _MergeAppendFiles:
        return _MergeAppendFiles(
            operation=Operation.APPEND,
            transaction=self._transaction,
            io=self._io,
            branch=self._branch,
            snapshot_properties=self._snapshot_properties,
        )

    def overwrite(self, commit_uuid: uuid.UUID | None = None) -> _OverwriteFiles:
        return _OverwriteFiles(
            commit_uuid=commit_uuid,
            operation=Operation.OVERWRITE
            if self._transaction.table_metadata.snapshot_by_name(name=self._branch) is not None
            else Operation.APPEND,
            transaction=self._transaction,
            io=self._io,
            branch=self._branch,
            snapshot_properties=self._snapshot_properties,
        )

    def delete(self) -> _DeleteFiles:
        return _DeleteFiles(
            operation=Operation.DELETE,
            transaction=self._transaction,
            io=self._io,
            branch=self._branch,
            snapshot_properties=self._snapshot_properties,
        )


class _ManifestMergeManager(Generic[U]):
    _target_size_bytes: int
    _min_count_to_merge: int
    _merge_enabled: bool
    _snapshot_producer: _SnapshotProducer[U]

    def __init__(
        self, target_size_bytes: int, min_count_to_merge: int, merge_enabled: bool, snapshot_producer: _SnapshotProducer[U]
    ) -> None:
        self._target_size_bytes = target_size_bytes
        self._min_count_to_merge = min_count_to_merge
        self._merge_enabled = merge_enabled
        self._snapshot_producer = snapshot_producer

    def _group_by_spec(self, manifests: list[ManifestFile]) -> dict[int, list[ManifestFile]]:
        groups = defaultdict(list)
        for manifest in manifests:
            groups[manifest.partition_spec_id].append(manifest)
        return groups

    def _create_manifest(self, spec_id: int, manifest_bin: list[ManifestFile]) -> ManifestFile:
        with self._snapshot_producer.new_manifest_writer(spec=self._snapshot_producer.spec(spec_id)) as writer:
            for manifest in manifest_bin:
                for entry in self._snapshot_producer.fetch_manifest_entry(manifest=manifest, discard_deleted=False):
                    if entry.status == ManifestEntryStatus.DELETED and entry.snapshot_id == self._snapshot_producer.snapshot_id:
                        #  only files deleted by this snapshot should be added to the new manifest
                        writer.delete(entry)
                    elif entry.status == ManifestEntryStatus.ADDED and entry.snapshot_id == self._snapshot_producer.snapshot_id:
                        # added entries from this snapshot are still added, otherwise they should be existing
                        writer.add(entry)
                    elif entry.status != ManifestEntryStatus.DELETED:
                        # add all non-deleted files from the old manifest as existing files
                        writer.existing(entry)

        return writer.to_manifest_file()

    def _merge_group(self, first_manifest: ManifestFile, spec_id: int, manifests: list[ManifestFile]) -> list[ManifestFile]:
        packer: ListPacker[ManifestFile] = ListPacker(target_weight=self._target_size_bytes, lookback=1, largest_bin_first=False)
        bins: list[list[ManifestFile]] = packer.pack_end(manifests, lambda m: m.manifest_length)

        def merge_bin(manifest_bin: list[ManifestFile]) -> list[ManifestFile]:
            output_manifests = []
            if len(manifest_bin) == 1:
                output_manifests.append(manifest_bin[0])
            elif first_manifest in manifest_bin and len(manifest_bin) < self._min_count_to_merge:
                #  if the bin has the first manifest (the new data files or an appended manifest file) then only
                #  merge it if the number of manifests is above the minimum count. this is applied only to bins
                #  with an in-memory manifest so that large manifests don't prevent merging older groups.
                output_manifests.extend(manifest_bin)
            else:
                output_manifests.append(self._create_manifest(spec_id, manifest_bin))

            return output_manifests

        executor = ExecutorFactory.get_or_create()
        futures = [executor.submit(merge_bin, b) for b in bins]
        bin_results: list[list[ManifestFile]] = [r for f in futures if (r := f.result())]

        return [manifest for bin_result in bin_results for manifest in bin_result]

    def merge_manifests(self, manifests: list[ManifestFile]) -> list[ManifestFile]:
        if not self._merge_enabled or len(manifests) == 0:
            return manifests

        first_manifest = manifests[0]
        groups = self._group_by_spec(manifests)

        merged_manifests = []
        for spec_id in reversed(groups.keys()):
            merged_manifests.extend(self._merge_group(first_manifest, spec_id, groups[spec_id]))

        return merged_manifests


class ManageSnapshots(UpdateTableMetadata["ManageSnapshots"]):
    """
    Run snapshot management operations using APIs.

    APIs include create branch, create tag, etc.

    Use table.manage_snapshots().<operation>().commit() to run a specific operation.
    Use table.manage_snapshots().<operation-one>().<operation-two>().commit() to run multiple operations.
    Pending changes are applied on commit.

    We can also use context managers to make more changes. For example,

    with table.manage_snapshots() as ms:
       ms.create_tag(snapshot_id1, "Tag_A").create_tag(snapshot_id2, "Tag_B")
    """

    _updates: tuple[TableUpdate, ...]
    _requirements: tuple[TableRequirement, ...]

    def __init__(self, transaction: Transaction) -> None:
        super().__init__(transaction)
        self._updates = ()
        self._requirements = ()

    def _commit(self) -> UpdatesAndRequirements:
        """Apply the pending changes and commit."""
        return self._updates, self._requirements

    def _commit_if_ref_updates_exist(self) -> None:
        """Stage any pending ref updates to the transaction state."""
        if self._updates:
            self._transaction._stage(*self._commit())
            self._updates = ()
            self._requirements = ()

    def _remove_ref_snapshot(self, ref_name: str) -> ManageSnapshots:
        """Remove a snapshot ref.

        Args:
            ref_name: branch / tag name to remove
        Stages the updates and requirements for the remove-snapshot-ref.
        Returns
            This method for chaining
        """
        updates = (RemoveSnapshotRefUpdate(ref_name=ref_name),)
        requirements = (
            AssertRefSnapshotId(
                snapshot_id=self._transaction.table_metadata.refs[ref_name].snapshot_id
                if ref_name in self._transaction.table_metadata.refs
                else None,
                ref=ref_name,
            ),
        )
        self._updates += updates
        self._requirements += requirements
        return self

    def create_tag(self, snapshot_id: int, tag_name: str, max_ref_age_ms: int | None = None) -> ManageSnapshots:
        """
        Create a new tag pointing to the given snapshot id.

        Args:
            snapshot_id (int): snapshot id of the existing snapshot to tag
            tag_name (str): name of the tag
            max_ref_age_ms (Optional[int]): max ref age in milliseconds

        Returns:
            This for method chaining
        """
        update, requirement = self._transaction._set_ref_snapshot(
            snapshot_id=snapshot_id,
            ref_name=tag_name,
            type=SnapshotRefType.TAG,
            max_ref_age_ms=max_ref_age_ms,
        )
        self._updates += update
        self._requirements += requirement
        return self

    def remove_tag(self, tag_name: str) -> ManageSnapshots:
        """
        Remove a tag.

        Args:
            tag_name (str): name of tag to remove
        Returns:
            This for method chaining
        """
        return self._remove_ref_snapshot(ref_name=tag_name)

    def create_branch(
        self,
        snapshot_id: int,
        branch_name: str,
        max_ref_age_ms: int | None = None,
        max_snapshot_age_ms: int | None = None,
        min_snapshots_to_keep: int | None = None,
    ) -> ManageSnapshots:
        """
        Create a new branch pointing to the given snapshot id.

        Args:
            snapshot_id (int): snapshot id of existing snapshot at which the branch is created.
            branch_name (str): name of the new branch
            max_ref_age_ms (Optional[int]): max ref age in milliseconds
            max_snapshot_age_ms (Optional[int]): max age of snapshots to keep in milliseconds
            min_snapshots_to_keep (Optional[int]): min number of snapshots to keep for the branch
        Returns:
            This for method chaining
        """
        update, requirement = self._transaction._set_ref_snapshot(
            snapshot_id=snapshot_id,
            ref_name=branch_name,
            type=SnapshotRefType.BRANCH,
            max_ref_age_ms=max_ref_age_ms,
            max_snapshot_age_ms=max_snapshot_age_ms,
            min_snapshots_to_keep=min_snapshots_to_keep,
        )
        self._updates += update
        self._requirements += requirement
        return self

    def remove_branch(self, branch_name: str) -> ManageSnapshots:
        """
        Remove a branch.

        Args:
            branch_name (str): name of branch to remove
        Returns:
            This for method chaining
        """
        return self._remove_ref_snapshot(ref_name=branch_name)

    def set_current_snapshot(self, snapshot_id: int | None = None, ref_name: str | None = None) -> ManageSnapshots:
        """Set the current snapshot to a specific snapshot ID or ref.

        Args:
            snapshot_id: The ID of the snapshot to set as current.
            ref_name: The snapshot reference (branch or tag) to set as current.

        Returns:
            This for method chaining.

        Raises:
            ValueError: If neither or both arguments are provided, or if the snapshot/ref does not exist.
        """
        self._commit_if_ref_updates_exist()

        if (snapshot_id is None) == (ref_name is None):
            raise ValueError("Either snapshot_id or ref_name must be provided, not both")

        target_snapshot_id: int
        if snapshot_id is not None:
            target_snapshot_id = snapshot_id
        else:
            if ref_name not in self._transaction.table_metadata.refs:
                raise ValueError(f"Cannot find matching snapshot ID for ref: {ref_name}")
            target_snapshot_id = self._transaction.table_metadata.refs[ref_name].snapshot_id

        if self._transaction.table_metadata.snapshot_by_id(target_snapshot_id) is None:
            raise ValueError(f"Cannot set current snapshot to unknown snapshot id: {target_snapshot_id}")

        update, requirement = self._transaction._set_ref_snapshot(
            snapshot_id=target_snapshot_id,
            ref_name=MAIN_BRANCH,
            type=SnapshotRefType.BRANCH,
        )
        self._transaction._stage(update, requirement)
        return self

    def rollback_to_snapshot(self, snapshot_id: int) -> ManageSnapshots:
        """Rollback the table to the given snapshot id.

        The snapshot needs to be an ancestor of the current table state.

        Args:
            snapshot_id (int): rollback to this snapshot_id that used to be current.

        Returns:
            This for method chaining

        Raises:
            ValueError: If the snapshot does not exist or is not an ancestor of the current table state.
        """
        if not self._transaction.table_metadata.snapshot_by_id(snapshot_id):
            raise ValueError(f"Cannot roll back to unknown snapshot id: {snapshot_id}")

        if not self._is_current_ancestor(snapshot_id):
            raise ValueError(f"Cannot roll back to snapshot, not an ancestor of the current state: {snapshot_id}")

        return self.set_current_snapshot(snapshot_id=snapshot_id)

    def rollback_to_timestamp(self, timestamp_ms: int) -> ManageSnapshots:
        """Rollback the table to the latest snapshot before the given timestamp.

        Finds the latest ancestor snapshot whose timestamp is before the given timestamp and rolls back to it.

        Args:
            timestamp_ms: Rollback to the latest snapshot before this timestamp in milliseconds.

        Returns:
            This for method chaining

        Raises:
            ValueError: If no valid snapshot exists older than the given timestamp.
        """
        snapshot = latest_ancestor_before_timestamp(self._transaction.table_metadata, timestamp_ms)
        if snapshot is None:
            raise ValueError(f"Cannot roll back, no valid snapshot older than: {timestamp_ms}")

        return self.set_current_snapshot(snapshot_id=snapshot.snapshot_id)

    def _is_current_ancestor(self, snapshot_id: int) -> bool:
        return snapshot_id in self._current_ancestors()

    def _current_ancestors(self) -> set[int]:
        return {
            a.snapshot_id
            for a in ancestors_of(
                self._transaction.table_metadata.current_snapshot(),
                self._transaction.table_metadata,
            )
        }


class ExpireSnapshots(UpdateTableMetadata["ExpireSnapshots"]):
    """Expire snapshots by ID.

    Use table.expire_snapshots().<operation>().commit() to run a specific operation.
    Use table.expire_snapshots().<operation-one>().<operation-two>().commit() to run multiple operations.
    Pending changes are applied on commit.
    """

    _updates: tuple[TableUpdate, ...]
    _requirements: tuple[TableRequirement, ...]
    _snapshot_ids_to_expire: set[int]

    def __init__(self, transaction: Transaction) -> None:
        super().__init__(transaction)
        self._updates = ()
        self._requirements = ()
        self._snapshot_ids_to_expire = set()

    def _commit(self) -> UpdatesAndRequirements:
        """
        Commit the staged updates and requirements.

        This will remove the snapshots with the given IDs, but will always skip protected snapshots (branch/tag heads).

        Returns:
            Tuple of updates and requirements to be committed,
            as required by the calling parent apply functions.
        """
        # Remove any protected snapshot IDs from the set to expire, just in case
        protected_ids = self._get_protected_snapshot_ids()
        self._snapshot_ids_to_expire -= protected_ids
        update = RemoveSnapshotsUpdate(snapshot_ids=self._snapshot_ids_to_expire)
        self._updates += (update,)
        return self._updates, self._requirements

    def _get_protected_snapshot_ids(self) -> set[int]:
        """
        Get the IDs of protected snapshots.

        These are the HEAD snapshots of all branches and all tagged snapshots.  These ids are to be excluded from expiration.

        Returns:
            Set of protected snapshot IDs to exclude from expiration.
        """
        return {
            ref.snapshot_id
            for ref in self._transaction.table_metadata.refs.values()
            if ref.snapshot_ref_type in [SnapshotRefType.TAG, SnapshotRefType.BRANCH]
        }

    def by_id(self, snapshot_id: int) -> ExpireSnapshots:
        """
        Expire a snapshot by its ID.

        This will mark the snapshot for expiration.

        Args:
            snapshot_id (int): The ID of the snapshot to expire.
        Returns:
            This for method chaining.
        """
        if self._transaction.table_metadata.snapshot_by_id(snapshot_id) is None:
            raise ValueError(f"Snapshot with ID {snapshot_id} does not exist.")

        if snapshot_id in self._get_protected_snapshot_ids():
            raise ValueError(f"Snapshot with ID {snapshot_id} is protected and cannot be expired.")

        self._snapshot_ids_to_expire.add(snapshot_id)

        return self

    def by_ids(self, snapshot_ids: list[int]) -> ExpireSnapshots:
        """
        Expire multiple snapshots by their IDs.

        This will mark the snapshots for expiration.

        Args:
            snapshot_ids (List[int]): List of snapshot IDs to expire.
        Returns:
            This for method chaining.
        """
        for snapshot_id in snapshot_ids:
            self.by_id(snapshot_id)
        return self

    def older_than(self, dt: datetime) -> ExpireSnapshots:
        """
        Expire all unprotected snapshots with a timestamp older than a given value.

        Args:
            dt (datetime): Only snapshots with datetime < this value will be expired.

        Returns:
            This for method chaining.
        """
        protected_ids = self._get_protected_snapshot_ids()
        expire_from = datetime_to_millis(dt)
        for snapshot in self._transaction.table_metadata.snapshots:
            if snapshot.timestamp_ms < expire_from and snapshot.snapshot_id not in protected_ids:
                self._snapshot_ids_to_expire.add(snapshot.snapshot_id)
        return self
