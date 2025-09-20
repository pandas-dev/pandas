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
import os
import uuid
import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cached_property
from itertools import chain
from types import TracebackType
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Set,
    Tuple,
    Type,
    TypeVar,
    Union,
)

from pydantic import Field
from sortedcontainers import SortedList

import pyiceberg.expressions.parser as parser
from pyiceberg.expressions import (
    AlwaysFalse,
    AlwaysTrue,
    And,
    BooleanExpression,
    EqualTo,
    IsNull,
    Or,
    Reference,
)
from pyiceberg.expressions.visitors import (
    ResidualEvaluator,
    _InclusiveMetricsEvaluator,
    bind,
    expression_evaluator,
    inclusive_projection,
    manifest_evaluator,
)
from pyiceberg.io import FileIO, load_file_io
from pyiceberg.manifest import (
    POSITIONAL_DELETE_SCHEMA,
    DataFile,
    DataFileContent,
    ManifestContent,
    ManifestEntry,
    ManifestFile,
)
from pyiceberg.partitioning import (
    PARTITION_FIELD_ID_START,
    UNPARTITIONED_PARTITION_SPEC,
    PartitionKey,
    PartitionSpec,
)
from pyiceberg.schema import Schema
from pyiceberg.table.inspect import InspectTable
from pyiceberg.table.locations import LocationProvider, load_location_provider
from pyiceberg.table.maintenance import MaintenanceTable
from pyiceberg.table.metadata import (
    INITIAL_SEQUENCE_NUMBER,
    TableMetadata,
)
from pyiceberg.table.name_mapping import (
    NameMapping,
)
from pyiceberg.table.refs import MAIN_BRANCH, SnapshotRef
from pyiceberg.table.snapshots import (
    Snapshot,
    SnapshotLogEntry,
)
from pyiceberg.table.sorting import UNSORTED_SORT_ORDER, SortOrder
from pyiceberg.table.update import (
    AddPartitionSpecUpdate,
    AddSchemaUpdate,
    AddSortOrderUpdate,
    AssertCreate,
    AssertRefSnapshotId,
    AssertTableUUID,
    AssignUUIDUpdate,
    RemovePropertiesUpdate,
    SetCurrentSchemaUpdate,
    SetDefaultSortOrderUpdate,
    SetDefaultSpecUpdate,
    SetLocationUpdate,
    SetPropertiesUpdate,
    SetSnapshotRefUpdate,
    TableRequirement,
    TableUpdate,
    UpdatesAndRequirements,
    UpgradeFormatVersionUpdate,
    update_table_metadata,
)
from pyiceberg.table.update.schema import UpdateSchema
from pyiceberg.table.update.snapshot import ManageSnapshots, UpdateSnapshot, _FastAppendFiles
from pyiceberg.table.update.spec import UpdateSpec
from pyiceberg.table.update.statistics import UpdateStatistics
from pyiceberg.transforms import IdentityTransform
from pyiceberg.typedef import (
    EMPTY_DICT,
    IcebergBaseModel,
    IcebergRootModel,
    Identifier,
    KeyDefaultDict,
    Properties,
    Record,
    TableVersion,
)
from pyiceberg.types import (
    strtobool,
)
from pyiceberg.utils.concurrent import ExecutorFactory
from pyiceberg.utils.config import Config
from pyiceberg.utils.properties import property_as_bool

if TYPE_CHECKING:
    import bodo.pandas as bd
    import daft
    import pandas as pd
    import polars as pl
    import pyarrow as pa
    import ray
    from duckdb import DuckDBPyConnection
    from pyiceberg_core.datafusion import IcebergDataFusionTable

    from pyiceberg.catalog import Catalog

ALWAYS_TRUE = AlwaysTrue()
DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE = "downcast-ns-timestamp-to-us-on-write"


@dataclass()
class UpsertResult:
    """Summary the upsert operation."""

    rows_updated: int = 0
    rows_inserted: int = 0


class TableProperties:
    PARQUET_ROW_GROUP_SIZE_BYTES = "write.parquet.row-group-size-bytes"
    PARQUET_ROW_GROUP_SIZE_BYTES_DEFAULT = 128 * 1024 * 1024  # 128 MB

    PARQUET_ROW_GROUP_LIMIT = "write.parquet.row-group-limit"
    PARQUET_ROW_GROUP_LIMIT_DEFAULT = 1048576

    PARQUET_PAGE_SIZE_BYTES = "write.parquet.page-size-bytes"
    PARQUET_PAGE_SIZE_BYTES_DEFAULT = 1024 * 1024  # 1 MB

    PARQUET_PAGE_ROW_LIMIT = "write.parquet.page-row-limit"
    PARQUET_PAGE_ROW_LIMIT_DEFAULT = 20000

    PARQUET_DICT_SIZE_BYTES = "write.parquet.dict-size-bytes"
    PARQUET_DICT_SIZE_BYTES_DEFAULT = 2 * 1024 * 1024  # 2 MB

    PARQUET_COMPRESSION = "write.parquet.compression-codec"
    PARQUET_COMPRESSION_DEFAULT = "zstd"

    PARQUET_COMPRESSION_LEVEL = "write.parquet.compression-level"
    PARQUET_COMPRESSION_LEVEL_DEFAULT = None

    PARQUET_BLOOM_FILTER_MAX_BYTES = "write.parquet.bloom-filter-max-bytes"
    PARQUET_BLOOM_FILTER_MAX_BYTES_DEFAULT = 1024 * 1024

    PARQUET_BLOOM_FILTER_COLUMN_ENABLED_PREFIX = "write.parquet.bloom-filter-enabled.column"

    WRITE_TARGET_FILE_SIZE_BYTES = "write.target-file-size-bytes"
    WRITE_TARGET_FILE_SIZE_BYTES_DEFAULT = 512 * 1024 * 1024  # 512 MB

    WRITE_AVRO_COMPRESSION = "write.avro.compression-codec"
    WRITE_AVRO_COMPRESSION_DEFAULT = "gzip"

    DEFAULT_WRITE_METRICS_MODE = "write.metadata.metrics.default"
    DEFAULT_WRITE_METRICS_MODE_DEFAULT = "truncate(16)"

    METRICS_MODE_COLUMN_CONF_PREFIX = "write.metadata.metrics.column"

    WRITE_PARTITION_SUMMARY_LIMIT = "write.summary.partition-limit"
    WRITE_PARTITION_SUMMARY_LIMIT_DEFAULT = 0

    WRITE_PY_LOCATION_PROVIDER_IMPL = "write.py-location-provider.impl"

    OBJECT_STORE_ENABLED = "write.object-storage.enabled"
    OBJECT_STORE_ENABLED_DEFAULT = False

    WRITE_OBJECT_STORE_PARTITIONED_PATHS = "write.object-storage.partitioned-paths"
    WRITE_OBJECT_STORE_PARTITIONED_PATHS_DEFAULT = True

    WRITE_DATA_PATH = "write.data.path"
    WRITE_METADATA_PATH = "write.metadata.path"

    DELETE_MODE = "write.delete.mode"
    DELETE_MODE_COPY_ON_WRITE = "copy-on-write"
    DELETE_MODE_MERGE_ON_READ = "merge-on-read"
    DELETE_MODE_DEFAULT = DELETE_MODE_COPY_ON_WRITE

    DEFAULT_NAME_MAPPING = "schema.name-mapping.default"
    FORMAT_VERSION = "format-version"
    DEFAULT_FORMAT_VERSION: TableVersion = 2

    MANIFEST_TARGET_SIZE_BYTES = "commit.manifest.target-size-bytes"
    MANIFEST_TARGET_SIZE_BYTES_DEFAULT = 8 * 1024 * 1024  # 8 MB

    MANIFEST_MIN_MERGE_COUNT = "commit.manifest.min-count-to-merge"
    MANIFEST_MIN_MERGE_COUNT_DEFAULT = 100

    MANIFEST_MERGE_ENABLED = "commit.manifest-merge.enabled"
    MANIFEST_MERGE_ENABLED_DEFAULT = False

    METADATA_PREVIOUS_VERSIONS_MAX = "write.metadata.previous-versions-max"
    METADATA_PREVIOUS_VERSIONS_MAX_DEFAULT = 100

    METADATA_DELETE_AFTER_COMMIT_ENABLED = "write.metadata.delete-after-commit.enabled"
    METADATA_DELETE_AFTER_COMMIT_ENABLED_DEFAULT = False

    MAX_SNAPSHOT_AGE_MS = "history.expire.max-snapshot-age-ms"
    MAX_SNAPSHOT_AGE_MS_DEFAULT = 5 * 24 * 60 * 60 * 1000  # 5 days

    MIN_SNAPSHOTS_TO_KEEP = "history.expire.min-snapshots-to-keep"
    MIN_SNAPSHOTS_TO_KEEP_DEFAULT = 1


class Transaction:
    _table: Table
    _autocommit: bool
    _updates: Tuple[TableUpdate, ...]
    _requirements: Tuple[TableRequirement, ...]

    def __init__(self, table: Table, autocommit: bool = False):
        """Open a transaction to stage and commit changes to a table.

        Args:
            table: The table that will be altered.
            autocommit: Option to automatically commit the changes when they are staged.
        """
        self._table = table
        self._autocommit = autocommit
        self._updates = ()
        self._requirements = ()

    @property
    def table_metadata(self) -> TableMetadata:
        return update_table_metadata(self._table.metadata, self._updates)

    def __enter__(self) -> Transaction:
        """Start a transaction to update the table."""
        return self

    def __exit__(
        self, exctype: Optional[Type[BaseException]], excinst: Optional[BaseException], exctb: Optional[TracebackType]
    ) -> None:
        """Close and commit the transaction if no exceptions have been raised."""
        if exctype is None and excinst is None and exctb is None:
            self.commit_transaction()

    def _apply(self, updates: Tuple[TableUpdate, ...], requirements: Tuple[TableRequirement, ...] = ()) -> Transaction:
        """Check if the requirements are met, and applies the updates to the metadata."""
        for requirement in requirements:
            requirement.validate(self.table_metadata)

        self._updates += updates

        # For the requirements, it does not make sense to add a requirement more than once
        # For example, you cannot assert that the current schema has two different IDs
        existing_requirements = {type(requirement) for requirement in self._requirements}
        for new_requirement in requirements:
            if type(new_requirement) not in existing_requirements:
                self._requirements = self._requirements + (new_requirement,)

        if self._autocommit:
            self.commit_transaction()

        return self

    def _scan(self, row_filter: Union[str, BooleanExpression] = ALWAYS_TRUE, case_sensitive: bool = True) -> DataScan:
        """Minimal data scan of the table with the current state of the transaction."""
        return DataScan(
            table_metadata=self.table_metadata, io=self._table.io, row_filter=row_filter, case_sensitive=case_sensitive
        )

    def upgrade_table_version(self, format_version: TableVersion) -> Transaction:
        """Set the table to a certain version.

        Args:
            format_version: The newly set version.

        Returns:
            The alter table builder.
        """
        if format_version not in {1, 2}:
            raise ValueError(f"Unsupported table format version: {format_version}")

        if format_version < self.table_metadata.format_version:
            raise ValueError(f"Cannot downgrade v{self.table_metadata.format_version} table to v{format_version}")

        if format_version > self.table_metadata.format_version:
            return self._apply((UpgradeFormatVersionUpdate(format_version=format_version),))

        return self

    def set_properties(self, properties: Properties = EMPTY_DICT, **kwargs: Any) -> Transaction:
        """Set properties.

        When a property is already set, it will be overwritten.

        Args:
            properties: The properties set on the table.
            kwargs: properties can also be pass as kwargs.

        Returns:
            The alter table builder.
        """
        if properties and kwargs:
            raise ValueError("Cannot pass both properties and kwargs")
        updates = properties or kwargs
        return self._apply((SetPropertiesUpdate(updates=updates),))

    def _set_ref_snapshot(
        self,
        snapshot_id: int,
        ref_name: str,
        type: str,
        max_ref_age_ms: Optional[int] = None,
        max_snapshot_age_ms: Optional[int] = None,
        min_snapshots_to_keep: Optional[int] = None,
    ) -> UpdatesAndRequirements:
        """Update a ref to a snapshot.

        Returns:
            The updates and requirements for the set-snapshot-ref staged
        """
        updates = (
            SetSnapshotRefUpdate(
                snapshot_id=snapshot_id,
                ref_name=ref_name,
                type=type,
                max_ref_age_ms=max_ref_age_ms,
                max_snapshot_age_ms=max_snapshot_age_ms,
                min_snapshots_to_keep=min_snapshots_to_keep,
            ),
        )
        requirements = (
            AssertRefSnapshotId(
                snapshot_id=self.table_metadata.refs[ref_name].snapshot_id if ref_name in self.table_metadata.refs else None,
                ref=ref_name,
            ),
        )

        return updates, requirements

    def _build_partition_predicate(self, partition_records: Set[Record]) -> BooleanExpression:
        """Build a filter predicate matching any of the input partition records.

        Args:
            partition_records: A set of partition records to match
        Returns:
            A predicate matching any of the input partition records.
        """
        partition_spec = self.table_metadata.spec()
        schema = self.table_metadata.schema()
        partition_fields = [schema.find_field(field.source_id).name for field in partition_spec.fields]

        expr: BooleanExpression = AlwaysFalse()
        for partition_record in partition_records:
            match_partition_expression: BooleanExpression = AlwaysTrue()

            for pos, partition_field in enumerate(partition_fields):
                predicate = (
                    EqualTo(Reference(partition_field), partition_record[pos])
                    if partition_record[pos] is not None
                    else IsNull(Reference(partition_field))
                )
                match_partition_expression = And(match_partition_expression, predicate)
            expr = Or(expr, match_partition_expression)
        return expr

    def _append_snapshot_producer(
        self, snapshot_properties: Dict[str, str], branch: Optional[str] = MAIN_BRANCH
    ) -> _FastAppendFiles:
        """Determine the append type based on table properties.

        Args:
            snapshot_properties: Custom properties to be added to the snapshot summary
        Returns:
            Either a fast-append or a merge-append snapshot producer.
        """
        manifest_merge_enabled = property_as_bool(
            self.table_metadata.properties,
            TableProperties.MANIFEST_MERGE_ENABLED,
            TableProperties.MANIFEST_MERGE_ENABLED_DEFAULT,
        )
        update_snapshot = self.update_snapshot(snapshot_properties=snapshot_properties, branch=branch)
        return update_snapshot.merge_append() if manifest_merge_enabled else update_snapshot.fast_append()

    def update_schema(self, allow_incompatible_changes: bool = False, case_sensitive: bool = True) -> UpdateSchema:
        """Create a new UpdateSchema to alter the columns of this table.

        Args:
            allow_incompatible_changes: If changes are allowed that might break downstream consumers.
            case_sensitive: If field names are case-sensitive.

        Returns:
            A new UpdateSchema.
        """
        return UpdateSchema(
            self,
            allow_incompatible_changes=allow_incompatible_changes,
            case_sensitive=case_sensitive,
            name_mapping=self.table_metadata.name_mapping(),
        )

    def update_snapshot(
        self, snapshot_properties: Dict[str, str] = EMPTY_DICT, branch: Optional[str] = MAIN_BRANCH
    ) -> UpdateSnapshot:
        """Create a new UpdateSnapshot to produce a new snapshot for the table.

        Returns:
            A new UpdateSnapshot
        """
        return UpdateSnapshot(self, io=self._table.io, branch=branch, snapshot_properties=snapshot_properties)

    def update_statistics(self) -> UpdateStatistics:
        """
        Create a new UpdateStatistics to update the statistics of the table.

        Returns:
            A new UpdateStatistics
        """
        return UpdateStatistics(transaction=self)

    def append(self, df: pa.Table, snapshot_properties: Dict[str, str] = EMPTY_DICT, branch: Optional[str] = MAIN_BRANCH) -> None:
        """
        Shorthand API for appending a PyArrow table to a table transaction.

        Args:
            df: The Arrow dataframe that will be appended to overwrite the table
            snapshot_properties: Custom properties to be added to the snapshot summary
            branch: Branch Reference to run the append operation
        """
        try:
            import pyarrow as pa
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError("For writes PyArrow needs to be installed") from e

        from pyiceberg.io.pyarrow import _check_pyarrow_schema_compatible, _dataframe_to_data_files

        if not isinstance(df, pa.Table):
            raise ValueError(f"Expected PyArrow table, got: {df}")

        downcast_ns_timestamp_to_us = Config().get_bool(DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE) or False
        _check_pyarrow_schema_compatible(
            self.table_metadata.schema(),
            provided_schema=df.schema,
            downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us,
            format_version=self.table_metadata.format_version,
        )

        with self._append_snapshot_producer(snapshot_properties, branch=branch) as append_files:
            # skip writing data files if the dataframe is empty
            if df.shape[0] > 0:
                data_files = list(
                    _dataframe_to_data_files(
                        table_metadata=self.table_metadata, write_uuid=append_files.commit_uuid, df=df, io=self._table.io
                    )
                )
                for data_file in data_files:
                    append_files.append_data_file(data_file)

    def dynamic_partition_overwrite(
        self, df: pa.Table, snapshot_properties: Dict[str, str] = EMPTY_DICT, branch: Optional[str] = MAIN_BRANCH
    ) -> None:
        """
        Shorthand for overwriting existing partitions with a PyArrow table.

        The function detects partition values in the provided arrow table using the current
        partition spec, and deletes existing partitions matching these values. Finally, the
        data in the table is appended to the table.

        Args:
            df: The Arrow dataframe that will be used to overwrite the table
            snapshot_properties: Custom properties to be added to the snapshot summary
            branch: Branch Reference to run the dynamic partition overwrite operation
        """
        try:
            import pyarrow as pa
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError("For writes PyArrow needs to be installed") from e

        from pyiceberg.io.pyarrow import _check_pyarrow_schema_compatible, _dataframe_to_data_files

        if not isinstance(df, pa.Table):
            raise ValueError(f"Expected PyArrow table, got: {df}")

        if self.table_metadata.spec().is_unpartitioned():
            raise ValueError("Cannot apply dynamic overwrite on an unpartitioned table.")

        for field in self.table_metadata.spec().fields:
            if not isinstance(field.transform, IdentityTransform):
                raise ValueError(
                    f"For now dynamic overwrite does not support a table with non-identity-transform field in the latest partition spec: {field}"
                )

        downcast_ns_timestamp_to_us = Config().get_bool(DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE) or False
        _check_pyarrow_schema_compatible(
            self.table_metadata.schema(),
            provided_schema=df.schema,
            downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us,
            format_version=self.table_metadata.format_version,
        )

        # If dataframe does not have data, there is no need to overwrite
        if df.shape[0] == 0:
            return

        append_snapshot_commit_uuid = uuid.uuid4()
        data_files: List[DataFile] = list(
            _dataframe_to_data_files(
                table_metadata=self._table.metadata, write_uuid=append_snapshot_commit_uuid, df=df, io=self._table.io
            )
        )

        partitions_to_overwrite = {data_file.partition for data_file in data_files}
        delete_filter = self._build_partition_predicate(partition_records=partitions_to_overwrite)
        self.delete(delete_filter=delete_filter, snapshot_properties=snapshot_properties, branch=branch)

        with self._append_snapshot_producer(snapshot_properties, branch=branch) as append_files:
            append_files.commit_uuid = append_snapshot_commit_uuid
            for data_file in data_files:
                append_files.append_data_file(data_file)

    def overwrite(
        self,
        df: pa.Table,
        overwrite_filter: Union[BooleanExpression, str] = ALWAYS_TRUE,
        snapshot_properties: Dict[str, str] = EMPTY_DICT,
        case_sensitive: bool = True,
        branch: Optional[str] = MAIN_BRANCH,
    ) -> None:
        """
        Shorthand for adding a table overwrite with a PyArrow table to the transaction.

        An overwrite may produce zero or more snapshots based on the operation:

            - DELETE: In case existing Parquet files can be dropped completely.
            - REPLACE: In case existing Parquet files need to be rewritten.
            - APPEND: In case new data is being inserted into the table.

        Args:
            df: The Arrow dataframe that will be used to overwrite the table
            overwrite_filter: ALWAYS_TRUE when you overwrite all the data,
                              or a boolean expression in case of a partial overwrite
            snapshot_properties: Custom properties to be added to the snapshot summary
            case_sensitive: A bool determine if the provided `overwrite_filter` is case-sensitive
            branch: Branch Reference to run the overwrite operation
        """
        try:
            import pyarrow as pa
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError("For writes PyArrow needs to be installed") from e

        from pyiceberg.io.pyarrow import _check_pyarrow_schema_compatible, _dataframe_to_data_files

        if not isinstance(df, pa.Table):
            raise ValueError(f"Expected PyArrow table, got: {df}")

        downcast_ns_timestamp_to_us = Config().get_bool(DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE) or False
        _check_pyarrow_schema_compatible(
            self.table_metadata.schema(),
            provided_schema=df.schema,
            downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us,
            format_version=self.table_metadata.format_version,
        )

        if overwrite_filter != AlwaysFalse():
            # Only delete when the filter is != AlwaysFalse
            self.delete(
                delete_filter=overwrite_filter,
                case_sensitive=case_sensitive,
                snapshot_properties=snapshot_properties,
                branch=branch,
            )

        with self._append_snapshot_producer(snapshot_properties, branch=branch) as append_files:
            # skip writing data files if the dataframe is empty
            if df.shape[0] > 0:
                data_files = _dataframe_to_data_files(
                    table_metadata=self.table_metadata, write_uuid=append_files.commit_uuid, df=df, io=self._table.io
                )
                for data_file in data_files:
                    append_files.append_data_file(data_file)

    def delete(
        self,
        delete_filter: Union[str, BooleanExpression],
        snapshot_properties: Dict[str, str] = EMPTY_DICT,
        case_sensitive: bool = True,
        branch: Optional[str] = MAIN_BRANCH,
    ) -> None:
        """
        Shorthand for deleting record from a table.

        A delete may produce zero or more snapshots based on the operation:

            - DELETE: In case existing Parquet files can be dropped completely.
            - REPLACE: In case existing Parquet files need to be rewritten

        Args:
            delete_filter: A boolean expression to delete rows from a table
            snapshot_properties: Custom properties to be added to the snapshot summary
            case_sensitive: A bool determine if the provided `delete_filter` is case-sensitive
            branch: Branch Reference to run the delete operation
        """
        from pyiceberg.io.pyarrow import (
            ArrowScan,
            _dataframe_to_data_files,
            _expression_to_complementary_pyarrow,
        )

        if (
            self.table_metadata.properties.get(TableProperties.DELETE_MODE, TableProperties.DELETE_MODE_DEFAULT)
            == TableProperties.DELETE_MODE_MERGE_ON_READ
        ):
            warnings.warn("Merge on read is not yet supported, falling back to copy-on-write")

        if isinstance(delete_filter, str):
            delete_filter = _parse_row_filter(delete_filter)

        with self.update_snapshot(snapshot_properties=snapshot_properties, branch=branch).delete() as delete_snapshot:
            delete_snapshot.delete_by_predicate(delete_filter, case_sensitive)

        # Check if there are any files that require an actual rewrite of a data file
        if delete_snapshot.rewrites_needed is True:
            bound_delete_filter = bind(self.table_metadata.schema(), delete_filter, case_sensitive)
            preserve_row_filter = _expression_to_complementary_pyarrow(bound_delete_filter)

            file_scan = self._scan(row_filter=delete_filter, case_sensitive=case_sensitive)
            if branch is not None:
                file_scan = file_scan.use_ref(branch)
            files = file_scan.plan_files()

            commit_uuid = uuid.uuid4()
            counter = itertools.count(0)

            replaced_files: List[Tuple[DataFile, List[DataFile]]] = []
            # This will load the Parquet file into memory, including:
            #   - Filter out the rows based on the delete filter
            #   - Projecting it to the current schema
            #   - Applying the positional deletes if they are there
            # When writing
            #   - Apply the latest partition-spec
            #   - And sort order when added
            for original_file in files:
                df = ArrowScan(
                    table_metadata=self.table_metadata,
                    io=self._table.io,
                    projected_schema=self.table_metadata.schema(),
                    row_filter=AlwaysTrue(),
                ).to_table(tasks=[original_file])
                filtered_df = df.filter(preserve_row_filter)

                # Only rewrite if there are records being deleted
                if len(filtered_df) == 0:
                    replaced_files.append((original_file.file, []))
                elif len(df) != len(filtered_df):
                    replaced_files.append(
                        (
                            original_file.file,
                            list(
                                _dataframe_to_data_files(
                                    io=self._table.io,
                                    df=filtered_df,
                                    table_metadata=self.table_metadata,
                                    write_uuid=commit_uuid,
                                    counter=counter,
                                )
                            ),
                        )
                    )

            if len(replaced_files) > 0:
                with self.update_snapshot(
                    snapshot_properties=snapshot_properties, branch=branch
                ).overwrite() as overwrite_snapshot:
                    overwrite_snapshot.commit_uuid = commit_uuid
                    for original_data_file, replaced_data_files in replaced_files:
                        overwrite_snapshot.delete_data_file(original_data_file)
                        for replaced_data_file in replaced_data_files:
                            overwrite_snapshot.append_data_file(replaced_data_file)

        if not delete_snapshot.files_affected and not delete_snapshot.rewrites_needed:
            warnings.warn("Delete operation did not match any records")

    def upsert(
        self,
        df: pa.Table,
        join_cols: Optional[List[str]] = None,
        when_matched_update_all: bool = True,
        when_not_matched_insert_all: bool = True,
        case_sensitive: bool = True,
        branch: Optional[str] = MAIN_BRANCH,
    ) -> UpsertResult:
        """Shorthand API for performing an upsert to an iceberg table.

        Args:

            df: The input dataframe to upsert with the table's data.
            join_cols: Columns to join on, if not provided, it will use the identifier-field-ids.
            when_matched_update_all: Bool indicating to update rows that are matched but require an update due to a value in a non-key column changing
            when_not_matched_insert_all: Bool indicating new rows to be inserted that do not match any existing rows in the table
            case_sensitive: Bool indicating if the match should be case-sensitive
            branch: Branch Reference to run the upsert operation

            To learn more about the identifier-field-ids: https://iceberg.apache.org/spec/#identifier-field-ids

                Example Use Cases:
                    Case 1: Both Parameters = True (Full Upsert)
                    Existing row found → Update it
                    New row found → Insert it

                    Case 2: when_matched_update_all = False, when_not_matched_insert_all = True
                    Existing row found → Do nothing (no updates)
                    New row found → Insert it

                    Case 3: when_matched_update_all = True, when_not_matched_insert_all = False
                    Existing row found → Update it
                    New row found → Do nothing (no inserts)

                    Case 4: Both Parameters = False (No Merge Effect)
                    Existing row found → Do nothing
                    New row found → Do nothing
                    (Function effectively does nothing)


        Returns:
            An UpsertResult class (contains details of rows updated and inserted)
        """
        try:
            import pyarrow as pa  # noqa: F401
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError("For writes PyArrow needs to be installed") from e

        from pyiceberg.io.pyarrow import expression_to_pyarrow
        from pyiceberg.table import upsert_util

        if join_cols is None:
            join_cols = []
            for field_id in self.table_metadata.schema().identifier_field_ids:
                col = self.table_metadata.schema().find_column_name(field_id)
                if col is not None:
                    join_cols.append(col)
                else:
                    raise ValueError(f"Field-ID could not be found: {join_cols}")

        if len(join_cols) == 0:
            raise ValueError("Join columns could not be found, please set identifier-field-ids or pass in explicitly.")

        if not when_matched_update_all and not when_not_matched_insert_all:
            raise ValueError("no upsert options selected...exiting")

        if upsert_util.has_duplicate_rows(df, join_cols):
            raise ValueError("Duplicate rows found in source dataset based on the key columns. No upsert executed")

        from pyiceberg.io.pyarrow import _check_pyarrow_schema_compatible

        downcast_ns_timestamp_to_us = Config().get_bool(DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE) or False
        _check_pyarrow_schema_compatible(
            self.table_metadata.schema(),
            provided_schema=df.schema,
            downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us,
            format_version=self.table_metadata.format_version,
        )

        # get list of rows that exist so we don't have to load the entire target table
        matched_predicate = upsert_util.create_match_filter(df, join_cols)

        # We must use Transaction.table_metadata for the scan. This includes all uncommitted - but relevant - changes.

        matched_iceberg_record_batches_scan = DataScan(
            table_metadata=self.table_metadata,
            io=self._table.io,
            row_filter=matched_predicate,
            case_sensitive=case_sensitive,
        )

        if branch in self.table_metadata.refs:
            matched_iceberg_record_batches_scan = matched_iceberg_record_batches_scan.use_ref(branch)

        matched_iceberg_record_batches = matched_iceberg_record_batches_scan.to_arrow_batch_reader()

        batches_to_overwrite = []
        overwrite_predicates = []
        rows_to_insert = df

        for batch in matched_iceberg_record_batches:
            rows = pa.Table.from_batches([batch])

            if when_matched_update_all:
                # function get_rows_to_update is doing a check on non-key columns to see if any of the values have actually changed
                # we don't want to do just a blanket overwrite for matched rows if the actual non-key column data hasn't changed
                # this extra step avoids unnecessary IO and writes
                rows_to_update = upsert_util.get_rows_to_update(df, rows, join_cols)

                if len(rows_to_update) > 0:
                    # build the match predicate filter
                    overwrite_mask_predicate = upsert_util.create_match_filter(rows_to_update, join_cols)

                    batches_to_overwrite.append(rows_to_update)
                    overwrite_predicates.append(overwrite_mask_predicate)

            if when_not_matched_insert_all:
                expr_match = upsert_util.create_match_filter(rows, join_cols)
                expr_match_bound = bind(self.table_metadata.schema(), expr_match, case_sensitive=case_sensitive)
                expr_match_arrow = expression_to_pyarrow(expr_match_bound)

                # Filter rows per batch.
                rows_to_insert = rows_to_insert.filter(~expr_match_arrow)

        update_row_cnt = 0
        insert_row_cnt = 0

        if batches_to_overwrite:
            rows_to_update = pa.concat_tables(batches_to_overwrite)
            update_row_cnt = len(rows_to_update)
            self.overwrite(
                rows_to_update,
                overwrite_filter=Or(*overwrite_predicates) if len(overwrite_predicates) > 1 else overwrite_predicates[0],
                branch=branch,
            )

        if when_not_matched_insert_all:
            insert_row_cnt = len(rows_to_insert)
            if rows_to_insert:
                self.append(rows_to_insert, branch=branch)

        return UpsertResult(rows_updated=update_row_cnt, rows_inserted=insert_row_cnt)

    def add_files(
        self, file_paths: List[str], snapshot_properties: Dict[str, str] = EMPTY_DICT, check_duplicate_files: bool = True
    ) -> None:
        """
        Shorthand API for adding files as data files to the table transaction.

        Args:
            file_paths: The list of full file paths to be added as data files to the table

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: Raises a ValueError given file_paths contains duplicate files
            ValueError: Raises a ValueError given file_paths already referenced by table
        """
        if len(file_paths) != len(set(file_paths)):
            raise ValueError("File paths must be unique")

        if check_duplicate_files:
            import pyarrow.compute as pc

            expr = pc.field("file_path").isin(file_paths)
            referenced_files = [file["file_path"] for file in self._table.inspect.data_files().filter(expr).to_pylist()]

            if referenced_files:
                raise ValueError(f"Cannot add files that are already referenced by table, files: {', '.join(referenced_files)}")

        if self.table_metadata.name_mapping() is None:
            self.set_properties(
                **{TableProperties.DEFAULT_NAME_MAPPING: self.table_metadata.schema().name_mapping.model_dump_json()}
            )
        with self.update_snapshot(snapshot_properties=snapshot_properties).fast_append() as update_snapshot:
            data_files = _parquet_files_to_data_files(
                table_metadata=self.table_metadata, file_paths=file_paths, io=self._table.io
            )
            for data_file in data_files:
                update_snapshot.append_data_file(data_file)

    def update_spec(self) -> UpdateSpec:
        """Create a new UpdateSpec to update the partitioning of the table.

        Returns:
            A new UpdateSpec.
        """
        return UpdateSpec(self)

    def remove_properties(self, *removals: str) -> Transaction:
        """Remove properties.

        Args:
            removals: Properties to be removed.

        Returns:
            The alter table builder.
        """
        return self._apply((RemovePropertiesUpdate(removals=removals),))

    def update_location(self, location: str) -> Transaction:
        """Set the new table location.

        Args:
            location: The new location of the table.

        Returns:
            The alter table builder.
        """
        raise NotImplementedError("Not yet implemented")

    def commit_transaction(self) -> Table:
        """Commit the changes to the catalog.

        Returns:
            The table with the updates applied.
        """
        if len(self._updates) > 0:
            self._requirements += (AssertTableUUID(uuid=self.table_metadata.table_uuid),)
            self._table._do_commit(  # pylint: disable=W0212
                updates=self._updates,
                requirements=self._requirements,
            )

        self._updates = ()
        self._requirements = ()

        return self._table


class CreateTableTransaction(Transaction):
    """A transaction that involves the creation of a new table."""

    def _initial_changes(self, table_metadata: TableMetadata) -> None:
        """Set the initial changes that can reconstruct the initial table metadata when creating the CreateTableTransaction."""
        self._updates += (
            AssignUUIDUpdate(uuid=table_metadata.table_uuid),
            UpgradeFormatVersionUpdate(format_version=table_metadata.format_version),
        )

        schema: Schema = table_metadata.schema()
        self._updates += (
            AddSchemaUpdate(schema_=schema),
            SetCurrentSchemaUpdate(schema_id=-1),
        )

        spec: PartitionSpec = table_metadata.spec()
        if spec.is_unpartitioned():
            self._updates += (AddPartitionSpecUpdate(spec=UNPARTITIONED_PARTITION_SPEC),)
        else:
            self._updates += (AddPartitionSpecUpdate(spec=spec),)
        self._updates += (SetDefaultSpecUpdate(spec_id=-1),)

        sort_order: Optional[SortOrder] = table_metadata.sort_order_by_id(table_metadata.default_sort_order_id)
        if sort_order is None or sort_order.is_unsorted:
            self._updates += (AddSortOrderUpdate(sort_order=UNSORTED_SORT_ORDER),)
        else:
            self._updates += (AddSortOrderUpdate(sort_order=sort_order),)
        self._updates += (SetDefaultSortOrderUpdate(sort_order_id=-1),)

        self._updates += (
            SetLocationUpdate(location=table_metadata.location),
            SetPropertiesUpdate(updates=table_metadata.properties),
        )

    def __init__(self, table: StagedTable):
        super().__init__(table, autocommit=False)
        self._initial_changes(table.metadata)

    def commit_transaction(self) -> Table:
        """Commit the changes to the catalog.

        In the case of a CreateTableTransaction, the only requirement is AssertCreate.
        Returns:
            The table with the updates applied.
        """
        if len(self._updates) > 0:
            self._table._do_commit(  # pylint: disable=W0212
                updates=self._updates,
                requirements=(AssertCreate(),),
            )

        self._updates = ()
        self._requirements = ()

        return self._table


class Namespace(IcebergRootModel[List[str]]):
    """Reference to one or more levels of a namespace."""

    root: List[str] = Field(
        ...,
        description="Reference to one or more levels of a namespace",
    )


class TableIdentifier(IcebergBaseModel):
    """Fully Qualified identifier to a table."""

    namespace: Namespace
    name: str


class CommitTableRequest(IcebergBaseModel):
    """A pydantic BaseModel for a table commit request."""

    identifier: TableIdentifier = Field()
    requirements: Tuple[TableRequirement, ...] = Field(default_factory=tuple)
    updates: Tuple[TableUpdate, ...] = Field(default_factory=tuple)


class CommitTableResponse(IcebergBaseModel):
    """A pydantic BaseModel for a table commit response."""

    metadata: TableMetadata
    metadata_location: str = Field(alias="metadata-location")


class Table:
    """An Iceberg table."""

    _identifier: Identifier = Field()
    metadata: TableMetadata
    metadata_location: str = Field()
    io: FileIO
    catalog: Catalog
    config: Dict[str, str]

    def __init__(
        self,
        identifier: Identifier,
        metadata: TableMetadata,
        metadata_location: str,
        io: FileIO,
        catalog: Catalog,
        config: Dict[str, str] = EMPTY_DICT,
    ) -> None:
        self._identifier = identifier
        self.metadata = metadata
        self.metadata_location = metadata_location
        self.io = io
        self.catalog = catalog
        self.config = config

    def transaction(self) -> Transaction:
        """Create a new transaction object to first stage the changes, and then commit them to the catalog.

        Returns:
            The transaction object
        """
        return Transaction(self)

    @property
    def inspect(self) -> InspectTable:
        """Return the InspectTable object to browse the table metadata.

        Returns:
            InspectTable object based on this Table.
        """
        return InspectTable(self)

    @property
    def maintenance(self) -> MaintenanceTable:
        """Return the MaintenanceTable object for maintenance.

        Returns:
            MaintenanceTable object based on this Table.
        """
        return MaintenanceTable(self)

    def refresh(self) -> Table:
        """Refresh the current table metadata.

        Returns:
            An updated instance of the same Iceberg table
        """
        fresh = self.catalog.load_table(self._identifier)
        self.metadata = fresh.metadata
        self.io = fresh.io
        self.metadata_location = fresh.metadata_location
        return self

    def name(self) -> Identifier:
        """Return the identifier of this table.

        Returns:
            An Identifier tuple of the table name
        """
        return self._identifier

    def scan(
        self,
        row_filter: Union[str, BooleanExpression] = ALWAYS_TRUE,
        selected_fields: Tuple[str, ...] = ("*",),
        case_sensitive: bool = True,
        snapshot_id: Optional[int] = None,
        options: Properties = EMPTY_DICT,
        limit: Optional[int] = None,
    ) -> DataScan:
        """Fetch a DataScan based on the table's current metadata.

            The data scan can be used to project the table's data
            that matches the provided row_filter onto the table's
            current schema.

        Args:
            row_filter:
                A string or BooleanExpression that describes the
                desired rows
            selected_fields:
                A tuple of strings representing the column names
                to return in the output dataframe.
            case_sensitive:
                If True column matching is case sensitive
            snapshot_id:
                Optional Snapshot ID to time travel to. If None,
                scans the table as of the current snapshot ID.
            options:
                Additional Table properties as a dictionary of
                string key value pairs to use for this scan.
            limit:
                An integer representing the number of rows to
                return in the scan result. If None, fetches all
                matching rows.

        Returns:
            A DataScan based on the table's current metadata.
        """
        return DataScan(
            table_metadata=self.metadata,
            io=self.io,
            row_filter=row_filter,
            selected_fields=selected_fields,
            case_sensitive=case_sensitive,
            snapshot_id=snapshot_id,
            options=options,
            limit=limit,
        )

    @property
    def format_version(self) -> TableVersion:
        return self.metadata.format_version

    def schema(self) -> Schema:
        """Return the schema for this table."""
        return next(schema for schema in self.metadata.schemas if schema.schema_id == self.metadata.current_schema_id)

    def schemas(self) -> Dict[int, Schema]:
        """Return a dict of the schema of this table."""
        return {schema.schema_id: schema for schema in self.metadata.schemas}

    def spec(self) -> PartitionSpec:
        """Return the partition spec of this table."""
        return next(spec for spec in self.metadata.partition_specs if spec.spec_id == self.metadata.default_spec_id)

    def specs(self) -> Dict[int, PartitionSpec]:
        """Return a dict the partition specs this table."""
        return {spec.spec_id: spec for spec in self.metadata.partition_specs}

    def sort_order(self) -> SortOrder:
        """Return the sort order of this table."""
        return next(
            sort_order for sort_order in self.metadata.sort_orders if sort_order.order_id == self.metadata.default_sort_order_id
        )

    def sort_orders(self) -> Dict[int, SortOrder]:
        """Return a dict of the sort orders of this table."""
        return {sort_order.order_id: sort_order for sort_order in self.metadata.sort_orders}

    def last_partition_id(self) -> int:
        """Return the highest assigned partition field ID across all specs or 999 if only the unpartitioned spec exists."""
        if self.metadata.last_partition_id:
            return self.metadata.last_partition_id
        return PARTITION_FIELD_ID_START - 1

    @property
    def properties(self) -> Dict[str, str]:
        """Properties of the table."""
        return self.metadata.properties

    def location(self) -> str:
        """Return the table's base location."""
        return self.metadata.location

    def location_provider(self) -> LocationProvider:
        """Return the table's location provider."""
        return load_location_provider(table_location=self.metadata.location, table_properties=self.metadata.properties)

    @property
    def last_sequence_number(self) -> int:
        return self.metadata.last_sequence_number

    def current_snapshot(self) -> Optional[Snapshot]:
        """Get the current snapshot for this table, or None if there is no current snapshot."""
        if self.metadata.current_snapshot_id is not None:
            return self.snapshot_by_id(self.metadata.current_snapshot_id)
        return None

    def snapshots(self) -> List[Snapshot]:
        return self.metadata.snapshots

    def snapshot_by_id(self, snapshot_id: int) -> Optional[Snapshot]:
        """Get the snapshot of this table with the given id, or None if there is no matching snapshot."""
        return self.metadata.snapshot_by_id(snapshot_id)

    def snapshot_by_name(self, name: str) -> Optional[Snapshot]:
        """Return the snapshot referenced by the given name or null if no such reference exists."""
        if ref := self.metadata.refs.get(name):
            return self.snapshot_by_id(ref.snapshot_id)
        return None

    def snapshot_as_of_timestamp(self, timestamp_ms: int, inclusive: bool = True) -> Optional[Snapshot]:
        """Get the snapshot that was current as of or right before the given timestamp, or None if there is no matching snapshot.

        Args:
            timestamp_ms: Find snapshot that was current at/before this timestamp
            inclusive: Includes timestamp_ms in search when True. Excludes timestamp_ms when False
        """
        for log_entry in reversed(self.history()):
            if (inclusive and log_entry.timestamp_ms <= timestamp_ms) or log_entry.timestamp_ms < timestamp_ms:
                return self.snapshot_by_id(log_entry.snapshot_id)
        return None

    def history(self) -> List[SnapshotLogEntry]:
        """Get the snapshot history of this table."""
        return self.metadata.snapshot_log

    def manage_snapshots(self) -> ManageSnapshots:
        """
        Shorthand to run snapshot management operations like create branch, create tag, etc.

        Use table.manage_snapshots().<operation>().commit() to run a specific operation.
        Use table.manage_snapshots().<operation-one>().<operation-two>().commit() to run multiple operations.
        Pending changes are applied on commit.

        We can also use context managers to make more changes. For example,

        with table.manage_snapshots() as ms:
           ms.create_tag(snapshot_id1, "Tag_A").create_tag(snapshot_id2, "Tag_B")
        """
        return ManageSnapshots(transaction=Transaction(self, autocommit=True))

    def update_statistics(self) -> UpdateStatistics:
        """
        Shorthand to run statistics management operations like add statistics and remove statistics.

        Use table.update_statistics().<operation>().commit() to run a specific operation.
        Use table.update_statistics().<operation-one>().<operation-two>().commit() to run multiple operations.

        Pending changes are applied on commit.

        We can also use context managers to make more changes. For example:

        with table.update_statistics() as update:
            update.set_statistics(statistics_file=statistics_file)
            update.remove_statistics(snapshot_id=2)
        """
        return UpdateStatistics(transaction=Transaction(self, autocommit=True))

    def update_schema(self, allow_incompatible_changes: bool = False, case_sensitive: bool = True) -> UpdateSchema:
        """Create a new UpdateSchema to alter the columns of this table.

        Args:
            allow_incompatible_changes: If changes are allowed that might break downstream consumers.
            case_sensitive: If field names are case-sensitive.

        Returns:
            A new UpdateSchema.
        """
        return UpdateSchema(
            transaction=Transaction(self, autocommit=True),
            allow_incompatible_changes=allow_incompatible_changes,
            case_sensitive=case_sensitive,
            name_mapping=self.name_mapping(),
        )

    def name_mapping(self) -> Optional[NameMapping]:
        """Return the table's field-id NameMapping."""
        return self.metadata.name_mapping()

    def upsert(
        self,
        df: pa.Table,
        join_cols: Optional[List[str]] = None,
        when_matched_update_all: bool = True,
        when_not_matched_insert_all: bool = True,
        case_sensitive: bool = True,
        branch: Optional[str] = MAIN_BRANCH,
    ) -> UpsertResult:
        """Shorthand API for performing an upsert to an iceberg table.

        Args:

            df: The input dataframe to upsert with the table's data.
            join_cols: Columns to join on, if not provided, it will use the identifier-field-ids.
            when_matched_update_all: Bool indicating to update rows that are matched but require an update due to a value in a non-key column changing
            when_not_matched_insert_all: Bool indicating new rows to be inserted that do not match any existing rows in the table
            case_sensitive: Bool indicating if the match should be case-sensitive
            branch: Branch Reference to run the upsert operation

            To learn more about the identifier-field-ids: https://iceberg.apache.org/spec/#identifier-field-ids

                Example Use Cases:
                    Case 1: Both Parameters = True (Full Upsert)
                    Existing row found → Update it
                    New row found → Insert it

                    Case 2: when_matched_update_all = False, when_not_matched_insert_all = True
                    Existing row found → Do nothing (no updates)
                    New row found → Insert it

                    Case 3: when_matched_update_all = True, when_not_matched_insert_all = False
                    Existing row found → Update it
                    New row found → Do nothing (no inserts)

                    Case 4: Both Parameters = False (No Merge Effect)
                    Existing row found → Do nothing
                    New row found → Do nothing
                    (Function effectively does nothing)


        Returns:
            An UpsertResult class (contains details of rows updated and inserted)
        """
        with self.transaction() as tx:
            return tx.upsert(
                df=df,
                join_cols=join_cols,
                when_matched_update_all=when_matched_update_all,
                when_not_matched_insert_all=when_not_matched_insert_all,
                case_sensitive=case_sensitive,
                branch=branch,
            )

    def append(self, df: pa.Table, snapshot_properties: Dict[str, str] = EMPTY_DICT, branch: Optional[str] = MAIN_BRANCH) -> None:
        """
        Shorthand API for appending a PyArrow table to the table.

        Args:
            df: The Arrow dataframe that will be appended to overwrite the table
            snapshot_properties: Custom properties to be added to the snapshot summary
            branch: Branch Reference to run the append operation
        """
        with self.transaction() as tx:
            tx.append(df=df, snapshot_properties=snapshot_properties, branch=branch)

    def dynamic_partition_overwrite(
        self, df: pa.Table, snapshot_properties: Dict[str, str] = EMPTY_DICT, branch: Optional[str] = MAIN_BRANCH
    ) -> None:
        """Shorthand for dynamic overwriting the table with a PyArrow table.

        Old partitions are auto detected and replaced with data files created for input arrow table.
        Args:
            df: The Arrow dataframe that will be used to overwrite the table
            snapshot_properties: Custom properties to be added to the snapshot summary
            branch: Branch Reference to run the dynamic partition overwrite operation
        """
        with self.transaction() as tx:
            tx.dynamic_partition_overwrite(df=df, snapshot_properties=snapshot_properties, branch=branch)

    def overwrite(
        self,
        df: pa.Table,
        overwrite_filter: Union[BooleanExpression, str] = ALWAYS_TRUE,
        snapshot_properties: Dict[str, str] = EMPTY_DICT,
        case_sensitive: bool = True,
        branch: Optional[str] = MAIN_BRANCH,
    ) -> None:
        """
        Shorthand for overwriting the table with a PyArrow table.

        An overwrite may produce zero or more snapshots based on the operation:

            - DELETE: In case existing Parquet files can be dropped completely.
            - REPLACE: In case existing Parquet files need to be rewritten.
            - APPEND: In case new data is being inserted into the table.

        Args:
            df: The Arrow dataframe that will be used to overwrite the table
            overwrite_filter: ALWAYS_TRUE when you overwrite all the data,
                              or a boolean expression in case of a partial overwrite
            snapshot_properties: Custom properties to be added to the snapshot summary
            case_sensitive: A bool determine if the provided `overwrite_filter` is case-sensitive
            branch: Branch Reference to run the overwrite operation
        """
        with self.transaction() as tx:
            tx.overwrite(
                df=df,
                overwrite_filter=overwrite_filter,
                case_sensitive=case_sensitive,
                snapshot_properties=snapshot_properties,
                branch=branch,
            )

    def delete(
        self,
        delete_filter: Union[BooleanExpression, str] = ALWAYS_TRUE,
        snapshot_properties: Dict[str, str] = EMPTY_DICT,
        case_sensitive: bool = True,
        branch: Optional[str] = MAIN_BRANCH,
    ) -> None:
        """
        Shorthand for deleting rows from the table.

        Args:
            delete_filter: The predicate that used to remove rows
            snapshot_properties: Custom properties to be added to the snapshot summary
            case_sensitive: A bool determine if the provided `delete_filter` is case-sensitive
            branch: Branch Reference to run the delete operation
        """
        with self.transaction() as tx:
            tx.delete(
                delete_filter=delete_filter, case_sensitive=case_sensitive, snapshot_properties=snapshot_properties, branch=branch
            )

    def add_files(
        self, file_paths: List[str], snapshot_properties: Dict[str, str] = EMPTY_DICT, check_duplicate_files: bool = True
    ) -> None:
        """
        Shorthand API for adding files as data files to the table.

        Args:
            file_paths: The list of full file paths to be added as data files to the table

        Raises:
            FileNotFoundError: If the file does not exist.
        """
        with self.transaction() as tx:
            tx.add_files(
                file_paths=file_paths, snapshot_properties=snapshot_properties, check_duplicate_files=check_duplicate_files
            )

    def update_spec(self, case_sensitive: bool = True) -> UpdateSpec:
        return UpdateSpec(Transaction(self, autocommit=True), case_sensitive=case_sensitive)

    def refs(self) -> Dict[str, SnapshotRef]:
        """Return the snapshot references in the table."""
        return self.metadata.refs

    def _do_commit(self, updates: Tuple[TableUpdate, ...], requirements: Tuple[TableRequirement, ...]) -> None:
        response = self.catalog.commit_table(self, requirements, updates)

        # https://github.com/apache/iceberg/blob/f6faa58/core/src/main/java/org/apache/iceberg/CatalogUtil.java#L527
        # delete old metadata if METADATA_DELETE_AFTER_COMMIT_ENABLED is set to true and uses
        # TableProperties.METADATA_PREVIOUS_VERSIONS_MAX to determine how many previous versions to keep -
        # everything else will be removed.
        try:
            self.catalog._delete_old_metadata(self.io, self.metadata, response.metadata)
        except Exception as e:
            warnings.warn(f"Failed to delete old metadata after commit: {e}")

        self.metadata = response.metadata
        self.metadata_location = response.metadata_location

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the Table class."""
        return (
            self.name() == other.name() and self.metadata == other.metadata and self.metadata_location == other.metadata_location
            if isinstance(other, Table)
            else False
        )

    def __repr__(self) -> str:
        """Return the string representation of the Table class."""
        table_name = self.catalog.table_name_from(self._identifier)
        schema_str = ",\n  ".join(str(column) for column in self.schema().columns if self.schema())
        partition_str = f"partition by: [{', '.join(field.name for field in self.spec().fields if self.spec())}]"
        sort_order_str = f"sort order: [{', '.join(str(field) for field in self.sort_order().fields if self.sort_order())}]"
        snapshot_str = f"snapshot: {str(self.current_snapshot()) if self.current_snapshot() else 'null'}"
        result_str = f"{table_name}(\n  {schema_str}\n),\n{partition_str},\n{sort_order_str},\n{snapshot_str}"
        return result_str

    def to_daft(self) -> daft.DataFrame:
        """Read a Daft DataFrame lazily from this Iceberg table.

        Returns:
            daft.DataFrame: Unmaterialized Daft Dataframe created from the Iceberg table
        """
        import daft

        return daft.read_iceberg(self)

    def to_bodo(self) -> bd.DataFrame:
        """Read a bodo DataFrame lazily from this Iceberg table.

        Returns:
            bd.DataFrame: Unmaterialized Bodo Dataframe created from the Iceberg table
        """
        import bodo.pandas as bd

        return bd.read_iceberg_table(self)

    def to_polars(self) -> pl.LazyFrame:
        """Lazily read from this Apache Iceberg table.

        Returns:
            pl.LazyFrame: Unmaterialized Polars LazyFrame created from the Iceberg table
        """
        import polars as pl

        return pl.scan_iceberg(self)

    def __datafusion_table_provider__(self) -> "IcebergDataFusionTable":
        """Return the DataFusion table provider PyCapsule interface.

        To support DataFusion features such as push down filtering, this function will return a PyCapsule
        interface that conforms to the FFI Table Provider required by DataFusion. From an end user perspective
        you should not need to call this function directly. Instead you can use ``register_table_provider`` in
        the DataFusion SessionContext.

        Returns:
            A PyCapsule DataFusion TableProvider interface.

        Example:
            ```python
            from datafusion import SessionContext
            from pyiceberg.catalog import load_catalog
            import pyarrow as pa
            catalog = load_catalog("catalog", type="in-memory")
            catalog.create_namespace_if_not_exists("default")
            data = pa.table({"x": [1, 2, 3], "y": [4, 5, 6]})
            iceberg_table = catalog.create_table("default.test", schema=data.schema)
            iceberg_table.append(data)
            ctx = SessionContext()
            ctx.register_table_provider("test", iceberg_table)
            ctx.table("test").show()
            ```
            Results in
            ```
            DataFrame()
            +---+---+
            | x | y |
            +---+---+
            | 1 | 4 |
            | 2 | 5 |
            | 3 | 6 |
            +---+---+
            ```
        """
        from pyiceberg_core.datafusion import IcebergDataFusionTable

        return IcebergDataFusionTable(
            identifier=self.name(),
            metadata_location=self.metadata_location,
            file_io_properties=self.io.properties,
        ).__datafusion_table_provider__()


class StaticTable(Table):
    """Load a table directly from a metadata file (i.e., without using a catalog)."""

    def refresh(self) -> Table:
        """Refresh the current table metadata."""
        raise NotImplementedError("To be implemented")

    @classmethod
    def _metadata_location_from_version_hint(cls, metadata_location: str, properties: Properties = EMPTY_DICT) -> str:
        version_hint_location = os.path.join(metadata_location, "metadata", "version-hint.text")
        io = load_file_io(properties=properties, location=version_hint_location)
        file = io.new_input(version_hint_location)

        with file.open() as stream:
            content = stream.read().decode("utf-8")

        if content.endswith(".metadata.json"):
            return os.path.join(metadata_location, "metadata", content)
        elif content.isnumeric():
            return os.path.join(metadata_location, "metadata", "v%s.metadata.json").format(content)
        else:
            return os.path.join(metadata_location, "metadata", "%s.metadata.json").format(content)

    @classmethod
    def from_metadata(cls, metadata_location: str, properties: Properties = EMPTY_DICT) -> StaticTable:
        if not metadata_location.endswith(".metadata.json"):
            metadata_location = StaticTable._metadata_location_from_version_hint(metadata_location, properties)

        io = load_file_io(properties=properties, location=metadata_location)
        file = io.new_input(metadata_location)

        from pyiceberg.serializers import FromInputFile

        metadata = FromInputFile.table_metadata(file)

        from pyiceberg.catalog.noop import NoopCatalog

        return cls(
            identifier=("static-table", metadata_location),
            metadata_location=metadata_location,
            metadata=metadata,
            io=load_file_io({**properties, **metadata.properties}),
            catalog=NoopCatalog("static-table"),
        )


class StagedTable(Table):
    def refresh(self) -> Table:
        raise ValueError("Cannot refresh a staged table")

    def scan(
        self,
        row_filter: Union[str, BooleanExpression] = ALWAYS_TRUE,
        selected_fields: Tuple[str, ...] = ("*",),
        case_sensitive: bool = True,
        snapshot_id: Optional[int] = None,
        options: Properties = EMPTY_DICT,
        limit: Optional[int] = None,
    ) -> DataScan:
        raise ValueError("Cannot scan a staged table")

    def to_daft(self) -> daft.DataFrame:
        raise ValueError("Cannot convert a staged table to a Daft DataFrame")


def _parse_row_filter(expr: Union[str, BooleanExpression]) -> BooleanExpression:
    """Accept an expression in the form of a BooleanExpression or a string.

    In the case of a string, it will be converted into a unbound BooleanExpression.

    Args:
        expr: Expression as a BooleanExpression or a string.

    Returns: An unbound BooleanExpression.
    """
    return parser.parse(expr) if isinstance(expr, str) else expr


S = TypeVar("S", bound="TableScan", covariant=True)


class TableScan(ABC):
    table_metadata: TableMetadata
    io: FileIO
    row_filter: BooleanExpression
    selected_fields: Tuple[str, ...]
    case_sensitive: bool
    snapshot_id: Optional[int]
    options: Properties
    limit: Optional[int]

    def __init__(
        self,
        table_metadata: TableMetadata,
        io: FileIO,
        row_filter: Union[str, BooleanExpression] = ALWAYS_TRUE,
        selected_fields: Tuple[str, ...] = ("*",),
        case_sensitive: bool = True,
        snapshot_id: Optional[int] = None,
        options: Properties = EMPTY_DICT,
        limit: Optional[int] = None,
    ):
        self.table_metadata = table_metadata
        self.io = io
        self.row_filter = _parse_row_filter(row_filter)
        self.selected_fields = selected_fields
        self.case_sensitive = case_sensitive
        self.snapshot_id = snapshot_id
        self.options = options
        self.limit = limit

    def snapshot(self) -> Optional[Snapshot]:
        if self.snapshot_id:
            return self.table_metadata.snapshot_by_id(self.snapshot_id)
        return self.table_metadata.current_snapshot()

    def projection(self) -> Schema:
        current_schema = self.table_metadata.schema()
        if self.snapshot_id is not None:
            snapshot = self.table_metadata.snapshot_by_id(self.snapshot_id)
            if snapshot is not None:
                if snapshot.schema_id is not None:
                    try:
                        current_schema = next(
                            schema for schema in self.table_metadata.schemas if schema.schema_id == snapshot.schema_id
                        )
                    except StopIteration:
                        warnings.warn(f"Metadata does not contain schema with id: {snapshot.schema_id}")
            else:
                raise ValueError(f"Snapshot not found: {self.snapshot_id}")

        if "*" in self.selected_fields:
            return current_schema

        return current_schema.select(*self.selected_fields, case_sensitive=self.case_sensitive)

    @abstractmethod
    def plan_files(self) -> Iterable[ScanTask]: ...

    @abstractmethod
    def to_arrow(self) -> pa.Table: ...

    @abstractmethod
    def to_pandas(self, **kwargs: Any) -> pd.DataFrame: ...

    @abstractmethod
    def to_polars(self) -> pl.DataFrame: ...

    def update(self: S, **overrides: Any) -> S:
        """Create a copy of this table scan with updated fields."""
        from inspect import signature

        # Extract those attributes that are constructor parameters. We don't use self.__dict__ as the kwargs to the
        # constructors because it may contain additional attributes that are not part of the constructor signature.
        params = signature(type(self).__init__).parameters.keys() - {"self"}  # Skip "self" parameter
        kwargs = {param: getattr(self, param) for param in params}  # Assume parameters are attributes

        return type(self)(**{**kwargs, **overrides})

    def use_ref(self: S, name: str) -> S:
        if self.snapshot_id:
            raise ValueError(f"Cannot override ref, already set snapshot id={self.snapshot_id}")
        if snapshot := self.table_metadata.snapshot_by_name(name):
            return self.update(snapshot_id=snapshot.snapshot_id)

        raise ValueError(f"Cannot scan unknown ref={name}")

    def select(self: S, *field_names: str) -> S:
        if "*" in self.selected_fields:
            return self.update(selected_fields=field_names)
        return self.update(selected_fields=tuple(set(self.selected_fields).intersection(set(field_names))))

    def filter(self: S, expr: Union[str, BooleanExpression]) -> S:
        return self.update(row_filter=And(self.row_filter, _parse_row_filter(expr)))

    def with_case_sensitive(self: S, case_sensitive: bool = True) -> S:
        return self.update(case_sensitive=case_sensitive)

    @abstractmethod
    def count(self) -> int: ...


class ScanTask(ABC):
    pass


@dataclass(init=False)
class FileScanTask(ScanTask):
    """Task representing a data file and its corresponding delete files."""

    file: DataFile
    delete_files: Set[DataFile]
    start: int
    length: int
    residual: BooleanExpression

    def __init__(
        self,
        data_file: DataFile,
        delete_files: Optional[Set[DataFile]] = None,
        start: Optional[int] = None,
        length: Optional[int] = None,
        residual: BooleanExpression = ALWAYS_TRUE,
    ) -> None:
        self.file = data_file
        self.delete_files = delete_files or set()
        self.start = start or 0
        self.length = length or data_file.file_size_in_bytes
        self.residual = residual


def _open_manifest(
    io: FileIO,
    manifest: ManifestFile,
    partition_filter: Callable[[DataFile], bool],
    metrics_evaluator: Callable[[DataFile], bool],
) -> List[ManifestEntry]:
    """Open a manifest file and return matching manifest entries.

    Returns:
        A list of ManifestEntry that matches the provided filters.
    """
    return [
        manifest_entry
        for manifest_entry in manifest.fetch_manifest_entry(io, discard_deleted=True)
        if partition_filter(manifest_entry.data_file) and metrics_evaluator(manifest_entry.data_file)
    ]


def _min_sequence_number(manifests: List[ManifestFile]) -> int:
    try:
        return min(
            manifest.min_sequence_number or INITIAL_SEQUENCE_NUMBER
            for manifest in manifests
            if manifest.content == ManifestContent.DATA
        )
    except ValueError:
        # In case of an empty iterator
        return INITIAL_SEQUENCE_NUMBER


def _match_deletes_to_data_file(data_entry: ManifestEntry, positional_delete_entries: SortedList[ManifestEntry]) -> Set[DataFile]:
    """Check if the delete file is relevant for the data file.

    Using the column metrics to see if the filename is in the lower and upper bound.

    Args:
        data_entry (ManifestEntry): The manifest entry path of the datafile.
        positional_delete_entries (List[ManifestEntry]): All the candidate positional deletes manifest entries.

    Returns:
        A set of files that are relevant for the data file.
    """
    relevant_entries = positional_delete_entries[positional_delete_entries.bisect_right(data_entry) :]

    if len(relevant_entries) > 0:
        evaluator = _InclusiveMetricsEvaluator(POSITIONAL_DELETE_SCHEMA, EqualTo("file_path", data_entry.data_file.file_path))
        return {
            positional_delete_entry.data_file
            for positional_delete_entry in relevant_entries
            if evaluator.eval(positional_delete_entry.data_file)
        }
    else:
        return set()


class DataScan(TableScan):
    def _build_partition_projection(self, spec_id: int) -> BooleanExpression:
        project = inclusive_projection(self.table_metadata.schema(), self.table_metadata.specs()[spec_id], self.case_sensitive)
        return project(self.row_filter)

    @cached_property
    def partition_filters(self) -> KeyDefaultDict[int, BooleanExpression]:
        return KeyDefaultDict(self._build_partition_projection)

    def _build_manifest_evaluator(self, spec_id: int) -> Callable[[ManifestFile], bool]:
        spec = self.table_metadata.specs()[spec_id]
        return manifest_evaluator(spec, self.table_metadata.schema(), self.partition_filters[spec_id], self.case_sensitive)

    def _build_partition_evaluator(self, spec_id: int) -> Callable[[DataFile], bool]:
        spec = self.table_metadata.specs()[spec_id]
        partition_type = spec.partition_type(self.table_metadata.schema())
        partition_schema = Schema(*partition_type.fields)
        partition_expr = self.partition_filters[spec_id]

        # The lambda created here is run in multiple threads.
        # So we avoid creating _EvaluatorExpression methods bound to a single
        # shared instance across multiple threads.
        return lambda data_file: expression_evaluator(partition_schema, partition_expr, self.case_sensitive)(data_file.partition)

    def _build_metrics_evaluator(self) -> Callable[[DataFile], bool]:
        schema = self.table_metadata.schema()
        include_empty_files = strtobool(self.options.get("include_empty_files", "false"))

        # The lambda created here is run in multiple threads.
        # So we avoid creating _InclusiveMetricsEvaluator methods bound to a single
        # shared instance across multiple threads.
        return lambda data_file: _InclusiveMetricsEvaluator(
            schema,
            self.row_filter,
            self.case_sensitive,
            include_empty_files,
        ).eval(data_file)

    def _build_residual_evaluator(self, spec_id: int) -> Callable[[DataFile], ResidualEvaluator]:
        spec = self.table_metadata.specs()[spec_id]

        from pyiceberg.expressions.visitors import residual_evaluator_of

        # The lambda created here is run in multiple threads.
        # So we avoid creating _EvaluatorExpression methods bound to a single
        # shared instance across multiple threads.
        return lambda datafile: (
            residual_evaluator_of(
                spec=spec,
                expr=self.row_filter,
                case_sensitive=self.case_sensitive,
                schema=self.table_metadata.schema(),
            )
        )

    @staticmethod
    def _check_sequence_number(min_sequence_number: int, manifest: ManifestFile) -> bool:
        """Ensure that no manifests are loaded that contain deletes that are older than the data.

        Args:
            min_sequence_number (int): The minimal sequence number.
            manifest (ManifestFile): A ManifestFile that can be either data or deletes.

        Returns:
            Boolean indicating if it is either a data file, or a relevant delete file.
        """
        return manifest.content == ManifestContent.DATA or (
            # Not interested in deletes that are older than the data
            manifest.content == ManifestContent.DELETES
            and (manifest.sequence_number or INITIAL_SEQUENCE_NUMBER) >= min_sequence_number
        )

    def plan_files(self) -> Iterable[FileScanTask]:
        """Plans the relevant files by filtering on the PartitionSpecs.

        Returns:
            List of FileScanTasks that contain both data and delete files.
        """
        snapshot = self.snapshot()
        if not snapshot:
            return iter([])

        # step 1: filter manifests using partition summaries
        # the filter depends on the partition spec used to write the manifest file, so create a cache of filters for each spec id

        manifest_evaluators: Dict[int, Callable[[ManifestFile], bool]] = KeyDefaultDict(self._build_manifest_evaluator)

        residual_evaluators: Dict[int, Callable[[DataFile], ResidualEvaluator]] = KeyDefaultDict(self._build_residual_evaluator)

        manifests = [
            manifest_file
            for manifest_file in snapshot.manifests(self.io)
            if manifest_evaluators[manifest_file.partition_spec_id](manifest_file)
        ]

        # step 2: filter the data files in each manifest
        # this filter depends on the partition spec used to write the manifest file

        partition_evaluators: Dict[int, Callable[[DataFile], bool]] = KeyDefaultDict(self._build_partition_evaluator)

        min_sequence_number = _min_sequence_number(manifests)

        data_entries: List[ManifestEntry] = []
        positional_delete_entries = SortedList(key=lambda entry: entry.sequence_number or INITIAL_SEQUENCE_NUMBER)

        executor = ExecutorFactory.get_or_create()
        for manifest_entry in chain(
            *executor.map(
                lambda args: _open_manifest(*args),
                [
                    (
                        self.io,
                        manifest,
                        partition_evaluators[manifest.partition_spec_id],
                        self._build_metrics_evaluator(),
                    )
                    for manifest in manifests
                    if self._check_sequence_number(min_sequence_number, manifest)
                ],
            )
        ):
            data_file = manifest_entry.data_file
            if data_file.content == DataFileContent.DATA:
                data_entries.append(manifest_entry)
            elif data_file.content == DataFileContent.POSITION_DELETES:
                positional_delete_entries.add(manifest_entry)
            elif data_file.content == DataFileContent.EQUALITY_DELETES:
                raise ValueError("PyIceberg does not yet support equality deletes: https://github.com/apache/iceberg/issues/6568")
            else:
                raise ValueError(f"Unknown DataFileContent ({data_file.content}): {manifest_entry}")

        return [
            FileScanTask(
                data_entry.data_file,
                delete_files=_match_deletes_to_data_file(
                    data_entry,
                    positional_delete_entries,
                ),
                residual=residual_evaluators[data_entry.data_file.spec_id](data_entry.data_file).residual_for(
                    data_entry.data_file.partition
                ),
            )
            for data_entry in data_entries
        ]

    def to_arrow(self) -> pa.Table:
        """Read an Arrow table eagerly from this DataScan.

        All rows will be loaded into memory at once.

        Returns:
            pa.Table: Materialized Arrow Table from the Iceberg table's DataScan
        """
        from pyiceberg.io.pyarrow import ArrowScan

        return ArrowScan(
            self.table_metadata, self.io, self.projection(), self.row_filter, self.case_sensitive, self.limit
        ).to_table(self.plan_files())

    def to_arrow_batch_reader(self) -> pa.RecordBatchReader:
        """Return an Arrow RecordBatchReader from this DataScan.

        For large results, using a RecordBatchReader requires less memory than
        loading an Arrow Table for the same DataScan, because a RecordBatch
        is read one at a time.

        Returns:
            pa.RecordBatchReader: Arrow RecordBatchReader from the Iceberg table's DataScan
                which can be used to read a stream of record batches one by one.
        """
        import pyarrow as pa

        from pyiceberg.io.pyarrow import ArrowScan, schema_to_pyarrow

        target_schema = schema_to_pyarrow(self.projection())
        batches = ArrowScan(
            self.table_metadata, self.io, self.projection(), self.row_filter, self.case_sensitive, self.limit
        ).to_record_batches(self.plan_files())

        return pa.RecordBatchReader.from_batches(
            target_schema,
            batches,
        ).cast(target_schema)

    def to_pandas(self, **kwargs: Any) -> pd.DataFrame:
        """Read a Pandas DataFrame eagerly from this Iceberg table.

        Returns:
            pd.DataFrame: Materialized Pandas Dataframe from the Iceberg table
        """
        return self.to_arrow().to_pandas(**kwargs)

    def to_duckdb(self, table_name: str, connection: Optional[DuckDBPyConnection] = None) -> DuckDBPyConnection:
        """Shorthand for loading the Iceberg Table in DuckDB.

        Returns:
            DuckDBPyConnection: In memory DuckDB connection with the Iceberg table.
        """
        import duckdb

        con = connection or duckdb.connect(database=":memory:")
        con.register(table_name, self.to_arrow())

        return con

    def to_ray(self) -> ray.data.dataset.Dataset:
        """Read a Ray Dataset eagerly from this Iceberg table.

        Returns:
            ray.data.dataset.Dataset: Materialized Ray Dataset from the Iceberg table
        """
        import ray

        return ray.data.from_arrow(self.to_arrow())

    def to_polars(self) -> pl.DataFrame:
        """Read a Polars DataFrame from this Iceberg table.

        Returns:
            pl.DataFrame: Materialized Polars Dataframe from the Iceberg table
        """
        import polars as pl

        result = pl.from_arrow(self.to_arrow())
        if isinstance(result, pl.Series):
            result = result.to_frame()

        return result

    def count(self) -> int:
        from pyiceberg.io.pyarrow import ArrowScan

        # Usage: Calculates the total number of records in a Scan that haven't had positional deletes.
        res = 0
        # every task is a FileScanTask
        tasks = self.plan_files()

        for task in tasks:
            # task.residual is a Boolean Expression if the filter condition is fully satisfied by the
            # partition value and task.delete_files represents that positional delete haven't been merged yet
            # hence those files have to read as a pyarrow table applying the filter and deletes
            if task.residual == AlwaysTrue() and len(task.delete_files) == 0:
                # Every File has a metadata stat that stores the file record count
                res += task.file.record_count
            else:
                arrow_scan = ArrowScan(
                    table_metadata=self.table_metadata,
                    io=self.io,
                    projected_schema=self.projection(),
                    row_filter=self.row_filter,
                    case_sensitive=self.case_sensitive,
                )
                tbl = arrow_scan.to_table([task])
                res += len(tbl)
        return res


@dataclass(frozen=True)
class WriteTask:
    """Task with the parameters for writing a DataFile."""

    write_uuid: uuid.UUID
    task_id: int
    schema: Schema
    record_batches: List[pa.RecordBatch]
    sort_order_id: Optional[int] = None
    partition_key: Optional[PartitionKey] = None

    def generate_data_file_filename(self, extension: str) -> str:
        # Mimics the behavior in the Java API:
        # https://github.com/apache/iceberg/blob/a582968975dd30ff4917fbbe999f1be903efac02/core/src/main/java/org/apache/iceberg/io/OutputFileFactory.java#L92-L101
        return f"00000-{self.task_id}-{self.write_uuid}.{extension}"


def _parquet_files_to_data_files(table_metadata: TableMetadata, file_paths: List[str], io: FileIO) -> Iterable[DataFile]:
    """Convert a list files into DataFiles.

    Returns:
        An iterable that supplies DataFiles that describe the parquet files.
    """
    from pyiceberg.io.pyarrow import parquet_file_to_data_file

    executor = ExecutorFactory.get_or_create()
    futures = [executor.submit(parquet_file_to_data_file, io, table_metadata, file_path) for file_path in file_paths]

    return [f.result() for f in futures if f.result()]
