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

from datetime import date, datetime, time
from decimal import Decimal
from typing import Annotated, Generic, Literal, TypeAlias, TypeVar
from uuid import UUID

from pydantic import Field, model_validator

from pyiceberg.catalog.rest.response import ErrorResponseMessage
from pyiceberg.expressions import BooleanExpression, SerializableBooleanExpression
from pyiceberg.manifest import FileFormat
from pyiceberg.typedef import IcebergBaseModel

# Primitive types that can appear in partition values and bounds
PrimitiveTypeValue: TypeAlias = bool | int | float | str | Decimal | UUID | date | time | datetime | bytes

V = TypeVar("V")


class KeyValueMap(IcebergBaseModel, Generic[V]):
    """Map serialized as parallel key/value arrays for column statistics."""

    keys: list[int] = Field(default_factory=list)
    values: list[V] = Field(default_factory=list)

    @model_validator(mode="after")
    def _validate_lengths_match(self) -> KeyValueMap[V]:
        if len(self.keys) != len(self.values):
            raise ValueError(f"keys and values must have same length: {len(self.keys)} != {len(self.values)}")
        return self

    def to_dict(self) -> dict[int, V]:
        """Convert to dictionary mapping field ID to value."""
        return dict(zip(self.keys, self.values, strict=True))


class CountMap(KeyValueMap[int]):
    """Map of field IDs to counts."""


class ValueMap(KeyValueMap[PrimitiveTypeValue]):
    """Map of field IDs to primitive values (for lower/upper bounds)."""


class StorageCredential(IcebergBaseModel):
    """Storage credential for accessing content files."""

    prefix: str = Field(description="Storage location prefix this credential applies to")
    config: dict[str, str] = Field(default_factory=dict)


class RESTContentFile(IcebergBaseModel):
    """Base model for data and delete files from REST API."""

    spec_id: int = Field(alias="spec-id")
    partition: list[PrimitiveTypeValue] = Field(default_factory=list)
    content: Literal["data", "position-deletes", "equality-deletes"]
    file_path: str = Field(alias="file-path")
    file_format: FileFormat = Field(alias="file-format")
    file_size_in_bytes: int = Field(alias="file-size-in-bytes")
    record_count: int = Field(alias="record-count")
    key_metadata: str | None = Field(alias="key-metadata", default=None)
    split_offsets: list[int] | None = Field(alias="split-offsets", default=None)
    sort_order_id: int | None = Field(alias="sort-order-id", default=None)


class RESTDataFile(RESTContentFile):
    """Data file from REST API."""

    content: Literal["data"] = Field(default="data")
    first_row_id: int | None = Field(alias="first-row-id", default=None)
    column_sizes: CountMap | None = Field(alias="column-sizes", default=None)
    value_counts: CountMap | None = Field(alias="value-counts", default=None)
    null_value_counts: CountMap | None = Field(alias="null-value-counts", default=None)
    nan_value_counts: CountMap | None = Field(alias="nan-value-counts", default=None)
    lower_bounds: ValueMap | None = Field(alias="lower-bounds", default=None)
    upper_bounds: ValueMap | None = Field(alias="upper-bounds", default=None)


class RESTPositionDeleteFile(RESTContentFile):
    """Position delete file from REST API."""

    content: Literal["position-deletes"] = Field(default="position-deletes")
    referenced_data_file: str | None = Field(alias="referenced-data-file", default=None)
    content_offset: int | None = Field(alias="content-offset", default=None)
    content_size_in_bytes: int | None = Field(alias="content-size-in-bytes", default=None)


class RESTEqualityDeleteFile(RESTContentFile):
    """Equality delete file from REST API."""

    content: Literal["equality-deletes"] = Field(default="equality-deletes")
    equality_ids: list[int] | None = Field(alias="equality-ids", default=None)


# Discriminated union for delete files
RESTDeleteFile = Annotated[
    RESTPositionDeleteFile | RESTEqualityDeleteFile,
    Field(discriminator="content"),
]


class RESTFileScanTask(IcebergBaseModel):
    """A file scan task from the REST server."""

    data_file: RESTDataFile = Field(alias="data-file")
    delete_file_references: list[int] | None = Field(alias="delete-file-references", default=None)
    residual_filter: BooleanExpression | None = Field(alias="residual-filter", default=None)


class ScanTasks(IcebergBaseModel):
    """Container for scan tasks returned by the server."""

    delete_files: list[RESTDeleteFile] = Field(alias="delete-files", default_factory=list)
    file_scan_tasks: list[RESTFileScanTask] = Field(alias="file-scan-tasks", default_factory=list)
    plan_tasks: list[str] = Field(alias="plan-tasks", default_factory=list)

    @model_validator(mode="after")
    def _validate_delete_file_references(self) -> ScanTasks:
        # validate delete file references are in bounds
        max_idx = len(self.delete_files) - 1
        for task in self.file_scan_tasks:
            for idx in task.delete_file_references or []:
                if idx < 0 or idx > max_idx:
                    raise ValueError(f"Invalid delete file reference: {idx} (valid range: 0-{max_idx})")

        if self.delete_files and not self.file_scan_tasks:
            raise ValueError("Invalid response: deleteFiles should only be returned with fileScanTasks that reference them")

        return self


class PlanCompleted(ScanTasks):
    """Completed scan plan result."""

    status: Literal["completed"] = "completed"
    plan_id: str | None = Field(alias="plan-id", default=None)
    storage_credentials: list[StorageCredential] | None = Field(alias="storage-credentials", default=None)


class PlanSubmitted(IcebergBaseModel):
    """Scan plan submitted, poll for completion."""

    status: Literal["submitted"] = "submitted"
    plan_id: str | None = Field(alias="plan-id", default=None)


class PlanCancelled(IcebergBaseModel):
    """Planning was cancelled."""

    status: Literal["cancelled"] = "cancelled"


class PlanFailed(IcebergBaseModel):
    """Planning failed with error."""

    status: Literal["failed"] = "failed"
    error: ErrorResponseMessage


PlanningResponse = Annotated[
    PlanCompleted | PlanSubmitted | PlanCancelled | PlanFailed,
    Field(discriminator="status"),
]


class PlanTableScanRequest(IcebergBaseModel):
    """Request body for planning a REST scan."""

    snapshot_id: int | None = Field(alias="snapshot-id", default=None)
    select: list[str] | None = Field(default=None)
    filter: SerializableBooleanExpression | None = Field(default=None)
    case_sensitive: bool = Field(alias="case-sensitive", default=True)
    use_snapshot_schema: bool = Field(alias="use-snapshot-schema", default=False)
    start_snapshot_id: int | None = Field(alias="start-snapshot-id", default=None)
    end_snapshot_id: int | None = Field(alias="end-snapshot-id", default=None)
    stats_fields: list[str] | None = Field(alias="stats-fields", default=None)
    min_rows_requested: int | None = Field(alias="min-rows-requested", default=None)

    @model_validator(mode="after")
    def _validate_snapshot_fields(self) -> PlanTableScanRequest:
        if self.start_snapshot_id is not None and self.end_snapshot_id is None:
            raise ValueError("end-snapshot-id is required when start-snapshot-id is specified")
        if self.snapshot_id is not None and self.start_snapshot_id is not None:
            raise ValueError("Cannot specify both snapshot-id and start-snapshot-id")
        return self


class FetchScanTasksRequest(IcebergBaseModel):
    """Request body for fetching scan tasks endpoint."""

    plan_task: str = Field(alias="plan-task")
