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
from enum import Enum
from typing import Annotated

from pydantic import Field, model_validator

from pyiceberg.exceptions import ValidationError
from pyiceberg.typedef import IcebergBaseModel

MAIN_BRANCH = "main"


class SnapshotRefType(str, Enum):
    BRANCH = "branch"
    TAG = "tag"

    def __repr__(self) -> str:
        """Return the string representation of the SnapshotRefType class."""
        return f"SnapshotRefType.{self.name}"

    def __str__(self) -> str:
        """Return the string representation of the SnapshotRefType class."""
        return self.value


class SnapshotRef(IcebergBaseModel):
    snapshot_id: int = Field(alias="snapshot-id")
    snapshot_ref_type: SnapshotRefType = Field(alias="type")
    min_snapshots_to_keep: Annotated[int | None, Field(alias="min-snapshots-to-keep", default=None, gt=0)]
    max_snapshot_age_ms: Annotated[int | None, Field(alias="max-snapshot-age-ms", default=None, gt=0)]
    max_ref_age_ms: Annotated[int | None, Field(alias="max-ref-age-ms", default=None, gt=0)]

    @model_validator(mode="after")
    def check_min_snapshots_to_keep(self) -> "SnapshotRef":
        if self.min_snapshots_to_keep is not None and self.snapshot_ref_type == SnapshotRefType.TAG:
            raise ValidationError("Tags do not support setting minSnapshotsToKeep")
        return self

    @model_validator(mode="after")
    def check_max_snapshot_age_ms(self) -> "SnapshotRef":
        if self.max_snapshot_age_ms is not None and self.snapshot_ref_type == SnapshotRefType.TAG:
            raise ValidationError("Tags do not support setting maxSnapshotAgeMs")
        return self
