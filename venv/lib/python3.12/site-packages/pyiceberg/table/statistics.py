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
from typing import Literal

from pydantic import Field

from pyiceberg.typedef import IcebergBaseModel


class BlobMetadata(IcebergBaseModel):
    type: Literal["apache-datasketches-theta-v1", "deletion-vector-v1"]
    snapshot_id: int = Field(alias="snapshot-id")
    sequence_number: int = Field(alias="sequence-number")
    fields: list[int]
    properties: dict[str, str] | None = None


class StatisticsCommonFields(IcebergBaseModel):
    """Common fields between table and partition statistics structs found on metadata."""

    snapshot_id: int = Field(alias="snapshot-id")
    statistics_path: str = Field(alias="statistics-path")
    file_size_in_bytes: int = Field(alias="file-size-in-bytes")


class StatisticsFile(StatisticsCommonFields):
    file_footer_size_in_bytes: int = Field(alias="file-footer-size-in-bytes")
    key_metadata: str | None = Field(alias="key-metadata", default=None)
    blob_metadata: list[BlobMetadata] = Field(alias="blob-metadata")


class PartitionStatisticsFile(StatisticsCommonFields):
    pass


def filter_statistics_by_snapshot_id(
    statistics: list[StatisticsFile | PartitionStatisticsFile],
    reject_snapshot_id: int,
) -> list[StatisticsFile | PartitionStatisticsFile]:
    return [stat for stat in statistics if stat.snapshot_id != reject_snapshot_id]
