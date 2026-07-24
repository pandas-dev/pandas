# -*- coding: utf-8 -*-
# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from __future__ import annotations

from typing import MutableMapping, MutableSequence

import google.protobuf.duration_pb2 as duration_pb2  # type: ignore
import google.protobuf.field_mask_pb2 as field_mask_pb2  # type: ignore
import google.protobuf.timestamp_pb2 as timestamp_pb2  # type: ignore
import google.rpc.status_pb2 as status_pb2  # type: ignore
import google.type.interval_pb2 as interval_pb2  # type: ignore
import proto  # type: ignore

__protobuf__ = proto.module(
    package="google.storage.control.v2",
    manifest={
        "FindingType",
        "FindingCategory",
        "FindingSeverity",
        "PendingRenameInfo",
        "Folder",
        "GetFolderRequest",
        "CreateFolderRequest",
        "DeleteFolderRequest",
        "ListFoldersRequest",
        "ListFoldersResponse",
        "RenameFolderRequest",
        "DeleteFolderRecursiveRequest",
        "CommonLongRunningOperationMetadata",
        "RenameFolderMetadata",
        "DeleteFolderRecursiveMetadata",
        "StorageLayout",
        "GetStorageLayoutRequest",
        "ManagedFolder",
        "GetManagedFolderRequest",
        "CreateManagedFolderRequest",
        "DeleteManagedFolderRequest",
        "ListManagedFoldersRequest",
        "ListManagedFoldersResponse",
        "CreateAnywhereCacheMetadata",
        "UpdateAnywhereCacheMetadata",
        "AnywhereCache",
        "CreateAnywhereCacheRequest",
        "UpdateAnywhereCacheRequest",
        "DisableAnywhereCacheRequest",
        "PauseAnywhereCacheRequest",
        "ResumeAnywhereCacheRequest",
        "GetAnywhereCacheRequest",
        "ListAnywhereCachesRequest",
        "ListAnywhereCachesResponse",
        "IntelligenceConfig",
        "UpdateOrganizationIntelligenceConfigRequest",
        "UpdateFolderIntelligenceConfigRequest",
        "UpdateProjectIntelligenceConfigRequest",
        "GetOrganizationIntelligenceConfigRequest",
        "GetFolderIntelligenceConfigRequest",
        "GetProjectIntelligenceConfigRequest",
        "IntelligenceFinding",
        "IntelligenceFindingRevision",
        "GetIntelligenceFindingRequest",
        "ListIntelligenceFindingsRequest",
        "ListIntelligenceFindingsResponse",
        "SummarizeIntelligenceFindingsRequest",
        "SummarizeIntelligenceFindingsResponse",
        "GetIntelligenceFindingRevisionRequest",
        "ListIntelligenceFindingRevisionsRequest",
        "ListIntelligenceFindingRevisionsResponse",
        "FindingSummary",
    },
)


class FindingType(proto.Enum):
    r"""List the finding types.

    Values:
        FINDING_TYPE_UNSPECIFIED (0):
            Finding type is unspecified.
        FINDING_TYPE_COLDLINE_AND_ARCHIVAL_STORAGE_OPERATIONS_SPIKE (1):
            Finding is about a spike in Class A/B
            operations on Coldline or Archive Cloud Storage
            objects.
        FINDING_TYPE_THROTTLED_REQUEST_SPIKE (2):
            Finding is about a spike in throttled
            requests (429 errors) within a project.
        FINDING_TYPE_CROSS_REGION_EGRESS_SPIKE (3):
            Finding is about a spike in cross region
            egress in Cloud Storage.
        FINDING_TYPE_STORAGE_GROWTH_ABOVE_TREND (4):
            Finding is about growth in storage above the
            expected trend.
    """

    FINDING_TYPE_UNSPECIFIED = 0
    FINDING_TYPE_COLDLINE_AND_ARCHIVAL_STORAGE_OPERATIONS_SPIKE = 1
    FINDING_TYPE_THROTTLED_REQUEST_SPIKE = 2
    FINDING_TYPE_CROSS_REGION_EGRESS_SPIKE = 3
    FINDING_TYPE_STORAGE_GROWTH_ABOVE_TREND = 4


class FindingCategory(proto.Enum):
    r"""List of categories a finding falls under.

    Values:
        FINDING_CATEGORY_UNSPECIFIED (0):
            Category is unspecified.
        FINDING_CATEGORY_DATA_MANAGEMENT (1):
            Category is 'Data Management'.
        FINDING_CATEGORY_PERFORMANCE (2):
            Category is 'Performance'.
    """

    FINDING_CATEGORY_UNSPECIFIED = 0
    FINDING_CATEGORY_DATA_MANAGEMENT = 1
    FINDING_CATEGORY_PERFORMANCE = 2


class FindingSeverity(proto.Enum):
    r"""Severity of the ``IntelligenceFinding`` resource.

    Values:
        FINDING_SEVERITY_UNSPECIFIED (0):
            Severity is unspecified.
        FINDING_SEVERITY_CRITICAL (1):
            Severity is critical.
    """

    FINDING_SEVERITY_UNSPECIFIED = 0
    FINDING_SEVERITY_CRITICAL = 1


class PendingRenameInfo(proto.Message):
    r"""Contains information about a pending rename operation.

    Attributes:
        operation (str):
            Output only. The name of the rename
            operation.
    """

    operation: str = proto.Field(
        proto.STRING,
        number=1,
    )


class Folder(proto.Message):
    r"""A folder resource. This resource can only exist in a
    hierarchical namespace enabled bucket.

    Attributes:
        name (str):
            Identifier. The name of this folder. Format:
            ``projects/{project}/buckets/{bucket}/folders/{folder}``
        metageneration (int):
            Output only. The version of the metadata for
            this folder. Used for preconditions and for
            detecting changes in metadata.
        create_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The creation time of the folder.
        update_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The modification time of the
            folder.
        pending_rename_info (google.cloud.storage_control_v2.types.PendingRenameInfo):
            Output only. Only present if the folder is
            part of an ongoing RenameFolder operation.
            Contains information which can be used to query
            the operation status. The presence of this field
            also indicates all write operations are blocked
            for this folder, including folder, managed
            folder, and object operations.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    metageneration: int = proto.Field(
        proto.INT64,
        number=3,
    )
    create_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=4,
        message=timestamp_pb2.Timestamp,
    )
    update_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=5,
        message=timestamp_pb2.Timestamp,
    )
    pending_rename_info: "PendingRenameInfo" = proto.Field(
        proto.MESSAGE,
        number=7,
        message="PendingRenameInfo",
    )


class GetFolderRequest(proto.Message):
    r"""Request message for GetFolder. This operation is only
    applicable to a hierarchical namespace enabled bucket.


    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        name (str):
            Required. Name of the folder. Format:
            ``projects/{project}/buckets/{bucket}/folders/{folder}``
        if_metageneration_match (int):
            Makes the operation only succeed conditional
            on whether the folder's current metageneration
            matches the given value.

            This field is a member of `oneof`_ ``_if_metageneration_match``.
        if_metageneration_not_match (int):
            Makes the operation only succeed conditional
            on whether the folder's current metageneration
            does not match the given value.

            This field is a member of `oneof`_ ``_if_metageneration_not_match``.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    name: str = proto.Field(
        proto.STRING,
        number=6,
    )
    if_metageneration_match: int = proto.Field(
        proto.INT64,
        number=3,
        optional=True,
    )
    if_metageneration_not_match: int = proto.Field(
        proto.INT64,
        number=4,
        optional=True,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=5,
    )


class CreateFolderRequest(proto.Message):
    r"""Request message for CreateFolder. This operation is only
    applicable to a hierarchical namespace enabled bucket.

    Attributes:
        parent (str):
            Required. Name of the bucket in which the
            folder will reside. The bucket must be a
            hierarchical namespace enabled bucket.
        folder (google.cloud.storage_control_v2.types.Folder):
            Required. Properties of the new folder being created. The
            bucket and name of the folder are specified in the parent
            and folder_id fields, respectively. Populating those fields
            in ``folder`` will result in an error.
        folder_id (str):
            Required. The full name of a folder, including all its
            parent folders. Folders use single '/' characters as a
            delimiter. The folder_id must end with a slash. For example,
            the folder_id of "books/biographies/" would create a new
            "biographies/" folder under the "books/" folder.
        recursive (bool):
            Optional. If true, parent folder doesn't have
            to be present and all missing ancestor folders
            will be created atomically.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    folder: "Folder" = proto.Field(
        proto.MESSAGE,
        number=2,
        message="Folder",
    )
    folder_id: str = proto.Field(
        proto.STRING,
        number=3,
    )
    recursive: bool = proto.Field(
        proto.BOOL,
        number=4,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=5,
    )


class DeleteFolderRequest(proto.Message):
    r"""Request message for DeleteFolder. This operation is only
    applicable to a hierarchical namespace enabled bucket.


    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        name (str):
            Required. Name of the folder. Format:
            ``projects/{project}/buckets/{bucket}/folders/{folder}``
        if_metageneration_match (int):
            Makes the operation only succeed conditional
            on whether the folder's current metageneration
            matches the given value.

            This field is a member of `oneof`_ ``_if_metageneration_match``.
        if_metageneration_not_match (int):
            Makes the operation only succeed conditional
            on whether the folder's current metageneration
            does not match the given value.

            This field is a member of `oneof`_ ``_if_metageneration_not_match``.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    name: str = proto.Field(
        proto.STRING,
        number=6,
    )
    if_metageneration_match: int = proto.Field(
        proto.INT64,
        number=3,
        optional=True,
    )
    if_metageneration_not_match: int = proto.Field(
        proto.INT64,
        number=4,
        optional=True,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=5,
    )


class ListFoldersRequest(proto.Message):
    r"""Request message for ListFolders. This operation is only
    applicable to a hierarchical namespace enabled bucket.

    Attributes:
        parent (str):
            Required. Name of the bucket in which to look
            for folders. The bucket must be a hierarchical
            namespace enabled bucket.
        page_size (int):
            Optional. Maximum number of folders to return
            in a single response. The service will use this
            parameter or 1,000 items, whichever is smaller.
        page_token (str):
            Optional. A previously-returned page token
            representing part of the larger set of results
            to view.
        prefix (str):
            Optional. Filter results to folders whose
            names begin with this prefix. If set, the value
            must either be an empty string or end with a
            '/'.
        delimiter (str):
            Optional. If set, returns results in a
            directory-like mode. The results will only
            include folders that either exactly match the
            above prefix, or are one level below the prefix.
            The only supported value is '/'.
        lexicographic_start (str):
            Optional. Filter results to folders whose names are
            lexicographically equal to or after lexicographic_start. If
            lexicographic_end is also set, the folders listed have names
            between lexicographic_start (inclusive) and
            lexicographic_end (exclusive).
        lexicographic_end (str):
            Optional. Filter results to folders whose names are
            lexicographically before lexicographic_end. If
            lexicographic_start is also set, the folders listed have
            names between lexicographic_start (inclusive) and
            lexicographic_end (exclusive).
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    page_size: int = proto.Field(
        proto.INT32,
        number=2,
    )
    page_token: str = proto.Field(
        proto.STRING,
        number=3,
    )
    prefix: str = proto.Field(
        proto.STRING,
        number=4,
    )
    delimiter: str = proto.Field(
        proto.STRING,
        number=8,
    )
    lexicographic_start: str = proto.Field(
        proto.STRING,
        number=6,
    )
    lexicographic_end: str = proto.Field(
        proto.STRING,
        number=7,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=9,
    )


class ListFoldersResponse(proto.Message):
    r"""Response message for ListFolders.

    Attributes:
        folders (MutableSequence[google.cloud.storage_control_v2.types.Folder]):
            The list of child folders
        next_page_token (str):
            The continuation token, used to page through
            large result sets. Provide this value in a
            subsequent request to return the next page of
            results.
    """

    @property
    def raw_page(self):
        return self

    folders: MutableSequence["Folder"] = proto.RepeatedField(
        proto.MESSAGE,
        number=1,
        message="Folder",
    )
    next_page_token: str = proto.Field(
        proto.STRING,
        number=2,
    )


class RenameFolderRequest(proto.Message):
    r"""Request message for RenameFolder. This operation is only
    applicable to a hierarchical namespace enabled bucket.


    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        name (str):
            Required. Name of the source folder being renamed. Format:
            ``projects/{project}/buckets/{bucket}/folders/{folder}``
        destination_folder_id (str):
            Required. The destination folder ID, e.g. ``foo/bar/``.
        if_metageneration_match (int):
            Makes the operation only succeed conditional
            on whether the source folder's current
            metageneration matches the given value.

            This field is a member of `oneof`_ ``_if_metageneration_match``.
        if_metageneration_not_match (int):
            Makes the operation only succeed conditional
            on whether the source folder's current
            metageneration does not match the given value.

            This field is a member of `oneof`_ ``_if_metageneration_not_match``.
        request_id (str):
            Optional. A unique identifier for this request. UUID is the
            recommended format, but other formats are still accepted.
            This request is only idempotent if a ``request_id`` is
            provided.
    """

    name: str = proto.Field(
        proto.STRING,
        number=7,
    )
    destination_folder_id: str = proto.Field(
        proto.STRING,
        number=8,
    )
    if_metageneration_match: int = proto.Field(
        proto.INT64,
        number=4,
        optional=True,
    )
    if_metageneration_not_match: int = proto.Field(
        proto.INT64,
        number=5,
        optional=True,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=6,
    )


class DeleteFolderRecursiveRequest(proto.Message):
    r"""Request message for DeleteFolderRecursive.

    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        name (str):
            Required. Name of the folder being deleted, however all of
            its contents will be deleted too. Format:
            ``projects/{project}/buckets/{bucket}/folders/{folder}``
        if_metageneration_match (int):
            Optional. Makes the operation only succeed
            conditional on whether the root folder's current
            metageneration matches the given value.

            This field is a member of `oneof`_ ``_if_metageneration_match``.
        if_metageneration_not_match (int):
            Optional. Makes the operation only succeed
            conditional on whether the root folder's current
            metageneration does not match the given value.

            This field is a member of `oneof`_ ``_if_metageneration_not_match``.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    if_metageneration_match: int = proto.Field(
        proto.INT64,
        number=2,
        optional=True,
    )
    if_metageneration_not_match: int = proto.Field(
        proto.INT64,
        number=3,
        optional=True,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=4,
    )


class CommonLongRunningOperationMetadata(proto.Message):
    r"""The message contains metadata that is common to all Storage Control
    long-running operations, present in its
    ``google.longrunning.Operation`` messages, and accessible via
    ``metadata.common_metadata``.

    Attributes:
        create_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The time the operation was
            created.
        end_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The time the operation finished
            running.
        update_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The time the operation was last
            modified.
        type_ (str):
            Output only. The type of operation invoked.
        requested_cancellation (bool):
            Output only. Identifies whether the user has
            requested cancellation.
        progress_percent (int):
            Output only. The estimated progress of the operation in
            percentage [0, 100]. The value -1 means the progress is
            unknown.
    """

    create_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=1,
        message=timestamp_pb2.Timestamp,
    )
    end_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=2,
        message=timestamp_pb2.Timestamp,
    )
    update_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=3,
        message=timestamp_pb2.Timestamp,
    )
    type_: str = proto.Field(
        proto.STRING,
        number=4,
    )
    requested_cancellation: bool = proto.Field(
        proto.BOOL,
        number=5,
    )
    progress_percent: int = proto.Field(
        proto.INT32,
        number=6,
    )


class RenameFolderMetadata(proto.Message):
    r"""Message returned in the metadata field of the Operation
    resource for RenameFolder operations.

    Attributes:
        common_metadata (google.cloud.storage_control_v2.types.CommonLongRunningOperationMetadata):
            Generic metadata for the long running
            operation.
        source_folder_id (str):
            The path of the source folder.
        destination_folder_id (str):
            The path of the destination folder.
    """

    common_metadata: "CommonLongRunningOperationMetadata" = proto.Field(
        proto.MESSAGE,
        number=1,
        message="CommonLongRunningOperationMetadata",
    )
    source_folder_id: str = proto.Field(
        proto.STRING,
        number=2,
    )
    destination_folder_id: str = proto.Field(
        proto.STRING,
        number=3,
    )


class DeleteFolderRecursiveMetadata(proto.Message):
    r"""Message returned in the metadata field of the Operation
    resource for DeleteFolderRecursive operations.

    Attributes:
        common_metadata (google.cloud.storage_control_v2.types.CommonLongRunningOperationMetadata):
            Generic metadata for the long running
            operation.
        folder_id (str):
            The path of the folder recursively deleted.
    """

    common_metadata: "CommonLongRunningOperationMetadata" = proto.Field(
        proto.MESSAGE,
        number=1,
        message="CommonLongRunningOperationMetadata",
    )
    folder_id: str = proto.Field(
        proto.STRING,
        number=2,
    )


class StorageLayout(proto.Message):
    r"""The storage layout configuration of a bucket.

    Attributes:
        name (str):
            Output only. The name of the StorageLayout resource. Format:
            ``projects/{project}/buckets/{bucket}/storageLayout``
        location (str):
            Output only. The location of the bucket.
        location_type (str):
            Output only. The location type of the bucket
            (region, dual-region, multi-region, etc).
        custom_placement_config (google.cloud.storage_control_v2.types.StorageLayout.CustomPlacementConfig):
            Output only. The data placement configuration
            for custom dual region. If there is no
            configuration, this is not a custom dual region
            bucket.
        hierarchical_namespace (google.cloud.storage_control_v2.types.StorageLayout.HierarchicalNamespace):
            Output only. The bucket's hierarchical
            namespace configuration. If there is no
            configuration, the hierarchical namespace is
            disabled.
    """

    class CustomPlacementConfig(proto.Message):
        r"""Configuration for Custom Dual Regions. It should specify precisely
        two eligible regions within the same Multiregion. More information
        on regions may be found
        `here <https://cloud.google.com/storage/docs/locations>`__.

        Attributes:
            data_locations (MutableSequence[str]):
                List of locations to use for data placement.
        """

        data_locations: MutableSequence[str] = proto.RepeatedField(
            proto.STRING,
            number=1,
        )

    class HierarchicalNamespace(proto.Message):
        r"""Configuration for a bucket's hierarchical namespace feature.

        Attributes:
            enabled (bool):
                Enables the hierarchical namespace feature.
        """

        enabled: bool = proto.Field(
            proto.BOOL,
            number=1,
        )

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    location: str = proto.Field(
        proto.STRING,
        number=2,
    )
    location_type: str = proto.Field(
        proto.STRING,
        number=3,
    )
    custom_placement_config: CustomPlacementConfig = proto.Field(
        proto.MESSAGE,
        number=4,
        message=CustomPlacementConfig,
    )
    hierarchical_namespace: HierarchicalNamespace = proto.Field(
        proto.MESSAGE,
        number=5,
        message=HierarchicalNamespace,
    )


class GetStorageLayoutRequest(proto.Message):
    r"""Request message for GetStorageLayout.

    Attributes:
        name (str):
            Required. The name of the StorageLayout resource. Format:
            ``projects/{project}/buckets/{bucket}/storageLayout``
        prefix (str):
            An optional prefix used for permission check.
            It is useful when the caller only has limited
            permissions under a specific prefix.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    prefix: str = proto.Field(
        proto.STRING,
        number=2,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=3,
    )


class ManagedFolder(proto.Message):
    r"""A managed folder.

    Attributes:
        name (str):
            Identifier. The name of this managed folder. Format:
            ``projects/{project}/buckets/{bucket}/managedFolders/{managedFolder}``
        metageneration (int):
            Output only. The metadata version of this
            managed folder. It increases whenever the
            metadata is updated. Used for preconditions and
            for detecting changes in metadata. Managed
            folders don't have a generation number.
        create_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The creation time of the managed
            folder.
        update_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The modification time of the
            managed folder.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    metageneration: int = proto.Field(
        proto.INT64,
        number=3,
    )
    create_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=4,
        message=timestamp_pb2.Timestamp,
    )
    update_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=5,
        message=timestamp_pb2.Timestamp,
    )


class GetManagedFolderRequest(proto.Message):
    r"""Request message for GetManagedFolder.

    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        name (str):
            Required. Name of the managed folder. Format:
            ``projects/{project}/buckets/{bucket}/managedFolders/{managedFolder}``
        if_metageneration_match (int):
            The operation succeeds conditional on the
            managed folder's current metageneration matching
            the value here specified.

            This field is a member of `oneof`_ ``_if_metageneration_match``.
        if_metageneration_not_match (int):
            The operation succeeds conditional on the
            managed folder's current metageneration NOT
            matching the value here specified.

            This field is a member of `oneof`_ ``_if_metageneration_not_match``.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    name: str = proto.Field(
        proto.STRING,
        number=6,
    )
    if_metageneration_match: int = proto.Field(
        proto.INT64,
        number=3,
        optional=True,
    )
    if_metageneration_not_match: int = proto.Field(
        proto.INT64,
        number=4,
        optional=True,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=5,
    )


class CreateManagedFolderRequest(proto.Message):
    r"""Request message for CreateManagedFolder.

    Attributes:
        parent (str):
            Required. Name of the bucket this managed
            folder belongs to.
        managed_folder (google.cloud.storage_control_v2.types.ManagedFolder):
            Required. Properties of the managed folder being created.
            The bucket and managed folder names are specified in the
            ``parent`` and ``managed_folder_id`` fields. Populating
            these fields in ``managed_folder`` will result in an error.
        managed_folder_id (str):
            Required. The name of the managed folder. It uses a single
            ``/`` as delimiter and leading and trailing ``/`` are
            allowed.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    managed_folder: "ManagedFolder" = proto.Field(
        proto.MESSAGE,
        number=2,
        message="ManagedFolder",
    )
    managed_folder_id: str = proto.Field(
        proto.STRING,
        number=3,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=4,
    )


class DeleteManagedFolderRequest(proto.Message):
    r"""DeleteManagedFolder RPC request message.

    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        name (str):
            Required. Name of the managed folder. Format:
            ``projects/{project}/buckets/{bucket}/managedFolders/{managedFolder}``
        if_metageneration_match (int):
            The operation succeeds conditional on the
            managed folder's current metageneration matching
            the value here specified.

            This field is a member of `oneof`_ ``_if_metageneration_match``.
        if_metageneration_not_match (int):
            The operation succeeds conditional on the
            managed folder's current metageneration NOT
            matching the value here specified.

            This field is a member of `oneof`_ ``_if_metageneration_not_match``.
        allow_non_empty (bool):
            Allows deletion of a managed folder even if
            it is not empty. A managed folder is empty if it
            manages no child managed folders or objects.
            Caller must have permission for
            storage.managedFolders.setIamPolicy.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    name: str = proto.Field(
        proto.STRING,
        number=7,
    )
    if_metageneration_match: int = proto.Field(
        proto.INT64,
        number=3,
        optional=True,
    )
    if_metageneration_not_match: int = proto.Field(
        proto.INT64,
        number=4,
        optional=True,
    )
    allow_non_empty: bool = proto.Field(
        proto.BOOL,
        number=5,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=6,
    )


class ListManagedFoldersRequest(proto.Message):
    r"""Request message for ListManagedFolders.

    Attributes:
        parent (str):
            Required. Name of the bucket this managed
            folder belongs to.
        page_size (int):
            Optional. Maximum number of managed folders
            to return in a single response. The service will
            use this parameter or 1,000 items, whichever is
            smaller.
        page_token (str):
            Optional. A previously-returned page token
            representing part of the larger set of results
            to view.
        prefix (str):
            Optional. Filter results to match managed
            folders with name starting with this prefix.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    page_size: int = proto.Field(
        proto.INT32,
        number=2,
    )
    page_token: str = proto.Field(
        proto.STRING,
        number=3,
    )
    prefix: str = proto.Field(
        proto.STRING,
        number=4,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=5,
    )


class ListManagedFoldersResponse(proto.Message):
    r"""Response message for ListManagedFolders.

    Attributes:
        managed_folders (MutableSequence[google.cloud.storage_control_v2.types.ManagedFolder]):
            The list of matching managed folders
        next_page_token (str):
            The continuation token, used to page through
            large result sets. Provide this value in a
            subsequent request to return the next page of
            results.
    """

    @property
    def raw_page(self):
        return self

    managed_folders: MutableSequence["ManagedFolder"] = proto.RepeatedField(
        proto.MESSAGE,
        number=1,
        message="ManagedFolder",
    )
    next_page_token: str = proto.Field(
        proto.STRING,
        number=2,
    )


class CreateAnywhereCacheMetadata(proto.Message):
    r"""Message returned in the metadata field of the Operation
    resource for CreateAnywhereCache operations.


    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        common_metadata (google.cloud.storage_control_v2.types.CommonLongRunningOperationMetadata):
            Generic metadata for the long running
            operation.
        anywhere_cache_id (str):
            Anywhere Cache ID.

            This field is a member of `oneof`_ ``_anywhere_cache_id``.
        zone (str):
            The zone in which the cache instance is
            running. For example, us-central1-a.

            This field is a member of `oneof`_ ``_zone``.
        ttl (google.protobuf.duration_pb2.Duration):
            Anywhere Cache entry's TTL. A cache-level
            config that is applied to all new cache entries
            on admission. Default ttl value (24hrs) is
            applied if not specified in the create request.

            This field is a member of `oneof`_ ``_ttl``.
        admission_policy (str):
            Anywhere Cache entry Admission Policy in
            kebab-case (e.g., "admit-on-first-miss").
            Default admission policy (admit-on-first-miss)
            is applied if not specified in the create
            request.

            This field is a member of `oneof`_ ``_admission_policy``.
    """

    common_metadata: "CommonLongRunningOperationMetadata" = proto.Field(
        proto.MESSAGE,
        number=1,
        message="CommonLongRunningOperationMetadata",
    )
    anywhere_cache_id: str = proto.Field(
        proto.STRING,
        number=2,
        optional=True,
    )
    zone: str = proto.Field(
        proto.STRING,
        number=6,
        optional=True,
    )
    ttl: duration_pb2.Duration = proto.Field(
        proto.MESSAGE,
        number=3,
        optional=True,
        message=duration_pb2.Duration,
    )
    admission_policy: str = proto.Field(
        proto.STRING,
        number=5,
        optional=True,
    )


class UpdateAnywhereCacheMetadata(proto.Message):
    r"""Message returned in the metadata field of the Operation
    resource for UpdateAnywhereCache operation.


    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        common_metadata (google.cloud.storage_control_v2.types.CommonLongRunningOperationMetadata):
            Generic metadata for the long running
            operation.
        anywhere_cache_id (str):
            Anywhere Cache ID.

            This field is a member of `oneof`_ ``_anywhere_cache_id``.
        zone (str):
            The zone in which the cache instance is
            running. For example, us-central1-a.

            This field is a member of `oneof`_ ``_zone``.
        ttl (google.protobuf.duration_pb2.Duration):
            Anywhere Cache entry's TTL between 1h and 7days. A
            cache-level config that is applied to all new cache entries
            on admission. If ``ttl`` is pending update, this field
            equals to the new value specified in the Update request.

            This field is a member of `oneof`_ ``_ttl``.
        admission_policy (str):
            L4 Cache entry Admission Policy in kebab-case (e.g.,
            "admit-on-first-miss"). If ``admission_policy`` is pending
            update, this field equals to the new value specified in the
            Update request.

            This field is a member of `oneof`_ ``_admission_policy``.
    """

    common_metadata: "CommonLongRunningOperationMetadata" = proto.Field(
        proto.MESSAGE,
        number=1,
        message="CommonLongRunningOperationMetadata",
    )
    anywhere_cache_id: str = proto.Field(
        proto.STRING,
        number=2,
        optional=True,
    )
    zone: str = proto.Field(
        proto.STRING,
        number=5,
        optional=True,
    )
    ttl: duration_pb2.Duration = proto.Field(
        proto.MESSAGE,
        number=3,
        optional=True,
        message=duration_pb2.Duration,
    )
    admission_policy: str = proto.Field(
        proto.STRING,
        number=4,
        optional=True,
    )


class AnywhereCache(proto.Message):
    r"""An Anywhere Cache Instance.

    Attributes:
        name (str):
            Immutable. The resource name of this AnywhereCache. Format:
            ``projects/{project}/buckets/{bucket}/anywhereCaches/{anywhere_cache}``
        zone (str):
            Immutable. The zone in which the cache
            instance is running. For example, us-central1-a.
        ttl (google.protobuf.duration_pb2.Duration):
            Cache entry TTL (ranges between 1h to 7d).
            This is a cache-level config that defines how
            long a cache entry can live. Default ttl value
            (24hrs) is applied if not specified in the
            create request. TTL must be in whole seconds.
        admission_policy (str):
            Cache admission policy. Valid policies includes:
            ``admit-on-first-miss`` and ``admit-on-second-miss``.
            Defaults to ``admit-on-first-miss``. Default value is
            applied if not specified in the create request.
        state (str):
            Output only. Cache state including RUNNING,
            CREATING, DISABLED and PAUSED.
        create_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. Time when Anywhere cache
            instance is allocated.
        update_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. Time when Anywhere cache
            instance is last updated, including creation.
        pending_update (bool):
            Output only. True if there is an active
            update operation against this cache instance.
            Subsequential update requests will be rejected
            if this field is true. Output only.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    zone: str = proto.Field(
        proto.STRING,
        number=10,
    )
    ttl: duration_pb2.Duration = proto.Field(
        proto.MESSAGE,
        number=3,
        message=duration_pb2.Duration,
    )
    admission_policy: str = proto.Field(
        proto.STRING,
        number=9,
    )
    state: str = proto.Field(
        proto.STRING,
        number=5,
    )
    create_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=6,
        message=timestamp_pb2.Timestamp,
    )
    update_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=7,
        message=timestamp_pb2.Timestamp,
    )
    pending_update: bool = proto.Field(
        proto.BOOL,
        number=8,
    )


class CreateAnywhereCacheRequest(proto.Message):
    r"""Request message for CreateAnywhereCache.

    Attributes:
        parent (str):
            Required. The bucket to which this cache belongs. Format:
            ``projects/{project}/buckets/{bucket}``
        anywhere_cache (google.cloud.storage_control_v2.types.AnywhereCache):
            Required. Properties of the Anywhere Cache instance being
            created. The parent bucket name is specified in the
            ``parent`` field. Server uses the default value of ``ttl``
            or ``admission_policy`` if not specified in request.
        request_id (str):
            Optional. A unique identifier for this request. UUID is the
            recommended format, but other formats are still accepted.
            This request is only idempotent if a ``request_id`` is
            provided.
    """

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    anywhere_cache: "AnywhereCache" = proto.Field(
        proto.MESSAGE,
        number=3,
        message="AnywhereCache",
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=4,
    )


class UpdateAnywhereCacheRequest(proto.Message):
    r"""Request message for UpdateAnywhereCache.

    Attributes:
        anywhere_cache (google.cloud.storage_control_v2.types.AnywhereCache):
            Required. The Anywhere Cache instance to be
            updated.
        update_mask (google.protobuf.field_mask_pb2.FieldMask):
            Required. List of fields to be updated. Mutable fields of
            AnywhereCache include ``ttl`` and ``admission_policy``.

            To specify ALL fields, specify a single field with the value
            ``*``. Note: We recommend against doing this. If a new field
            is introduced at a later time, an older client updating with
            the ``*`` may accidentally reset the new field's value.

            Not specifying any fields is an error.
        request_id (str):
            Optional. A unique identifier for this request. UUID is the
            recommended format, but other formats are still accepted.
            This request is only idempotent if a ``request_id`` is
            provided.
    """

    anywhere_cache: "AnywhereCache" = proto.Field(
        proto.MESSAGE,
        number=1,
        message="AnywhereCache",
    )
    update_mask: field_mask_pb2.FieldMask = proto.Field(
        proto.MESSAGE,
        number=2,
        message=field_mask_pb2.FieldMask,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=3,
    )


class DisableAnywhereCacheRequest(proto.Message):
    r"""Request message for DisableAnywhereCache.

    Attributes:
        name (str):
            Required. The name field in the request should be:
            ``projects/{project}/buckets/{bucket}/anywhereCaches/{anywhere_cache}``
        request_id (str):
            Optional. A unique identifier for this request. UUID is the
            recommended format, but other formats are still accepted.
            This request is only idempotent if a ``request_id`` is
            provided.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=2,
    )


class PauseAnywhereCacheRequest(proto.Message):
    r"""Request message for PauseAnywhereCache.

    Attributes:
        name (str):
            Required. The name field in the request should be:
            ``projects/{project}/buckets/{bucket}/anywhereCaches/{anywhere_cache}``
        request_id (str):
            Optional. A unique identifier for this request. UUID is the
            recommended format, but other formats are still accepted.
            This request is only idempotent if a ``request_id`` is
            provided.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=2,
    )


class ResumeAnywhereCacheRequest(proto.Message):
    r"""Request message for ResumeAnywhereCache.

    Attributes:
        name (str):
            Required. The name field in the request should be:
            ``projects/{project}/buckets/{bucket}/anywhereCaches/{anywhere_cache}``
        request_id (str):
            Optional. A unique identifier for this request. UUID is the
            recommended format, but other formats are still accepted.
            This request is only idempotent if a ``request_id`` is
            provided.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=2,
    )


class GetAnywhereCacheRequest(proto.Message):
    r"""Request message for GetAnywhereCache.

    Attributes:
        name (str):
            Required. The name field in the request should be:
            ``projects/{project}/buckets/{bucket}/anywhereCaches/{anywhere_cache}``
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=2,
    )


class ListAnywhereCachesRequest(proto.Message):
    r"""Request message for ListAnywhereCaches.

    Attributes:
        parent (str):
            Required. The bucket to which this cache
            belongs.
        page_size (int):
            Maximum number of caches to return in a
            single response. The service will use this
            parameter or 1,000 items, whichever is smaller.
        page_token (str):
            A previously-returned page token representing
            part of the larger set of results to view.
        request_id (str):
            Optional. A unique identifier for this
            request. UUID is the recommended format, but
            other formats are still accepted.
    """

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    page_size: int = proto.Field(
        proto.INT32,
        number=2,
    )
    page_token: str = proto.Field(
        proto.STRING,
        number=3,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=4,
    )


class ListAnywhereCachesResponse(proto.Message):
    r"""Response message for ListAnywhereCaches.

    Attributes:
        anywhere_caches (MutableSequence[google.cloud.storage_control_v2.types.AnywhereCache]):
            The list of items.
        next_page_token (str):
            A token, which can be sent as ``page_token`` to retrieve the
            next page. If this field is omitted, there are no subsequent
            pages.
    """

    @property
    def raw_page(self):
        return self

    anywhere_caches: MutableSequence["AnywhereCache"] = proto.RepeatedField(
        proto.MESSAGE,
        number=1,
        message="AnywhereCache",
    )
    next_page_token: str = proto.Field(
        proto.STRING,
        number=2,
    )


class IntelligenceConfig(proto.Message):
    r"""The ``IntelligenceConfig`` resource associated with your
    organization, folder, or project.

    Attributes:
        name (str):
            Identifier. The name of the ``IntelligenceConfig`` resource
            associated with your organization, folder, or project.

            The name format varies based on the GCP resource hierarchy
            as follows:

            - For project:
              ``projects/{project_number}/locations/global/intelligenceConfig``
            - For organization:
              ``organizations/{org_id}/locations/global/intelligenceConfig``
            - For folder:
              ``folders/{folder_id}/locations/global/intelligenceConfig``
        edition_config (google.cloud.storage_control_v2.types.IntelligenceConfig.EditionConfig):
            Optional. The edition configuration of the
            ``IntelligenceConfig`` resource.
        update_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The time at which the ``IntelligenceConfig``
            resource is last updated.
        filter (google.cloud.storage_control_v2.types.IntelligenceConfig.Filter):
            Optional. Filter over location and bucket.
        effective_intelligence_config (google.cloud.storage_control_v2.types.IntelligenceConfig.EffectiveIntelligenceConfig):
            Output only. The ``IntelligenceConfig`` resource that is
            applicable for the resource.
        trial_config (google.cloud.storage_control_v2.types.IntelligenceConfig.TrialConfig):
            The trial configuration of the ``IntelligenceConfig``
            resource.
    """

    class EditionConfig(proto.Enum):
        r"""The edition configuration of the ``IntelligenceConfig`` resource.
        This signifies the edition used for configuring the
        ``IntelligenceConfig`` resource and can only take the following
        values: ``EDITION_CONFIG_UNSPECIFIED``, ``INHERIT``, ``DISABLED``,
        ``STANDARD`` and ``TRIAL``.

        Values:
            EDITION_CONFIG_UNSPECIFIED (0):
                This is an unknown edition of the resource.
            INHERIT (1):
                The inherited edition from the parent and filters. This is
                the default edition when there is no ``IntelligenceConfig``
                setup for a GCP resource.
            DISABLED (2):
                The edition configuration is disabled for the
                ``IntelligenceConfig`` resource and its children. Filters
                are not applicable.
            STANDARD (3):
                The ``IntelligenceConfig`` resource is of STANDARD edition.
            TRIAL (5):
                The ``IntelligenceConfig`` resource is available in
                ``TRIAL`` edition. During the trial period, Cloud Storage
                does not charge for Storage Intelligence usage. You can
                specify the buckets to include in the trial period by using
                filters. At the end of the trial period, the
                ``IntelligenceConfig`` resource is upgraded to ``STANDARD``
                edition.
        """

        EDITION_CONFIG_UNSPECIFIED = 0
        INHERIT = 1
        DISABLED = 2
        STANDARD = 3
        TRIAL = 5

    class Filter(proto.Message):
        r"""Filter over location and bucket using include or exclude
        semantics. Resources that match the include or exclude filter
        are exclusively included or excluded from the Storage
        Intelligence plan.

        This message has `oneof`_ fields (mutually exclusive fields).
        For each oneof, at most one member field can be set at the same time.
        Setting any member of the oneof automatically clears all other
        members.

        .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

        Attributes:
            included_cloud_storage_locations (google.cloud.storage_control_v2.types.IntelligenceConfig.Filter.CloudStorageLocations):
                Bucket locations to include.

                This field is a member of `oneof`_ ``cloud_storage_locations``.
            excluded_cloud_storage_locations (google.cloud.storage_control_v2.types.IntelligenceConfig.Filter.CloudStorageLocations):
                Bucket locations to exclude.

                This field is a member of `oneof`_ ``cloud_storage_locations``.
            included_cloud_storage_buckets (google.cloud.storage_control_v2.types.IntelligenceConfig.Filter.CloudStorageBuckets):
                Buckets to include.

                This field is a member of `oneof`_ ``cloud_storage_buckets``.
            excluded_cloud_storage_buckets (google.cloud.storage_control_v2.types.IntelligenceConfig.Filter.CloudStorageBuckets):
                Buckets to exclude.

                This field is a member of `oneof`_ ``cloud_storage_buckets``.
        """

        class CloudStorageLocations(proto.Message):
            r"""Collection of bucket locations.

            Attributes:
                locations (MutableSequence[str]):
                    Optional. Bucket locations. Location can be any of the Cloud
                    Storage regions specified in lower case format. For example,
                    ``us-east1``, ``us-west1``.
            """

            locations: MutableSequence[str] = proto.RepeatedField(
                proto.STRING,
                number=1,
            )

        class CloudStorageBuckets(proto.Message):
            r"""Collection of buckets.

            Attributes:
                bucket_id_regexes (MutableSequence[str]):
                    Optional. A regex pattern for matching bucket names. Regex
                    should follow the syntax specified in
                    `google/re2 <https://github.com/google/re2>`__. For example,
                    ``^sample_.*`` matches all buckets of the form
                    ``gs://sample_bucket-1``, ``gs://sample_bucket-2``,
                    ``gs://sample_bucket-n`` but not
                    ``gs://test_sample_bucket``. If you want to match a single
                    bucket, say ``gs://sample_bucket``, use ``sample_bucket``.
            """

            bucket_id_regexes: MutableSequence[str] = proto.RepeatedField(
                proto.STRING,
                number=1,
            )

        included_cloud_storage_locations: "IntelligenceConfig.Filter.CloudStorageLocations" = proto.Field(
            proto.MESSAGE,
            number=1,
            oneof="cloud_storage_locations",
            message="IntelligenceConfig.Filter.CloudStorageLocations",
        )
        excluded_cloud_storage_locations: "IntelligenceConfig.Filter.CloudStorageLocations" = proto.Field(
            proto.MESSAGE,
            number=2,
            oneof="cloud_storage_locations",
            message="IntelligenceConfig.Filter.CloudStorageLocations",
        )
        included_cloud_storage_buckets: "IntelligenceConfig.Filter.CloudStorageBuckets" = proto.Field(
            proto.MESSAGE,
            number=3,
            oneof="cloud_storage_buckets",
            message="IntelligenceConfig.Filter.CloudStorageBuckets",
        )
        excluded_cloud_storage_buckets: "IntelligenceConfig.Filter.CloudStorageBuckets" = proto.Field(
            proto.MESSAGE,
            number=4,
            oneof="cloud_storage_buckets",
            message="IntelligenceConfig.Filter.CloudStorageBuckets",
        )

    class EffectiveIntelligenceConfig(proto.Message):
        r"""The effective ``IntelligenceConfig`` for the resource.

        Attributes:
            effective_edition (google.cloud.storage_control_v2.types.IntelligenceConfig.EffectiveIntelligenceConfig.EffectiveEdition):
                Output only. The ``IntelligenceConfig`` edition that is
                applicable for the resource.
            intelligence_config (str):
                Output only. The ``IntelligenceConfig`` resource that is
                applied for the target resource. Format:
                ``{organizations|folders|projects}/{id}/locations/{location}/intelligenceConfig``
        """

        class EffectiveEdition(proto.Enum):
            r"""The effective edition of the ``IntelligenceConfig`` resource.

            Values:
                EFFECTIVE_EDITION_UNSPECIFIED (0):
                    This is an unknown edition of the resource.
                NONE (1):
                    No edition.
                STANDARD (2):
                    The ``IntelligenceConfig`` resource is of STANDARD edition.
            """

            EFFECTIVE_EDITION_UNSPECIFIED = 0
            NONE = 1
            STANDARD = 2

        effective_edition: "IntelligenceConfig.EffectiveIntelligenceConfig.EffectiveEdition" = proto.Field(
            proto.ENUM,
            number=1,
            enum="IntelligenceConfig.EffectiveIntelligenceConfig.EffectiveEdition",
        )
        intelligence_config: str = proto.Field(
            proto.STRING,
            number=2,
        )

    class TrialConfig(proto.Message):
        r"""The trial configuration of the ``IntelligenceConfig`` resource.

        Attributes:
            expire_time (google.protobuf.timestamp_pb2.Timestamp):
                Output only. The time at which the trial
                expires.
        """

        expire_time: timestamp_pb2.Timestamp = proto.Field(
            proto.MESSAGE,
            number=3,
            message=timestamp_pb2.Timestamp,
        )

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    edition_config: EditionConfig = proto.Field(
        proto.ENUM,
        number=2,
        enum=EditionConfig,
    )
    update_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=3,
        message=timestamp_pb2.Timestamp,
    )
    filter: Filter = proto.Field(
        proto.MESSAGE,
        number=4,
        message=Filter,
    )
    effective_intelligence_config: EffectiveIntelligenceConfig = proto.Field(
        proto.MESSAGE,
        number=5,
        message=EffectiveIntelligenceConfig,
    )
    trial_config: TrialConfig = proto.Field(
        proto.MESSAGE,
        number=7,
        message=TrialConfig,
    )


class UpdateOrganizationIntelligenceConfigRequest(proto.Message):
    r"""Request message to update the ``IntelligenceConfig`` resource
    associated with your organization.

    **IAM Permissions**:

    Requires ``storage.intelligenceConfigs.update``
    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
    permission on the organization.

    Attributes:
        intelligence_config (google.cloud.storage_control_v2.types.IntelligenceConfig):
            Required. The ``IntelligenceConfig`` resource to be updated.
        update_mask (google.protobuf.field_mask_pb2.FieldMask):
            Required. The ``update_mask`` that specifies the fields
            within the ``IntelligenceConfig`` resource that should be
            modified by this update. Only the listed fields are updated.
        request_id (str):
            Optional. The ID that uniquely identifies the
            request, preventing duplicate processing.
    """

    intelligence_config: "IntelligenceConfig" = proto.Field(
        proto.MESSAGE,
        number=1,
        message="IntelligenceConfig",
    )
    update_mask: field_mask_pb2.FieldMask = proto.Field(
        proto.MESSAGE,
        number=2,
        message=field_mask_pb2.FieldMask,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=3,
    )


class UpdateFolderIntelligenceConfigRequest(proto.Message):
    r"""Request message to update the ``IntelligenceConfig`` resource
    associated with your folder.

    **IAM Permissions**:

    Requires ``storage.intelligenceConfigs.update``
    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
    permission on the folder.

    Attributes:
        intelligence_config (google.cloud.storage_control_v2.types.IntelligenceConfig):
            Required. The ``IntelligenceConfig`` resource to be updated.
        update_mask (google.protobuf.field_mask_pb2.FieldMask):
            Required. The ``update_mask`` that specifies the fields
            within the ``IntelligenceConfig`` resource that should be
            modified by this update. Only the listed fields are updated.
        request_id (str):
            Optional. The ID that uniquely identifies the
            request, preventing duplicate processing.
    """

    intelligence_config: "IntelligenceConfig" = proto.Field(
        proto.MESSAGE,
        number=1,
        message="IntelligenceConfig",
    )
    update_mask: field_mask_pb2.FieldMask = proto.Field(
        proto.MESSAGE,
        number=2,
        message=field_mask_pb2.FieldMask,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=3,
    )


class UpdateProjectIntelligenceConfigRequest(proto.Message):
    r"""Request message to update the ``IntelligenceConfig`` resource
    associated with your project.

    **IAM Permissions**:

    Requires ``storage.intelligenceConfigs.update``
    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
    permission on the folder.

    Attributes:
        intelligence_config (google.cloud.storage_control_v2.types.IntelligenceConfig):
            Required. The ``IntelligenceConfig`` resource to be updated.
        update_mask (google.protobuf.field_mask_pb2.FieldMask):
            Required. The ``update_mask`` that specifies the fields
            within the ``IntelligenceConfig`` resource that should be
            modified by this update. Only the listed fields are updated.
        request_id (str):
            Optional. The ID that uniquely identifies the
            request, preventing duplicate processing.
    """

    intelligence_config: "IntelligenceConfig" = proto.Field(
        proto.MESSAGE,
        number=1,
        message="IntelligenceConfig",
    )
    update_mask: field_mask_pb2.FieldMask = proto.Field(
        proto.MESSAGE,
        number=2,
        message=field_mask_pb2.FieldMask,
    )
    request_id: str = proto.Field(
        proto.STRING,
        number=3,
    )


class GetOrganizationIntelligenceConfigRequest(proto.Message):
    r"""Request message to get the ``IntelligenceConfig`` resource
    associated with your organization.

    **IAM Permissions**

    Requires ``storage.intelligenceConfigs.get``
    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
    permission on the organization.

    Attributes:
        name (str):
            Required. The name of the ``IntelligenceConfig`` resource
            associated with your organization.

            Format:
            ``organizations/{org_id}/locations/global/intelligenceConfig``
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )


class GetFolderIntelligenceConfigRequest(proto.Message):
    r"""Request message to get the ``IntelligenceConfig`` resource
    associated with your folder.

    **IAM Permissions**

    Requires ``storage.intelligenceConfigs.get``
    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
    permission on the folder.

    Attributes:
        name (str):
            Required. The name of the ``IntelligenceConfig`` resource
            associated with your folder.

            Format: ``folders/{id}/locations/global/intelligenceConfig``
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )


class GetProjectIntelligenceConfigRequest(proto.Message):
    r"""Request message to get the ``IntelligenceConfig`` resource
    associated with your project.

    **IAM Permissions**:

    Requires ``storage.intelligenceConfigs.get``
    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
    permission on the project.

    Attributes:
        name (str):
            Required. The name of the ``IntelligenceConfig`` resource
            associated with your project.

            Format:
            ``projects/{id}/locations/global/intelligenceConfig``
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )


class IntelligenceFinding(proto.Message):
    r"""The ``IntelligenceFinding`` resource that represents a security,
    performance, or cost-related finding about a project or bucket.

    This message has `oneof`_ fields (mutually exclusive fields).
    For each oneof, at most one member field can be set at the same time.
    Setting any member of the oneof automatically clears all other
    members.

    .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

    Attributes:
        name (str):
            Identifier. The resource name of ``IntelligenceFinding``.
            Format:
            ``projects/{project}/locations/{location}/intelligenceFindings/{intelligence_finding}``
        description (str):
            Output only. A short description about the
            finding.
        type_ (google.cloud.storage_control_v2.types.FindingType):
            Output only. Type of this finding.
        category (google.cloud.storage_control_v2.types.FindingCategory):
            Output only. Category of this finding.
        severity (google.cloud.storage_control_v2.types.FindingSeverity):
            Output only. Severity of the finding.
        create_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The time at which the finding
            was created.
        update_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The time at which the finding
            was last updated.
        target_resource (str):
            Output only. The fully qualified resource name of the
            resource that this ``IntelligenceFinding`` applies to. eg:

            - ``storage.googleapis.com/projects/_/buckets/b1``
            - ``cloudresourecemanager.googleapis.com/projects/p1``
        associated_resources (MutableSequence[str]):
            Output only. Contains GCP resource names that are relevant
            to this ``IntelligenceFinding``. The ``target_resource`` is
            also added as part of ``associated_resources``. eg:

            - ``storage.googleapis.com/projects/_/buckets/b1``
            - ``cloudresourecemanager.googleapis.com/projects/p1``
        observation_period (google.type.interval_pb2.Interval):
            Output only. The time interval during which the underlying
            data was used to generate this ``IntelligenceFinding``.
        coldline_and_archival_storage_operations_spike (google.cloud.storage_control_v2.types.IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike):
            Output only. ``IntelligenceFinding`` about a spike in Class
            A/B operations on Coldline or Archive Cloud Storage objects.

            This field is a member of `oneof`_ ``intelligence_finding_details``.
        throttled_requests_spike (google.cloud.storage_control_v2.types.IntelligenceFinding.ThrottledRequestSpike):
            Output only. ``IntelligenceFinding`` about a spike in
            throttled requests (429 errors) within a project.

            This field is a member of `oneof`_ ``intelligence_finding_details``.
        cross_region_egress_spike (google.cloud.storage_control_v2.types.IntelligenceFinding.CrossRegionEgressSpike):
            Output only. ``IntelligenceFinding`` about a spike in
            cross-region egress.

            This field is a member of `oneof`_ ``intelligence_finding_details``.
        storage_growth_above_trend (google.cloud.storage_control_v2.types.IntelligenceFinding.StorageGrowthAboveTrend):
            Output only. ``IntelligenceFinding`` about growth in storage
            above the expected trend.

            This field is a member of `oneof`_ ``intelligence_finding_details``.
    """

    class ColdlineAndArchivalStorageOperationsSpike(proto.Message):
        r"""Represents a finding about a spike in Class A/B operations on
        Coldline or Archive Cloud Storage objects. This corresponds to the
        ``COLD_AND_ARCHIVAL_STORAGE_OPERATIONS_SPIKE`` finding type.

        Attributes:
            percentage_increase (float):
                Output only. The percentage increase in
                operations across the project.
            total_operations_count (int):
                Output only. The total count of operations
                across the project.
            top_buckets (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution]):
                Output only. A list of the top buckets
                driving the increase in operations.
        """

        class BucketContribution(proto.Message):
            r"""Represents the operation spike details for a bucket.

            This message has `oneof`_ fields (mutually exclusive fields).
            For each oneof, at most one member field can be set at the same time.
            Setting any member of the oneof automatically clears all other
            members.

            .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

            Attributes:
                bucket (str):
                    Output only. The name of the bucket.
                percentage_increase (float):
                    Output only. The percentage increase in
                    operations for the bucket.
                total_operations_count (int):
                    Output only. The total count of operations
                    for the bucket.
                contribution (google.cloud.storage_control_v2.types.IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution.Contribution):
                    Output only. The details about the
                    contribution of the bucket.

                    This field is a member of `oneof`_ ``details``.
                error (google.rpc.status_pb2.Status):
                    Output only. The error related to accessing
                    the details about the contribution of the
                    bucket.

                    This field is a member of `oneof`_ ``details``.
            """

            class Contribution(proto.Message):
                r"""Represents the contribution of the bucket towards the
                ``IntelligenceFinding``.

                Attributes:
                    top_prefixes (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution.Contribution.PrefixContribution]):
                        Output only. A list of the top object
                        prefixes driving the increase in operations.
                """

                class PrefixContribution(proto.Message):
                    r"""Represents the operation spike details for an object prefix.

                    Attributes:
                        prefix (str):
                            Output only. The object prefix. Format: ``a/b/c``, 'a/b/d',
                            etc.
                        percentage_increase (float):
                            Output only. The percentage increase in
                            operations for the object prefix.
                        total_operations_count (int):
                            Output only. The total count of operations
                            for the object prefix.
                    """

                    prefix: str = proto.Field(
                        proto.STRING,
                        number=1,
                    )
                    percentage_increase: float = proto.Field(
                        proto.DOUBLE,
                        number=2,
                    )
                    total_operations_count: int = proto.Field(
                        proto.INT64,
                        number=3,
                    )

                top_prefixes: MutableSequence[
                    "IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution.Contribution.PrefixContribution"
                ] = proto.RepeatedField(
                    proto.MESSAGE,
                    number=1,
                    message="IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution.Contribution.PrefixContribution",
                )

            bucket: str = proto.Field(
                proto.STRING,
                number=1,
            )
            percentage_increase: float = proto.Field(
                proto.DOUBLE,
                number=2,
            )
            total_operations_count: int = proto.Field(
                proto.INT64,
                number=3,
            )
            contribution: "IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution.Contribution" = proto.Field(
                proto.MESSAGE,
                number=4,
                oneof="details",
                message="IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution.Contribution",
            )
            error: status_pb2.Status = proto.Field(
                proto.MESSAGE,
                number=5,
                oneof="details",
                message=status_pb2.Status,
            )

        percentage_increase: float = proto.Field(
            proto.DOUBLE,
            number=1,
        )
        total_operations_count: int = proto.Field(
            proto.INT64,
            number=2,
        )
        top_buckets: MutableSequence[
            "IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution"
        ] = proto.RepeatedField(
            proto.MESSAGE,
            number=3,
            message="IntelligenceFinding.ColdlineAndArchivalStorageOperationsSpike.BucketContribution",
        )

    class CrossRegionEgressSpike(proto.Message):
        r"""Represents a finding about a spike in cross-region egress from Cloud
        Storage. This corresponds to the ``CROSS_REGION_EGRESS_SPIKE``
        finding type.

        Attributes:
            total_egress_bytes (int):
                Output only. The total cross-region egress
                volume in bytes across the project.
            percentage_increase (float):
                Output only. The percentage increase in
                cross-region egress across the project.
            top_buckets (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFinding.CrossRegionEgressSpike.BucketContribution]):
                Output only. A list of top buckets driving
                the increase in cross-region egress.
        """

        class BucketContribution(proto.Message):
            r"""Represents the cross-region egress spike details for a
            bucket.

            This message has `oneof`_ fields (mutually exclusive fields).
            For each oneof, at most one member field can be set at the same time.
            Setting any member of the oneof automatically clears all other
            members.

            .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

            Attributes:
                bucket (str):
                    Output only. The name of the bucket.
                total_egress_bytes (int):
                    Output only. The total cross-region egress
                    volume in bytes for the bucket.
                percentage_increase (float):
                    Output only. The percentage increase in
                    cross-region egress for the bucket.
                contribution (google.cloud.storage_control_v2.types.IntelligenceFinding.CrossRegionEgressSpike.BucketContribution.Contribution):
                    Output only. The details about the
                    contribution of the bucket.

                    This field is a member of `oneof`_ ``details``.
                error (google.rpc.status_pb2.Status):
                    Output only. The error related to accessing
                    the details about the contribution of the
                    bucket.

                    This field is a member of `oneof`_ ``details``.
            """

            class Contribution(proto.Message):
                r"""Represents the contribution of the bucket towards the
                ``IntelligenceFinding``.

                Attributes:
                    top_prefixes (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFinding.CrossRegionEgressSpike.BucketContribution.Contribution.PrefixContribution]):
                        Output only. A list of the top object
                        prefixes driving the increase in cross-region
                        egress.
                """

                class PrefixContribution(proto.Message):
                    r"""Represents the cross-region egress spike details for an
                    object prefix.

                    Attributes:
                        prefix (str):
                            Output only. The object prefix. Format: ``a/b/c``, 'a/b/d',
                            etc.
                        total_egress_bytes (int):
                            Output only. The total cross-region egress
                            volume in bytes from the object prefix.
                        percentage_increase (float):
                            Output only. The percentage increase in
                            cross-region egress for the object prefix.
                    """

                    prefix: str = proto.Field(
                        proto.STRING,
                        number=1,
                    )
                    total_egress_bytes: int = proto.Field(
                        proto.INT64,
                        number=2,
                    )
                    percentage_increase: float = proto.Field(
                        proto.DOUBLE,
                        number=3,
                    )

                top_prefixes: MutableSequence[
                    "IntelligenceFinding.CrossRegionEgressSpike.BucketContribution.Contribution.PrefixContribution"
                ] = proto.RepeatedField(
                    proto.MESSAGE,
                    number=1,
                    message="IntelligenceFinding.CrossRegionEgressSpike.BucketContribution.Contribution.PrefixContribution",
                )

            bucket: str = proto.Field(
                proto.STRING,
                number=1,
            )
            total_egress_bytes: int = proto.Field(
                proto.INT64,
                number=2,
            )
            percentage_increase: float = proto.Field(
                proto.DOUBLE,
                number=3,
            )
            contribution: "IntelligenceFinding.CrossRegionEgressSpike.BucketContribution.Contribution" = proto.Field(
                proto.MESSAGE,
                number=4,
                oneof="details",
                message="IntelligenceFinding.CrossRegionEgressSpike.BucketContribution.Contribution",
            )
            error: status_pb2.Status = proto.Field(
                proto.MESSAGE,
                number=5,
                oneof="details",
                message=status_pb2.Status,
            )

        total_egress_bytes: int = proto.Field(
            proto.INT64,
            number=1,
        )
        percentage_increase: float = proto.Field(
            proto.DOUBLE,
            number=2,
        )
        top_buckets: MutableSequence[
            "IntelligenceFinding.CrossRegionEgressSpike.BucketContribution"
        ] = proto.RepeatedField(
            proto.MESSAGE,
            number=3,
            message="IntelligenceFinding.CrossRegionEgressSpike.BucketContribution",
        )

    class ThrottledRequestSpike(proto.Message):
        r"""Represents a finding about a spike in throttled requests (429
        errors) within a project. This corresponds to the
        ``THROTTLED_REQUEST_SPIKE`` finding type.

        Attributes:
            throttled_requests (int):
                Output only. The count of throttled requests
                across the project.
            percentage_increase (float):
                Output only. The percentage increase in
                throttled requests across the project.
            top_buckets (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFinding.ThrottledRequestSpike.BucketContribution]):
                Output only. A list of top buckets driving
                the increase in throttled requests.
        """

        class BucketContribution(proto.Message):
            r"""Represents the throttled requests details for a bucket.

            This message has `oneof`_ fields (mutually exclusive fields).
            For each oneof, at most one member field can be set at the same time.
            Setting any member of the oneof automatically clears all other
            members.

            .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

            Attributes:
                bucket (str):
                    Output only. The name of the bucket.
                throttled_requests (int):
                    Output only. The count of throttled requests
                    for the bucket.
                percentage_increase (float):
                    Output only. The percentage increase in
                    throttled requests for the bucket.
                contribution (google.cloud.storage_control_v2.types.IntelligenceFinding.ThrottledRequestSpike.BucketContribution.Contribution):
                    Output only. The details about the
                    contribution of the bucket.

                    This field is a member of `oneof`_ ``details``.
                error (google.rpc.status_pb2.Status):
                    Output only. The error related to accessing
                    the details about the contribution of the
                    bucket.

                    This field is a member of `oneof`_ ``details``.
            """

            class Contribution(proto.Message):
                r"""Represents the contribution of the bucket towards the
                ``IntelligenceFinding``.

                Attributes:
                    top_prefixes (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFinding.ThrottledRequestSpike.BucketContribution.Contribution.PrefixContribution]):
                        Output only. A list of top object prefixes
                        driving the increase in throttled requests.
                """

                class PrefixContribution(proto.Message):
                    r"""Represents throttled requests details for an object prefix.

                    Attributes:
                        prefix (str):
                            Output only. The object prefix. Format: ``a/b/c``, 'a/b/d',
                            etc.
                        throttled_requests (int):
                            Output only. The count of throttled requests
                            for the object prefix.
                        percentage_increase (float):
                            Output only. The percentage increase in
                            throttled requests for the object prefix.
                    """

                    prefix: str = proto.Field(
                        proto.STRING,
                        number=1,
                    )
                    throttled_requests: int = proto.Field(
                        proto.INT64,
                        number=2,
                    )
                    percentage_increase: float = proto.Field(
                        proto.DOUBLE,
                        number=3,
                    )

                top_prefixes: MutableSequence[
                    "IntelligenceFinding.ThrottledRequestSpike.BucketContribution.Contribution.PrefixContribution"
                ] = proto.RepeatedField(
                    proto.MESSAGE,
                    number=1,
                    message="IntelligenceFinding.ThrottledRequestSpike.BucketContribution.Contribution.PrefixContribution",
                )

            bucket: str = proto.Field(
                proto.STRING,
                number=1,
            )
            throttled_requests: int = proto.Field(
                proto.INT64,
                number=2,
            )
            percentage_increase: float = proto.Field(
                proto.DOUBLE,
                number=3,
            )
            contribution: "IntelligenceFinding.ThrottledRequestSpike.BucketContribution.Contribution" = proto.Field(
                proto.MESSAGE,
                number=4,
                oneof="details",
                message="IntelligenceFinding.ThrottledRequestSpike.BucketContribution.Contribution",
            )
            error: status_pb2.Status = proto.Field(
                proto.MESSAGE,
                number=5,
                oneof="details",
                message=status_pb2.Status,
            )

        throttled_requests: int = proto.Field(
            proto.INT64,
            number=1,
        )
        percentage_increase: float = proto.Field(
            proto.DOUBLE,
            number=2,
        )
        top_buckets: MutableSequence[
            "IntelligenceFinding.ThrottledRequestSpike.BucketContribution"
        ] = proto.RepeatedField(
            proto.MESSAGE,
            number=3,
            message="IntelligenceFinding.ThrottledRequestSpike.BucketContribution",
        )

    class StorageGrowthAboveTrend(proto.Message):
        r"""Represents a finding about a storage growth above the expected
        trend. This corresponds to the ``STORAGE_GROWTH_ABOVE_TREND``
        finding type.

        Attributes:
            total_storage_growth_bytes (int):
                Output only. The total storage growth in
                bytes.
            percentage_increase (float):
                Output only. The percentage increase in
                storage growth.
            top_buckets (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFinding.StorageGrowthAboveTrend.BucketContribution]):
                Output only. A list of top buckets driving
                the increase in storage growth.
        """

        class BucketContribution(proto.Message):
            r"""Represents the storage growth details for a bucket.

            .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

            Attributes:
                bucket (str):
                    Output only. The name of the bucket.
                total_storage_growth_bytes (int):
                    Output only. The total storage growth in
                    bytes for the bucket.
                percentage_increase (float):
                    Output only. The percentage increase in
                    storage growth for the bucket.
                error (google.rpc.status_pb2.Status):
                    Output only. The error related to accessing
                    the details about the contribution of the
                    bucket.

                    This field is a member of `oneof`_ ``details``.
            """

            bucket: str = proto.Field(
                proto.STRING,
                number=1,
            )
            total_storage_growth_bytes: int = proto.Field(
                proto.INT64,
                number=2,
            )
            percentage_increase: float = proto.Field(
                proto.DOUBLE,
                number=3,
            )
            error: status_pb2.Status = proto.Field(
                proto.MESSAGE,
                number=5,
                oneof="details",
                message=status_pb2.Status,
            )

        total_storage_growth_bytes: int = proto.Field(
            proto.INT64,
            number=1,
        )
        percentage_increase: float = proto.Field(
            proto.DOUBLE,
            number=2,
        )
        top_buckets: MutableSequence[
            "IntelligenceFinding.StorageGrowthAboveTrend.BucketContribution"
        ] = proto.RepeatedField(
            proto.MESSAGE,
            number=3,
            message="IntelligenceFinding.StorageGrowthAboveTrend.BucketContribution",
        )

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    description: str = proto.Field(
        proto.STRING,
        number=2,
    )
    type_: "FindingType" = proto.Field(
        proto.ENUM,
        number=3,
        enum="FindingType",
    )
    category: "FindingCategory" = proto.Field(
        proto.ENUM,
        number=4,
        enum="FindingCategory",
    )
    severity: "FindingSeverity" = proto.Field(
        proto.ENUM,
        number=5,
        enum="FindingSeverity",
    )
    create_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=6,
        message=timestamp_pb2.Timestamp,
    )
    update_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=7,
        message=timestamp_pb2.Timestamp,
    )
    target_resource: str = proto.Field(
        proto.STRING,
        number=8,
    )
    associated_resources: MutableSequence[str] = proto.RepeatedField(
        proto.STRING,
        number=9,
    )
    observation_period: interval_pb2.Interval = proto.Field(
        proto.MESSAGE,
        number=10,
        message=interval_pb2.Interval,
    )
    coldline_and_archival_storage_operations_spike: ColdlineAndArchivalStorageOperationsSpike = proto.Field(
        proto.MESSAGE,
        number=11,
        oneof="intelligence_finding_details",
        message=ColdlineAndArchivalStorageOperationsSpike,
    )
    throttled_requests_spike: ThrottledRequestSpike = proto.Field(
        proto.MESSAGE,
        number=12,
        oneof="intelligence_finding_details",
        message=ThrottledRequestSpike,
    )
    cross_region_egress_spike: CrossRegionEgressSpike = proto.Field(
        proto.MESSAGE,
        number=13,
        oneof="intelligence_finding_details",
        message=CrossRegionEgressSpike,
    )
    storage_growth_above_trend: StorageGrowthAboveTrend = proto.Field(
        proto.MESSAGE,
        number=14,
        oneof="intelligence_finding_details",
        message=StorageGrowthAboveTrend,
    )


class IntelligenceFindingRevision(proto.Message):
    r"""An ``IntelligenceFindingRevision`` represents a specific revision of
    an ``IntelligenceFinding`` resource.

    Attributes:
        name (str):
            Identifier. The resource name of
            ``IntelligenceFindingRevision``. Format:
            ``projects/{project}/locations/{location}/intelligenceFindings/{intelligence_finding}/revisions/{revision}``
        snapshot (google.cloud.storage_control_v2.types.IntelligenceFinding):
            Output only. The snapshot of the ``IntelligenceFinding`` at
            the time the revision was created. This field contains the
            full finding details as they existed for the revision.
        create_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The timestamp when the revision
            was created.
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )
    snapshot: "IntelligenceFinding" = proto.Field(
        proto.MESSAGE,
        number=2,
        message="IntelligenceFinding",
    )
    create_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=3,
        message=timestamp_pb2.Timestamp,
    )


class GetIntelligenceFindingRequest(proto.Message):
    r"""Request message to get the ``IntelligenceFinding`` resource
    associated with a project.

    Attributes:
        name (str):
            Required. The name of the ``IntelligenceFinding`` resource.

            Format:
            ``projects/{project}/locations/{location}/intelligenceFindings/{intelligence_finding}``
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )


class ListIntelligenceFindingsRequest(proto.Message):
    r"""Request message to list ``IntelligenceFinding`` resources associated
    with a project.

    Attributes:
        parent (str):
            Required. The parent of the ``IntelligenceFinding``
            resource.

            Format: ``projects/{project}/locations/{location}``
        filter (str):
            Optional. The filter expression to be applied. Supports
            filtering by ``type`` and ``associated_resources``.
        page_size (int):
            Optional. The maximum number of ``IntelligenceFinding``
            resources to return.

            The maximum value is ``100``; values above ``100`` will be
            coerced to ``100``. The default value is ``100``.
        page_token (str):
            Optional. A page token, received from a previous
            ``ListIntelligenceFindings`` call. Provide this to retrieve
            the subsequent page.

            When paginating, all other parameters provided to
            ``ListIntelligenceFindings`` must match the call that
            provided the page token.
    """

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    filter: str = proto.Field(
        proto.STRING,
        number=2,
    )
    page_size: int = proto.Field(
        proto.INT32,
        number=3,
    )
    page_token: str = proto.Field(
        proto.STRING,
        number=4,
    )


class ListIntelligenceFindingsResponse(proto.Message):
    r"""Response message to list the ``IntelligenceFinding`` resources
    associated with a project.

    Attributes:
        intelligence_findings (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFinding]):
            The ``IntelligenceFinding`` resources from the specified
            project.
        next_page_token (str):
            A token to retrieve the next page of results. Pass this
            value in the ``page_token`` field in the subsequent call.
    """

    @property
    def raw_page(self):
        return self

    intelligence_findings: MutableSequence["IntelligenceFinding"] = proto.RepeatedField(
        proto.MESSAGE,
        number=1,
        message="IntelligenceFinding",
    )
    next_page_token: str = proto.Field(
        proto.STRING,
        number=2,
    )


class SummarizeIntelligenceFindingsRequest(proto.Message):
    r"""Request message to summarize the intelligence findings for
    the specified scope(org, folder or project).

    Attributes:
        parent (str):
            Required. The scope to summarize the findings for. Format:

            - ``organizations/{organization}/locations/{location}``
            - ``folders/{folder}/locations/{location}``
            - ``projects/{project}/locations/{location}``
        resource_scope (google.cloud.storage_control_v2.types.SummarizeIntelligenceFindingsRequest.ResourceScope):
            Optional. Determines the granularity of the findings when
            the ``parent`` is an organization or folder.

            - ``PARENT`` (or not set): A single summary is returned for
              each insight type, aggregated across the entire ``parent``
              scope.
            - ``PROJECT``: A separate summary is returned for each
              insight type for every project within the ``parent``
              scope.

            The only supported values are ``PARENT`` and ``PROJECT``. If
            no value is specified, the API behaviour defaults to the
            ``PARENT``.
        filter (str):
            Optional. The filter expression, following
            AIP-160. Supports filtering by FindingType.
        page_size (int):
            Optional. The maximum number of findings to return.

            The maximum value is ``100``; values above ``100`` will be
            coerced to ``100``. The default value is ``100``.
        page_token (str):
            Optional. A page token, received from a previous
            ``SummarizeIntelligenceFindings`` call. Provide this to
            retrieve the subsequent page.

            When paginating, all other parameters provided to
            ``SummarizeIntelligenceFindings`` must match the call that
            provided the page token.
    """

    class ResourceScope(proto.Enum):
        r"""The list of resource scopes.

        Values:
            RESOURCE_SCOPE_UNSPECIFIED (0):
                The default behavior. Falls back to PARENT
                behaviour
            PARENT (1):
                Summaries are aggregated at the level of the ``parent``
                resource.
            PROJECT (2):
                Summaries are broken down by each project within the
                ``parent`` scope.
        """

        RESOURCE_SCOPE_UNSPECIFIED = 0
        PARENT = 1
        PROJECT = 2

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    resource_scope: ResourceScope = proto.Field(
        proto.ENUM,
        number=2,
        enum=ResourceScope,
    )
    filter: str = proto.Field(
        proto.STRING,
        number=3,
    )
    page_size: int = proto.Field(
        proto.INT32,
        number=4,
    )
    page_token: str = proto.Field(
        proto.STRING,
        number=5,
    )


class SummarizeIntelligenceFindingsResponse(proto.Message):
    r"""Response message to summarize the intelligence findings for a
    specified scope(org, folder or project).

    Attributes:
        finding_summaries (MutableSequence[google.cloud.storage_control_v2.types.FindingSummary]):
            The list of ``FindingSummary`` summaries.
        next_page_token (str):
            A token to retrieve the next page of results. Pass this
            value in the ``page_token`` field in the subsequent call to
            ``SummarizeIntelligenceFindings`` to retrieve the next page
            of results.
    """

    @property
    def raw_page(self):
        return self

    finding_summaries: MutableSequence["FindingSummary"] = proto.RepeatedField(
        proto.MESSAGE,
        number=1,
        message="FindingSummary",
    )
    next_page_token: str = proto.Field(
        proto.STRING,
        number=2,
    )


class GetIntelligenceFindingRevisionRequest(proto.Message):
    r"""Request message to get the ``IntelligenceFindingRevision`` resource
    associated with a project.

    Attributes:
        name (str):
            Required. The name of the ``IntelligenceFindingRevision``
            resource.

            Format:
            -------

            ``projects/{project}/locations/{location}/intelligenceFindings/{intelligence_finding}/revisions/{revision}``
    """

    name: str = proto.Field(
        proto.STRING,
        number=1,
    )


class ListIntelligenceFindingRevisionsRequest(proto.Message):
    r"""Request message to list ``IntelligenceFindingRevision`` resources
    associated with a project.

    Attributes:
        parent (str):
            Required. The parent of the ``IntelligenceFindingRevision``
            resource.

            Format:
            -------

            ``projects/{project}/locations/{location}/intelligenceFindings/{intelligence_finding}``
        page_size (int):
            Optional. The maximum number of
            ``IntelligenceFindingRevision`` resources to return.

            The maximum value is ``100``; values above ``100`` will be
            coerced to ``100``. The default value is ``100``.
        page_token (str):
            Optional. A page token, received from a previous
            ``ListIntelligenceFindingRevisions`` call. Provide this to
            retrieve the subsequent page.
    """

    parent: str = proto.Field(
        proto.STRING,
        number=1,
    )
    page_size: int = proto.Field(
        proto.INT32,
        number=2,
    )
    page_token: str = proto.Field(
        proto.STRING,
        number=3,
    )


class ListIntelligenceFindingRevisionsResponse(proto.Message):
    r"""Response message to list ``IntelligenceFindingRevision`` resources
    associated with a project.

    Attributes:
        intelligence_finding_revisions (MutableSequence[google.cloud.storage_control_v2.types.IntelligenceFindingRevision]):
            The ``IntelligenceFindingRevision`` resources from the
            specified project.
        next_page_token (str):
            A token that can be sent as ``page_token`` to retrieve the
            next page.
    """

    @property
    def raw_page(self):
        return self

    intelligence_finding_revisions: MutableSequence["IntelligenceFindingRevision"] = (
        proto.RepeatedField(
            proto.MESSAGE,
            number=1,
            message="IntelligenceFindingRevision",
        )
    )
    next_page_token: str = proto.Field(
        proto.STRING,
        number=2,
    )


class FindingSummary(proto.Message):
    r"""A summary of findings generated for an organization, a
    folder, or a project.

    Attributes:
        type_ (google.cloud.storage_control_v2.types.FindingType):
            Output only. The type of the finding.
        category (google.cloud.storage_control_v2.types.FindingCategory):
            Output only. The category of finding.
        target_resource (str):
            Output only. The fully qualified Cloud resource name for
            which this summary was generated. eg:
            ``//cloudresourcemanager.googleapis.com/projects/p1``
        create_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The creation time of the
            earliest finding that this summary is based on.
        update_time (google.protobuf.timestamp_pb2.Timestamp):
            Output only. The time of the most recent
            update among all the findings that this summary
            is based on.
        severity (google.cloud.storage_control_v2.types.FindingSeverity):
            Severity of the finding.
        summary_details (MutableSequence[google.cloud.storage_control_v2.types.FindingSummary.SummaryDetails]):
            Output only. List of ``SummaryDetails``.
    """

    class SummaryDetails(proto.Message):
        r"""Details about the ``FindingSummary`` resource.

        This message has `oneof`_ fields (mutually exclusive fields).
        For each oneof, at most one member field can be set at the same time.
        Setting any member of the oneof automatically clears all other
        members.

        .. _oneof: https://proto-plus-python.readthedocs.io/en/stable/fields.html#oneofs-mutually-exclusive-fields

        Attributes:
            count (int):
                The count of impacted resources.

                This field is a member of `oneof`_ ``magnitude``.
            percentage (float):
                The percentage of impacted resources.

                This field is a member of `oneof`_ ``magnitude``.
            resource_type (google.cloud.storage_control_v2.types.FindingSummary.SummaryDetails.ResourceType):
                Output only. The type of Cloud resource this
                summary detail applies to.
            description (str):
                Output only. A short description about the
                FindingSummary
        """

        class ResourceType(proto.Enum):
            r"""The list of resource types.

            Values:
                RESOURCE_TYPE_UNSPECIFIED (0):
                    Resource type is unspecified.
                PROJECT (1):
                    Resource type is project.
                BUCKET (2):
                    Resource type is bucket.
            """

            RESOURCE_TYPE_UNSPECIFIED = 0
            PROJECT = 1
            BUCKET = 2

        count: int = proto.Field(
            proto.INT64,
            number=1,
            oneof="magnitude",
        )
        percentage: float = proto.Field(
            proto.FLOAT,
            number=2,
            oneof="magnitude",
        )
        resource_type: "FindingSummary.SummaryDetails.ResourceType" = proto.Field(
            proto.ENUM,
            number=3,
            enum="FindingSummary.SummaryDetails.ResourceType",
        )
        description: str = proto.Field(
            proto.STRING,
            number=4,
        )

    type_: "FindingType" = proto.Field(
        proto.ENUM,
        number=1,
        enum="FindingType",
    )
    category: "FindingCategory" = proto.Field(
        proto.ENUM,
        number=2,
        enum="FindingCategory",
    )
    target_resource: str = proto.Field(
        proto.STRING,
        number=4,
    )
    create_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=5,
        message=timestamp_pb2.Timestamp,
    )
    update_time: timestamp_pb2.Timestamp = proto.Field(
        proto.MESSAGE,
        number=6,
        message=timestamp_pb2.Timestamp,
    )
    severity: "FindingSeverity" = proto.Field(
        proto.ENUM,
        number=7,
        enum="FindingSeverity",
    )
    summary_details: MutableSequence[SummaryDetails] = proto.RepeatedField(
        proto.MESSAGE,
        number=8,
        message=SummaryDetails,
    )


__all__ = tuple(sorted(__protobuf__.manifest))
