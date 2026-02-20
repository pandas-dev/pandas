# -*- coding: utf-8 -*-
# Copyright 2025 Google LLC
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
import proto  # type: ignore

__protobuf__ = proto.module(
    package="google.storage.control.v2",
    manifest={
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
    },
)


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


__all__ = tuple(sorted(__protobuf__.manifest))
