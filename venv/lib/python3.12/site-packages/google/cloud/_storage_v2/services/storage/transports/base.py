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
import abc
from typing import Awaitable, Callable, Dict, Optional, Sequence, Union

from google.cloud._storage_v2 import gapic_version as package_version

import google.auth  # type: ignore
import google.api_core
from google.api_core import exceptions as core_exceptions
from google.api_core import gapic_v1
from google.api_core import retry as retries
from google.auth import credentials as ga_credentials  # type: ignore
from google.oauth2 import service_account  # type: ignore
import google.protobuf

from google.cloud._storage_v2.types import storage
from google.iam.v1 import iam_policy_pb2  # type: ignore
from google.iam.v1 import policy_pb2  # type: ignore
from google.longrunning import operations_pb2  # type: ignore
from google.protobuf import empty_pb2  # type: ignore

DEFAULT_CLIENT_INFO = gapic_v1.client_info.ClientInfo(
    gapic_version=package_version.__version__
)

if hasattr(DEFAULT_CLIENT_INFO, "protobuf_runtime_version"):  # pragma: NO COVER
    DEFAULT_CLIENT_INFO.protobuf_runtime_version = google.protobuf.__version__


class StorageTransport(abc.ABC):
    """Abstract transport class for Storage."""

    AUTH_SCOPES = (
        "https://www.googleapis.com/auth/cloud-platform",
        "https://www.googleapis.com/auth/cloud-platform.read-only",
        "https://www.googleapis.com/auth/devstorage.full_control",
        "https://www.googleapis.com/auth/devstorage.read_only",
        "https://www.googleapis.com/auth/devstorage.read_write",
    )

    DEFAULT_HOST: str = "storage.googleapis.com"

    def __init__(
        self,
        *,
        host: str = DEFAULT_HOST,
        credentials: Optional[ga_credentials.Credentials] = None,
        credentials_file: Optional[str] = None,
        scopes: Optional[Sequence[str]] = None,
        quota_project_id: Optional[str] = None,
        client_info: gapic_v1.client_info.ClientInfo = DEFAULT_CLIENT_INFO,
        always_use_jwt_access: Optional[bool] = False,
        api_audience: Optional[str] = None,
        **kwargs,
    ) -> None:
        """Instantiate the transport.

        Args:
            host (Optional[str]):
                 The hostname to connect to (default: 'storage.googleapis.com').
            credentials (Optional[google.auth.credentials.Credentials]): The
                authorization credentials to attach to requests. These
                credentials identify the application to the service; if none
                are specified, the client will attempt to ascertain the
                credentials from the environment.
            credentials_file (Optional[str]): Deprecated. A file with credentials that can
                be loaded with :func:`google.auth.load_credentials_from_file`.
                This argument is mutually exclusive with credentials. This argument will be
                removed in the next major version of this library.
            scopes (Optional[Sequence[str]]): A list of scopes.
            quota_project_id (Optional[str]): An optional project to use for billing
                and quota.
            client_info (google.api_core.gapic_v1.client_info.ClientInfo):
                The client info used to send a user-agent string along with
                API requests. If ``None``, then default info will be used.
                Generally, you only need to set this if you're developing
                your own client library.
            always_use_jwt_access (Optional[bool]): Whether self signed JWT should
                be used for service account credentials.
        """

        scopes_kwargs = {"scopes": scopes, "default_scopes": self.AUTH_SCOPES}

        # Save the scopes.
        self._scopes = scopes
        if not hasattr(self, "_ignore_credentials"):
            self._ignore_credentials: bool = False

        # If no credentials are provided, then determine the appropriate
        # defaults.
        if credentials and credentials_file:
            raise core_exceptions.DuplicateCredentialArgs(
                "'credentials_file' and 'credentials' are mutually exclusive"
            )

        if credentials_file is not None:
            credentials, _ = google.auth.load_credentials_from_file(
                credentials_file, **scopes_kwargs, quota_project_id=quota_project_id
            )
        elif credentials is None and not self._ignore_credentials:
            credentials, _ = google.auth.default(
                **scopes_kwargs, quota_project_id=quota_project_id
            )
            # Don't apply audience if the credentials file passed from user.
            if hasattr(credentials, "with_gdch_audience"):
                credentials = credentials.with_gdch_audience(
                    api_audience if api_audience else host
                )

        # If the credentials are service account credentials, then always try to use self signed JWT.
        if (
            always_use_jwt_access
            and isinstance(credentials, service_account.Credentials)
            and hasattr(service_account.Credentials, "with_always_use_jwt_access")
        ):
            credentials = credentials.with_always_use_jwt_access(True)

        # Save the credentials.
        self._credentials = credentials

        # Save the hostname. Default to port 443 (HTTPS) if none is specified.
        if ":" not in host:
            host += ":443"
        self._host = host

    @property
    def host(self):
        return self._host

    def _prep_wrapped_messages(self, client_info):
        # Precompute the wrapped methods.
        self._wrapped_methods = {
            self.delete_bucket: gapic_v1.method.wrap_method(
                self.delete_bucket,
                default_timeout=None,
                client_info=client_info,
            ),
            self.get_bucket: gapic_v1.method.wrap_method(
                self.get_bucket,
                default_timeout=None,
                client_info=client_info,
            ),
            self.create_bucket: gapic_v1.method.wrap_method(
                self.create_bucket,
                default_timeout=None,
                client_info=client_info,
            ),
            self.list_buckets: gapic_v1.method.wrap_method(
                self.list_buckets,
                default_timeout=None,
                client_info=client_info,
            ),
            self.lock_bucket_retention_policy: gapic_v1.method.wrap_method(
                self.lock_bucket_retention_policy,
                default_timeout=None,
                client_info=client_info,
            ),
            self.get_iam_policy: gapic_v1.method.wrap_method(
                self.get_iam_policy,
                default_timeout=None,
                client_info=client_info,
            ),
            self.set_iam_policy: gapic_v1.method.wrap_method(
                self.set_iam_policy,
                default_timeout=None,
                client_info=client_info,
            ),
            self.test_iam_permissions: gapic_v1.method.wrap_method(
                self.test_iam_permissions,
                default_timeout=None,
                client_info=client_info,
            ),
            self.update_bucket: gapic_v1.method.wrap_method(
                self.update_bucket,
                default_timeout=None,
                client_info=client_info,
            ),
            self.compose_object: gapic_v1.method.wrap_method(
                self.compose_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.delete_object: gapic_v1.method.wrap_method(
                self.delete_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.restore_object: gapic_v1.method.wrap_method(
                self.restore_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.cancel_resumable_write: gapic_v1.method.wrap_method(
                self.cancel_resumable_write,
                default_timeout=None,
                client_info=client_info,
            ),
            self.get_object: gapic_v1.method.wrap_method(
                self.get_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.read_object: gapic_v1.method.wrap_method(
                self.read_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.bidi_read_object: gapic_v1.method.wrap_method(
                self.bidi_read_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.update_object: gapic_v1.method.wrap_method(
                self.update_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.write_object: gapic_v1.method.wrap_method(
                self.write_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.bidi_write_object: gapic_v1.method.wrap_method(
                self.bidi_write_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.list_objects: gapic_v1.method.wrap_method(
                self.list_objects,
                default_timeout=None,
                client_info=client_info,
            ),
            self.rewrite_object: gapic_v1.method.wrap_method(
                self.rewrite_object,
                default_timeout=None,
                client_info=client_info,
            ),
            self.start_resumable_write: gapic_v1.method.wrap_method(
                self.start_resumable_write,
                default_timeout=None,
                client_info=client_info,
            ),
            self.query_write_status: gapic_v1.method.wrap_method(
                self.query_write_status,
                default_timeout=None,
                client_info=client_info,
            ),
            self.move_object: gapic_v1.method.wrap_method(
                self.move_object,
                default_timeout=None,
                client_info=client_info,
            ),
        }

    def close(self):
        """Closes resources associated with the transport.

        .. warning::
             Only call this method if the transport is NOT shared
             with other clients - this may cause errors in other clients!
        """
        raise NotImplementedError()

    @property
    def delete_bucket(
        self,
    ) -> Callable[
        [storage.DeleteBucketRequest],
        Union[empty_pb2.Empty, Awaitable[empty_pb2.Empty]],
    ]:
        raise NotImplementedError()

    @property
    def get_bucket(
        self,
    ) -> Callable[
        [storage.GetBucketRequest], Union[storage.Bucket, Awaitable[storage.Bucket]]
    ]:
        raise NotImplementedError()

    @property
    def create_bucket(
        self,
    ) -> Callable[
        [storage.CreateBucketRequest], Union[storage.Bucket, Awaitable[storage.Bucket]]
    ]:
        raise NotImplementedError()

    @property
    def list_buckets(
        self,
    ) -> Callable[
        [storage.ListBucketsRequest],
        Union[storage.ListBucketsResponse, Awaitable[storage.ListBucketsResponse]],
    ]:
        raise NotImplementedError()

    @property
    def lock_bucket_retention_policy(
        self,
    ) -> Callable[
        [storage.LockBucketRetentionPolicyRequest],
        Union[storage.Bucket, Awaitable[storage.Bucket]],
    ]:
        raise NotImplementedError()

    @property
    def get_iam_policy(
        self,
    ) -> Callable[
        [iam_policy_pb2.GetIamPolicyRequest],
        Union[policy_pb2.Policy, Awaitable[policy_pb2.Policy]],
    ]:
        raise NotImplementedError()

    @property
    def set_iam_policy(
        self,
    ) -> Callable[
        [iam_policy_pb2.SetIamPolicyRequest],
        Union[policy_pb2.Policy, Awaitable[policy_pb2.Policy]],
    ]:
        raise NotImplementedError()

    @property
    def test_iam_permissions(
        self,
    ) -> Callable[
        [iam_policy_pb2.TestIamPermissionsRequest],
        Union[
            iam_policy_pb2.TestIamPermissionsResponse,
            Awaitable[iam_policy_pb2.TestIamPermissionsResponse],
        ],
    ]:
        raise NotImplementedError()

    @property
    def update_bucket(
        self,
    ) -> Callable[
        [storage.UpdateBucketRequest], Union[storage.Bucket, Awaitable[storage.Bucket]]
    ]:
        raise NotImplementedError()

    @property
    def compose_object(
        self,
    ) -> Callable[
        [storage.ComposeObjectRequest], Union[storage.Object, Awaitable[storage.Object]]
    ]:
        raise NotImplementedError()

    @property
    def delete_object(
        self,
    ) -> Callable[
        [storage.DeleteObjectRequest],
        Union[empty_pb2.Empty, Awaitable[empty_pb2.Empty]],
    ]:
        raise NotImplementedError()

    @property
    def restore_object(
        self,
    ) -> Callable[
        [storage.RestoreObjectRequest], Union[storage.Object, Awaitable[storage.Object]]
    ]:
        raise NotImplementedError()

    @property
    def cancel_resumable_write(
        self,
    ) -> Callable[
        [storage.CancelResumableWriteRequest],
        Union[
            storage.CancelResumableWriteResponse,
            Awaitable[storage.CancelResumableWriteResponse],
        ],
    ]:
        raise NotImplementedError()

    @property
    def get_object(
        self,
    ) -> Callable[
        [storage.GetObjectRequest], Union[storage.Object, Awaitable[storage.Object]]
    ]:
        raise NotImplementedError()

    @property
    def read_object(
        self,
    ) -> Callable[
        [storage.ReadObjectRequest],
        Union[storage.ReadObjectResponse, Awaitable[storage.ReadObjectResponse]],
    ]:
        raise NotImplementedError()

    @property
    def bidi_read_object(
        self,
    ) -> Callable[
        [storage.BidiReadObjectRequest],
        Union[
            storage.BidiReadObjectResponse, Awaitable[storage.BidiReadObjectResponse]
        ],
    ]:
        raise NotImplementedError()

    @property
    def update_object(
        self,
    ) -> Callable[
        [storage.UpdateObjectRequest], Union[storage.Object, Awaitable[storage.Object]]
    ]:
        raise NotImplementedError()

    @property
    def write_object(
        self,
    ) -> Callable[
        [storage.WriteObjectRequest],
        Union[storage.WriteObjectResponse, Awaitable[storage.WriteObjectResponse]],
    ]:
        raise NotImplementedError()

    @property
    def bidi_write_object(
        self,
    ) -> Callable[
        [storage.BidiWriteObjectRequest],
        Union[
            storage.BidiWriteObjectResponse, Awaitable[storage.BidiWriteObjectResponse]
        ],
    ]:
        raise NotImplementedError()

    @property
    def list_objects(
        self,
    ) -> Callable[
        [storage.ListObjectsRequest],
        Union[storage.ListObjectsResponse, Awaitable[storage.ListObjectsResponse]],
    ]:
        raise NotImplementedError()

    @property
    def rewrite_object(
        self,
    ) -> Callable[
        [storage.RewriteObjectRequest],
        Union[storage.RewriteResponse, Awaitable[storage.RewriteResponse]],
    ]:
        raise NotImplementedError()

    @property
    def start_resumable_write(
        self,
    ) -> Callable[
        [storage.StartResumableWriteRequest],
        Union[
            storage.StartResumableWriteResponse,
            Awaitable[storage.StartResumableWriteResponse],
        ],
    ]:
        raise NotImplementedError()

    @property
    def query_write_status(
        self,
    ) -> Callable[
        [storage.QueryWriteStatusRequest],
        Union[
            storage.QueryWriteStatusResponse,
            Awaitable[storage.QueryWriteStatusResponse],
        ],
    ]:
        raise NotImplementedError()

    @property
    def move_object(
        self,
    ) -> Callable[
        [storage.MoveObjectRequest], Union[storage.Object, Awaitable[storage.Object]]
    ]:
        raise NotImplementedError()

    @property
    def kind(self) -> str:
        raise NotImplementedError()


__all__ = ("StorageTransport",)
