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

import google.api_core
from google.api_core import exceptions as core_exceptions
from google.api_core import gapic_v1, operations_v1
from google.api_core import retry as retries
import google.auth  # type: ignore
from google.auth import credentials as ga_credentials  # type: ignore
import google.iam.v1.iam_policy_pb2 as iam_policy_pb2  # type: ignore
import google.iam.v1.policy_pb2 as policy_pb2  # type: ignore
from google.longrunning import operations_pb2  # type: ignore
from google.oauth2 import service_account  # type: ignore
import google.protobuf
import google.protobuf.empty_pb2 as empty_pb2  # type: ignore

from google.cloud.storage_control_v2 import gapic_version as package_version
from google.cloud.storage_control_v2.types import storage_control

DEFAULT_CLIENT_INFO = gapic_v1.client_info.ClientInfo(
    gapic_version=package_version.__version__
)

if hasattr(DEFAULT_CLIENT_INFO, "protobuf_runtime_version"):  # pragma: NO COVER
    DEFAULT_CLIENT_INFO.protobuf_runtime_version = google.protobuf.__version__


class StorageControlTransport(abc.ABC):
    """Abstract transport class for StorageControl."""

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
                credentials_file,
                scopes=scopes,
                quota_project_id=quota_project_id,
                default_scopes=self.AUTH_SCOPES,
            )
        elif credentials is None and not self._ignore_credentials:
            credentials, _ = google.auth.default(
                scopes=scopes,
                quota_project_id=quota_project_id,
                default_scopes=self.AUTH_SCOPES,
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
            self.create_folder: gapic_v1.method.wrap_method(
                self.create_folder,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.delete_folder: gapic_v1.method.wrap_method(
                self.delete_folder,
                default_timeout=None,
                client_info=client_info,
            ),
            self.get_folder: gapic_v1.method.wrap_method(
                self.get_folder,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.list_folders: gapic_v1.method.wrap_method(
                self.list_folders,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.rename_folder: gapic_v1.method.wrap_method(
                self.rename_folder,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.delete_folder_recursive: gapic_v1.method.wrap_method(
                self.delete_folder_recursive,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.get_storage_layout: gapic_v1.method.wrap_method(
                self.get_storage_layout,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.create_managed_folder: gapic_v1.method.wrap_method(
                self.create_managed_folder,
                default_timeout=None,
                client_info=client_info,
            ),
            self.delete_managed_folder: gapic_v1.method.wrap_method(
                self.delete_managed_folder,
                default_timeout=None,
                client_info=client_info,
            ),
            self.get_managed_folder: gapic_v1.method.wrap_method(
                self.get_managed_folder,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.list_managed_folders: gapic_v1.method.wrap_method(
                self.list_managed_folders,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.create_anywhere_cache: gapic_v1.method.wrap_method(
                self.create_anywhere_cache,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.update_anywhere_cache: gapic_v1.method.wrap_method(
                self.update_anywhere_cache,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.disable_anywhere_cache: gapic_v1.method.wrap_method(
                self.disable_anywhere_cache,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.pause_anywhere_cache: gapic_v1.method.wrap_method(
                self.pause_anywhere_cache,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.resume_anywhere_cache: gapic_v1.method.wrap_method(
                self.resume_anywhere_cache,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.get_anywhere_cache: gapic_v1.method.wrap_method(
                self.get_anywhere_cache,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.list_anywhere_caches: gapic_v1.method.wrap_method(
                self.list_anywhere_caches,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.get_project_intelligence_config: gapic_v1.method.wrap_method(
                self.get_project_intelligence_config,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.update_project_intelligence_config: gapic_v1.method.wrap_method(
                self.update_project_intelligence_config,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.get_folder_intelligence_config: gapic_v1.method.wrap_method(
                self.get_folder_intelligence_config,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.update_folder_intelligence_config: gapic_v1.method.wrap_method(
                self.update_folder_intelligence_config,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.get_organization_intelligence_config: gapic_v1.method.wrap_method(
                self.get_organization_intelligence_config,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
                client_info=client_info,
            ),
            self.update_organization_intelligence_config: gapic_v1.method.wrap_method(
                self.update_organization_intelligence_config,
                default_retry=retries.Retry(
                    initial=1.0,
                    maximum=60.0,
                    multiplier=2,
                    predicate=retries.if_exception_type(
                        core_exceptions.DeadlineExceeded,
                        core_exceptions.InternalServerError,
                        core_exceptions.ResourceExhausted,
                        core_exceptions.ServiceUnavailable,
                        core_exceptions.Unknown,
                    ),
                    deadline=60.0,
                ),
                default_timeout=60.0,
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
        }

    def close(self):
        """Closes resources associated with the transport.

        .. warning::
             Only call this method if the transport is NOT shared
             with other clients - this may cause errors in other clients!
        """
        raise NotImplementedError()

    @property
    def operations_client(self):
        """Return the client designed to process long-running operations."""
        raise NotImplementedError()

    @property
    def create_folder(
        self,
    ) -> Callable[
        [storage_control.CreateFolderRequest],
        Union[storage_control.Folder, Awaitable[storage_control.Folder]],
    ]:
        raise NotImplementedError()

    @property
    def delete_folder(
        self,
    ) -> Callable[
        [storage_control.DeleteFolderRequest],
        Union[empty_pb2.Empty, Awaitable[empty_pb2.Empty]],
    ]:
        raise NotImplementedError()

    @property
    def get_folder(
        self,
    ) -> Callable[
        [storage_control.GetFolderRequest],
        Union[storage_control.Folder, Awaitable[storage_control.Folder]],
    ]:
        raise NotImplementedError()

    @property
    def list_folders(
        self,
    ) -> Callable[
        [storage_control.ListFoldersRequest],
        Union[
            storage_control.ListFoldersResponse,
            Awaitable[storage_control.ListFoldersResponse],
        ],
    ]:
        raise NotImplementedError()

    @property
    def rename_folder(
        self,
    ) -> Callable[
        [storage_control.RenameFolderRequest],
        Union[operations_pb2.Operation, Awaitable[operations_pb2.Operation]],
    ]:
        raise NotImplementedError()

    @property
    def delete_folder_recursive(
        self,
    ) -> Callable[
        [storage_control.DeleteFolderRecursiveRequest],
        Union[operations_pb2.Operation, Awaitable[operations_pb2.Operation]],
    ]:
        raise NotImplementedError()

    @property
    def get_storage_layout(
        self,
    ) -> Callable[
        [storage_control.GetStorageLayoutRequest],
        Union[storage_control.StorageLayout, Awaitable[storage_control.StorageLayout]],
    ]:
        raise NotImplementedError()

    @property
    def create_managed_folder(
        self,
    ) -> Callable[
        [storage_control.CreateManagedFolderRequest],
        Union[storage_control.ManagedFolder, Awaitable[storage_control.ManagedFolder]],
    ]:
        raise NotImplementedError()

    @property
    def delete_managed_folder(
        self,
    ) -> Callable[
        [storage_control.DeleteManagedFolderRequest],
        Union[empty_pb2.Empty, Awaitable[empty_pb2.Empty]],
    ]:
        raise NotImplementedError()

    @property
    def get_managed_folder(
        self,
    ) -> Callable[
        [storage_control.GetManagedFolderRequest],
        Union[storage_control.ManagedFolder, Awaitable[storage_control.ManagedFolder]],
    ]:
        raise NotImplementedError()

    @property
    def list_managed_folders(
        self,
    ) -> Callable[
        [storage_control.ListManagedFoldersRequest],
        Union[
            storage_control.ListManagedFoldersResponse,
            Awaitable[storage_control.ListManagedFoldersResponse],
        ],
    ]:
        raise NotImplementedError()

    @property
    def create_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.CreateAnywhereCacheRequest],
        Union[operations_pb2.Operation, Awaitable[operations_pb2.Operation]],
    ]:
        raise NotImplementedError()

    @property
    def update_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.UpdateAnywhereCacheRequest],
        Union[operations_pb2.Operation, Awaitable[operations_pb2.Operation]],
    ]:
        raise NotImplementedError()

    @property
    def disable_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.DisableAnywhereCacheRequest],
        Union[storage_control.AnywhereCache, Awaitable[storage_control.AnywhereCache]],
    ]:
        raise NotImplementedError()

    @property
    def pause_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.PauseAnywhereCacheRequest],
        Union[storage_control.AnywhereCache, Awaitable[storage_control.AnywhereCache]],
    ]:
        raise NotImplementedError()

    @property
    def resume_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.ResumeAnywhereCacheRequest],
        Union[storage_control.AnywhereCache, Awaitable[storage_control.AnywhereCache]],
    ]:
        raise NotImplementedError()

    @property
    def get_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.GetAnywhereCacheRequest],
        Union[storage_control.AnywhereCache, Awaitable[storage_control.AnywhereCache]],
    ]:
        raise NotImplementedError()

    @property
    def list_anywhere_caches(
        self,
    ) -> Callable[
        [storage_control.ListAnywhereCachesRequest],
        Union[
            storage_control.ListAnywhereCachesResponse,
            Awaitable[storage_control.ListAnywhereCachesResponse],
        ],
    ]:
        raise NotImplementedError()

    @property
    def get_project_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetProjectIntelligenceConfigRequest],
        Union[
            storage_control.IntelligenceConfig,
            Awaitable[storage_control.IntelligenceConfig],
        ],
    ]:
        raise NotImplementedError()

    @property
    def update_project_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateProjectIntelligenceConfigRequest],
        Union[
            storage_control.IntelligenceConfig,
            Awaitable[storage_control.IntelligenceConfig],
        ],
    ]:
        raise NotImplementedError()

    @property
    def get_folder_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetFolderIntelligenceConfigRequest],
        Union[
            storage_control.IntelligenceConfig,
            Awaitable[storage_control.IntelligenceConfig],
        ],
    ]:
        raise NotImplementedError()

    @property
    def update_folder_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateFolderIntelligenceConfigRequest],
        Union[
            storage_control.IntelligenceConfig,
            Awaitable[storage_control.IntelligenceConfig],
        ],
    ]:
        raise NotImplementedError()

    @property
    def get_organization_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetOrganizationIntelligenceConfigRequest],
        Union[
            storage_control.IntelligenceConfig,
            Awaitable[storage_control.IntelligenceConfig],
        ],
    ]:
        raise NotImplementedError()

    @property
    def update_organization_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateOrganizationIntelligenceConfigRequest],
        Union[
            storage_control.IntelligenceConfig,
            Awaitable[storage_control.IntelligenceConfig],
        ],
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
    def kind(self) -> str:
        raise NotImplementedError()


__all__ = ("StorageControlTransport",)
