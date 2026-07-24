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
import dataclasses
import json  # type: ignore
import logging
import warnings
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import google.iam.v1.iam_policy_pb2 as iam_policy_pb2  # type: ignore
import google.iam.v1.policy_pb2 as policy_pb2  # type: ignore
import google.protobuf
import google.protobuf.empty_pb2 as empty_pb2  # type: ignore
from google.api_core import exceptions as core_exceptions
from google.api_core import gapic_v1, operations_v1, rest_helpers, rest_streaming
from google.api_core import retry as retries
from google.auth import credentials as ga_credentials  # type: ignore
from google.auth.transport.requests import AuthorizedSession  # type: ignore
from google.longrunning import operations_pb2  # type: ignore
from google.protobuf import json_format
from requests import __version__ as requests_version

from google.cloud.storage_control_v2.types import storage_control

from .base import DEFAULT_CLIENT_INFO as BASE_DEFAULT_CLIENT_INFO
from .rest_base import _BaseStorageControlRestTransport

try:
    OptionalRetry = Union[retries.Retry, gapic_v1.method._MethodDefault, None]
except AttributeError:  # pragma: NO COVER
    OptionalRetry = Union[retries.Retry, object, None]  # type: ignore

try:
    from google.api_core import client_logging  # type: ignore

    CLIENT_LOGGING_SUPPORTED = True  # pragma: NO COVER
except ImportError:  # pragma: NO COVER
    CLIENT_LOGGING_SUPPORTED = False

_LOGGER = logging.getLogger(__name__)

DEFAULT_CLIENT_INFO = gapic_v1.client_info.ClientInfo(
    gapic_version=BASE_DEFAULT_CLIENT_INFO.gapic_version,
    grpc_version=None,
    rest_version=f"requests@{requests_version}",
)

if hasattr(DEFAULT_CLIENT_INFO, "protobuf_runtime_version"):  # pragma: NO COVER
    DEFAULT_CLIENT_INFO.protobuf_runtime_version = google.protobuf.__version__


class StorageControlRestInterceptor:
    """Interceptor for StorageControl.

    Interceptors are used to manipulate requests, request metadata, and responses
    in arbitrary ways.
    Example use cases include:
    * Logging
    * Verifying requests according to service or custom semantics
    * Stripping extraneous information from responses

    These use cases and more can be enabled by injecting an
    instance of a custom subclass when constructing the StorageControlRestTransport.

    .. code-block:: python
        class MyCustomStorageControlInterceptor(StorageControlRestInterceptor):
            def pre_create_anywhere_cache(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_create_anywhere_cache(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_create_folder(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_create_folder(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_create_managed_folder(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_create_managed_folder(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_delete_folder(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def pre_delete_folder_recursive(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_delete_folder_recursive(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_delete_managed_folder(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def pre_disable_anywhere_cache(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_disable_anywhere_cache(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_anywhere_cache(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_anywhere_cache(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_folder(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_folder(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_folder_intelligence_config(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_folder_intelligence_config(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_iam_policy(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_iam_policy(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_intelligence_finding(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_intelligence_finding(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_intelligence_finding_revision(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_intelligence_finding_revision(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_managed_folder(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_managed_folder(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_organization_intelligence_config(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_organization_intelligence_config(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_project_intelligence_config(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_project_intelligence_config(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_get_storage_layout(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_get_storage_layout(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_list_anywhere_caches(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_list_anywhere_caches(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_list_folders(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_list_folders(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_list_intelligence_finding_revisions(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_list_intelligence_finding_revisions(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_list_intelligence_findings(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_list_intelligence_findings(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_list_managed_folders(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_list_managed_folders(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_pause_anywhere_cache(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_pause_anywhere_cache(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_rename_folder(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_rename_folder(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_resume_anywhere_cache(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_resume_anywhere_cache(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_set_iam_policy(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_set_iam_policy(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_summarize_intelligence_findings(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_summarize_intelligence_findings(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_test_iam_permissions(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_test_iam_permissions(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_update_anywhere_cache(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_update_anywhere_cache(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_update_folder_intelligence_config(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_update_folder_intelligence_config(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_update_organization_intelligence_config(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_update_organization_intelligence_config(self, response):
                logging.log(f"Received response: {response}")
                return response

            def pre_update_project_intelligence_config(self, request, metadata):
                logging.log(f"Received request: {request}")
                return request, metadata

            def post_update_project_intelligence_config(self, response):
                logging.log(f"Received response: {response}")
                return response

        transport = StorageControlRestTransport(interceptor=MyCustomStorageControlInterceptor())
        client = StorageControlClient(transport=transport)


    """

    def pre_get_folder_intelligence_config(
        self,
        request: storage_control.GetFolderIntelligenceConfigRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.GetFolderIntelligenceConfigRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for get_folder_intelligence_config

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_get_folder_intelligence_config(
        self, response: storage_control.IntelligenceConfig
    ) -> storage_control.IntelligenceConfig:
        """Post-rpc interceptor for get_folder_intelligence_config

        DEPRECATED. Please use the `post_get_folder_intelligence_config_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_get_folder_intelligence_config` interceptor runs
        before the `post_get_folder_intelligence_config_with_metadata` interceptor.
        """
        return response

    def post_get_folder_intelligence_config_with_metadata(
        self,
        response: storage_control.IntelligenceConfig,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.IntelligenceConfig, Sequence[Tuple[str, Union[str, bytes]]]
    ]:
        """Post-rpc interceptor for get_folder_intelligence_config

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_get_folder_intelligence_config_with_metadata`
        interceptor in new development instead of the `post_get_folder_intelligence_config` interceptor.
        When both interceptors are used, this `post_get_folder_intelligence_config_with_metadata` interceptor runs after the
        `post_get_folder_intelligence_config` interceptor. The (possibly modified) response returned by
        `post_get_folder_intelligence_config` will be passed to
        `post_get_folder_intelligence_config_with_metadata`.
        """
        return response, metadata

    def pre_get_intelligence_finding(
        self,
        request: storage_control.GetIntelligenceFindingRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.GetIntelligenceFindingRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for get_intelligence_finding

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_get_intelligence_finding(
        self, response: storage_control.IntelligenceFinding
    ) -> storage_control.IntelligenceFinding:
        """Post-rpc interceptor for get_intelligence_finding

        DEPRECATED. Please use the `post_get_intelligence_finding_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_get_intelligence_finding` interceptor runs
        before the `post_get_intelligence_finding_with_metadata` interceptor.
        """
        return response

    def post_get_intelligence_finding_with_metadata(
        self,
        response: storage_control.IntelligenceFinding,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.IntelligenceFinding, Sequence[Tuple[str, Union[str, bytes]]]
    ]:
        """Post-rpc interceptor for get_intelligence_finding

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_get_intelligence_finding_with_metadata`
        interceptor in new development instead of the `post_get_intelligence_finding` interceptor.
        When both interceptors are used, this `post_get_intelligence_finding_with_metadata` interceptor runs after the
        `post_get_intelligence_finding` interceptor. The (possibly modified) response returned by
        `post_get_intelligence_finding` will be passed to
        `post_get_intelligence_finding_with_metadata`.
        """
        return response, metadata

    def pre_get_intelligence_finding_revision(
        self,
        request: storage_control.GetIntelligenceFindingRevisionRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.GetIntelligenceFindingRevisionRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for get_intelligence_finding_revision

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_get_intelligence_finding_revision(
        self, response: storage_control.IntelligenceFindingRevision
    ) -> storage_control.IntelligenceFindingRevision:
        """Post-rpc interceptor for get_intelligence_finding_revision

        DEPRECATED. Please use the `post_get_intelligence_finding_revision_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_get_intelligence_finding_revision` interceptor runs
        before the `post_get_intelligence_finding_revision_with_metadata` interceptor.
        """
        return response

    def post_get_intelligence_finding_revision_with_metadata(
        self,
        response: storage_control.IntelligenceFindingRevision,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.IntelligenceFindingRevision,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Post-rpc interceptor for get_intelligence_finding_revision

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_get_intelligence_finding_revision_with_metadata`
        interceptor in new development instead of the `post_get_intelligence_finding_revision` interceptor.
        When both interceptors are used, this `post_get_intelligence_finding_revision_with_metadata` interceptor runs after the
        `post_get_intelligence_finding_revision` interceptor. The (possibly modified) response returned by
        `post_get_intelligence_finding_revision` will be passed to
        `post_get_intelligence_finding_revision_with_metadata`.
        """
        return response, metadata

    def pre_get_organization_intelligence_config(
        self,
        request: storage_control.GetOrganizationIntelligenceConfigRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.GetOrganizationIntelligenceConfigRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for get_organization_intelligence_config

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_get_organization_intelligence_config(
        self, response: storage_control.IntelligenceConfig
    ) -> storage_control.IntelligenceConfig:
        """Post-rpc interceptor for get_organization_intelligence_config

        DEPRECATED. Please use the `post_get_organization_intelligence_config_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_get_organization_intelligence_config` interceptor runs
        before the `post_get_organization_intelligence_config_with_metadata` interceptor.
        """
        return response

    def post_get_organization_intelligence_config_with_metadata(
        self,
        response: storage_control.IntelligenceConfig,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.IntelligenceConfig, Sequence[Tuple[str, Union[str, bytes]]]
    ]:
        """Post-rpc interceptor for get_organization_intelligence_config

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_get_organization_intelligence_config_with_metadata`
        interceptor in new development instead of the `post_get_organization_intelligence_config` interceptor.
        When both interceptors are used, this `post_get_organization_intelligence_config_with_metadata` interceptor runs after the
        `post_get_organization_intelligence_config` interceptor. The (possibly modified) response returned by
        `post_get_organization_intelligence_config` will be passed to
        `post_get_organization_intelligence_config_with_metadata`.
        """
        return response, metadata

    def pre_get_project_intelligence_config(
        self,
        request: storage_control.GetProjectIntelligenceConfigRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.GetProjectIntelligenceConfigRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for get_project_intelligence_config

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_get_project_intelligence_config(
        self, response: storage_control.IntelligenceConfig
    ) -> storage_control.IntelligenceConfig:
        """Post-rpc interceptor for get_project_intelligence_config

        DEPRECATED. Please use the `post_get_project_intelligence_config_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_get_project_intelligence_config` interceptor runs
        before the `post_get_project_intelligence_config_with_metadata` interceptor.
        """
        return response

    def post_get_project_intelligence_config_with_metadata(
        self,
        response: storage_control.IntelligenceConfig,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.IntelligenceConfig, Sequence[Tuple[str, Union[str, bytes]]]
    ]:
        """Post-rpc interceptor for get_project_intelligence_config

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_get_project_intelligence_config_with_metadata`
        interceptor in new development instead of the `post_get_project_intelligence_config` interceptor.
        When both interceptors are used, this `post_get_project_intelligence_config_with_metadata` interceptor runs after the
        `post_get_project_intelligence_config` interceptor. The (possibly modified) response returned by
        `post_get_project_intelligence_config` will be passed to
        `post_get_project_intelligence_config_with_metadata`.
        """
        return response, metadata

    def pre_list_intelligence_finding_revisions(
        self,
        request: storage_control.ListIntelligenceFindingRevisionsRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.ListIntelligenceFindingRevisionsRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for list_intelligence_finding_revisions

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_list_intelligence_finding_revisions(
        self, response: storage_control.ListIntelligenceFindingRevisionsResponse
    ) -> storage_control.ListIntelligenceFindingRevisionsResponse:
        """Post-rpc interceptor for list_intelligence_finding_revisions

        DEPRECATED. Please use the `post_list_intelligence_finding_revisions_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_list_intelligence_finding_revisions` interceptor runs
        before the `post_list_intelligence_finding_revisions_with_metadata` interceptor.
        """
        return response

    def post_list_intelligence_finding_revisions_with_metadata(
        self,
        response: storage_control.ListIntelligenceFindingRevisionsResponse,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.ListIntelligenceFindingRevisionsResponse,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Post-rpc interceptor for list_intelligence_finding_revisions

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_list_intelligence_finding_revisions_with_metadata`
        interceptor in new development instead of the `post_list_intelligence_finding_revisions` interceptor.
        When both interceptors are used, this `post_list_intelligence_finding_revisions_with_metadata` interceptor runs after the
        `post_list_intelligence_finding_revisions` interceptor. The (possibly modified) response returned by
        `post_list_intelligence_finding_revisions` will be passed to
        `post_list_intelligence_finding_revisions_with_metadata`.
        """
        return response, metadata

    def pre_list_intelligence_findings(
        self,
        request: storage_control.ListIntelligenceFindingsRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.ListIntelligenceFindingsRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for list_intelligence_findings

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_list_intelligence_findings(
        self, response: storage_control.ListIntelligenceFindingsResponse
    ) -> storage_control.ListIntelligenceFindingsResponse:
        """Post-rpc interceptor for list_intelligence_findings

        DEPRECATED. Please use the `post_list_intelligence_findings_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_list_intelligence_findings` interceptor runs
        before the `post_list_intelligence_findings_with_metadata` interceptor.
        """
        return response

    def post_list_intelligence_findings_with_metadata(
        self,
        response: storage_control.ListIntelligenceFindingsResponse,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.ListIntelligenceFindingsResponse,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Post-rpc interceptor for list_intelligence_findings

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_list_intelligence_findings_with_metadata`
        interceptor in new development instead of the `post_list_intelligence_findings` interceptor.
        When both interceptors are used, this `post_list_intelligence_findings_with_metadata` interceptor runs after the
        `post_list_intelligence_findings` interceptor. The (possibly modified) response returned by
        `post_list_intelligence_findings` will be passed to
        `post_list_intelligence_findings_with_metadata`.
        """
        return response, metadata

    def pre_summarize_intelligence_findings(
        self,
        request: storage_control.SummarizeIntelligenceFindingsRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.SummarizeIntelligenceFindingsRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for summarize_intelligence_findings

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_summarize_intelligence_findings(
        self, response: storage_control.SummarizeIntelligenceFindingsResponse
    ) -> storage_control.SummarizeIntelligenceFindingsResponse:
        """Post-rpc interceptor for summarize_intelligence_findings

        DEPRECATED. Please use the `post_summarize_intelligence_findings_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_summarize_intelligence_findings` interceptor runs
        before the `post_summarize_intelligence_findings_with_metadata` interceptor.
        """
        return response

    def post_summarize_intelligence_findings_with_metadata(
        self,
        response: storage_control.SummarizeIntelligenceFindingsResponse,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.SummarizeIntelligenceFindingsResponse,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Post-rpc interceptor for summarize_intelligence_findings

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_summarize_intelligence_findings_with_metadata`
        interceptor in new development instead of the `post_summarize_intelligence_findings` interceptor.
        When both interceptors are used, this `post_summarize_intelligence_findings_with_metadata` interceptor runs after the
        `post_summarize_intelligence_findings` interceptor. The (possibly modified) response returned by
        `post_summarize_intelligence_findings` will be passed to
        `post_summarize_intelligence_findings_with_metadata`.
        """
        return response, metadata

    def pre_update_folder_intelligence_config(
        self,
        request: storage_control.UpdateFolderIntelligenceConfigRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.UpdateFolderIntelligenceConfigRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for update_folder_intelligence_config

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_update_folder_intelligence_config(
        self, response: storage_control.IntelligenceConfig
    ) -> storage_control.IntelligenceConfig:
        """Post-rpc interceptor for update_folder_intelligence_config

        DEPRECATED. Please use the `post_update_folder_intelligence_config_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_update_folder_intelligence_config` interceptor runs
        before the `post_update_folder_intelligence_config_with_metadata` interceptor.
        """
        return response

    def post_update_folder_intelligence_config_with_metadata(
        self,
        response: storage_control.IntelligenceConfig,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.IntelligenceConfig, Sequence[Tuple[str, Union[str, bytes]]]
    ]:
        """Post-rpc interceptor for update_folder_intelligence_config

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_update_folder_intelligence_config_with_metadata`
        interceptor in new development instead of the `post_update_folder_intelligence_config` interceptor.
        When both interceptors are used, this `post_update_folder_intelligence_config_with_metadata` interceptor runs after the
        `post_update_folder_intelligence_config` interceptor. The (possibly modified) response returned by
        `post_update_folder_intelligence_config` will be passed to
        `post_update_folder_intelligence_config_with_metadata`.
        """
        return response, metadata

    def pre_update_organization_intelligence_config(
        self,
        request: storage_control.UpdateOrganizationIntelligenceConfigRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.UpdateOrganizationIntelligenceConfigRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for update_organization_intelligence_config

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_update_organization_intelligence_config(
        self, response: storage_control.IntelligenceConfig
    ) -> storage_control.IntelligenceConfig:
        """Post-rpc interceptor for update_organization_intelligence_config

        DEPRECATED. Please use the `post_update_organization_intelligence_config_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_update_organization_intelligence_config` interceptor runs
        before the `post_update_organization_intelligence_config_with_metadata` interceptor.
        """
        return response

    def post_update_organization_intelligence_config_with_metadata(
        self,
        response: storage_control.IntelligenceConfig,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.IntelligenceConfig, Sequence[Tuple[str, Union[str, bytes]]]
    ]:
        """Post-rpc interceptor for update_organization_intelligence_config

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_update_organization_intelligence_config_with_metadata`
        interceptor in new development instead of the `post_update_organization_intelligence_config` interceptor.
        When both interceptors are used, this `post_update_organization_intelligence_config_with_metadata` interceptor runs after the
        `post_update_organization_intelligence_config` interceptor. The (possibly modified) response returned by
        `post_update_organization_intelligence_config` will be passed to
        `post_update_organization_intelligence_config_with_metadata`.
        """
        return response, metadata

    def pre_update_project_intelligence_config(
        self,
        request: storage_control.UpdateProjectIntelligenceConfigRequest,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.UpdateProjectIntelligenceConfigRequest,
        Sequence[Tuple[str, Union[str, bytes]]],
    ]:
        """Pre-rpc interceptor for update_project_intelligence_config

        Override in a subclass to manipulate the request or metadata
        before they are sent to the StorageControl server.
        """
        return request, metadata

    def post_update_project_intelligence_config(
        self, response: storage_control.IntelligenceConfig
    ) -> storage_control.IntelligenceConfig:
        """Post-rpc interceptor for update_project_intelligence_config

        DEPRECATED. Please use the `post_update_project_intelligence_config_with_metadata`
        interceptor instead.

        Override in a subclass to read or manipulate the response
        after it is returned by the StorageControl server but before
        it is returned to user code. This `post_update_project_intelligence_config` interceptor runs
        before the `post_update_project_intelligence_config_with_metadata` interceptor.
        """
        return response

    def post_update_project_intelligence_config_with_metadata(
        self,
        response: storage_control.IntelligenceConfig,
        metadata: Sequence[Tuple[str, Union[str, bytes]]],
    ) -> Tuple[
        storage_control.IntelligenceConfig, Sequence[Tuple[str, Union[str, bytes]]]
    ]:
        """Post-rpc interceptor for update_project_intelligence_config

        Override in a subclass to read or manipulate the response or metadata after it
        is returned by the StorageControl server but before it is returned to user code.

        We recommend only using this `post_update_project_intelligence_config_with_metadata`
        interceptor in new development instead of the `post_update_project_intelligence_config` interceptor.
        When both interceptors are used, this `post_update_project_intelligence_config_with_metadata` interceptor runs after the
        `post_update_project_intelligence_config` interceptor. The (possibly modified) response returned by
        `post_update_project_intelligence_config` will be passed to
        `post_update_project_intelligence_config_with_metadata`.
        """
        return response, metadata


@dataclasses.dataclass
class StorageControlRestStub:
    _session: AuthorizedSession
    _host: str
    _interceptor: StorageControlRestInterceptor


class StorageControlRestTransport(_BaseStorageControlRestTransport):
    """REST backend synchronous transport for StorageControl.

    StorageControl service includes selected control plane
    operations.

    This class defines the same methods as the primary client, so the
    primary client can load the underlying transport implementation
    and call it.

    It sends JSON representations of protocol buffers over HTTP/1.1
    """

    def __init__(
        self,
        *,
        host: str = "storage.googleapis.com",
        credentials: Optional[ga_credentials.Credentials] = None,
        credentials_file: Optional[str] = None,
        scopes: Optional[Sequence[str]] = None,
        client_cert_source_for_mtls: Optional[Callable[[], Tuple[bytes, bytes]]] = None,
        quota_project_id: Optional[str] = None,
        client_info: gapic_v1.client_info.ClientInfo = DEFAULT_CLIENT_INFO,
        always_use_jwt_access: Optional[bool] = False,
        url_scheme: str = "https",
        interceptor: Optional[StorageControlRestInterceptor] = None,
        api_audience: Optional[str] = None,
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
                This argument is ignored if ``channel`` is provided. This argument will be
                removed in the next major version of this library.
            scopes (Optional(Sequence[str])): A list of scopes. This argument is
                ignored if ``channel`` is provided.
            client_cert_source_for_mtls (Callable[[], Tuple[bytes, bytes]]): Client
                certificate to configure mutual TLS HTTP channel. It is ignored
                if ``channel`` is provided.
            quota_project_id (Optional[str]): An optional project to use for billing
                and quota.
            client_info (google.api_core.gapic_v1.client_info.ClientInfo):
                The client info used to send a user-agent string along with
                API requests. If ``None``, then default info will be used.
                Generally, you only need to set this if you are developing
                your own client library.
            always_use_jwt_access (Optional[bool]): Whether self signed JWT should
                be used for service account credentials.
            url_scheme: the protocol scheme for the API endpoint.  Normally
                "https", but for testing or local servers,
                "http" can be specified.
            interceptor (Optional[StorageControlRestInterceptor]): Interceptor used
                to manipulate requests, request metadata, and responses.
            api_audience (Optional[str]): The intended audience for the API calls
                to the service that will be set when using certain 3rd party
                authentication flows. Audience is typically a resource identifier.
                If not set, the host value will be used as a default.
        """
        # Run the base constructor
        # TODO(yon-mg): resolve other ctor params i.e. scopes, quota, etc.
        # TODO: When custom host (api_endpoint) is set, `scopes` must *also* be set on the
        # credentials object
        super().__init__(
            host=host,
            credentials=credentials,
            client_info=client_info,
            always_use_jwt_access=always_use_jwt_access,
            url_scheme=url_scheme,
            api_audience=api_audience,
        )
        self._session = AuthorizedSession(
            self._credentials, default_host=self.DEFAULT_HOST
        )
        self._operations_client: Optional[operations_v1.AbstractOperationsClient] = None
        if client_cert_source_for_mtls:
            self._session.configure_mtls_channel(client_cert_source_for_mtls)
        self._interceptor = interceptor or StorageControlRestInterceptor()
        self._prep_wrapped_messages(client_info)

    @property
    def operations_client(self) -> operations_v1.AbstractOperationsClient:
        """Create the client designed to process long-running operations.

        This property caches on the instance; repeated calls return the same
        client.
        """
        # Only create a new client if we do not already have one.
        if self._operations_client is None:
            http_options: Dict[str, List[Dict[str, str]]] = {}

            rest_transport = operations_v1.OperationsRestTransport(
                host=self._host,
                # use the credentials which are saved
                credentials=self._credentials,
                scopes=self._scopes,
                http_options=http_options,
                path_prefix="v2",
            )

            self._operations_client = operations_v1.AbstractOperationsClient(
                transport=rest_transport
            )

        # Return the client from cache.
        return self._operations_client

    class _CreateAnywhereCache(
        _BaseStorageControlRestTransport._BaseCreateAnywhereCache,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.CreateAnywhereCache")

        def __call__(
            self,
            request: storage_control.CreateAnywhereCacheRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> operations_pb2.Operation:
            raise NotImplementedError(
                "Method CreateAnywhereCache is not available over REST transport"
            )

    class _CreateFolder(
        _BaseStorageControlRestTransport._BaseCreateFolder, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.CreateFolder")

        def __call__(
            self,
            request: storage_control.CreateFolderRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.Folder:
            raise NotImplementedError(
                "Method CreateFolder is not available over REST transport"
            )

    class _CreateManagedFolder(
        _BaseStorageControlRestTransport._BaseCreateManagedFolder,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.CreateManagedFolder")

        def __call__(
            self,
            request: storage_control.CreateManagedFolderRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.ManagedFolder:
            raise NotImplementedError(
                "Method CreateManagedFolder is not available over REST transport"
            )

    class _DeleteFolder(
        _BaseStorageControlRestTransport._BaseDeleteFolder, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.DeleteFolder")

        def __call__(
            self,
            request: storage_control.DeleteFolderRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ):
            raise NotImplementedError(
                "Method DeleteFolder is not available over REST transport"
            )

    class _DeleteFolderRecursive(
        _BaseStorageControlRestTransport._BaseDeleteFolderRecursive,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.DeleteFolderRecursive")

        def __call__(
            self,
            request: storage_control.DeleteFolderRecursiveRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> operations_pb2.Operation:
            raise NotImplementedError(
                "Method DeleteFolderRecursive is not available over REST transport"
            )

    class _DeleteManagedFolder(
        _BaseStorageControlRestTransport._BaseDeleteManagedFolder,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.DeleteManagedFolder")

        def __call__(
            self,
            request: storage_control.DeleteManagedFolderRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ):
            raise NotImplementedError(
                "Method DeleteManagedFolder is not available over REST transport"
            )

    class _DisableAnywhereCache(
        _BaseStorageControlRestTransport._BaseDisableAnywhereCache,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.DisableAnywhereCache")

        def __call__(
            self,
            request: storage_control.DisableAnywhereCacheRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.AnywhereCache:
            raise NotImplementedError(
                "Method DisableAnywhereCache is not available over REST transport"
            )

    class _GetAnywhereCache(
        _BaseStorageControlRestTransport._BaseGetAnywhereCache, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetAnywhereCache")

        def __call__(
            self,
            request: storage_control.GetAnywhereCacheRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.AnywhereCache:
            raise NotImplementedError(
                "Method GetAnywhereCache is not available over REST transport"
            )

    class _GetFolder(
        _BaseStorageControlRestTransport._BaseGetFolder, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetFolder")

        def __call__(
            self,
            request: storage_control.GetFolderRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.Folder:
            raise NotImplementedError(
                "Method GetFolder is not available over REST transport"
            )

    class _GetFolderIntelligenceConfig(
        _BaseStorageControlRestTransport._BaseGetFolderIntelligenceConfig,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetFolderIntelligenceConfig")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
            )
            return response

        def __call__(
            self,
            request: storage_control.GetFolderIntelligenceConfigRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.IntelligenceConfig:
            r"""Call the get folder intelligence
            config method over HTTP.

                Args:
                    request (~.storage_control.GetFolderIntelligenceConfigRequest):
                        The request object. Request message to get the ``IntelligenceConfig``
                    resource associated with your folder.

                    **IAM Permissions**

                    Requires ``storage.intelligenceConfigs.get``
                    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
                    permission on the folder.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.IntelligenceConfig:
                        The ``IntelligenceConfig`` resource associated with your
                    organization, folder, or project.

            """

            http_options = _BaseStorageControlRestTransport._BaseGetFolderIntelligenceConfig._get_http_options()

            request, metadata = self._interceptor.pre_get_folder_intelligence_config(
                request, metadata
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseGetFolderIntelligenceConfig._get_transcoded_request(
                http_options, request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseGetFolderIntelligenceConfig._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.GetFolderIntelligenceConfig",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetFolderIntelligenceConfig",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = (
                StorageControlRestTransport._GetFolderIntelligenceConfig._get_response(
                    self._host,
                    metadata,
                    query_params,
                    self._session,
                    timeout,
                    transcoded_request,
                )
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.IntelligenceConfig()
            pb_resp = storage_control.IntelligenceConfig.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_get_folder_intelligence_config(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_get_folder_intelligence_config_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = storage_control.IntelligenceConfig.to_json(
                        response
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.get_folder_intelligence_config",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetFolderIntelligenceConfig",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _GetIamPolicy(
        _BaseStorageControlRestTransport._BaseGetIamPolicy, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetIamPolicy")

        def __call__(
            self,
            request: iam_policy_pb2.GetIamPolicyRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> policy_pb2.Policy:
            raise NotImplementedError(
                "Method GetIamPolicy is not available over REST transport"
            )

    class _GetIntelligenceFinding(
        _BaseStorageControlRestTransport._BaseGetIntelligenceFinding,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetIntelligenceFinding")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
            )
            return response

        def __call__(
            self,
            request: storage_control.GetIntelligenceFindingRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.IntelligenceFinding:
            r"""Call the get intelligence finding method over HTTP.

            Args:
                request (~.storage_control.GetIntelligenceFindingRequest):
                    The request object. Request message to get the ``IntelligenceFinding``
                resource associated with a project.
                retry (google.api_core.retry.Retry): Designation of what errors, if any,
                    should be retried.
                timeout (float): The timeout for this request.
                metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                    sent along with the request as metadata. Normally, each value must be of type `str`,
                    but for metadata keys ending with the suffix `-bin`, the corresponding values must
                    be of type `bytes`.

            Returns:
                ~.storage_control.IntelligenceFinding:
                    The ``IntelligenceFinding`` resource that represents a
                security, performance, or cost-related finding about a
                project or bucket.

            """

            http_options = _BaseStorageControlRestTransport._BaseGetIntelligenceFinding._get_http_options()

            request, metadata = self._interceptor.pre_get_intelligence_finding(
                request, metadata
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseGetIntelligenceFinding._get_transcoded_request(
                http_options, request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseGetIntelligenceFinding._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.GetIntelligenceFinding",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetIntelligenceFinding",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = (
                StorageControlRestTransport._GetIntelligenceFinding._get_response(
                    self._host,
                    metadata,
                    query_params,
                    self._session,
                    timeout,
                    transcoded_request,
                )
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.IntelligenceFinding()
            pb_resp = storage_control.IntelligenceFinding.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_get_intelligence_finding(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = self._interceptor.post_get_intelligence_finding_with_metadata(
                resp, response_metadata
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = storage_control.IntelligenceFinding.to_json(
                        response
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.get_intelligence_finding",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetIntelligenceFinding",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _GetIntelligenceFindingRevision(
        _BaseStorageControlRestTransport._BaseGetIntelligenceFindingRevision,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetIntelligenceFindingRevision")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
            )
            return response

        def __call__(
            self,
            request: storage_control.GetIntelligenceFindingRevisionRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.IntelligenceFindingRevision:
            r"""Call the get intelligence finding
            revision method over HTTP.

                Args:
                    request (~.storage_control.GetIntelligenceFindingRevisionRequest):
                        The request object. Request message to get the
                    ``IntelligenceFindingRevision`` resource associated with
                    a project.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.IntelligenceFindingRevision:
                        An ``IntelligenceFindingRevision`` represents a specific
                    revision of an ``IntelligenceFinding`` resource.

            """

            http_options = _BaseStorageControlRestTransport._BaseGetIntelligenceFindingRevision._get_http_options()

            request, metadata = self._interceptor.pre_get_intelligence_finding_revision(
                request, metadata
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseGetIntelligenceFindingRevision._get_transcoded_request(
                http_options, request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseGetIntelligenceFindingRevision._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.GetIntelligenceFindingRevision",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetIntelligenceFindingRevision",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = StorageControlRestTransport._GetIntelligenceFindingRevision._get_response(
                self._host,
                metadata,
                query_params,
                self._session,
                timeout,
                transcoded_request,
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.IntelligenceFindingRevision()
            pb_resp = storage_control.IntelligenceFindingRevision.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_get_intelligence_finding_revision(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_get_intelligence_finding_revision_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = (
                        storage_control.IntelligenceFindingRevision.to_json(response)
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.get_intelligence_finding_revision",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetIntelligenceFindingRevision",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _GetManagedFolder(
        _BaseStorageControlRestTransport._BaseGetManagedFolder, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetManagedFolder")

        def __call__(
            self,
            request: storage_control.GetManagedFolderRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.ManagedFolder:
            raise NotImplementedError(
                "Method GetManagedFolder is not available over REST transport"
            )

    class _GetOrganizationIntelligenceConfig(
        _BaseStorageControlRestTransport._BaseGetOrganizationIntelligenceConfig,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetOrganizationIntelligenceConfig")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
            )
            return response

        def __call__(
            self,
            request: storage_control.GetOrganizationIntelligenceConfigRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.IntelligenceConfig:
            r"""Call the get organization
            intelligence config method over HTTP.

                Args:
                    request (~.storage_control.GetOrganizationIntelligenceConfigRequest):
                        The request object. Request message to get the ``IntelligenceConfig``
                    resource associated with your organization.

                    **IAM Permissions**

                    Requires ``storage.intelligenceConfigs.get``
                    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
                    permission on the organization.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.IntelligenceConfig:
                        The ``IntelligenceConfig`` resource associated with your
                    organization, folder, or project.

            """

            http_options = _BaseStorageControlRestTransport._BaseGetOrganizationIntelligenceConfig._get_http_options()

            request, metadata = (
                self._interceptor.pre_get_organization_intelligence_config(
                    request, metadata
                )
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseGetOrganizationIntelligenceConfig._get_transcoded_request(
                http_options, request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseGetOrganizationIntelligenceConfig._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.GetOrganizationIntelligenceConfig",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetOrganizationIntelligenceConfig",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = StorageControlRestTransport._GetOrganizationIntelligenceConfig._get_response(
                self._host,
                metadata,
                query_params,
                self._session,
                timeout,
                transcoded_request,
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.IntelligenceConfig()
            pb_resp = storage_control.IntelligenceConfig.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_get_organization_intelligence_config(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_get_organization_intelligence_config_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = storage_control.IntelligenceConfig.to_json(
                        response
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.get_organization_intelligence_config",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetOrganizationIntelligenceConfig",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _GetProjectIntelligenceConfig(
        _BaseStorageControlRestTransport._BaseGetProjectIntelligenceConfig,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetProjectIntelligenceConfig")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
            )
            return response

        def __call__(
            self,
            request: storage_control.GetProjectIntelligenceConfigRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.IntelligenceConfig:
            r"""Call the get project intelligence
            config method over HTTP.

                Args:
                    request (~.storage_control.GetProjectIntelligenceConfigRequest):
                        The request object. Request message to get the ``IntelligenceConfig``
                    resource associated with your project.

                    **IAM Permissions**:

                    Requires ``storage.intelligenceConfigs.get``
                    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
                    permission on the project.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.IntelligenceConfig:
                        The ``IntelligenceConfig`` resource associated with your
                    organization, folder, or project.

            """

            http_options = _BaseStorageControlRestTransport._BaseGetProjectIntelligenceConfig._get_http_options()

            request, metadata = self._interceptor.pre_get_project_intelligence_config(
                request, metadata
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseGetProjectIntelligenceConfig._get_transcoded_request(
                http_options, request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseGetProjectIntelligenceConfig._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.GetProjectIntelligenceConfig",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetProjectIntelligenceConfig",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = (
                StorageControlRestTransport._GetProjectIntelligenceConfig._get_response(
                    self._host,
                    metadata,
                    query_params,
                    self._session,
                    timeout,
                    transcoded_request,
                )
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.IntelligenceConfig()
            pb_resp = storage_control.IntelligenceConfig.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_get_project_intelligence_config(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_get_project_intelligence_config_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = storage_control.IntelligenceConfig.to_json(
                        response
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.get_project_intelligence_config",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "GetProjectIntelligenceConfig",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _GetStorageLayout(
        _BaseStorageControlRestTransport._BaseGetStorageLayout, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.GetStorageLayout")

        def __call__(
            self,
            request: storage_control.GetStorageLayoutRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.StorageLayout:
            raise NotImplementedError(
                "Method GetStorageLayout is not available over REST transport"
            )

    class _ListAnywhereCaches(
        _BaseStorageControlRestTransport._BaseListAnywhereCaches, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.ListAnywhereCaches")

        def __call__(
            self,
            request: storage_control.ListAnywhereCachesRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.ListAnywhereCachesResponse:
            raise NotImplementedError(
                "Method ListAnywhereCaches is not available over REST transport"
            )

    class _ListFolders(
        _BaseStorageControlRestTransport._BaseListFolders, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.ListFolders")

        def __call__(
            self,
            request: storage_control.ListFoldersRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.ListFoldersResponse:
            raise NotImplementedError(
                "Method ListFolders is not available over REST transport"
            )

    class _ListIntelligenceFindingRevisions(
        _BaseStorageControlRestTransport._BaseListIntelligenceFindingRevisions,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.ListIntelligenceFindingRevisions")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
            )
            return response

        def __call__(
            self,
            request: storage_control.ListIntelligenceFindingRevisionsRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.ListIntelligenceFindingRevisionsResponse:
            r"""Call the list intelligence finding
            revisions method over HTTP.

                Args:
                    request (~.storage_control.ListIntelligenceFindingRevisionsRequest):
                        The request object. Request message to list ``IntelligenceFindingRevision``
                    resources associated with a project.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.ListIntelligenceFindingRevisionsResponse:
                        Response message to list ``IntelligenceFindingRevision``
                    resources associated with a project.

            """

            http_options = _BaseStorageControlRestTransport._BaseListIntelligenceFindingRevisions._get_http_options()

            request, metadata = (
                self._interceptor.pre_list_intelligence_finding_revisions(
                    request, metadata
                )
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseListIntelligenceFindingRevisions._get_transcoded_request(
                http_options, request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseListIntelligenceFindingRevisions._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.ListIntelligenceFindingRevisions",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "ListIntelligenceFindingRevisions",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = StorageControlRestTransport._ListIntelligenceFindingRevisions._get_response(
                self._host,
                metadata,
                query_params,
                self._session,
                timeout,
                transcoded_request,
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.ListIntelligenceFindingRevisionsResponse()
            pb_resp = storage_control.ListIntelligenceFindingRevisionsResponse.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_list_intelligence_finding_revisions(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_list_intelligence_finding_revisions_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = storage_control.ListIntelligenceFindingRevisionsResponse.to_json(
                        response
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.list_intelligence_finding_revisions",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "ListIntelligenceFindingRevisions",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _ListIntelligenceFindings(
        _BaseStorageControlRestTransport._BaseListIntelligenceFindings,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.ListIntelligenceFindings")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
            )
            return response

        def __call__(
            self,
            request: storage_control.ListIntelligenceFindingsRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.ListIntelligenceFindingsResponse:
            r"""Call the list intelligence
            findings method over HTTP.

                Args:
                    request (~.storage_control.ListIntelligenceFindingsRequest):
                        The request object. Request message to list ``IntelligenceFinding``
                    resources associated with a project.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.ListIntelligenceFindingsResponse:
                        Response message to list the ``IntelligenceFinding``
                    resources associated with a project.

            """

            http_options = _BaseStorageControlRestTransport._BaseListIntelligenceFindings._get_http_options()

            request, metadata = self._interceptor.pre_list_intelligence_findings(
                request, metadata
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseListIntelligenceFindings._get_transcoded_request(
                http_options, request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseListIntelligenceFindings._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.ListIntelligenceFindings",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "ListIntelligenceFindings",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = (
                StorageControlRestTransport._ListIntelligenceFindings._get_response(
                    self._host,
                    metadata,
                    query_params,
                    self._session,
                    timeout,
                    transcoded_request,
                )
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.ListIntelligenceFindingsResponse()
            pb_resp = storage_control.ListIntelligenceFindingsResponse.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_list_intelligence_findings(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = self._interceptor.post_list_intelligence_findings_with_metadata(
                resp, response_metadata
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = (
                        storage_control.ListIntelligenceFindingsResponse.to_json(
                            response
                        )
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.list_intelligence_findings",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "ListIntelligenceFindings",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _ListManagedFolders(
        _BaseStorageControlRestTransport._BaseListManagedFolders, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.ListManagedFolders")

        def __call__(
            self,
            request: storage_control.ListManagedFoldersRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.ListManagedFoldersResponse:
            raise NotImplementedError(
                "Method ListManagedFolders is not available over REST transport"
            )

    class _PauseAnywhereCache(
        _BaseStorageControlRestTransport._BasePauseAnywhereCache, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.PauseAnywhereCache")

        def __call__(
            self,
            request: storage_control.PauseAnywhereCacheRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.AnywhereCache:
            raise NotImplementedError(
                "Method PauseAnywhereCache is not available over REST transport"
            )

    class _RenameFolder(
        _BaseStorageControlRestTransport._BaseRenameFolder, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.RenameFolder")

        def __call__(
            self,
            request: storage_control.RenameFolderRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> operations_pb2.Operation:
            raise NotImplementedError(
                "Method RenameFolder is not available over REST transport"
            )

    class _ResumeAnywhereCache(
        _BaseStorageControlRestTransport._BaseResumeAnywhereCache,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.ResumeAnywhereCache")

        def __call__(
            self,
            request: storage_control.ResumeAnywhereCacheRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.AnywhereCache:
            raise NotImplementedError(
                "Method ResumeAnywhereCache is not available over REST transport"
            )

    class _SetIamPolicy(
        _BaseStorageControlRestTransport._BaseSetIamPolicy, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.SetIamPolicy")

        def __call__(
            self,
            request: iam_policy_pb2.SetIamPolicyRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> policy_pb2.Policy:
            raise NotImplementedError(
                "Method SetIamPolicy is not available over REST transport"
            )

    class _SummarizeIntelligenceFindings(
        _BaseStorageControlRestTransport._BaseSummarizeIntelligenceFindings,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.SummarizeIntelligenceFindings")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
            )
            return response

        def __call__(
            self,
            request: storage_control.SummarizeIntelligenceFindingsRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.SummarizeIntelligenceFindingsResponse:
            r"""Call the summarize intelligence
            findings method over HTTP.

                Args:
                    request (~.storage_control.SummarizeIntelligenceFindingsRequest):
                        The request object. Request message to summarize the
                    intelligence findings for the specified
                    scope(org, folder or project).
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.SummarizeIntelligenceFindingsResponse:
                        Response message to summarize the
                    intelligence findings for a specified
                    scope(org, folder or project).

            """

            http_options = _BaseStorageControlRestTransport._BaseSummarizeIntelligenceFindings._get_http_options()

            request, metadata = self._interceptor.pre_summarize_intelligence_findings(
                request, metadata
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseSummarizeIntelligenceFindings._get_transcoded_request(
                http_options, request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseSummarizeIntelligenceFindings._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.SummarizeIntelligenceFindings",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "SummarizeIntelligenceFindings",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = StorageControlRestTransport._SummarizeIntelligenceFindings._get_response(
                self._host,
                metadata,
                query_params,
                self._session,
                timeout,
                transcoded_request,
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.SummarizeIntelligenceFindingsResponse()
            pb_resp = storage_control.SummarizeIntelligenceFindingsResponse.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_summarize_intelligence_findings(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_summarize_intelligence_findings_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = (
                        storage_control.SummarizeIntelligenceFindingsResponse.to_json(
                            response
                        )
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.summarize_intelligence_findings",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "SummarizeIntelligenceFindings",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _TestIamPermissions(
        _BaseStorageControlRestTransport._BaseTestIamPermissions, StorageControlRestStub
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.TestIamPermissions")

        def __call__(
            self,
            request: iam_policy_pb2.TestIamPermissionsRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> iam_policy_pb2.TestIamPermissionsResponse:
            raise NotImplementedError(
                "Method TestIamPermissions is not available over REST transport"
            )

    class _UpdateAnywhereCache(
        _BaseStorageControlRestTransport._BaseUpdateAnywhereCache,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.UpdateAnywhereCache")

        def __call__(
            self,
            request: storage_control.UpdateAnywhereCacheRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> operations_pb2.Operation:
            raise NotImplementedError(
                "Method UpdateAnywhereCache is not available over REST transport"
            )

    class _UpdateFolderIntelligenceConfig(
        _BaseStorageControlRestTransport._BaseUpdateFolderIntelligenceConfig,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.UpdateFolderIntelligenceConfig")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
                data=body,
            )
            return response

        def __call__(
            self,
            request: storage_control.UpdateFolderIntelligenceConfigRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.IntelligenceConfig:
            r"""Call the update folder
            intelligence config method over HTTP.

                Args:
                    request (~.storage_control.UpdateFolderIntelligenceConfigRequest):
                        The request object. Request message to update the ``IntelligenceConfig``
                    resource associated with your folder.

                    **IAM Permissions**:

                    Requires ``storage.intelligenceConfigs.update``
                    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
                    permission on the folder.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.IntelligenceConfig:
                        The ``IntelligenceConfig`` resource associated with your
                    organization, folder, or project.

            """

            http_options = _BaseStorageControlRestTransport._BaseUpdateFolderIntelligenceConfig._get_http_options()

            request, metadata = self._interceptor.pre_update_folder_intelligence_config(
                request, metadata
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseUpdateFolderIntelligenceConfig._get_transcoded_request(
                http_options, request
            )

            body = _BaseStorageControlRestTransport._BaseUpdateFolderIntelligenceConfig._get_request_body_json(
                transcoded_request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseUpdateFolderIntelligenceConfig._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.UpdateFolderIntelligenceConfig",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "UpdateFolderIntelligenceConfig",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = StorageControlRestTransport._UpdateFolderIntelligenceConfig._get_response(
                self._host,
                metadata,
                query_params,
                self._session,
                timeout,
                transcoded_request,
                body,
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.IntelligenceConfig()
            pb_resp = storage_control.IntelligenceConfig.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_update_folder_intelligence_config(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_update_folder_intelligence_config_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = storage_control.IntelligenceConfig.to_json(
                        response
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.update_folder_intelligence_config",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "UpdateFolderIntelligenceConfig",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _UpdateOrganizationIntelligenceConfig(
        _BaseStorageControlRestTransport._BaseUpdateOrganizationIntelligenceConfig,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash(
                "StorageControlRestTransport.UpdateOrganizationIntelligenceConfig"
            )

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
                data=body,
            )
            return response

        def __call__(
            self,
            request: storage_control.UpdateOrganizationIntelligenceConfigRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.IntelligenceConfig:
            r"""Call the update organization
            intelligence config method over HTTP.

                Args:
                    request (~.storage_control.UpdateOrganizationIntelligenceConfigRequest):
                        The request object. Request message to update the ``IntelligenceConfig``
                    resource associated with your organization.

                    **IAM Permissions**:

                    Requires ``storage.intelligenceConfigs.update``
                    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
                    permission on the organization.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.IntelligenceConfig:
                        The ``IntelligenceConfig`` resource associated with your
                    organization, folder, or project.

            """

            http_options = _BaseStorageControlRestTransport._BaseUpdateOrganizationIntelligenceConfig._get_http_options()

            request, metadata = (
                self._interceptor.pre_update_organization_intelligence_config(
                    request, metadata
                )
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseUpdateOrganizationIntelligenceConfig._get_transcoded_request(
                http_options, request
            )

            body = _BaseStorageControlRestTransport._BaseUpdateOrganizationIntelligenceConfig._get_request_body_json(
                transcoded_request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseUpdateOrganizationIntelligenceConfig._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.UpdateOrganizationIntelligenceConfig",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "UpdateOrganizationIntelligenceConfig",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = StorageControlRestTransport._UpdateOrganizationIntelligenceConfig._get_response(
                self._host,
                metadata,
                query_params,
                self._session,
                timeout,
                transcoded_request,
                body,
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.IntelligenceConfig()
            pb_resp = storage_control.IntelligenceConfig.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_update_organization_intelligence_config(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_update_organization_intelligence_config_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = storage_control.IntelligenceConfig.to_json(
                        response
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.update_organization_intelligence_config",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "UpdateOrganizationIntelligenceConfig",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    class _UpdateProjectIntelligenceConfig(
        _BaseStorageControlRestTransport._BaseUpdateProjectIntelligenceConfig,
        StorageControlRestStub,
    ):
        def __hash__(self):
            return hash("StorageControlRestTransport.UpdateProjectIntelligenceConfig")

        @staticmethod
        def _get_response(
            host,
            metadata,
            query_params,
            session,
            timeout,
            transcoded_request,
            body=None,
        ):
            uri = transcoded_request["uri"]
            method = transcoded_request["method"]
            headers = dict(metadata)
            headers["Content-Type"] = "application/json"
            response = getattr(session, method)(
                "{host}{uri}".format(host=host, uri=uri),
                timeout=timeout,
                headers=headers,
                params=rest_helpers.flatten_query_params(query_params, strict=True),
                data=body,
            )
            return response

        def __call__(
            self,
            request: storage_control.UpdateProjectIntelligenceConfigRequest,
            *,
            retry: OptionalRetry = gapic_v1.method.DEFAULT,
            timeout: Optional[float] = None,
            metadata: Sequence[Tuple[str, Union[str, bytes]]] = (),
        ) -> storage_control.IntelligenceConfig:
            r"""Call the update project
            intelligence config method over HTTP.

                Args:
                    request (~.storage_control.UpdateProjectIntelligenceConfigRequest):
                        The request object. Request message to update the ``IntelligenceConfig``
                    resource associated with your project.

                    **IAM Permissions**:

                    Requires ``storage.intelligenceConfigs.update``
                    `IAM <https://cloud.google.com/iam/docs/overview#permissions>`__
                    permission on the folder.
                    retry (google.api_core.retry.Retry): Designation of what errors, if any,
                        should be retried.
                    timeout (float): The timeout for this request.
                    metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                        sent along with the request as metadata. Normally, each value must be of type `str`,
                        but for metadata keys ending with the suffix `-bin`, the corresponding values must
                        be of type `bytes`.

                Returns:
                    ~.storage_control.IntelligenceConfig:
                        The ``IntelligenceConfig`` resource associated with your
                    organization, folder, or project.

            """

            http_options = _BaseStorageControlRestTransport._BaseUpdateProjectIntelligenceConfig._get_http_options()

            request, metadata = (
                self._interceptor.pre_update_project_intelligence_config(
                    request, metadata
                )
            )
            transcoded_request = _BaseStorageControlRestTransport._BaseUpdateProjectIntelligenceConfig._get_transcoded_request(
                http_options, request
            )

            body = _BaseStorageControlRestTransport._BaseUpdateProjectIntelligenceConfig._get_request_body_json(
                transcoded_request
            )

            # Jsonify the query params
            query_params = _BaseStorageControlRestTransport._BaseUpdateProjectIntelligenceConfig._get_query_params_json(
                transcoded_request
            )

            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                request_url = "{host}{uri}".format(
                    host=self._host, uri=transcoded_request["uri"]
                )
                method = transcoded_request["method"]
                try:
                    request_payload = type(request).to_json(request)
                except:
                    request_payload = None
                http_request = {
                    "payload": request_payload,
                    "requestMethod": method,
                    "requestUrl": request_url,
                    "headers": dict(metadata),
                }
                _LOGGER.debug(
                    f"Sending request for google.storage.control_v2.StorageControlClient.UpdateProjectIntelligenceConfig",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "UpdateProjectIntelligenceConfig",
                        "httpRequest": http_request,
                        "metadata": http_request["headers"],
                    },
                )

            # Send the request
            response = StorageControlRestTransport._UpdateProjectIntelligenceConfig._get_response(
                self._host,
                metadata,
                query_params,
                self._session,
                timeout,
                transcoded_request,
                body,
            )

            # In case of error, raise the appropriate core_exceptions.GoogleAPICallError exception
            # subclass.
            if response.status_code >= 400:
                raise core_exceptions.from_http_response(response)

            # Return the response
            resp = storage_control.IntelligenceConfig()
            pb_resp = storage_control.IntelligenceConfig.pb(resp)

            json_format.Parse(response.content, pb_resp, ignore_unknown_fields=True)

            resp = self._interceptor.post_update_project_intelligence_config(resp)
            response_metadata = [(k, str(v)) for k, v in response.headers.items()]
            resp, _ = (
                self._interceptor.post_update_project_intelligence_config_with_metadata(
                    resp, response_metadata
                )
            )
            if CLIENT_LOGGING_SUPPORTED and _LOGGER.isEnabledFor(
                logging.DEBUG
            ):  # pragma: NO COVER
                try:
                    response_payload = storage_control.IntelligenceConfig.to_json(
                        response
                    )
                except:
                    response_payload = None
                http_response = {
                    "payload": response_payload,
                    "headers": dict(response.headers),
                    "status": response.status_code,
                }
                _LOGGER.debug(
                    "Received response for google.storage.control_v2.StorageControlClient.update_project_intelligence_config",
                    extra={
                        "serviceName": "google.storage.control.v2.StorageControl",
                        "rpcName": "UpdateProjectIntelligenceConfig",
                        "metadata": http_response["headers"],
                        "httpResponse": http_response,
                    },
                )
            return resp

    @property
    def create_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.CreateAnywhereCacheRequest], operations_pb2.Operation
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._CreateAnywhereCache(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def create_folder(
        self,
    ) -> Callable[[storage_control.CreateFolderRequest], storage_control.Folder]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._CreateFolder(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def create_managed_folder(
        self,
    ) -> Callable[
        [storage_control.CreateManagedFolderRequest], storage_control.ManagedFolder
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._CreateManagedFolder(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def delete_folder(
        self,
    ) -> Callable[[storage_control.DeleteFolderRequest], empty_pb2.Empty]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._DeleteFolder(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def delete_folder_recursive(
        self,
    ) -> Callable[
        [storage_control.DeleteFolderRecursiveRequest], operations_pb2.Operation
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._DeleteFolderRecursive(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def delete_managed_folder(
        self,
    ) -> Callable[[storage_control.DeleteManagedFolderRequest], empty_pb2.Empty]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._DeleteManagedFolder(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def disable_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.DisableAnywhereCacheRequest], storage_control.AnywhereCache
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._DisableAnywhereCache(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def get_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.GetAnywhereCacheRequest], storage_control.AnywhereCache
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetAnywhereCache(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def get_folder(
        self,
    ) -> Callable[[storage_control.GetFolderRequest], storage_control.Folder]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetFolder(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def get_folder_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetFolderIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetFolderIntelligenceConfig(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def get_iam_policy(
        self,
    ) -> Callable[[iam_policy_pb2.GetIamPolicyRequest], policy_pb2.Policy]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetIamPolicy(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def get_intelligence_finding(
        self,
    ) -> Callable[
        [storage_control.GetIntelligenceFindingRequest],
        storage_control.IntelligenceFinding,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetIntelligenceFinding(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def get_intelligence_finding_revision(
        self,
    ) -> Callable[
        [storage_control.GetIntelligenceFindingRevisionRequest],
        storage_control.IntelligenceFindingRevision,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetIntelligenceFindingRevision(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def get_managed_folder(
        self,
    ) -> Callable[
        [storage_control.GetManagedFolderRequest], storage_control.ManagedFolder
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetManagedFolder(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def get_organization_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetOrganizationIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetOrganizationIntelligenceConfig(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def get_project_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.GetProjectIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetProjectIntelligenceConfig(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def get_storage_layout(
        self,
    ) -> Callable[
        [storage_control.GetStorageLayoutRequest], storage_control.StorageLayout
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._GetStorageLayout(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def list_anywhere_caches(
        self,
    ) -> Callable[
        [storage_control.ListAnywhereCachesRequest],
        storage_control.ListAnywhereCachesResponse,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._ListAnywhereCaches(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def list_folders(
        self,
    ) -> Callable[
        [storage_control.ListFoldersRequest], storage_control.ListFoldersResponse
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._ListFolders(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def list_intelligence_finding_revisions(
        self,
    ) -> Callable[
        [storage_control.ListIntelligenceFindingRevisionsRequest],
        storage_control.ListIntelligenceFindingRevisionsResponse,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._ListIntelligenceFindingRevisions(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def list_intelligence_findings(
        self,
    ) -> Callable[
        [storage_control.ListIntelligenceFindingsRequest],
        storage_control.ListIntelligenceFindingsResponse,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._ListIntelligenceFindings(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def list_managed_folders(
        self,
    ) -> Callable[
        [storage_control.ListManagedFoldersRequest],
        storage_control.ListManagedFoldersResponse,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._ListManagedFolders(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def pause_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.PauseAnywhereCacheRequest], storage_control.AnywhereCache
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._PauseAnywhereCache(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def rename_folder(
        self,
    ) -> Callable[[storage_control.RenameFolderRequest], operations_pb2.Operation]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._RenameFolder(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def resume_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.ResumeAnywhereCacheRequest], storage_control.AnywhereCache
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._ResumeAnywhereCache(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def set_iam_policy(
        self,
    ) -> Callable[[iam_policy_pb2.SetIamPolicyRequest], policy_pb2.Policy]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._SetIamPolicy(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def summarize_intelligence_findings(
        self,
    ) -> Callable[
        [storage_control.SummarizeIntelligenceFindingsRequest],
        storage_control.SummarizeIntelligenceFindingsResponse,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._SummarizeIntelligenceFindings(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def test_iam_permissions(
        self,
    ) -> Callable[
        [iam_policy_pb2.TestIamPermissionsRequest],
        iam_policy_pb2.TestIamPermissionsResponse,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._TestIamPermissions(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def update_anywhere_cache(
        self,
    ) -> Callable[
        [storage_control.UpdateAnywhereCacheRequest], operations_pb2.Operation
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._UpdateAnywhereCache(self._session, self._host, self._interceptor)  # type: ignore

    @property
    def update_folder_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateFolderIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._UpdateFolderIntelligenceConfig(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def update_organization_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateOrganizationIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._UpdateOrganizationIntelligenceConfig(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def update_project_intelligence_config(
        self,
    ) -> Callable[
        [storage_control.UpdateProjectIntelligenceConfigRequest],
        storage_control.IntelligenceConfig,
    ]:
        # The return type is fine, but mypy isn't sophisticated enough to determine what's going on here.
        # In C++ this would require a dynamic_cast
        return self._UpdateProjectIntelligenceConfig(
            self._session, self._host, self._interceptor
        )  # type: ignore

    @property
    def kind(self) -> str:
        return "rest"

    def close(self):
        self._session.close()


__all__ = ("StorageControlRestTransport",)
