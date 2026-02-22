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
import sys

import google.api_core as api_core

from google.cloud.storage_control_v2 import gapic_version as package_version

__version__ = package_version.__version__

if sys.version_info >= (3, 8):  # pragma: NO COVER
    from importlib import metadata
else:  # pragma: NO COVER
    # TODO(https://github.com/googleapis/python-api-core/issues/835): Remove
    # this code path once we drop support for Python 3.7
    import importlib_metadata as metadata

from .services.storage_control import StorageControlAsyncClient, StorageControlClient
from .types.storage_control import (
    AnywhereCache,
    CommonLongRunningOperationMetadata,
    CreateAnywhereCacheMetadata,
    CreateAnywhereCacheRequest,
    CreateFolderRequest,
    CreateManagedFolderRequest,
    DeleteFolderRecursiveMetadata,
    DeleteFolderRecursiveRequest,
    DeleteFolderRequest,
    DeleteManagedFolderRequest,
    DisableAnywhereCacheRequest,
    Folder,
    GetAnywhereCacheRequest,
    GetFolderIntelligenceConfigRequest,
    GetFolderRequest,
    GetManagedFolderRequest,
    GetOrganizationIntelligenceConfigRequest,
    GetProjectIntelligenceConfigRequest,
    GetStorageLayoutRequest,
    IntelligenceConfig,
    ListAnywhereCachesRequest,
    ListAnywhereCachesResponse,
    ListFoldersRequest,
    ListFoldersResponse,
    ListManagedFoldersRequest,
    ListManagedFoldersResponse,
    ManagedFolder,
    PauseAnywhereCacheRequest,
    PendingRenameInfo,
    RenameFolderMetadata,
    RenameFolderRequest,
    ResumeAnywhereCacheRequest,
    StorageLayout,
    UpdateAnywhereCacheMetadata,
    UpdateAnywhereCacheRequest,
    UpdateFolderIntelligenceConfigRequest,
    UpdateOrganizationIntelligenceConfigRequest,
    UpdateProjectIntelligenceConfigRequest,
)

if hasattr(api_core, "check_python_version") and hasattr(
    api_core, "check_dependency_versions"
):  # pragma: NO COVER
    api_core.check_python_version("google.cloud.storage_control_v2")  # type: ignore
    api_core.check_dependency_versions("google.cloud.storage_control_v2")  # type: ignore
else:  # pragma: NO COVER
    # An older version of api_core is installed which does not define the
    # functions above. We do equivalent checks manually.
    try:
        import sys
        import warnings

        _py_version_str = sys.version.split()[0]
        _package_label = "google.cloud.storage_control_v2"
        if sys.version_info < (3, 9):
            warnings.warn(
                "You are using a non-supported Python version "
                + f"({_py_version_str}).  Google will not post any further "
                + f"updates to {_package_label} supporting this Python version. "
                + "Please upgrade to the latest Python version, or at "
                + f"least to Python 3.9, and then update {_package_label}.",
                FutureWarning,
            )
        if sys.version_info[:2] == (3, 9):
            warnings.warn(
                f"You are using a Python version ({_py_version_str}) "
                + f"which Google will stop supporting in {_package_label} in "
                + "January 2026. Please "
                + "upgrade to the latest Python version, or at "
                + "least to Python 3.10, before then, and "
                + f"then update {_package_label}.",
                FutureWarning,
            )

        def parse_version_to_tuple(version_string: str):
            """Safely converts a semantic version string to a comparable tuple of integers.
            Example: "4.25.8" -> (4, 25, 8)
            Ignores non-numeric parts and handles common version formats.
            Args:
                version_string: Version string in the format "x.y.z" or "x.y.z<suffix>"
            Returns:
                Tuple of integers for the parsed version string.
            """
            parts = []
            for part in version_string.split("."):
                try:
                    parts.append(int(part))
                except ValueError:
                    # If it's a non-numeric part (e.g., '1.0.0b1' -> 'b1'), stop here.
                    # This is a simplification compared to 'packaging.parse_version', but sufficient
                    # for comparing strictly numeric semantic versions.
                    break
            return tuple(parts)

        def _get_version(dependency_name):
            try:
                version_string: str = metadata.version(dependency_name)
                parsed_version = parse_version_to_tuple(version_string)
                return (parsed_version, version_string)
            except Exception:
                # Catch exceptions from metadata.version() (e.g., PackageNotFoundError)
                # or errors during parse_version_to_tuple
                return (None, "--")

        _dependency_package = "google.protobuf"
        _next_supported_version = "4.25.8"
        _next_supported_version_tuple = (4, 25, 8)
        _recommendation = " (we recommend 6.x)"
        (_version_used, _version_used_string) = _get_version(_dependency_package)
        if _version_used and _version_used < _next_supported_version_tuple:
            warnings.warn(
                f"Package {_package_label} depends on "
                + f"{_dependency_package}, currently installed at version "
                + f"{_version_used_string}. Future updates to "
                + f"{_package_label} will require {_dependency_package} at "
                + f"version {_next_supported_version} or higher{_recommendation}."
                + " Please ensure "
                + "that either (a) your Python environment doesn't pin the "
                + f"version of {_dependency_package}, so that updates to "
                + f"{_package_label} can require the higher version, or "
                + "(b) you manually update your Python environment to use at "
                + f"least version {_next_supported_version} of "
                + f"{_dependency_package}.",
                FutureWarning,
            )
    except Exception:
        warnings.warn(
            "Could not determine the version of Python "
            + "currently being used. To continue receiving "
            + "updates for {_package_label}, ensure you are "
            + "using a supported version of Python; see "
            + "https://devguide.python.org/versions/"
        )

__all__ = (
    "StorageControlAsyncClient",
    "AnywhereCache",
    "CommonLongRunningOperationMetadata",
    "CreateAnywhereCacheMetadata",
    "CreateAnywhereCacheRequest",
    "CreateFolderRequest",
    "CreateManagedFolderRequest",
    "DeleteFolderRecursiveMetadata",
    "DeleteFolderRecursiveRequest",
    "DeleteFolderRequest",
    "DeleteManagedFolderRequest",
    "DisableAnywhereCacheRequest",
    "Folder",
    "GetAnywhereCacheRequest",
    "GetFolderIntelligenceConfigRequest",
    "GetFolderRequest",
    "GetManagedFolderRequest",
    "GetOrganizationIntelligenceConfigRequest",
    "GetProjectIntelligenceConfigRequest",
    "GetStorageLayoutRequest",
    "IntelligenceConfig",
    "ListAnywhereCachesRequest",
    "ListAnywhereCachesResponse",
    "ListFoldersRequest",
    "ListFoldersResponse",
    "ListManagedFoldersRequest",
    "ListManagedFoldersResponse",
    "ManagedFolder",
    "PauseAnywhereCacheRequest",
    "PendingRenameInfo",
    "RenameFolderMetadata",
    "RenameFolderRequest",
    "ResumeAnywhereCacheRequest",
    "StorageControlClient",
    "StorageLayout",
    "UpdateAnywhereCacheMetadata",
    "UpdateAnywhereCacheRequest",
    "UpdateFolderIntelligenceConfigRequest",
    "UpdateOrganizationIntelligenceConfigRequest",
    "UpdateProjectIntelligenceConfigRequest",
)
