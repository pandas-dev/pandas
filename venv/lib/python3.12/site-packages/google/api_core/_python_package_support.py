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

"""Code to check versions of dependencies used by Google Cloud Client Libraries."""

import warnings
import sys
from typing import Optional, Tuple

from collections import namedtuple

from ._python_version_support import (
    _flatten_message,
    _get_distribution_and_import_packages,
)

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    # TODO(https://github.com/googleapis/python-api-core/issues/835): Remove
    # this code path once we drop support for Python 3.7
    import importlib_metadata as metadata

ParsedVersion = Tuple[int, ...]

# Here we list all the packages for which we want to issue warnings
# about deprecated and unsupported versions.
DependencyConstraint = namedtuple(
    "DependencyConstraint",
    ["package_name", "minimum_fully_supported_version", "recommended_version"],
)
_PACKAGE_DEPENDENCY_WARNINGS = [
    DependencyConstraint(
        "google.protobuf",
        minimum_fully_supported_version="4.25.8",
        recommended_version="6.x",
    )
]


DependencyVersion = namedtuple("DependencyVersion", ["version", "version_string"])
# Version string we provide in a DependencyVersion when we can't determine the version of a
# package.
UNKNOWN_VERSION_STRING = "--"


def parse_version_to_tuple(version_string: str) -> ParsedVersion:
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


def get_dependency_version(
    dependency_name: str,
) -> DependencyVersion:
    """Get the parsed version of an installed package dependency.

    This function checks for an installed package and returns its version
    as a comparable tuple of integers object for safe comparison. It handles
    both modern (Python 3.8+) and legacy (Python 3.7) environments.

    Args:
        dependency_name: The distribution name of the package (e.g., 'requests').

    Returns:
        A DependencyVersion namedtuple with `version`  (a tuple of integers) and
        `version_string` attributes, or `DependencyVersion(None,
        UNKNOWN_VERSION_STRING)` if the package is not found or
        another error occurs during version discovery.

    """
    try:
        version_string: str = metadata.version(dependency_name)
        parsed_version = parse_version_to_tuple(version_string)
        return DependencyVersion(parsed_version, version_string)
    except Exception:
        # Catch exceptions from metadata.version() (e.g., PackageNotFoundError)
        # or errors during parse_version_to_tuple
        return DependencyVersion(None, UNKNOWN_VERSION_STRING)


def warn_deprecation_for_versions_less_than(
    consumer_import_package: str,
    dependency_import_package: str,
    minimum_fully_supported_version: str,
    recommended_version: Optional[str] = None,
    message_template: Optional[str] = None,
):
    """Issue any needed deprecation warnings for `dependency_import_package`.

    If `dependency_import_package` is installed at a version less than
    `minimum_fully_supported_version`, this issues a warning using either a
    default `message_template` or one provided by the user. The
    default `message_template` informs the user that they will not receive
    future updates for `consumer_import_package` if
    `dependency_import_package` is somehow pinned to a version lower
    than `minimum_fully_supported_version`.

    Args:
      consumer_import_package: The import name of the package that
        needs `dependency_import_package`.
      dependency_import_package: The import name of the dependency to check.
      minimum_fully_supported_version: The dependency_import_package version number
        below which a deprecation warning will be logged.
      recommended_version: If provided, the recommended next version, which
        could be higher than `minimum_fully_supported_version`.
      message_template: A custom default message template to replace
        the default. This `message_template` is treated as an
        f-string, where the following variables are defined:
        `dependency_import_package`, `consumer_import_package` and
        `dependency_distribution_package` and
        `consumer_distribution_package` and `dependency_package`,
        `consumer_package` , which contain the import packages, the
        distribution packages, and pretty string with both the
        distribution and import packages for the dependency and the
        consumer, respectively; and `minimum_fully_supported_version`,
        `version_used`, and `version_used_string`, which refer to supported
        and currently-used versions of the dependency.

    """
    if (
        not consumer_import_package
        or not dependency_import_package
        or not minimum_fully_supported_version
    ):  # pragma: NO COVER
        return

    dependency_version = get_dependency_version(dependency_import_package)
    if not dependency_version.version:
        return

    if dependency_version.version < parse_version_to_tuple(
        minimum_fully_supported_version
    ):
        (
            dependency_package,
            dependency_distribution_package,
        ) = _get_distribution_and_import_packages(dependency_import_package)
        (
            consumer_package,
            consumer_distribution_package,
        ) = _get_distribution_and_import_packages(consumer_import_package)

        recommendation = (
            " (we recommend {recommended_version})" if recommended_version else ""
        )
        message_template = message_template or _flatten_message(
            """
            DEPRECATION: Package {consumer_package} depends on
            {dependency_package}, currently installed at version
            {version_used_string}. Future updates to
            {consumer_package} will require {dependency_package} at
            version {minimum_fully_supported_version} or
            higher{recommendation}. Please ensure that either (a) your
            Python environment doesn't pin the version of
            {dependency_package}, so that updates to
            {consumer_package} can require the higher version, or (b)
            you manually update your Python environment to use at
            least version {minimum_fully_supported_version} of
            {dependency_package}.
            """
        )
        warnings.warn(
            message_template.format(
                consumer_import_package=consumer_import_package,
                dependency_import_package=dependency_import_package,
                consumer_distribution_package=consumer_distribution_package,
                dependency_distribution_package=dependency_distribution_package,
                dependency_package=dependency_package,
                consumer_package=consumer_package,
                minimum_fully_supported_version=minimum_fully_supported_version,
                recommendation=recommendation,
                version_used=dependency_version.version,
                version_used_string=dependency_version.version_string,
            ),
            FutureWarning,
        )


def check_dependency_versions(
    consumer_import_package: str, *package_dependency_warnings: DependencyConstraint
):
    """Bundle checks for all package dependencies.

    This function can be called by all consumers of google.api_core,
    to emit needed deprecation warnings for any of their
    dependencies. The dependencies to check can be passed as arguments, or if
    none are provided, it will default to the list in
    `_PACKAGE_DEPENDENCY_WARNINGS`.

    Args:
      consumer_import_package: The distribution name of the calling package, whose
        dependencies we're checking.
      *package_dependency_warnings: A variable number of DependencyConstraint
        objects, each specifying a dependency to check.
    """
    if not package_dependency_warnings:
        package_dependency_warnings = tuple(_PACKAGE_DEPENDENCY_WARNINGS)
    for package_info in package_dependency_warnings:
        warn_deprecation_for_versions_less_than(
            consumer_import_package,
            package_info.package_name,
            package_info.minimum_fully_supported_version,
            recommended_version=package_info.recommended_version,
        )
