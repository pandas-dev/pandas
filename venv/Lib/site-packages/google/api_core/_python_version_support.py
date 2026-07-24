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

"""Code to check Python versions supported by Google Cloud Client Libraries."""

import datetime
import enum
import functools
from importlib import metadata
import logging
import warnings
import sys
import textwrap
from typing import Any, List, NamedTuple, Optional, Dict, Tuple


_LOGGER = logging.getLogger(__name__)


class PythonVersionStatus(enum.Enum):
    """Support status of a Python version in this client library artifact release.

    "Support", in this context, means that this release of a client library
    artifact is configured to run on the currently configured version of
    Python.
    """

    PYTHON_VERSION_STATUS_UNSPECIFIED = "PYTHON_VERSION_STATUS_UNSPECIFIED"

    PYTHON_VERSION_SUPPORTED = "PYTHON_VERSION_SUPPORTED"
    """This Python version is fully supported, so the artifact running on this
    version will have all features and bug fixes."""

    PYTHON_VERSION_DEPRECATED = "PYTHON_VERSION_DEPRECATED"
    """This Python version is still supported, but support will end within a
    year. At that time, there will be no more releases for this artifact
    running under this Python version."""

    PYTHON_VERSION_EOL = "PYTHON_VERSION_EOL"
    """This Python version has reached its end of life in the Python community
    (see https://devguide.python.org/versions/), and this artifact will cease
    supporting this Python version within the next few releases."""

    PYTHON_VERSION_UNSUPPORTED = "PYTHON_VERSION_UNSUPPORTED"
    """This release of the client library artifact may not be the latest, since
    current releases no longer support this Python version."""


class VersionInfo(NamedTuple):
    """Hold release and support date information for a Python version."""

    version: str
    python_beta: Optional[datetime.date]
    python_start: datetime.date
    python_eol: datetime.date
    gapic_start: Optional[datetime.date] = None  # unused
    gapic_deprecation: Optional[datetime.date] = None
    gapic_end: Optional[datetime.date] = None
    dep_unpatchable_cve: Optional[datetime.date] = None  # unused


PYTHON_VERSIONS: List[VersionInfo] = [
    # Refer to https://devguide.python.org/versions/ and the PEPs linked therefrom.
    VersionInfo(
        version="3.10",
        python_beta=datetime.date(2021, 5, 3),
        python_start=datetime.date(2021, 10, 4),
        python_eol=datetime.date(2026, 10, 4),  # TODO: specify day when announced
    ),
    VersionInfo(
        version="3.11",
        python_beta=datetime.date(2022, 5, 8),
        python_start=datetime.date(2022, 10, 24),
        python_eol=datetime.date(2027, 10, 24),  # TODO: specify day when announced
    ),
    VersionInfo(
        version="3.12",
        python_beta=datetime.date(2023, 5, 22),
        python_start=datetime.date(2023, 10, 2),
        python_eol=datetime.date(2028, 10, 2),  # TODO: specify day when announced
    ),
    VersionInfo(
        version="3.13",
        python_beta=datetime.date(2024, 5, 8),
        python_start=datetime.date(2024, 10, 7),
        python_eol=datetime.date(2029, 10, 7),  # TODO: specify day when announced
    ),
    VersionInfo(
        version="3.14",
        python_beta=datetime.date(2025, 5, 7),
        python_start=datetime.date(2025, 10, 7),
        python_eol=datetime.date(2030, 10, 7),  # TODO: specify day when announced
    ),
]

PYTHON_VERSION_INFO: Dict[Tuple[int, int], VersionInfo] = {}
for info in PYTHON_VERSIONS:
    major, minor = map(int, info.version.split("."))
    PYTHON_VERSION_INFO[(major, minor)] = info


LOWEST_TRACKED_VERSION = min(PYTHON_VERSION_INFO.keys())
_FAKE_PAST_DATE = datetime.date.min + datetime.timedelta(days=900)
_FAKE_PAST_VERSION = VersionInfo(
    version="0.0",
    python_beta=_FAKE_PAST_DATE,
    python_start=_FAKE_PAST_DATE,
    python_eol=_FAKE_PAST_DATE,
)
_FAKE_FUTURE_DATE = datetime.date.max - datetime.timedelta(days=900)
_FAKE_FUTURE_VERSION = VersionInfo(
    version="999.0",
    python_beta=_FAKE_FUTURE_DATE,
    python_start=_FAKE_FUTURE_DATE,
    python_eol=_FAKE_FUTURE_DATE,
)
DEPRECATION_WARNING_PERIOD = datetime.timedelta(days=365)
EOL_GRACE_PERIOD = datetime.timedelta(weeks=1)


def _flatten_message(text: str) -> str:
    """Dedent a multi-line string and flatten it into a single line."""
    return " ".join(textwrap.dedent(text).strip().split())


@functools.cache
def _cached_packages_distributions():
    return metadata.packages_distributions()


def _get_pypi_package_name(module_name):
    """Determine the PyPI package name for a given module name."""
    try:
        module_to_distributions = _cached_packages_distributions()

        if module_name in module_to_distributions:  # pragma: NO COVER
            return module_to_distributions[module_name][0]
    except Exception as e:  # pragma: NO COVER
        _LOGGER.info(
            "An error occurred while determining PyPI package name for %s: %s",
            module_name,
            e,
        )

    return None


def _get_distribution_and_import_packages(import_package: str) -> Tuple[str, Any]:
    """Return a pretty string with distribution & import package names."""
    distribution_package = _get_pypi_package_name(import_package)
    dependency_distribution_and_import_packages = (
        f"package {distribution_package} ({import_package})"
        if distribution_package
        else import_package
    )
    return dependency_distribution_and_import_packages, distribution_package


def check_python_version(
    package: str = "this package", today: Optional[datetime.date] = None
) -> PythonVersionStatus:
    """Check the running Python version and issue a support warning if needed.

    Args:
        today: The date to check against. Defaults to the current date.

    Returns:
        The support status of the current Python version.
    """
    today = today or datetime.date.today()

    python_version = sys.version_info
    version_tuple = (python_version.major, python_version.minor)
    py_version_str = sys.version.split()[0]

    version_info = PYTHON_VERSION_INFO.get(version_tuple)

    if not version_info:
        if version_tuple < LOWEST_TRACKED_VERSION:
            version_info = _FAKE_PAST_VERSION
        else:
            version_info = _FAKE_FUTURE_VERSION

    gapic_deprecation = version_info.gapic_deprecation or (
        version_info.python_eol - DEPRECATION_WARNING_PERIOD
    )
    gapic_end = version_info.gapic_end or (version_info.python_eol + EOL_GRACE_PERIOD)

    def min_python(date: datetime.date) -> str:
        """Find the minimum supported Python version for a given date."""
        for version, info in sorted(PYTHON_VERSION_INFO.items()):
            if info.python_start <= date < info.python_eol:
                return f"{version[0]}.{version[1]}"
        return "at a currently supported version [https://devguide.python.org/versions]"

    # Resolve the pretty package label lazily so we avoid any work on
    # the happy path (supported Python version, no warning needed).
    def get_package_label():
        label, _ = _get_distribution_and_import_packages(package)
        return label

    if gapic_end < today:
        package_label = get_package_label()
        message = _flatten_message(
            f"""
            You are using a non-supported Python version ({py_version_str}).
            Google will not post any further updates to {package_label}
            supporting this Python version. Please upgrade to the latest Python
            version, or at least Python {min_python(today)}, and then update
            {package_label}.
            """
        )
        warnings.warn(message, FutureWarning)
        return PythonVersionStatus.PYTHON_VERSION_UNSUPPORTED

    eol_date = version_info.python_eol + EOL_GRACE_PERIOD
    if eol_date <= today <= gapic_end:
        package_label = get_package_label()
        message = _flatten_message(
            f"""
            You are using a Python version ({py_version_str})
            past its end of life. Google will update {package_label}
            with critical bug fixes on a best-effort basis, but not
            with any other fixes or features. Please upgrade
            to the latest Python version, or at least Python
            {min_python(today)}, and then update {package_label}.
            """
        )
        warnings.warn(message, FutureWarning)
        return PythonVersionStatus.PYTHON_VERSION_EOL

    if gapic_deprecation <= today <= gapic_end:
        package_label = get_package_label()
        message = _flatten_message(
            f"""
            You are using a Python version ({py_version_str}) which Google will
            stop supporting in new releases of {package_label} once it reaches
            its end of life ({version_info.python_eol}). Please upgrade to the
            latest Python version, or at least Python
            {min_python(version_info.python_eol)}, to continue receiving updates
            for {package_label} past that date.
            """
        )
        warnings.warn(message, FutureWarning)
        return PythonVersionStatus.PYTHON_VERSION_DEPRECATED

    return PythonVersionStatus.PYTHON_VERSION_SUPPORTED
