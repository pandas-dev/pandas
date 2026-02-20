# Copyright 2017 Google LLC
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

"""Google API Core.

This package contains common code and utilities used by Google client libraries.
"""

from google.api_core import _python_package_support
from google.api_core import _python_version_support
from google.api_core import version as api_core_version

__version__ = api_core_version.__version__

# NOTE: Until dependent artifacts require this version of
# google.api_core, the functionality below must be made available
# manually in those artifacts.

# expose dependency checks for external callers
check_python_version = _python_version_support.check_python_version
check_dependency_versions = _python_package_support.check_dependency_versions
parse_version_to_tuple = _python_package_support.parse_version_to_tuple
warn_deprecation_for_versions_less_than = (
    _python_package_support.warn_deprecation_for_versions_less_than
)
DependencyConstraint = _python_package_support.DependencyConstraint

# perform version checks against api_core, and emit warnings if needed
check_python_version(package="google.api_core")
check_dependency_versions("google.api_core")
