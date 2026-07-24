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

import importlib.util
from typing import Set

_has_grpc = importlib.util.find_spec("grpc") is not None

# PEP 0810: Explicit Lazy Imports
# Python 3.15+ natively intercepts and defers these imports.
# Developers can disable this behavior and force eager imports.
# For more information, see:
# https://docs.python.org/3.15/library/sys.html#sys.set_lazy_imports_filter
# Older Python versions safely ignore this variable.
__lazy_modules__: Set[str] = {
    "google.api_core.gapic_v1.client_info",
    "google.api_core.gapic_v1.routing_header",
}
__all__ = ["client_info", "routing_header"]

if _has_grpc:
    __lazy_modules__.update(
        {
            "google.api_core.gapic_v1.config",
            "google.api_core.gapic_v1.config_async",
            "google.api_core.gapic_v1.method",
            "google.api_core.gapic_v1.method_async",
        }
    )

from google.api_core.gapic_v1 import (  # noqa: E402
    client_info,
    routing_header,
)

if _has_grpc:
    from google.api_core.gapic_v1 import (  # noqa: F401
        config,
        config_async,
        method,
        method_async,
    )

    __all__.extend(["config", "config_async", "method", "method_async"])
