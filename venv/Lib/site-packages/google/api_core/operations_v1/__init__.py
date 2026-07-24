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

"""Package for interacting with the google.longrunning.operations meta-API."""

import importlib.util
from typing import Set

try:
    _has_async_rest = (
        importlib.util.find_spec("google.auth.aio.transport.sessions") is not None
    )
except ModuleNotFoundError:
    _has_async_rest = False

# PEP 0810: Explicit Lazy Imports
# Python 3.15+ natively intercepts and defers these imports.
# Developers can disable this behavior and force eager imports.
# For more information, see:
# https://docs.python.org/3.15/library/sys.html#sys.set_lazy_imports_filter
# Older Python versions safely ignore this variable.
# NOTE: We statically define all modules here (including async ones) to ensure
# static analysis tools (mypy, pyright, Ruff) can easily parse them. If async
# support is not present, the imports are ignored, making their presence safe.
__lazy_modules__: Set[str] = {
    "google.api_core.operations_v1.abstract_operations_client",
    "google.api_core.operations_v1.operations_async_client",
    "google.api_core.operations_v1.operations_client",
    "google.api_core.operations_v1.transports.rest",
    "google.api_core.operations_v1.transports.rest_asyncio",
    "google.api_core.operations_v1.operations_rest_client_async",
}

__all__ = [
    "AbstractOperationsClient",
    "OperationsAsyncClient",
    "OperationsClient",
    "OperationsRestTransport",
]


from google.api_core.operations_v1.abstract_operations_client import (  # noqa: E402
    AbstractOperationsClient,
)
from google.api_core.operations_v1.operations_async_client import (  # noqa: E402
    OperationsAsyncClient,
)
from google.api_core.operations_v1.operations_client import (  # noqa: E402
    OperationsClient,
)
from google.api_core.operations_v1.transports.rest import (  # noqa: E402
    OperationsRestTransport,
)

if _has_async_rest:
    try:
        # On Python 3.15+, PEP 0810 lazy loading means these imports will succeed
        # instantly (returning a lazy proxy). Any actual ImportErrors (e.g., due to
        # missing aiohttp/auth dependencies) are deferred until the proxies are accessed.
        from google.api_core.operations_v1.transports.rest_asyncio import (  # noqa: E402, F401
            AsyncOperationsRestTransport,
        )
        from google.api_core.operations_v1.operations_rest_client_async import (  # noqa: E402, F401
            AsyncOperationsRestClient,
        )

        __all__.extend(["AsyncOperationsRestClient", "AsyncOperationsRestTransport"])
    except ImportError:
        # Fallback for older python/environments when importlib find_spec succeeds but actual import fails
        pass
