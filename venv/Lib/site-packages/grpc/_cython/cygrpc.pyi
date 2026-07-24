# Copyright 2026 gRPC authors.
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

# This file serves as a type stub for the cygrpc Cython extension module.
# Since type checkers (like pyright/mypy) and IDEs cannot directly inspect the
# compiled C extension, this file provides necessary type information.
# The catch-all __getattr__ allows any attribute access to typecheck as Any.

from typing import Any

def __getattr__(name: str) -> Any: ...
