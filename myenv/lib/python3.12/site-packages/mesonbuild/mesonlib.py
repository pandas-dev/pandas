# SPDX-license-identifier: Apache-2.0
# Copyright 2012-2021 The Meson development team
# Copyright © 2021 Intel Corporation

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# pylint: skip-file

"""Helper functions and classes."""

import os

from .utils.core import *
from .utils.vsenv import *

from .utils.universal import *

# Here we import either the posix implementations, the windows implementations,
# or a generic no-op implementation
if os.name == 'posix':
    from .utils.posix import *
elif os.name == 'nt':
    from .utils.win32 import *
else:
    from .utils.platform import *
