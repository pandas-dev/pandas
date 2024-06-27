# Copyright 2018 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys

import pkg_resources


def modify_path():
    """Modify the module search path."""
    # Only modify path on Windows.
    if os.name != "nt":
        return

    path = os.environ.get("PATH")
    if path is None:
        return

    try:
        extra_dll_dir = pkg_resources.resource_filename("google_crc32c", "extra-dll")
        if os.path.isdir(extra_dll_dir):
            # Python 3.6, 3.7 use path
            os.environ["PATH"] = path + os.pathsep + extra_dll_dir
            # Python 3.8+ uses add_dll_directory.
            if sys.version_info[0] == 3 and sys.version_info[1] >= 8:
                os.add_dll_directory(extra_dll_dir)
    except ImportError:
        pass

modify_path()
