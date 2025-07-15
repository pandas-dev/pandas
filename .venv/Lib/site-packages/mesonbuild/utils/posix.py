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
from __future__ import annotations

"""Posix specific implementations of mesonlib functionality."""

import fcntl
import typing as T

from .universal import MesonException
from .platform import BuildDirLock as BuildDirLockBase

__all__ = ['BuildDirLock']

class BuildDirLock(BuildDirLockBase):

    def __enter__(self) -> None:
        self.lockfile = open(self.lockfilename, 'w', encoding='utf-8')
        try:
            fcntl.flock(self.lockfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except (BlockingIOError, PermissionError):
            self.lockfile.close()
            raise MesonException('Some other Meson process is already using this build directory. Exiting.')
        except OSError as e:
            self.lockfile.close()
            raise MesonException(f'Failed to lock the build directory: {e.strerror}')

    def __exit__(self, *args: T.Any) -> None:
        fcntl.flock(self.lockfile, fcntl.LOCK_UN)
        self.lockfile.close()
