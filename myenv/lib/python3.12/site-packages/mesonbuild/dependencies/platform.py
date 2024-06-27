# Copyright 2013-2017 The Meson development team

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This file contains the detection logic for external dependencies that are
# platform-specific (generally speaking).
from __future__ import annotations

from .base import DependencyTypeName, ExternalDependency, DependencyException
from .detect import packages
from ..mesonlib import MesonException
import typing as T

if T.TYPE_CHECKING:
    from ..environment import Environment

class AppleFrameworks(ExternalDependency):
    def __init__(self, env: 'Environment', kwargs: T.Dict[str, T.Any]) -> None:
        super().__init__(DependencyTypeName('appleframeworks'), env, kwargs)
        modules = kwargs.get('modules', [])
        if isinstance(modules, str):
            modules = [modules]
        if not modules:
            raise DependencyException("AppleFrameworks dependency requires at least one module.")
        self.frameworks = modules
        if not self.clib_compiler:
            raise DependencyException('No C-like compilers are available, cannot find the framework')
        self.is_found = True
        for f in self.frameworks:
            try:
                args = self.clib_compiler.find_framework(f, env, [])
            except MesonException as e:
                if 'non-clang' in str(e):
                    self.is_found = False
                    self.link_args = []
                    self.compile_args = []
                    return
                raise

            if args is not None:
                # No compile args are needed for system frameworks
                self.link_args += args
            else:
                self.is_found = False

    def log_info(self) -> str:
        return ', '.join(self.frameworks)

    @staticmethod
    def log_tried() -> str:
        return 'framework'

packages['appleframeworks'] = AppleFrameworks
