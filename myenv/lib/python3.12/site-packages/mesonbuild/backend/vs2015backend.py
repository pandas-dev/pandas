# Copyright 2014-2016 The Meson development team

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

import typing as T

from .vs2010backend import Vs2010Backend
from ..mesonlib import MesonException

if T.TYPE_CHECKING:
    from ..build import Build
    from ..interpreter import Interpreter

class Vs2015Backend(Vs2010Backend):

    name = 'vs2015'

    def __init__(self, build: T.Optional[Build], interpreter: T.Optional[Interpreter]):
        super().__init__(build, interpreter)
        self.vs_version = '2015'
        self.sln_file_version = '12.00'
        self.sln_version_comment = '14'
        if self.environment is not None:
            # TODO: we assume host == build
            comps = self.environment.coredata.compilers.host
            if comps and all(c.id == 'intel-cl' for c in comps.values()):
                c = list(comps.values())[0]
                if c.version.startswith('19'):
                    self.platform_toolset = 'Intel C++ Compiler 19.0'
                else:
                    # We don't have support for versions older than 2019 right now.
                    raise MesonException('There is currently no support for ICL before 19, patches welcome.')
            if self.platform_toolset is None:
                self.platform_toolset = 'v140'
