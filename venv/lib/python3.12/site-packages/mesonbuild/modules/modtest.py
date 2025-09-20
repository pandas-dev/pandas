# Copyright 2015 The Meson development team

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

from . import ExtensionModule, ModuleInfo
from ..interpreterbase import noKwargs, noPosargs

if T.TYPE_CHECKING:
    from . import ModuleState
    from ..interpreter.interpreter import Interpreter
    from ..interpreterbase.baseobjects import TYPE_kwargs, TYPE_var


class TestModule(ExtensionModule):

    INFO = ModuleInfo('modtest')

    def __init__(self, interpreter: Interpreter) -> None:
        super().__init__(interpreter)
        self.methods.update({
            'print_hello': self.print_hello,
        })

    @noKwargs
    @noPosargs
    def print_hello(self, state: ModuleState, args: T.List[TYPE_var], kwargs: TYPE_kwargs) -> None:
        print('Hello from a Meson module')


def initialize(interp: Interpreter) -> TestModule:
    return TestModule(interp)
