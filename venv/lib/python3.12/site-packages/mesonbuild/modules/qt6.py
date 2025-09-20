# Copyright 2020 The Meson development team

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

from .qt import QtBaseModule
from . import ModuleInfo

if T.TYPE_CHECKING:
    from ..interpreter import Interpreter

class Qt6Module(QtBaseModule):

    INFO = ModuleInfo('qt6', '0.57.0')

    def __init__(self, interpreter: Interpreter):
        QtBaseModule.__init__(self, interpreter, qt_version=6)


def initialize(interp: Interpreter) -> Qt6Module:
    return Qt6Module(interp)
