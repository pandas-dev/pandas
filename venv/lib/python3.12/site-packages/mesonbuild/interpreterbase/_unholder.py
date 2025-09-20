# Copyright 2013-2021 The Meson development team

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

from .baseobjects import InterpreterObject, MesonInterpreterObject, ObjectHolder, HoldableTypes
from .exceptions import InvalidArguments
from ..mesonlib import HoldableObject, MesonBugException

if T.TYPE_CHECKING:
    from .baseobjects import TYPE_var

def _unholder(obj: InterpreterObject) -> TYPE_var:
    if isinstance(obj, ObjectHolder):
        assert isinstance(obj.held_object, HoldableTypes)
        return obj.held_object
    elif isinstance(obj, MesonInterpreterObject):
        return obj
    elif isinstance(obj, HoldableObject):
        raise MesonBugException(f'Argument {obj} of type {type(obj).__name__} is not held by an ObjectHolder.')
    elif isinstance(obj, InterpreterObject):
        raise InvalidArguments(f'Argument {obj} of type {type(obj).__name__} cannot be passed to a method or function')
    raise MesonBugException(f'Unknown object {obj} of type {type(obj).__name__} in the parameters.')
