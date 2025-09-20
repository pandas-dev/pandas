#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.
"""
This is a singleton metaclass that can be used to cache and reuse existing objects.

In the Iceberg codebase we have a lot of objects that are stateless (for example Types such as StringType,
BooleanType etc). FixedTypes have arguments (eg. Fixed[22]) that we also make part of the key when caching
the newly created object.

The Singleton uses a metaclass which essentially defines a new type. When the Type gets created, it will first
evaluate the `__call__` method with all the arguments. If we already initialized a class earlier, we'll just
return it.

More information on metaclasses: https://docs.python.org/3/reference/datamodel.html#metaclasses
"""

from typing import Any, ClassVar, Dict


def _convert_to_hashable_type(element: Any) -> Any:
    if isinstance(element, dict):
        return tuple((_convert_to_hashable_type(k), _convert_to_hashable_type(v)) for k, v in element.items())
    elif isinstance(element, list):
        return tuple(map(_convert_to_hashable_type, element))
    return element


class Singleton:
    _instances: ClassVar[Dict] = {}  # type: ignore

    def __new__(cls, *args, **kwargs):  # type: ignore
        key = (cls, tuple(args), _convert_to_hashable_type(kwargs))
        if key not in cls._instances:
            cls._instances[key] = super().__new__(cls)
        return cls._instances[key]

    def __deepcopy__(self, memo: Dict[int, Any]) -> Any:
        """
        Prevent deep copy operations for singletons.

        The IcebergRootModel inherits from Pydantic RootModel,
        which has its own implementation of deepcopy. When deepcopy
        runs, it calls the RootModel __deepcopy__ method and ignores
        that it's a Singleton. To handle this, the order of inheritance
        is adjusted and a __deepcopy__ method is implemented for
        singletons that simply returns itself.
        """
        return self
