# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
from __future__ import annotations

from abc import abstractmethod
from datetime import date, datetime
from decimal import Decimal
from functools import lru_cache
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Generic,
    Literal,
    Optional,
    Protocol,
    Set,
    Tuple,
    TypeVar,
    Union,
    runtime_checkable,
)
from uuid import UUID

from pydantic import BaseModel, ConfigDict, RootModel
from typing_extensions import TypeAlias

if TYPE_CHECKING:
    from pyiceberg.types import StructType


class FrozenDict(Dict[Any, Any]):
    def __setitem__(self, instance: Any, value: Any) -> None:
        """Assign a value to a FrozenDict."""
        raise AttributeError("FrozenDict does not support assignment")

    def update(self, *args: Any, **kwargs: Any) -> None:
        raise AttributeError("FrozenDict does not support .update()")


UTF8 = "utf-8"

EMPTY_DICT = FrozenDict()

K = TypeVar("K")
V = TypeVar("V")


# from https://stackoverflow.com/questions/2912231/is-there-a-clever-way-to-pass-the-key-to-defaultdicts-default-factory
class KeyDefaultDict(Dict[K, V]):
    def __init__(self, default_factory: Callable[[K], V]):
        super().__init__()
        self.default_factory = default_factory

    def __missing__(self, key: K) -> V:
        """Define behavior if you access a non-existent key in a KeyDefaultDict."""
        val = self.default_factory(key)
        self[key] = val
        return val


Identifier = Tuple[str, ...]
"""A tuple of strings representing a table identifier.

Each string in the tuple represents a part of the table's unique path. For example,
a table in a namespace might be identified as:

    ("namespace", "table_name")

Examples:
    >>> identifier: Identifier = ("namespace", "table_name")
"""

Properties = Dict[str, Any]
"""A dictionary type for properties in PyIceberg."""


RecursiveDict = Dict[str, Union[str, "RecursiveDict"]]
"""A recursive dictionary type for nested structures in PyIceberg."""

# Represents the literal value
L = TypeVar("L", str, bool, int, float, bytes, UUID, Decimal, datetime, date, covariant=True)


@runtime_checkable
class StructProtocol(Protocol):  # pragma: no cover
    """A generic protocol used by accessors to get and set at positions of an object."""

    @abstractmethod
    def __getitem__(self, pos: int) -> Any:
        """Fetch a value from a StructProtocol."""

    @abstractmethod
    def __setitem__(self, pos: int, value: Any) -> None:
        """Assign a value to a StructProtocol."""


class IcebergBaseModel(BaseModel):
    """
    This class extends the Pydantic BaseModel to set default values by overriding them.

    This is because we always want to set by_alias to True. In Python, the dash can't
    be used in variable names, and this is used throughout the Iceberg spec.

    The same goes for exclude_none, if a field is None we want to omit it from
    serialization, for example, the doc attribute on the NestedField object.
    Default non-null values will be serialized.

    This is recommended by Pydantic:
    https://pydantic-docs.helpmanual.io/usage/model_config/#change-behaviour-globally
    """

    model_config = ConfigDict(populate_by_name=True, frozen=True)

    def _exclude_private_properties(self, exclude: Optional[Set[str]] = None) -> Set[str]:
        # A small trick to exclude private properties. Properties are serialized by pydantic,
        # regardless if they start with an underscore.
        # This will look at the dict, and find the fields and exclude them
        return set.union(
            {field for field in self.__dict__ if field.startswith("_") and not field == "__root__"}, exclude or set()
        )

    def model_dump(
        self, exclude_none: bool = True, exclude: Optional[Set[str]] = None, by_alias: bool = True, **kwargs: Any
    ) -> Dict[str, Any]:
        return super().model_dump(
            exclude_none=exclude_none, exclude=self._exclude_private_properties(exclude), by_alias=by_alias, **kwargs
        )

    def model_dump_json(
        self, exclude_none: bool = True, exclude: Optional[Set[str]] = None, by_alias: bool = True, **kwargs: Any
    ) -> str:
        return super().model_dump_json(
            exclude_none=exclude_none, exclude=self._exclude_private_properties(exclude), by_alias=by_alias, **kwargs
        )


T = TypeVar("T")


class IcebergRootModel(RootModel[T], Generic[T]):
    """
    This class extends the Pydantic BaseModel to set default values by overriding them.

    This is because we always want to set by_alias to True. In Python, the dash can't
    be used in variable names, and this is used throughout the Iceberg spec.

    The same goes for exclude_none, if a field is None we want to omit it from
    serialization, for example, the doc attribute on the NestedField object.
    Default non-null values will be serialized.

    This is recommended by Pydantic:
    https://pydantic-docs.helpmanual.io/usage/model_config/#change-behaviour-globally
    """

    model_config = ConfigDict(frozen=True)


@lru_cache
def _get_struct_fields(struct_type: StructType) -> Tuple[str, ...]:
    return tuple(field.name for field in struct_type.fields)


class Record(StructProtocol):
    __slots__ = ("_position_to_field_name",)
    _position_to_field_name: Tuple[str, ...]

    def __init__(self, *data: Any, struct: Optional[StructType] = None, **named_data: Any) -> None:
        if struct is not None:
            self._position_to_field_name = _get_struct_fields(struct)
        elif named_data:
            # Order of named_data is preserved (PEP 468) so this can be used to generate the position dict
            self._position_to_field_name = tuple(named_data.keys())
        else:
            self._position_to_field_name = tuple(f"field{idx + 1}" for idx in range(len(data)))

        for idx, d in enumerate(data):
            self[idx] = d

        for field_name, d in named_data.items():
            self.__setattr__(field_name, d)

    def __setitem__(self, pos: int, value: Any) -> None:
        """Assign a value to a Record."""
        self.__setattr__(self._position_to_field_name[pos], value)

    def __getitem__(self, pos: int) -> Any:
        """Fetch a value from a Record."""
        return self.__getattribute__(self._position_to_field_name[pos])

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the Record class."""
        if not isinstance(other, Record):
            return False
        return self.__dict__ == other.__dict__

    def __repr__(self) -> str:
        """Return the string representation of the Record class."""
        return f"{self.__class__.__name__}[{', '.join(f'{key}={repr(value)}' for key, value in self.__dict__.items() if not key.startswith('_'))}]"

    def __len__(self) -> int:
        """Return the number of fields in the Record class."""
        return len(self._position_to_field_name)

    def __hash__(self) -> int:
        """Return hash value of the Record class."""
        return hash(str(self))


TableVersion: TypeAlias = Literal[1, 2, 3]
