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
# pylint: disable=keyword-arg-before-vararg
from enum import Enum
from typing import Annotated, Any, Callable, Dict, List, Optional, Union

from pydantic import (
    BeforeValidator,
    Field,
    PlainSerializer,
    WithJsonSchema,
    model_validator,
)

from pyiceberg.schema import Schema
from pyiceberg.transforms import IdentityTransform, Transform, parse_transform
from pyiceberg.typedef import IcebergBaseModel
from pyiceberg.types import IcebergType


class SortDirection(Enum):
    ASC = "asc"
    DESC = "desc"

    def __str__(self) -> str:
        """Return the string representation of the SortDirection class."""
        return self.name

    def __repr__(self) -> str:
        """Return the string representation of the SortDirection class."""
        return f"SortDirection.{self.name}"


class NullOrder(Enum):
    NULLS_FIRST = "nulls-first"
    NULLS_LAST = "nulls-last"

    def __str__(self) -> str:
        """Return the string representation of the NullOrder class."""
        return self.name.replace("_", " ")

    def __repr__(self) -> str:
        """Return the string representation of the NullOrder class."""
        return f"NullOrder.{self.name}"


class SortField(IcebergBaseModel):
    """Sort order field.

    Args:
      source_id (int): Source column id from the table’s schema.
      transform (str): Transform that is used to produce values to be sorted on from the source column.
                       This is the same transform as described in partition transforms.
      direction (SortDirection): Sort direction, that can only be either asc or desc.
      null_order (NullOrder): Null order that describes the order of null values when sorted. Can only be either nulls-first or nulls-last.
    """

    def __init__(
        self,
        source_id: Optional[int] = None,
        transform: Optional[Union[Transform[Any, Any], Callable[[IcebergType], Transform[Any, Any]]]] = None,
        direction: Optional[SortDirection] = None,
        null_order: Optional[NullOrder] = None,
        **data: Any,
    ):
        if source_id is not None:
            data["source-id"] = source_id
        if transform is not None:
            data["transform"] = transform
        if direction is not None:
            data["direction"] = direction
        if null_order is not None:
            data["null-order"] = null_order
        super().__init__(**data)

    @model_validator(mode="before")
    def set_null_order(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        values["direction"] = values["direction"] if values.get("direction") else SortDirection.ASC
        if not values.get("null-order"):
            values["null-order"] = NullOrder.NULLS_FIRST if values["direction"] == SortDirection.ASC else NullOrder.NULLS_LAST
        return values

    @model_validator(mode="before")
    @classmethod
    def map_source_ids_onto_source_id(cls, data: Any) -> Any:
        if isinstance(data, dict):
            if "source-id" not in data and (source_ids := data["source-ids"]):
                if isinstance(source_ids, list):
                    if len(source_ids) == 0:
                        raise ValueError("Empty source-ids is not allowed")
                    if len(source_ids) > 1:
                        raise ValueError("Multi argument transforms are not yet supported")
                    data["source-id"] = source_ids[0]
        return data

    source_id: int = Field(alias="source-id")
    transform: Annotated[  # type: ignore
        Transform,
        BeforeValidator(parse_transform),
        PlainSerializer(lambda c: str(c), return_type=str),  # pylint: disable=W0108
        WithJsonSchema({"type": "string"}, mode="serialization"),
    ] = Field(default=IdentityTransform())
    direction: SortDirection = Field()
    null_order: NullOrder = Field(alias="null-order")

    def __str__(self) -> str:
        """Return the string representation of the SortField class."""
        if isinstance(self.transform, IdentityTransform):
            # In the case of an identity transform, we can omit the transform
            return f"{self.source_id} {self.direction} {self.null_order}"
        else:
            return f"{self.transform}({self.source_id}) {self.direction} {self.null_order}"


INITIAL_SORT_ORDER_ID = 1


class SortOrder(IcebergBaseModel):
    """Describes how the data is sorted within the table.

    Users can sort their data within partitions by columns to gain performance.

    The order of the sort fields within the list defines the order in which the sort is applied to the data.

    Args:
      fields (List[SortField]): The fields how the table is sorted.

    Keyword Args:
      order_id (int): An unique id of the sort-order of a table.
    """

    order_id: int = Field(alias="order-id", default=INITIAL_SORT_ORDER_ID)
    fields: List[SortField] = Field(default_factory=list)

    def __init__(self, *fields: SortField, **data: Any):
        if fields:
            data["fields"] = fields
        super().__init__(**data)

    @property
    def is_unsorted(self) -> bool:
        return len(self.fields) == 0

    def __str__(self) -> str:
        """Return the string representation of the SortOrder class."""
        result_str = "["
        if self.fields:
            result_str += "\n  " + "\n  ".join([str(field) for field in self.fields]) + "\n"
        result_str += "]"
        return result_str

    def __repr__(self) -> str:
        """Return the string representation of the SortOrder class."""
        fields = f"{', '.join(repr(column) for column in self.fields)}, " if self.fields else ""
        return f"SortOrder({fields}order_id={self.order_id})"


UNSORTED_SORT_ORDER_ID = 0
UNSORTED_SORT_ORDER = SortOrder(order_id=UNSORTED_SORT_ORDER_ID)


def assign_fresh_sort_order_ids(
    sort_order: SortOrder, old_schema: Schema, fresh_schema: Schema, sort_order_id: int = INITIAL_SORT_ORDER_ID
) -> SortOrder:
    if sort_order.is_unsorted:
        return UNSORTED_SORT_ORDER

    fresh_fields = []
    for field in sort_order.fields:
        original_field = old_schema.find_column_name(field.source_id)
        if original_field is None:
            raise ValueError(f"Could not find in old schema: {field}")
        fresh_field = fresh_schema.find_field(original_field)
        if fresh_field is None:
            raise ValueError(f"Could not find field in fresh schema: {original_field}")
        fresh_fields.append(
            SortField(
                source_id=fresh_field.field_id,
                transform=field.transform,
                direction=field.direction,
                null_order=field.null_order,
            )
        )

    return SortOrder(*fresh_fields, order_id=sort_order_id)
