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

import uuid
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import date, datetime, time
from functools import cached_property, singledispatch
from typing import Annotated, Any, Generic, TypeVar
from urllib.parse import quote_plus

from pydantic import (
    BeforeValidator,
    Field,
    PlainSerializer,
    WithJsonSchema,
    model_validator,
)

from pyiceberg.exceptions import ValidationError
from pyiceberg.schema import Schema
from pyiceberg.transforms import (
    BucketTransform,
    DayTransform,
    HourTransform,
    IdentityTransform,
    MonthTransform,
    Transform,
    TruncateTransform,
    UnknownTransform,
    VoidTransform,
    YearTransform,
    parse_transform,
)
from pyiceberg.typedef import IcebergBaseModel, Record
from pyiceberg.types import (
    DateType,
    IcebergType,
    NestedField,
    PrimitiveType,
    StructType,
    TimestampType,
    TimestamptzType,
    TimeType,
    UnknownType,
    UUIDType,
)
from pyiceberg.utils.datetime import date_to_days, datetime_to_micros, time_to_micros

INITIAL_PARTITION_SPEC_ID = 0
PARTITION_FIELD_ID_START: int = 1000


class PartitionField(IcebergBaseModel):
    """PartitionField represents how one partition value is derived from the source column via transformation.

    Attributes:
        source_id(int): The source column id of table's schema.
        field_id(int): The partition field id across all the table partition specs.
        transform(Transform): The transform used to produce partition values from source column.
        name(str): The name of this partition field.
    """

    source_id: int = Field(alias="source-id")
    field_id: int = Field(alias="field-id")
    transform: Annotated[  # type: ignore
        Transform,
        BeforeValidator(parse_transform),
        PlainSerializer(lambda c: str(c), return_type=str),  # pylint: disable=W0108
        WithJsonSchema({"type": "string"}, mode="serialization"),
    ] = Field()
    name: str = Field()

    def __init__(
        self,
        source_id: int | None = None,
        field_id: int | None = None,
        transform: Transform[Any, Any] | None = None,
        name: str | None = None,
        **data: Any,
    ):
        if source_id is not None:
            data["source-id"] = source_id
        if field_id is not None:
            data["field-id"] = field_id
        if transform is not None:
            data["transform"] = transform
        if name is not None:
            data["name"] = name

        super().__init__(**data)

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

    def __str__(self) -> str:
        """Return the string representation of the PartitionField class."""
        return f"{self.field_id}: {self.name}: {self.transform}({self.source_id})"


class PartitionSpec(IcebergBaseModel):
    """
    PartitionSpec captures the transformation from table data to partition values.

    Attributes:
        spec_id(int): any change to PartitionSpec will produce a new specId.
        fields(Tuple[PartitionField): list of partition fields to produce partition values.
    """

    spec_id: int = Field(alias="spec-id", default=INITIAL_PARTITION_SPEC_ID)
    fields: tuple[PartitionField, ...] = Field(default_factory=tuple)

    def __init__(
        self,
        *fields: PartitionField,
        **data: Any,
    ):
        if fields:
            data["fields"] = tuple(fields)
        super().__init__(**data)

    def __eq__(self, other: Any) -> bool:
        """
        Produce a boolean to return True if two objects are considered equal.

        Note:
            Equality of PartitionSpec is determined by spec_id and partition fields only.
        """
        if not isinstance(other, PartitionSpec):
            return False
        return self.spec_id == other.spec_id and self.fields == other.fields

    def __str__(self) -> str:
        """
        Produce a human-readable string representation of PartitionSpec.

        Note:
            Only include list of partition fields in the PartitionSpec's string representation.
        """
        result_str = "["
        if self.fields:
            result_str += "\n  " + "\n  ".join([str(field) for field in self.fields]) + "\n"
        result_str += "]"
        return result_str

    def __repr__(self) -> str:
        """Return the string representation of the PartitionSpec class."""
        fields = f"{', '.join(repr(column) for column in self.fields)}, " if self.fields else ""
        return f"PartitionSpec({fields}spec_id={self.spec_id})"

    def is_unpartitioned(self) -> bool:
        return not self.fields

    @property
    def last_assigned_field_id(self) -> int:
        if self.fields:
            return max(pf.field_id for pf in self.fields)
        return PARTITION_FIELD_ID_START - 1

    @cached_property
    def source_id_to_fields_map(self) -> dict[int, list[PartitionField]]:
        source_id_to_fields_map: dict[int, list[PartitionField]] = {}
        for partition_field in self.fields:
            existing = source_id_to_fields_map.get(partition_field.source_id, [])
            existing.append(partition_field)
            source_id_to_fields_map[partition_field.source_id] = existing
        return source_id_to_fields_map

    def fields_by_source_id(self, field_id: int) -> list[PartitionField]:
        return self.source_id_to_fields_map.get(field_id, [])

    def compatible_with(self, other: PartitionSpec) -> bool:
        """Produce a boolean to return True if two PartitionSpec are considered compatible."""
        if self == other:
            return True
        if len(self.fields) != len(other.fields):
            return False
        return all(
            this_field.source_id == that_field.source_id
            and this_field.transform == that_field.transform
            and this_field.name == that_field.name
            for this_field, that_field in zip(self.fields, other.fields, strict=True)
        )

    def partition_type(self, schema: Schema) -> StructType:
        """Produce a struct of the PartitionSpec.

        The partition fields should be optional:

        - All partition transforms are required to produce null if the input value is null, so it can
          happen when the source column is optional.
        - Partition fields may be added later, in which case not all files would have the result field,
          and it may be null.

        There is a case where we can guarantee that a partition field in the first and only partition spec
        that uses a required source column will never be null, but it doesn't seem worth tracking this case.

        :param schema: The schema to bind to.
        :return: A StructType that represents the PartitionSpec, with a NestedField for each PartitionField.
        """
        nested_fields = []
        schema_ids = schema._lazy_id_to_field
        for field in self.fields:
            if source_field := schema_ids.get(field.source_id):
                result_type = field.transform.result_type(source_field.field_type)
                nested_fields.append(NestedField(field.field_id, field.name, result_type, required=source_field.required))
            else:
                # Since the source field has been drop we cannot determine the type
                nested_fields.append(NestedField(field.field_id, field.name, UnknownType()))
        return StructType(*nested_fields)

    def partition_to_path(self, data: Record, schema: Schema) -> str:
        partition_type = self.partition_type(schema)
        field_types = partition_type.fields

        field_strs = []
        value_strs = []
        for pos in range(len(self.fields)):
            partition_field = self.fields[pos]
            value_str = partition_field.transform.to_human_string(field_types[pos].field_type, value=data[pos])

            value_strs.append(quote_plus(value_str, safe=""))
            field_strs.append(quote_plus(partition_field.name, safe=""))

        path = "/".join([field_str + "=" + value_str for field_str, value_str in zip(field_strs, value_strs, strict=True)])
        return path

    def check_compatible(self, schema: Schema, allow_missing_fields: bool = False) -> None:
        # if the underlying field is dropped, we cannot check they are compatible -- continue
        schema_fields = schema._lazy_id_to_field
        parents = schema._lazy_id_to_parent

        for field in self.fields:
            source_field = schema_fields.get(field.source_id)

            if allow_missing_fields and source_field is None:
                continue

            if isinstance(field.transform, VoidTransform):
                continue

            if not source_field:
                raise ValidationError(f"Cannot find source column for partition field: {field}")

            source_type = source_field.field_type
            if not source_type.is_primitive:
                raise ValidationError(f"Cannot partition by non-primitive source field: {source_field}")
            if not field.transform.can_transform(source_type):
                raise ValidationError(
                    f"Invalid source field {source_field.name} with type {source_type} " + f"for transform: {field.transform}"
                )

            # The only valid parent types for a PartitionField are StructTypes. This must be checked recursively
            parent_id = parents.get(field.source_id)
            while parent_id:
                parent_type = schema.find_type(parent_id)
                if not parent_type.is_struct:
                    raise ValidationError(f"Invalid partition field parent: {parent_type}")
                parent_id = parents.get(parent_id)


UNPARTITIONED_PARTITION_SPEC = PartitionSpec(spec_id=0)


def validate_partition_name(
    field_name: str,
    partition_transform: Transform[Any, Any],
    source_id: int,
    schema: Schema,
    partition_names: set[str],
) -> None:
    """Validate that a partition field name doesn't conflict with schema field names."""
    try:
        schema_field = schema.find_field(field_name)
    except ValueError:
        return  # No conflict if field doesn't exist in schema

    if isinstance(partition_transform, (IdentityTransform, VoidTransform)):
        # For identity and void transforms, allow conflict only if sourced from the same schema field
        if schema_field.field_id != source_id:
            raise ValueError(f"Cannot create identity partition sourced from different field in schema: {field_name}")
    else:
        raise ValueError(f"Cannot create partition with a name that exists in schema: {field_name}")
    if not field_name:
        raise ValueError("Undefined name")
    if field_name in partition_names:
        raise ValueError(f"Partition name has to be unique: {field_name}")


def assign_fresh_partition_spec_ids(spec: PartitionSpec, old_schema: Schema, fresh_schema: Schema) -> PartitionSpec:
    partition_fields = []
    for pos, field in enumerate(spec.fields):
        original_column_name = old_schema.find_column_name(field.source_id)
        if original_column_name is None:
            raise ValueError(f"Could not find in old schema: {field}")
        fresh_field = fresh_schema.find_field(original_column_name)
        if fresh_field is None:
            raise ValueError(f"Could not find field in fresh schema: {original_column_name}")

        validate_partition_name(field.name, field.transform, fresh_field.field_id, fresh_schema, set())

        partition_fields.append(
            PartitionField(
                name=field.name,
                source_id=fresh_field.field_id,
                field_id=PARTITION_FIELD_ID_START + pos,
                transform=field.transform,
            )
        )
    return PartitionSpec(*partition_fields, spec_id=INITIAL_PARTITION_SPEC_ID)


T = TypeVar("T")


class PartitionSpecVisitor(Generic[T], ABC):
    @abstractmethod
    def identity(self, field_id: int, source_name: str, source_id: int) -> T:
        """Visit identity partition field."""

    @abstractmethod
    def bucket(self, field_id: int, source_name: str, source_id: int, num_buckets: int) -> T:
        """Visit bucket partition field."""

    @abstractmethod
    def truncate(self, field_id: int, source_name: str, source_id: int, width: int) -> T:
        """Visit truncate partition field."""

    @abstractmethod
    def year(self, field_id: int, source_name: str, source_id: int) -> T:
        """Visit year partition field."""

    @abstractmethod
    def month(self, field_id: int, source_name: str, source_id: int) -> T:
        """Visit month partition field."""

    @abstractmethod
    def day(self, field_id: int, source_name: str, source_id: int) -> T:
        """Visit day partition field."""

    @abstractmethod
    def hour(self, field_id: int, source_name: str, source_id: int) -> T:
        """Visit hour partition field."""

    @abstractmethod
    def always_null(self, field_id: int, source_name: str, source_id: int) -> T:
        """Visit void partition field."""

    @abstractmethod
    def unknown(self, field_id: int, source_name: str, source_id: int, transform: str) -> T:
        """Visit unknown partition field."""
        raise ValueError(f"Unknown transform is not supported: {transform}")


class _PartitionNameGenerator(PartitionSpecVisitor[str]):
    def identity(self, field_id: int, source_name: str, source_id: int) -> str:
        return source_name

    def bucket(self, field_id: int, source_name: str, source_id: int, num_buckets: int) -> str:
        return f"{source_name}_bucket_{num_buckets}"

    def truncate(self, field_id: int, source_name: str, source_id: int, width: int) -> str:
        return source_name + "_trunc_" + str(width)

    def year(self, field_id: int, source_name: str, source_id: int) -> str:
        return source_name + "_year"

    def month(self, field_id: int, source_name: str, source_id: int) -> str:
        return source_name + "_month"

    def day(self, field_id: int, source_name: str, source_id: int) -> str:
        return source_name + "_day"

    def hour(self, field_id: int, source_name: str, source_id: int) -> str:
        return source_name + "_hour"

    def always_null(self, field_id: int, source_name: str, source_id: int) -> str:
        return source_name + "_null"

    def unknown(self, field_id: int, source_name: str, source_id: int, transform: str) -> str:
        return super().unknown(field_id, source_name, source_id, transform)


R = TypeVar("R")


@singledispatch
def _visit(spec: PartitionSpec, schema: Schema, visitor: PartitionSpecVisitor[R]) -> list[R]:
    return [_visit_partition_field(schema, field, visitor) for field in spec.fields]


def _visit_partition_field(schema: Schema, field: PartitionField, visitor: PartitionSpecVisitor[R]) -> R:
    source_name = schema.find_column_name(field.source_id)
    if not source_name:
        raise ValueError(f"Could not find field with id {field.source_id}")

    transform = field.transform
    if isinstance(transform, IdentityTransform):
        return visitor.identity(field.field_id, source_name, field.source_id)
    elif isinstance(transform, BucketTransform):
        return visitor.bucket(field.field_id, source_name, field.source_id, transform.num_buckets)
    elif isinstance(transform, TruncateTransform):
        return visitor.truncate(field.field_id, source_name, field.source_id, transform.width)
    elif isinstance(transform, DayTransform):
        return visitor.day(field.field_id, source_name, field.source_id)
    elif isinstance(transform, HourTransform):
        return visitor.hour(field.field_id, source_name, field.source_id)
    elif isinstance(transform, MonthTransform):
        return visitor.month(field.field_id, source_name, field.source_id)
    elif isinstance(transform, YearTransform):
        return visitor.year(field.field_id, source_name, field.source_id)
    elif isinstance(transform, VoidTransform):
        return visitor.always_null(field.field_id, source_name, field.source_id)
    elif isinstance(transform, UnknownTransform):
        return visitor.unknown(field.field_id, source_name, field.source_id, repr(transform))
    else:
        raise ValueError(f"Unknown transform {transform}")


@dataclass(frozen=True)
class PartitionFieldValue:
    field: PartitionField
    value: Any


@dataclass(frozen=True)
class PartitionKey:
    field_values: list[PartitionFieldValue]
    partition_spec: PartitionSpec
    schema: Schema

    @cached_property
    def partition(self) -> Record:  # partition key transformed with iceberg internal representation as input
        iceberg_typed_key_values = []
        for raw_partition_field_value in self.field_values:
            partition_fields = self.partition_spec.source_id_to_fields_map[raw_partition_field_value.field.source_id]
            if len(partition_fields) != 1:
                raise ValueError(f"Cannot have redundant partitions: {partition_fields}")
            partition_field = partition_fields[0]
            iceberg_typed_key_values.append(
                partition_record_value(
                    partition_field=partition_field,
                    value=raw_partition_field_value.value,
                    schema=self.schema,
                )
            )
        return Record(*iceberg_typed_key_values)

    def to_path(self) -> str:
        return self.partition_spec.partition_to_path(self.partition, self.schema)


def partition_record_value(partition_field: PartitionField, value: Any, schema: Schema) -> Any:
    """
    Return the Partition Record representation of the value.

    The value is first converted to internal partition representation.
    For example, UUID is converted to bytes[16], DateType to days since epoch, etc.

    Then the corresponding PartitionField's transform is applied to return
    the final partition record value.
    """
    iceberg_type = schema.find_field(name_or_id=partition_field.source_id).field_type
    return _to_partition_representation(iceberg_type, value)


@singledispatch
def _to_partition_representation(type: IcebergType, value: Any) -> Any:
    """Strip the logical type into the physical type.

    It can be that the value is already transformed into its physical type,
    in this case it will return the original value. Keep in mind that the
    bucket transform always will return an int, but an identity transform
    can return date that still needs to be transformed into an int (days
    since epoch).
    """
    return TypeError(f"Unsupported partition field type: {type}")


@_to_partition_representation.register(TimestampType)
@_to_partition_representation.register(TimestamptzType)
def _(type: IcebergType, value: int | datetime | None) -> int | None:
    if value is None:
        return None
    elif isinstance(value, int):
        return value
    elif isinstance(value, datetime):
        return datetime_to_micros(value)
    else:
        raise ValueError(f"Type not recognized: {value}")


@_to_partition_representation.register(DateType)
def _(type: IcebergType, value: int | date | None) -> int | None:
    if value is None:
        return None
    elif isinstance(value, int):
        return value
    elif isinstance(value, date):
        return date_to_days(value)
    else:
        raise ValueError(f"Type not recognized: {value}")


@_to_partition_representation.register(TimeType)
def _(type: IcebergType, value: time | None) -> int | None:
    return time_to_micros(value) if value is not None else None


@_to_partition_representation.register(UUIDType)
def _(type: IcebergType, value: uuid.UUID | int | bytes | None) -> bytes | int | None:
    if value is None:
        return None
    elif isinstance(value, bytes):
        return value  # IdentityTransform
    elif isinstance(value, uuid.UUID):
        return value.bytes  # IdentityTransform
    elif isinstance(value, int):
        return value  # BucketTransform
    else:
        raise ValueError(f"Type not recognized: {value}")


@_to_partition_representation.register(PrimitiveType)
def _(type: IcebergType, value: Any | None) -> Any | None:
    return value
