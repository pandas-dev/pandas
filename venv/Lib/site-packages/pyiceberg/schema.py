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
# pylint: disable=W0511
from __future__ import annotations

import builtins
import itertools
from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass
from functools import cached_property, partial, singledispatch
from typing import (
    TYPE_CHECKING,
    Any,
    Generic,
    Literal,
    TypeVar,
)

from pydantic import Field, PrivateAttr, model_validator

from pyiceberg.exceptions import ResolveError
from pyiceberg.typedef import EMPTY_DICT, IcebergBaseModel, StructProtocol
from pyiceberg.types import (
    BinaryType,
    BooleanType,
    DateType,
    DecimalType,
    DoubleType,
    FixedType,
    FloatType,
    IcebergType,
    IntegerType,
    ListType,
    LongType,
    MapType,
    NestedField,
    PrimitiveType,
    StringType,
    StructType,
    TimestampNanoType,
    TimestampType,
    TimestamptzNanoType,
    TimestamptzType,
    TimeType,
    UnknownType,
    UUIDType,
)

if TYPE_CHECKING:
    import pyarrow as pa

    from pyiceberg.table.name_mapping import (
        NameMapping,
    )

T = TypeVar("T")
P = TypeVar("P")

INITIAL_SCHEMA_ID = 0

FIELD_ID_PROP = "field-id"
ICEBERG_FIELD_NAME_PROP = "iceberg-field-name"


class Schema(IcebergBaseModel):
    """A table Schema.

    Example:
        >>> from pyiceberg import schema
        >>> from pyiceberg import types
    """

    type: Literal["struct"] = "struct"
    fields: tuple[NestedField, ...] = Field(default_factory=tuple)
    schema_id: int = Field(alias="schema-id", default=INITIAL_SCHEMA_ID)
    identifier_field_ids: list[int] = Field(alias="identifier-field-ids", default_factory=list)

    _name_to_id: dict[str, int] = PrivateAttr()

    def __init__(self, *fields: NestedField, **data: Any):
        if fields:
            data["fields"] = fields
        super().__init__(**data)
        self._name_to_id = index_by_name(self)

    def __str__(self) -> str:
        """Return the string representation of the Schema class."""
        return "table {\n" + "\n".join(["  " + str(field) for field in self.columns]) + "\n}"

    def __repr__(self) -> str:
        """Return the string representation of the Schema class."""
        columns_repr = ", ".join(repr(column) for column in self.columns)
        return f"Schema({columns_repr}, schema_id={self.schema_id}, identifier_field_ids={self.identifier_field_ids})"

    def __len__(self) -> int:
        """Return the length of an instance of the Literal class."""
        return len(self.fields)

    def __eq__(self, other: Any) -> bool:
        """Return the equality of two instances of the Schema class."""
        if not other:
            return False

        if not isinstance(other, Schema):
            return False

        if len(self.columns) != len(other.columns):
            return False

        identifier_field_ids_is_equal = self.identifier_field_ids == other.identifier_field_ids
        schema_is_equal = all(lhs == rhs for lhs, rhs in zip(self.columns, other.columns, strict=True))

        return identifier_field_ids_is_equal and schema_is_equal

    @model_validator(mode="after")
    def check_schema(self) -> Schema:
        if self.identifier_field_ids:
            for field_id in self.identifier_field_ids:
                self._validate_identifier_field(field_id)

        return self

    @property
    def columns(self) -> tuple[NestedField, ...]:
        """A tuple of the top-level fields."""
        return self.fields

    @cached_property
    def _lazy_id_to_field(self) -> dict[int, NestedField]:
        """Return an index of field ID to NestedField instance.

        This is calculated once when called for the first time. Subsequent calls to this method will use a cached index.
        """
        return index_by_id(self)

    @cached_property
    def _lazy_id_to_parent(self) -> dict[int, int]:
        """Returns an index of field ID to parent field IDs.

        This is calculated once when called for the first time. Subsequent calls to this method will use a cached index.
        """
        return _index_parents(self)

    @cached_property
    def _lazy_name_to_id_lower(self) -> dict[str, int]:
        """Return an index of lower-case field names to field IDs.

        This is calculated once when called for the first time. Subsequent calls to this method will use a cached index.
        """
        return {name.lower(): field_id for name, field_id in self._name_to_id.items()}

    @cached_property
    def _lazy_id_to_name(self) -> dict[int, str]:
        """Return an index of field ID to full name.

        This is calculated once when called for the first time. Subsequent calls to this method will use a cached index.
        """
        return index_name_by_id(self)

    @cached_property
    def _lazy_id_to_accessor(self) -> dict[int, Accessor]:
        """Return an index of field ID to accessor.

        This is calculated once when called for the first time. Subsequent calls to this method will use a cached index.
        """
        return build_position_accessors(self)

    def as_struct(self) -> StructType:
        """Return the schema as a struct."""
        return StructType(*self.fields)

    def as_arrow(self) -> pa.Schema:
        """Return the schema as an Arrow schema."""
        from pyiceberg.io.pyarrow import schema_to_pyarrow

        return schema_to_pyarrow(self)

    def find_field(self, name_or_id: str | int, case_sensitive: bool = True) -> NestedField:
        """Find a field using a field name or field ID.

        Args:
            name_or_id (Union[str, int]): Either a field name or a field ID.
            case_sensitive (bool, optional): Whether to perform a case-sensitive lookup using a field name. Defaults to True.

        Raises:
            ValueError: When the value cannot be found.

        Returns:
            NestedField: The matched NestedField.
        """
        if isinstance(name_or_id, int):
            if name_or_id not in self._lazy_id_to_field:
                raise ValueError(f"Could not find field with id: {name_or_id}")
            return self._lazy_id_to_field[name_or_id]

        if case_sensitive:
            field_id = self._name_to_id.get(name_or_id)
        else:
            field_id = self._lazy_name_to_id_lower.get(name_or_id.lower())

        if field_id is None:
            raise ValueError(f"Could not find field with name {name_or_id}, case_sensitive={case_sensitive}")

        return self._lazy_id_to_field[field_id]

    def find_type(self, name_or_id: str | int, case_sensitive: bool = True) -> IcebergType:
        """Find a field type using a field name or field ID.

        Args:
            name_or_id (Union[str, int]): Either a field name or a field ID.
            case_sensitive (bool, optional): Whether to perform a case-sensitive lookup using a field name. Defaults to True.

        Returns:
            NestedField: The type of the matched NestedField.
        """
        field = self.find_field(name_or_id=name_or_id, case_sensitive=case_sensitive)
        if not field:
            raise ValueError(f"Could not find field with name or id {name_or_id}, case_sensitive={case_sensitive}")
        return field.field_type

    @property
    def highest_field_id(self) -> int:
        return max(self._lazy_id_to_name.keys(), default=0)

    @cached_property
    def name_mapping(self) -> NameMapping:
        from pyiceberg.table.name_mapping import create_mapping_from_schema

        return create_mapping_from_schema(self)

    def find_column_name(self, column_id: int) -> str | None:
        """Find a column name given a column ID.

        Args:
            column_id (int): The ID of the column.

        Returns:
            str: The column name (or None if the column ID cannot be found).
        """
        return self._lazy_id_to_name.get(column_id)

    @property
    def column_names(self) -> list[str]:
        """
        Return a list of all the column names, including nested fields.

        Excludes short names.

        Returns:
            List[str]: The column names.
        """
        return list(self._lazy_id_to_name.values())

    def accessor_for_field(self, field_id: int) -> Accessor:
        """Find a schema position accessor given a field ID.

        Args:
            field_id (int): The ID of the field.

        Raises:
            ValueError: When the value cannot be found.

        Returns:
            Accessor: An accessor for the given field ID.
        """
        if field_id not in self._lazy_id_to_accessor:
            raise ValueError(f"Could not find accessor for field with id: {field_id}")

        return self._lazy_id_to_accessor[field_id]

    def identifier_field_names(self) -> set[str]:
        """Return the names of the identifier fields.

        Returns:
            Set of names of the identifier fields
        """
        ids = set()
        for field_id in self.identifier_field_ids:
            column_name = self.find_column_name(field_id)
            if column_name is None:
                raise ValueError(f"Could not find identifier column id: {field_id}")
            ids.add(column_name)

        return ids

    def select(self, *names: str, case_sensitive: bool = True) -> Schema:
        """Return a new schema instance pruned to a subset of columns.

        Args:
            names (List[str]): A list of column names.
            case_sensitive (bool, optional): Whether to perform a case-sensitive lookup for each column name. Defaults to True.

        Returns:
            Schema: A new schema with pruned columns.

        Raises:
            ValueError: If a column is selected that doesn't exist.
        """
        try:
            if case_sensitive:
                ids = {self._name_to_id[name] for name in names}
            else:
                ids = {self._lazy_name_to_id_lower[name.lower()] for name in names}
        except KeyError as e:
            raise ValueError(f"Could not find column: {e}") from e

        return prune_columns(self, ids)

    @property
    def field_ids(self) -> set[int]:
        """Return the IDs of the current schema."""
        return set(self._name_to_id.values())

    def _validate_identifier_field(self, field_id: int) -> None:
        """Validate that the field with the given ID is a valid identifier field.

        Args:
          field_id: The ID of the field to validate.

        Raises:
          ValueError: If the field is not valid.
        """
        field = self.find_field(field_id)
        if not field.field_type.is_primitive:
            raise ValueError(f"Identifier field {field_id} invalid: not a primitive type field")

        if not field.required:
            raise ValueError(f"Identifier field {field_id} invalid: not a required field")

        if isinstance(field.field_type, (DoubleType, FloatType)):
            raise ValueError(f"Identifier field {field_id} invalid: must not be float or double field")

        # Check whether the nested field is in a chain of required struct fields
        # Exploring from root for better error message for list and map types
        parent_id = self._lazy_id_to_parent.get(field.field_id)
        fields: list[int] = []
        while parent_id is not None:
            fields.append(parent_id)
            parent_id = self._lazy_id_to_parent.get(parent_id)

        while fields:
            parent = self.find_field(fields.pop())
            if not parent.field_type.is_struct:
                raise ValueError(f"Cannot add field {field.name} as an identifier field: must not be nested in {parent}")

            if not parent.required:
                raise ValueError(
                    f"Cannot add field {field.name} as an identifier field: must not be nested in an optional field {parent}"
                )

    def check_format_version_compatibility(self, format_version: int) -> None:
        """Check that the schema is compatible for the given table format version.

        Args:
          format_version: The Iceberg table format version.

        Raises:
          ValueError: If the schema is not compatible for the format version.
        """
        for field in self._lazy_id_to_field.values():
            if format_version < field.field_type.minimum_format_version():
                raise ValueError(
                    f"{field.field_type} is only supported in {field.field_type.minimum_format_version()} or higher. "
                    f"Current format version is: {format_version}"
                )


class SchemaVisitor(Generic[T], ABC):
    def before_field(self, field: NestedField) -> None:
        """Override this method to perform an action immediately before visiting a field."""

    def after_field(self, field: NestedField) -> None:
        """Override this method to perform an action immediately after visiting a field."""

    def before_list_element(self, element: NestedField) -> None:
        """Override this method to perform an action immediately before visiting an element within a ListType."""
        self.before_field(element)

    def after_list_element(self, element: NestedField) -> None:
        """Override this method to perform an action immediately after visiting an element within a ListType."""
        self.after_field(element)

    def before_map_key(self, key: NestedField) -> None:
        """Override this method to perform an action immediately before visiting a key within a MapType."""
        self.before_field(key)

    def after_map_key(self, key: NestedField) -> None:
        """Override this method to perform an action immediately after visiting a key within a MapType."""
        self.after_field(key)

    def before_map_value(self, value: NestedField) -> None:
        """Override this method to perform an action immediately before visiting a value within a MapType."""
        self.before_field(value)

    def after_map_value(self, value: NestedField) -> None:
        """Override this method to perform an action immediately after visiting a value within a MapType."""
        self.after_field(value)

    @abstractmethod
    def schema(self, schema: Schema, struct_result: T) -> T:
        """Visit a Schema."""

    @abstractmethod
    def struct(self, struct: StructType, field_results: builtins.list[T]) -> T:
        """Visit a StructType."""

    @abstractmethod
    def field(self, field: NestedField, field_result: T) -> T:
        """Visit a NestedField."""

    @abstractmethod
    def list(self, list_type: ListType, element_result: T) -> T:
        """Visit a ListType."""

    @abstractmethod
    def map(self, map_type: MapType, key_result: T, value_result: T) -> T:
        """Visit a MapType."""

    @abstractmethod
    def primitive(self, primitive: PrimitiveType) -> T:
        """Visit a PrimitiveType."""


class PreOrderSchemaVisitor(Generic[T], ABC):
    @abstractmethod
    def schema(self, schema: Schema, struct_result: Callable[[], T]) -> T:
        """Visit a Schema."""

    @abstractmethod
    def struct(self, struct: StructType, field_results: builtins.list[Callable[[], T]]) -> T:
        """Visit a StructType."""

    @abstractmethod
    def field(self, field: NestedField, field_result: Callable[[], T]) -> T:
        """Visit a NestedField."""

    @abstractmethod
    def list(self, list_type: ListType, element_result: Callable[[], T]) -> T:
        """Visit a ListType."""

    @abstractmethod
    def map(self, map_type: MapType, key_result: Callable[[], T], value_result: Callable[[], T]) -> T:
        """Visit a MapType."""

    @abstractmethod
    def primitive(self, primitive: PrimitiveType) -> T:
        """Visit a PrimitiveType."""


class SchemaWithPartnerVisitor(Generic[P, T], ABC):
    def before_field(self, field: NestedField, field_partner: P | None) -> None:
        """Override this method to perform an action immediately before visiting a field."""

    def after_field(self, field: NestedField, field_partner: P | None) -> None:
        """Override this method to perform an action immediately after visiting a field."""

    def before_list_element(self, element: NestedField, element_partner: P | None) -> None:
        """Override this method to perform an action immediately before visiting an element within a ListType."""
        self.before_field(element, element_partner)

    def after_list_element(self, element: NestedField, element_partner: P | None) -> None:
        """Override this method to perform an action immediately after visiting an element within a ListType."""
        self.after_field(element, element_partner)

    def before_map_key(self, key: NestedField, key_partner: P | None) -> None:
        """Override this method to perform an action immediately before visiting a key within a MapType."""
        self.before_field(key, key_partner)

    def after_map_key(self, key: NestedField, key_partner: P | None) -> None:
        """Override this method to perform an action immediately after visiting a key within a MapType."""
        self.after_field(key, key_partner)

    def before_map_value(self, value: NestedField, value_partner: P | None) -> None:
        """Override this method to perform an action immediately before visiting a value within a MapType."""
        self.before_field(value, value_partner)

    def after_map_value(self, value: NestedField, value_partner: P | None) -> None:
        """Override this method to perform an action immediately after visiting a value within a MapType."""
        self.after_field(value, value_partner)

    @abstractmethod
    def schema(self, schema: Schema, schema_partner: P | None, struct_result: T) -> T:
        """Visit a schema with a partner."""

    @abstractmethod
    def struct(self, struct: StructType, struct_partner: P | None, field_results: builtins.list[T]) -> T:
        """Visit a struct type with a partner."""

    @abstractmethod
    def field(self, field: NestedField, field_partner: P | None, field_result: T) -> T:
        """Visit a nested field with a partner."""

    @abstractmethod
    def list(self, list_type: ListType, list_partner: P | None, element_result: T) -> T:
        """Visit a list type with a partner."""

    @abstractmethod
    def map(self, map_type: MapType, map_partner: P | None, key_result: T, value_result: T) -> T:
        """Visit a map type with a partner."""

    @abstractmethod
    def primitive(self, primitive: PrimitiveType, primitive_partner: P | None) -> T:
        """Visit a primitive type with a partner."""


class PrimitiveWithPartnerVisitor(SchemaWithPartnerVisitor[P, T]):
    def primitive(self, primitive: PrimitiveType, primitive_partner: P | None) -> T:
        """Visit a PrimitiveType."""
        if isinstance(primitive, BooleanType):
            return self.visit_boolean(primitive, primitive_partner)
        elif isinstance(primitive, IntegerType):
            return self.visit_integer(primitive, primitive_partner)
        elif isinstance(primitive, LongType):
            return self.visit_long(primitive, primitive_partner)
        elif isinstance(primitive, FloatType):
            return self.visit_float(primitive, primitive_partner)
        elif isinstance(primitive, DoubleType):
            return self.visit_double(primitive, primitive_partner)
        elif isinstance(primitive, DecimalType):
            return self.visit_decimal(primitive, primitive_partner)
        elif isinstance(primitive, DateType):
            return self.visit_date(primitive, primitive_partner)
        elif isinstance(primitive, TimeType):
            return self.visit_time(primitive, primitive_partner)
        elif isinstance(primitive, TimestampType):
            return self.visit_timestamp(primitive, primitive_partner)
        elif isinstance(primitive, TimestampNanoType):
            return self.visit_timestamp_ns(primitive, primitive_partner)
        elif isinstance(primitive, TimestamptzType):
            return self.visit_timestamptz(primitive, primitive_partner)
        elif isinstance(primitive, TimestamptzNanoType):
            return self.visit_timestamptz_ns(primitive, primitive_partner)
        elif isinstance(primitive, StringType):
            return self.visit_string(primitive, primitive_partner)
        elif isinstance(primitive, UUIDType):
            return self.visit_uuid(primitive, primitive_partner)
        elif isinstance(primitive, FixedType):
            return self.visit_fixed(primitive, primitive_partner)
        elif isinstance(primitive, BinaryType):
            return self.visit_binary(primitive, primitive_partner)
        elif isinstance(primitive, UnknownType):
            return self.visit_unknown(primitive, primitive_partner)
        else:
            raise ValueError(f"Type not recognized: {primitive}")

    @abstractmethod
    def visit_boolean(self, boolean_type: BooleanType, partner: P | None) -> T:
        """Visit a BooleanType."""

    @abstractmethod
    def visit_integer(self, integer_type: IntegerType, partner: P | None) -> T:
        """Visit a IntegerType."""

    @abstractmethod
    def visit_long(self, long_type: LongType, partner: P | None) -> T:
        """Visit a LongType."""

    @abstractmethod
    def visit_float(self, float_type: FloatType, partner: P | None) -> T:
        """Visit a FloatType."""

    @abstractmethod
    def visit_double(self, double_type: DoubleType, partner: P | None) -> T:
        """Visit a DoubleType."""

    @abstractmethod
    def visit_decimal(self, decimal_type: DecimalType, partner: P | None) -> T:
        """Visit a DecimalType."""

    @abstractmethod
    def visit_date(self, date_type: DateType, partner: P | None) -> T:
        """Visit a DecimalType."""

    @abstractmethod
    def visit_time(self, time_type: TimeType, partner: P | None) -> T:
        """Visit a DecimalType."""

    @abstractmethod
    def visit_timestamp(self, timestamp_type: TimestampType, partner: P | None) -> T:
        """Visit a TimestampType."""

    @abstractmethod
    def visit_timestamp_ns(self, timestamp_ns_type: TimestampNanoType, partner: P | None) -> T:
        """Visit a TimestampNanoType."""

    @abstractmethod
    def visit_timestamptz(self, timestamptz_type: TimestamptzType, partner: P | None) -> T:
        """Visit a TimestamptzType."""

    @abstractmethod
    def visit_timestamptz_ns(self, timestamptz_ns_type: TimestamptzNanoType, partner: P | None) -> T:
        """Visit a TimestamptzNanoType."""

    @abstractmethod
    def visit_string(self, string_type: StringType, partner: P | None) -> T:
        """Visit a StringType."""

    @abstractmethod
    def visit_uuid(self, uuid_type: UUIDType, partner: P | None) -> T:
        """Visit a UUIDType."""

    @abstractmethod
    def visit_fixed(self, fixed_type: FixedType, partner: P | None) -> T:
        """Visit a FixedType."""

    @abstractmethod
    def visit_binary(self, binary_type: BinaryType, partner: P | None) -> T:
        """Visit a BinaryType."""

    @abstractmethod
    def visit_unknown(self, unknown_type: UnknownType, partner: P | None) -> T:
        """Visit a UnknownType."""


class PartnerAccessor(Generic[P], ABC):
    @abstractmethod
    def schema_partner(self, partner: P | None) -> P | None:
        """Return the equivalent of the schema as a struct."""

    @abstractmethod
    def field_partner(self, partner_struct: P | None, field_id: int, field_name: str) -> P | None:
        """Return the equivalent struct field by name or id in the partner struct."""

    @abstractmethod
    def list_element_partner(self, partner_list: P | None) -> P | None:
        """Return the equivalent list element in the partner list."""

    @abstractmethod
    def map_key_partner(self, partner_map: P | None) -> P | None:
        """Return the equivalent map key in the partner map."""

    @abstractmethod
    def map_value_partner(self, partner_map: P | None) -> P | None:
        """Return the equivalent map value in the partner map."""


@singledispatch
def visit_with_partner(
    schema_or_type: Schema | IcebergType, partner: P, visitor: SchemaWithPartnerVisitor[T, P], accessor: PartnerAccessor[P]
) -> T:
    raise ValueError(f"Unsupported type: {schema_or_type}")


@visit_with_partner.register(Schema)
def _(schema: Schema, partner: P, visitor: SchemaWithPartnerVisitor[P, T], accessor: PartnerAccessor[P]) -> T:
    struct_partner = accessor.schema_partner(partner)
    return visitor.schema(schema, partner, visit_with_partner(schema.as_struct(), struct_partner, visitor, accessor))  # type: ignore


@visit_with_partner.register(StructType)
def _(struct: StructType, partner: P, visitor: SchemaWithPartnerVisitor[P, T], accessor: PartnerAccessor[P]) -> T:
    field_results = []
    for field in struct.fields:
        field_partner = accessor.field_partner(partner, field.field_id, field.name)
        visitor.before_field(field, field_partner)
        try:
            field_result = visit_with_partner(field.field_type, field_partner, visitor, accessor)  # type: ignore
            field_results.append(visitor.field(field, field_partner, field_result))
        finally:
            visitor.after_field(field, field_partner)

    return visitor.struct(struct, partner, field_results)


@visit_with_partner.register(ListType)
def _(list_type: ListType, partner: P, visitor: SchemaWithPartnerVisitor[P, T], accessor: PartnerAccessor[P]) -> T:
    element_partner = accessor.list_element_partner(partner)
    visitor.before_list_element(list_type.element_field, element_partner)
    try:
        element_result = visit_with_partner(list_type.element_type, element_partner, visitor, accessor)  # type: ignore
    finally:
        visitor.after_list_element(list_type.element_field, element_partner)

    return visitor.list(list_type, partner, element_result)


@visit_with_partner.register(MapType)
def _(map_type: MapType, partner: P, visitor: SchemaWithPartnerVisitor[P, T], accessor: PartnerAccessor[P]) -> T:
    key_partner = accessor.map_key_partner(partner)
    visitor.before_map_key(map_type.key_field, key_partner)
    try:
        key_result = visit_with_partner(map_type.key_type, key_partner, visitor, accessor)  # type: ignore
    finally:
        visitor.after_map_key(map_type.key_field, key_partner)

    value_partner = accessor.map_value_partner(partner)
    visitor.before_map_value(map_type.value_field, value_partner)
    try:
        value_result = visit_with_partner(map_type.value_type, value_partner, visitor, accessor)  # type: ignore
    finally:
        visitor.after_map_value(map_type.value_field, value_partner)
    return visitor.map(map_type, partner, key_result, value_result)


@visit_with_partner.register(PrimitiveType)
def _(primitive: PrimitiveType, partner: P, visitor: SchemaWithPartnerVisitor[P, T], _: PartnerAccessor[P]) -> T:
    return visitor.primitive(primitive, partner)


class SchemaVisitorPerPrimitiveType(SchemaVisitor[T], ABC):
    def primitive(self, primitive: PrimitiveType) -> T:
        """Visit a PrimitiveType."""
        if isinstance(primitive, FixedType):
            return self.visit_fixed(primitive)
        elif isinstance(primitive, DecimalType):
            return self.visit_decimal(primitive)
        elif isinstance(primitive, BooleanType):
            return self.visit_boolean(primitive)
        elif isinstance(primitive, IntegerType):
            return self.visit_integer(primitive)
        elif isinstance(primitive, LongType):
            return self.visit_long(primitive)
        elif isinstance(primitive, FloatType):
            return self.visit_float(primitive)
        elif isinstance(primitive, DoubleType):
            return self.visit_double(primitive)
        elif isinstance(primitive, DateType):
            return self.visit_date(primitive)
        elif isinstance(primitive, TimeType):
            return self.visit_time(primitive)
        elif isinstance(primitive, TimestampType):
            return self.visit_timestamp(primitive)
        elif isinstance(primitive, TimestampNanoType):
            return self.visit_timestamp_ns(primitive)
        elif isinstance(primitive, TimestamptzType):
            return self.visit_timestamptz(primitive)
        elif isinstance(primitive, TimestamptzNanoType):
            return self.visit_timestamptz_ns(primitive)
        elif isinstance(primitive, StringType):
            return self.visit_string(primitive)
        elif isinstance(primitive, UUIDType):
            return self.visit_uuid(primitive)
        elif isinstance(primitive, BinaryType):
            return self.visit_binary(primitive)
        elif isinstance(primitive, UnknownType):
            return self.visit_unknown(primitive)
        else:
            raise ValueError(f"Type not recognized: {primitive}")

    @abstractmethod
    def visit_fixed(self, fixed_type: FixedType) -> T:
        """Visit a FixedType."""

    @abstractmethod
    def visit_decimal(self, decimal_type: DecimalType) -> T:
        """Visit a DecimalType."""

    @abstractmethod
    def visit_boolean(self, boolean_type: BooleanType) -> T:
        """Visit a BooleanType."""

    @abstractmethod
    def visit_integer(self, integer_type: IntegerType) -> T:
        """Visit a IntegerType."""

    @abstractmethod
    def visit_long(self, long_type: LongType) -> T:
        """Visit a LongType."""

    @abstractmethod
    def visit_float(self, float_type: FloatType) -> T:
        """Visit a FloatType."""

    @abstractmethod
    def visit_double(self, double_type: DoubleType) -> T:
        """Visit a DoubleType."""

    @abstractmethod
    def visit_date(self, date_type: DateType) -> T:
        """Visit a DecimalType."""

    @abstractmethod
    def visit_time(self, time_type: TimeType) -> T:
        """Visit a DecimalType."""

    @abstractmethod
    def visit_timestamp(self, timestamp_type: TimestampType) -> T:
        """Visit a TimestampType."""

    @abstractmethod
    def visit_timestamp_ns(self, timestamp_type: TimestampNanoType) -> T:
        """Visit a TimestampNanoType."""

    @abstractmethod
    def visit_timestamptz(self, timestamptz_type: TimestamptzType) -> T:
        """Visit a TimestamptzType."""

    @abstractmethod
    def visit_timestamptz_ns(self, timestamptz_ns_type: TimestamptzNanoType) -> T:
        """Visit a TimestamptzNanoType."""

    @abstractmethod
    def visit_string(self, string_type: StringType) -> T:
        """Visit a StringType."""

    @abstractmethod
    def visit_uuid(self, uuid_type: UUIDType) -> T:
        """Visit a UUIDType."""

    @abstractmethod
    def visit_binary(self, binary_type: BinaryType) -> T:
        """Visit a BinaryType."""

    @abstractmethod
    def visit_unknown(self, unknown_type: UnknownType) -> T:
        """Visit a UnknownType."""


@dataclass(init=True, eq=True, frozen=True)
class Accessor:
    """An accessor for a specific position in a container that implements the StructProtocol."""

    position: int
    inner: Accessor | None = None

    def __str__(self) -> str:
        """Return the string representation of the Accessor class."""
        return f"Accessor(position={self.position},inner={self.inner})"

    def __repr__(self) -> str:
        """Return the string representation of the Accessor class."""
        return self.__str__()

    def get(self, container: StructProtocol) -> Any:
        """Return the value at self.position in `container`.

        Args:
            container (StructProtocol): A container to access at position `self.position`.

        Returns:
            Any: The value at position `self.position` in the container.
        """
        pos = self.position
        val = container[pos]
        inner = self
        while inner.inner:
            inner = inner.inner
            val = val[inner.position]

        return val


@singledispatch
def visit(obj: Schema | IcebergType, visitor: SchemaVisitor[T]) -> T:
    """Apply a schema visitor to any point within a schema.

    The function traverses the schema in post-order fashion.

    Args:
        obj (Union[Schema, IcebergType]): An instance of a Schema or an IcebergType.
        visitor (SchemaVisitor[T]): An instance of an implementation of the generic SchemaVisitor base class.

    Raises:
        NotImplementedError: If attempting to visit an unrecognized object type.
    """
    raise NotImplementedError(f"Cannot visit non-type: {obj}")


@visit.register(Schema)
def _(obj: Schema, visitor: SchemaVisitor[T]) -> T:
    """Visit a Schema with a concrete SchemaVisitor."""
    return visitor.schema(obj, visit(obj.as_struct(), visitor))


@visit.register(StructType)
def _(obj: StructType, visitor: SchemaVisitor[T]) -> T:
    """Visit a StructType with a concrete SchemaVisitor."""
    results = []

    for field in obj.fields:
        visitor.before_field(field)
        result = visit(field.field_type, visitor)
        visitor.after_field(field)
        results.append(visitor.field(field, result))

    return visitor.struct(obj, results)


@visit.register(ListType)
def _(obj: ListType, visitor: SchemaVisitor[T]) -> T:
    """Visit a ListType with a concrete SchemaVisitor."""
    visitor.before_list_element(obj.element_field)
    result = visit(obj.element_type, visitor)
    visitor.after_list_element(obj.element_field)

    return visitor.list(obj, result)


@visit.register(MapType)
def _(obj: MapType, visitor: SchemaVisitor[T]) -> T:
    """Visit a MapType with a concrete SchemaVisitor."""
    visitor.before_map_key(obj.key_field)
    key_result = visit(obj.key_type, visitor)
    visitor.after_map_key(obj.key_field)

    visitor.before_map_value(obj.value_field)
    value_result = visit(obj.value_type, visitor)
    visitor.after_map_value(obj.value_field)

    return visitor.map(obj, key_result, value_result)


@visit.register(PrimitiveType)
def _(obj: PrimitiveType, visitor: SchemaVisitor[T]) -> T:
    """Visit a PrimitiveType with a concrete SchemaVisitor."""
    return visitor.primitive(obj)


@singledispatch
def pre_order_visit(obj: Schema | IcebergType, visitor: PreOrderSchemaVisitor[T]) -> T:
    """Apply a schema visitor to any point within a schema.

    The function traverses the schema in pre-order fashion. This is a slimmed down version
    compared to the post-order traversal (missing before and after methods), mostly
    because we don't use the pre-order traversal much.

    Args:
        obj (Union[Schema, IcebergType]): An instance of a Schema or an IcebergType.
        visitor (PreOrderSchemaVisitor[T]): An instance of an implementation of the generic PreOrderSchemaVisitor base class.

    Raises:
        NotImplementedError: If attempting to visit an unrecognized object type.
    """
    raise NotImplementedError(f"Cannot visit non-type: {obj}")


@pre_order_visit.register(Schema)
def _(obj: Schema, visitor: PreOrderSchemaVisitor[T]) -> T:
    """Visit a Schema with a concrete PreOrderSchemaVisitor."""
    return visitor.schema(obj, lambda: pre_order_visit(obj.as_struct(), visitor))


@pre_order_visit.register(StructType)
def _(obj: StructType, visitor: PreOrderSchemaVisitor[T]) -> T:
    """Visit a StructType with a concrete PreOrderSchemaVisitor."""
    return visitor.struct(
        obj,
        [
            partial(
                lambda field: visitor.field(field, partial(lambda field: pre_order_visit(field.field_type, visitor), field)),
                field,
            )
            for field in obj.fields
        ],
    )


@pre_order_visit.register(ListType)
def _(obj: ListType, visitor: PreOrderSchemaVisitor[T]) -> T:
    """Visit a ListType with a concrete PreOrderSchemaVisitor."""
    return visitor.list(obj, lambda: pre_order_visit(obj.element_type, visitor))


@pre_order_visit.register(MapType)
def _(obj: MapType, visitor: PreOrderSchemaVisitor[T]) -> T:
    """Visit a MapType with a concrete PreOrderSchemaVisitor."""
    return visitor.map(obj, lambda: pre_order_visit(obj.key_type, visitor), lambda: pre_order_visit(obj.value_type, visitor))


@pre_order_visit.register(PrimitiveType)
def _(obj: PrimitiveType, visitor: PreOrderSchemaVisitor[T]) -> T:
    """Visit a PrimitiveType with a concrete PreOrderSchemaVisitor."""
    return visitor.primitive(obj)


class _IndexById(SchemaVisitor[dict[int, NestedField]]):
    """A schema visitor for generating a field ID to NestedField index."""

    def __init__(self) -> None:
        self._index: dict[int, NestedField] = {}

    def schema(self, schema: Schema, struct_result: dict[int, NestedField]) -> dict[int, NestedField]:
        return self._index

    def struct(self, struct: StructType, field_results: builtins.list[dict[int, NestedField]]) -> dict[int, NestedField]:
        return self._index

    def field(self, field: NestedField, field_result: dict[int, NestedField]) -> dict[int, NestedField]:
        """Add the field ID to the index."""
        self._index[field.field_id] = field
        return self._index

    def list(self, list_type: ListType, element_result: dict[int, NestedField]) -> dict[int, NestedField]:
        """Add the list element ID to the index."""
        self._index[list_type.element_field.field_id] = list_type.element_field
        return self._index

    def map(
        self, map_type: MapType, key_result: dict[int, NestedField], value_result: dict[int, NestedField]
    ) -> dict[int, NestedField]:
        """Add the key ID and value ID as individual items in the index."""
        self._index[map_type.key_field.field_id] = map_type.key_field
        self._index[map_type.value_field.field_id] = map_type.value_field
        return self._index

    def primitive(self, primitive: PrimitiveType) -> dict[int, NestedField]:
        return self._index


def index_by_id(schema_or_type: Schema | IcebergType) -> dict[int, NestedField]:
    """Generate an index of field IDs to NestedField instances.

    Args:
        schema_or_type (Union[Schema, IcebergType]): A schema or type to index.

    Returns:
        Dict[int, NestedField]: An index of field IDs to NestedField instances.
    """
    return visit(schema_or_type, _IndexById())


class _IndexParents(SchemaVisitor[dict[int, int]]):
    def __init__(self) -> None:
        self.id_to_parent: dict[int, int] = {}
        self.id_stack: list[int] = []

    def before_field(self, field: NestedField) -> None:
        self.id_stack.append(field.field_id)

    def after_field(self, field: NestedField) -> None:
        self.id_stack.pop()

    def schema(self, schema: Schema, struct_result: dict[int, int]) -> dict[int, int]:
        return self.id_to_parent

    def struct(self, struct: StructType, field_results: builtins.list[dict[int, int]]) -> dict[int, int]:
        for field in struct.fields:
            parent_id = self.id_stack[-1] if self.id_stack else None
            if parent_id is not None:
                # fields in the root struct are not added
                self.id_to_parent[field.field_id] = parent_id

        return self.id_to_parent

    def field(self, field: NestedField, field_result: dict[int, int]) -> dict[int, int]:
        return self.id_to_parent

    def list(self, list_type: ListType, element_result: dict[int, int]) -> dict[int, int]:
        self.id_to_parent[list_type.element_id] = self.id_stack[-1]
        return self.id_to_parent

    def map(self, map_type: MapType, key_result: dict[int, int], value_result: dict[int, int]) -> dict[int, int]:
        self.id_to_parent[map_type.key_id] = self.id_stack[-1]
        self.id_to_parent[map_type.value_id] = self.id_stack[-1]
        return self.id_to_parent

    def primitive(self, primitive: PrimitiveType) -> dict[int, int]:
        return self.id_to_parent


def _index_parents(schema_or_type: Schema | IcebergType) -> dict[int, int]:
    """Generate an index of field IDs to their parent field IDs.

    Args:
        schema_or_type (Union[Schema, IcebergType]): A schema or type to index.

    Returns:
        Dict[int, int]: An index of field IDs to their parent field IDs.
    """
    return visit(schema_or_type, _IndexParents())


class _IndexByName(SchemaVisitor[dict[str, int]]):
    """A schema visitor for generating a field name to field ID index."""

    def __init__(self) -> None:
        self._index: dict[str, int] = {}
        self._short_name_to_id: dict[str, int] = {}
        self._combined_index: dict[str, int] = {}
        self._field_names: list[str] = []
        self._short_field_names: list[str] = []

    def before_map_value(self, value: NestedField) -> None:
        if not isinstance(value.field_type, StructType):
            self._short_field_names.append(value.name)
        self._field_names.append(value.name)

    def after_map_value(self, value: NestedField) -> None:
        if not isinstance(value.field_type, StructType):
            self._short_field_names.pop()
        self._field_names.pop()

    def before_list_element(self, element: NestedField) -> None:
        """Short field names omit element when the element is a StructType."""
        if not isinstance(element.field_type, StructType):
            self._short_field_names.append(element.name)
        self._field_names.append(element.name)

    def after_list_element(self, element: NestedField) -> None:
        if not isinstance(element.field_type, StructType):
            self._short_field_names.pop()
        self._field_names.pop()

    def before_field(self, field: NestedField) -> None:
        """Store the field name."""
        self._field_names.append(field.name)
        self._short_field_names.append(field.name)

    def after_field(self, field: NestedField) -> None:
        """Remove the last field name stored."""
        self._field_names.pop()
        self._short_field_names.pop()

    def schema(self, schema: Schema, struct_result: dict[str, int]) -> dict[str, int]:
        return self._index

    def struct(self, struct: StructType, field_results: builtins.list[dict[str, int]]) -> dict[str, int]:
        return self._index

    def field(self, field: NestedField, field_result: dict[str, int]) -> dict[str, int]:
        """Add the field name to the index."""
        self._add_field(field.name, field.field_id)
        return self._index

    def list(self, list_type: ListType, element_result: dict[str, int]) -> dict[str, int]:
        """Add the list element name to the index."""
        self._add_field(list_type.element_field.name, list_type.element_field.field_id)
        return self._index

    def map(self, map_type: MapType, key_result: dict[str, int], value_result: dict[str, int]) -> dict[str, int]:
        """Add the key name and value name as individual items in the index."""
        self._add_field(map_type.key_field.name, map_type.key_field.field_id)
        self._add_field(map_type.value_field.name, map_type.value_field.field_id)
        return self._index

    def _add_field(self, name: str, field_id: int) -> None:
        """Add a field name to the index, mapping its full name to its field ID.

        Args:
            name (str): The field name.
            field_id (int): The field ID.

        Raises:
            ValueError: If the field name is already contained in the index.
        """
        full_name = name

        if self._field_names:
            full_name = ".".join([".".join(self._field_names), name])

        if full_name in self._index:
            raise ValueError(f"Invalid schema, multiple fields for name {full_name}: {self._index[full_name]} and {field_id}")
        self._index[full_name] = field_id

        if self._short_field_names:
            short_name = ".".join([".".join(self._short_field_names), name])
            self._short_name_to_id[short_name] = field_id

    def primitive(self, primitive: PrimitiveType) -> dict[str, int]:
        return self._index

    def by_name(self) -> dict[str, int]:
        """Return an index of combined full and short names.

        Note: Only short names that do not conflict with full names are included.
        """
        combined_index = self._short_name_to_id.copy()
        combined_index.update(self._index)
        return combined_index

    def by_id(self) -> dict[int, str]:
        """Return an index of ID to full names."""
        id_to_full_name = {value: key for key, value in self._index.items()}
        return id_to_full_name


def index_by_name(schema_or_type: Schema | IcebergType) -> dict[str, int]:
    """Generate an index of field names to field IDs.

    Args:
        schema_or_type (Union[Schema, IcebergType]): A schema or type to index.

    Returns:
        Dict[str, int]: An index of field names to field IDs.
    """
    if len(schema_or_type.fields) > 0:
        indexer = _IndexByName()
        visit(schema_or_type, indexer)
        return indexer.by_name()
    else:
        return EMPTY_DICT


def index_name_by_id(schema_or_type: Schema | IcebergType) -> dict[int, str]:
    """Generate an index of field IDs full field names.

    Args:
        schema_or_type (Union[Schema, IcebergType]): A schema or type to index.

    Returns:
        Dict[str, int]: An index of field IDs to full names.
    """
    indexer = _IndexByName()
    visit(schema_or_type, indexer)
    return indexer.by_id()


Position = int


class _BuildPositionAccessors(SchemaVisitor[dict[Position, Accessor]]):
    """A schema visitor for generating a field ID to accessor index.

    Example:
        >>> from pyiceberg.schema import Schema
        >>> from pyiceberg.types import *
        >>> schema = Schema(
        ...     NestedField(field_id=2, name="id", field_type=IntegerType(), required=False),
        ...     NestedField(field_id=1, name="data", field_type=StringType(), required=True),
        ...     NestedField(
        ...         field_id=3,
        ...         name="location",
        ...         field_type=StructType(
        ...             NestedField(field_id=5, name="latitude", field_type=FloatType(), required=False),
        ...             NestedField(field_id=6, name="longitude", field_type=FloatType(), required=False),
        ...         ),
        ...         required=True,
        ...     ),
        ...     schema_id=1,
        ...     identifier_field_ids=[1],
        ... )
        >>> result = build_position_accessors(schema)
        >>> expected = {
        ...     2: Accessor(position=0, inner=None),
        ...     1: Accessor(position=1, inner=None),
        ...     5: Accessor(position=2, inner=Accessor(position=0, inner=None)),
        ...     6: Accessor(position=2, inner=Accessor(position=1, inner=None))
        ...     3: Accessor(position=2, inner=None),
        ... }
        >>> result == expected
        True
    """

    def schema(self, schema: Schema, struct_result: dict[Position, Accessor]) -> dict[Position, Accessor]:
        return struct_result

    def struct(self, struct: StructType, field_results: builtins.list[dict[Position, Accessor]]) -> dict[Position, Accessor]:
        result = {}

        for position, field in enumerate(struct.fields):
            if field_results[position]:
                for inner_field_id, acc in field_results[position].items():
                    result[inner_field_id] = Accessor(position, inner=acc)
            result[field.field_id] = Accessor(position)

        return result

    def field(self, field: NestedField, field_result: dict[Position, Accessor]) -> dict[Position, Accessor]:
        return field_result

    def list(self, list_type: ListType, element_result: dict[Position, Accessor]) -> dict[Position, Accessor]:
        return {}

    def map(
        self, map_type: MapType, key_result: dict[Position, Accessor], value_result: dict[Position, Accessor]
    ) -> dict[Position, Accessor]:
        return {}

    def primitive(self, primitive: PrimitiveType) -> dict[Position, Accessor]:
        return {}


def build_position_accessors(schema_or_type: Schema | IcebergType) -> dict[int, Accessor]:
    """Generate an index of field IDs to schema position accessors.

    Args:
        schema_or_type (Union[Schema, IcebergType]): A schema or type to index.

    Returns:
        Dict[int, Accessor]: An index of field IDs to accessors.
    """
    return visit(schema_or_type, _BuildPositionAccessors())


def assign_fresh_schema_ids(schema_or_type: Schema | IcebergType, next_id: Callable[[], int] | None = None) -> Schema:
    """Traverses the schema, and sets new IDs."""
    return pre_order_visit(schema_or_type, _SetFreshIDs(next_id_func=next_id))


class _SetFreshIDs(PreOrderSchemaVisitor[IcebergType]):
    """Traverses the schema and assigns monotonically increasing ids."""

    old_id_to_new_id: dict[int, int]

    def __init__(self, next_id_func: Callable[[], int] | None = None) -> None:
        self.old_id_to_new_id = {}
        counter = itertools.count(1)
        self.next_id_func = next_id_func if next_id_func is not None else lambda: next(counter)

    def _get_and_increment(self, current_id: int) -> int:
        new_id = self.next_id_func()
        self.old_id_to_new_id[current_id] = new_id
        return new_id

    def schema(self, schema: Schema, struct_result: Callable[[], StructType]) -> Schema:
        return Schema(
            *struct_result().fields,
            identifier_field_ids=[self.old_id_to_new_id[field_id] for field_id in schema.identifier_field_ids],
        )

    def struct(self, struct: StructType, field_results: builtins.list[Callable[[], IcebergType]]) -> StructType:
        new_ids = [self._get_and_increment(field.field_id) for field in struct.fields]
        new_fields = []
        for field_id, field, field_type in zip(new_ids, struct.fields, field_results, strict=True):
            new_fields.append(
                NestedField(
                    field_id=field_id,
                    name=field.name,
                    field_type=field_type(),
                    required=field.required,
                    doc=field.doc,
                )
            )
        return StructType(*new_fields)

    def field(self, field: NestedField, field_result: Callable[[], IcebergType]) -> IcebergType:
        return field_result()

    def list(self, list_type: ListType, element_result: Callable[[], IcebergType]) -> ListType:
        element_id = self._get_and_increment(list_type.element_id)
        return ListType(
            element_id=element_id,
            element=element_result(),
            element_required=list_type.element_required,
        )

    def map(self, map_type: MapType, key_result: Callable[[], IcebergType], value_result: Callable[[], IcebergType]) -> MapType:
        key_id = self._get_and_increment(map_type.key_id)
        value_id = self._get_and_increment(map_type.value_id)
        return MapType(
            key_id=key_id,
            key_type=key_result(),
            value_id=value_id,
            value_type=value_result(),
            value_required=map_type.value_required,
        )

    def primitive(self, primitive: PrimitiveType) -> PrimitiveType:
        return primitive


# Implementation copied from Apache Iceberg repo.
def make_compatible_name(name: str) -> str:
    """Make a field name compatible with Avro specification.

    This function sanitizes field names to comply with Avro naming rules:
    - Names must start with [A-Za-z_]
    - Subsequent characters must be [A-Za-z0-9_]

    Invalid characters are replaced with _xHHHH where HHHH is the hex code.
    Names starting with digits get a leading underscore.

    Args:
        name: The original field name

    Returns:
        A sanitized name that complies with Avro specification
    """
    if not _valid_avro_name(name):
        return _sanitize_name(name)
    return name


def _valid_avro_name(name: str) -> bool:
    if not len(name):
        raise ValueError("Can not validate empty avro name")
    first = name[0]
    if not (first.isalpha() or first == "_"):
        return False

    for character in name[1:]:
        if not (character.isalnum() or character == "_"):
            return False
    return True


def _sanitize_name(name: str) -> str:
    sb = []
    first = name[0]
    if not (first.isalpha() or first == "_"):
        sb.append(_sanitize_char(first))
    else:
        sb.append(first)

    for character in name[1:]:
        if not (character.isalnum() or character == "_"):
            sb.append(_sanitize_char(character))
        else:
            sb.append(character)
    return "".join(sb)


def _sanitize_char(character: str) -> str:
    if character.isdigit():
        return "_" + character
    return "_x" + hex(ord(character))[2:].upper()


def sanitize_column_names(schema: Schema) -> Schema:
    """Sanitize column names to make them compatible with Avro.

    The column name should be starting with '_' or digit followed by a string only contains '_', digit or alphabet,
    otherwise it will be sanitized to conform the avro naming convention.

    Args:
        schema: The schema to be sanitized.

    Returns:
        The sanitized schema.
    """
    result = visit(schema.as_struct(), _SanitizeColumnsVisitor())
    return Schema(
        *(result or StructType()).fields,
        schema_id=schema.schema_id,
        identifier_field_ids=schema.identifier_field_ids,
    )


class _SanitizeColumnsVisitor(SchemaVisitor[IcebergType | None]):
    def schema(self, schema: Schema, struct_result: IcebergType | None) -> IcebergType | None:
        return struct_result

    def field(self, field: NestedField, field_result: IcebergType | None) -> IcebergType | None:
        return NestedField(
            field_id=field.field_id,
            name=make_compatible_name(field.name),
            field_type=field_result,
            doc=field.doc,
            required=field.required,
        )

    def struct(self, struct: StructType, field_results: builtins.list[IcebergType | None]) -> IcebergType | None:
        return StructType(*[field for field in field_results if field is not None])

    def list(self, list_type: ListType, element_result: IcebergType | None) -> IcebergType | None:
        return ListType(element_id=list_type.element_id, element_type=element_result, element_required=list_type.element_required)

    def map(self, map_type: MapType, key_result: IcebergType | None, value_result: IcebergType | None) -> IcebergType | None:
        return MapType(
            key_id=map_type.key_id,
            value_id=map_type.value_id,
            key_type=key_result,
            value_type=value_result,
            value_required=map_type.value_required,
        )

    def primitive(self, primitive: PrimitiveType) -> IcebergType | None:
        return primitive


def prune_columns(schema: Schema, selected: set[int], select_full_types: bool = True) -> Schema:
    """Prunes a column by only selecting a set of field-ids.

    Args:
        schema: The schema to be pruned.
        selected: The field-ids to be included.
        select_full_types: Return the full struct when a subset is recorded

    Returns:
        The pruned schema.
    """
    result = visit(schema.as_struct(), _PruneColumnsVisitor(selected, select_full_types))
    return Schema(
        *(result or StructType()).fields,
        schema_id=schema.schema_id,
        identifier_field_ids=list(selected.intersection(schema.identifier_field_ids)),
    )


class _PruneColumnsVisitor(SchemaVisitor[IcebergType | None]):
    selected: set[int]
    select_full_types: bool

    def __init__(self, selected: set[int], select_full_types: bool):
        self.selected = selected
        self.select_full_types = select_full_types

    def schema(self, schema: Schema, struct_result: IcebergType | None) -> IcebergType | None:
        return struct_result

    def struct(self, struct: StructType, field_results: builtins.list[IcebergType | None]) -> IcebergType | None:
        fields = struct.fields
        selected_fields = []
        same_type = True

        for idx, projected_type in enumerate(field_results):
            field = fields[idx]
            if field.field_type == projected_type:
                selected_fields.append(field)
            elif projected_type is not None:
                same_type = False
                # Type has changed, create a new field with the projected type
                selected_fields.append(
                    NestedField(
                        field_id=field.field_id,
                        name=field.name,
                        field_type=projected_type,
                        doc=field.doc,
                        required=field.required,
                    )
                )

        if selected_fields:
            if len(selected_fields) == len(fields) and same_type is True:
                # Nothing has changed, and we can return the original struct
                return struct
            else:
                return StructType(*selected_fields)
        return None

    def field(self, field: NestedField, field_result: IcebergType | None) -> IcebergType | None:
        if field.field_id in self.selected:
            if self.select_full_types:
                return field.field_type
            elif field.field_type.is_struct:
                return self._project_selected_struct(field_result)
            else:
                if not field.field_type.is_primitive:
                    raise ValueError(
                        f"Cannot explicitly project List or Map types, "
                        f"{field.field_id}:{field.name} of type {field.field_type} was selected"
                    )
                # Selected non-struct field
                return field.field_type
        elif field_result is not None:
            # This field wasn't selected but a subfield was so include that
            return field_result
        else:
            return None

    def list(self, list_type: ListType, element_result: IcebergType | None) -> IcebergType | None:
        if list_type.element_id in self.selected:
            if self.select_full_types:
                return list_type
            elif list_type.element_type and list_type.element_type.is_struct:
                projected_struct = self._project_selected_struct(element_result)
                return self._project_list(list_type, projected_struct)
            else:
                if not list_type.element_type.is_primitive:
                    raise ValueError(
                        f"Cannot explicitly project List or Map types, "
                        f"{list_type.element_id} of type {list_type.element_type} was selected"
                    )
                return list_type
        elif element_result is not None:
            return self._project_list(list_type, element_result)
        else:
            return None

    def map(self, map_type: MapType, key_result: IcebergType | None, value_result: IcebergType | None) -> IcebergType | None:
        if map_type.value_id in self.selected:
            if self.select_full_types:
                return map_type
            elif map_type.value_type and map_type.value_type.is_struct:
                projected_struct = self._project_selected_struct(value_result)
                return self._project_map(map_type, projected_struct)
            if not map_type.value_type.is_primitive:
                raise ValueError(
                    f"Cannot explicitly project List or Map types, "
                    f"Map value {map_type.value_id} of type {map_type.value_type} was selected"
                )
            return map_type
        elif value_result is not None:
            return self._project_map(map_type, value_result)
        elif map_type.key_id in self.selected:
            return map_type
        return None

    def primitive(self, primitive: PrimitiveType) -> IcebergType | None:
        return None

    @staticmethod
    def _project_selected_struct(projected_field: IcebergType | None) -> StructType:
        if projected_field and not isinstance(projected_field, StructType):
            raise ValueError("Expected a struct")

        if projected_field is None:
            return StructType()
        else:
            return projected_field

    @staticmethod
    def _project_list(list_type: ListType, element_result: IcebergType) -> ListType:
        if list_type.element_type == element_result:
            return list_type
        else:
            return ListType(
                element_id=list_type.element_id, element_type=element_result, element_required=list_type.element_required
            )

    @staticmethod
    def _project_map(map_type: MapType, value_result: IcebergType) -> MapType:
        if map_type.value_type == value_result:
            return map_type
        else:
            return MapType(
                key_id=map_type.key_id,
                value_id=map_type.value_id,
                key_type=map_type.key_type,
                value_type=value_result,
                value_required=map_type.value_required,
            )


@singledispatch
def promote(file_type: IcebergType, read_type: IcebergType) -> IcebergType:
    """Promotes reading a file type to a read type.

    Args:
        file_type (IcebergType): The type of the Avro file.
        read_type (IcebergType): The requested read type.

    Raises:
        ResolveError: If attempting to resolve an unrecognized object type.
    """
    if file_type == read_type:
        return file_type
    else:
        raise ResolveError(f"Cannot promote {file_type} to {read_type}")


@promote.register(IntegerType)
def _(file_type: IntegerType, read_type: IcebergType) -> IcebergType:
    if isinstance(read_type, LongType):
        # Ints/Longs are binary compatible in Avro, so this is okay
        return read_type
    else:
        raise ResolveError(f"Cannot promote an int to {read_type}")


@promote.register(FloatType)
def _(file_type: FloatType, read_type: IcebergType) -> IcebergType:
    if isinstance(read_type, DoubleType):
        # A double type is wider
        return read_type
    else:
        raise ResolveError(f"Cannot promote an float to {read_type}")


@promote.register(StringType)
def _(file_type: StringType, read_type: IcebergType) -> IcebergType:
    if isinstance(read_type, BinaryType):
        return read_type
    else:
        raise ResolveError(f"Cannot promote an string to {read_type}")


@promote.register(BinaryType)
def _(file_type: BinaryType, read_type: IcebergType) -> IcebergType:
    if isinstance(read_type, StringType):
        return read_type
    else:
        raise ResolveError(f"Cannot promote an binary to {read_type}")


@promote.register(DecimalType)
def _(file_type: DecimalType, read_type: IcebergType) -> IcebergType:
    if isinstance(read_type, DecimalType):
        if file_type.precision <= read_type.precision and file_type.scale == file_type.scale:
            return read_type
        else:
            raise ResolveError(f"Cannot reduce precision from {file_type} to {read_type}")
    else:
        raise ResolveError(f"Cannot promote an decimal to {read_type}")


@promote.register(FixedType)
def _(file_type: FixedType, read_type: IcebergType) -> IcebergType:
    if isinstance(read_type, UUIDType) and len(file_type) == 16:
        # Since pyarrow reads parquet UUID as fixed 16-byte binary, the promotion is needed to ensure read compatibility
        return read_type
    else:
        raise ResolveError(f"Cannot promote {file_type} to {read_type}")


@promote.register(UnknownType)
def _(file_type: UnknownType, read_type: IcebergType) -> IcebergType:
    # Per V3 Spec, "Unknown" can be promoted to any Primitive type
    if isinstance(read_type, PrimitiveType):
        return read_type
    else:
        raise ResolveError(f"Cannot promote {file_type} to {read_type}")


def _check_schema_compatible(requested_schema: Schema, provided_schema: Schema) -> None:
    """
    Check if the `provided_schema` is compatible with `requested_schema`.

    Both Schemas must have valid IDs and share the same ID for the same field names.

    Two schemas are considered compatible when:
    1. All `required` fields in `requested_schema` are present and are also `required` in the `provided_schema`
    2. Field Types are consistent for fields that are present in both schemas. I.e. the field type
       in the `provided_schema` can be promoted to the field type of the same field ID in `requested_schema`

    Raises:
        ValueError: If the schemas are not compatible.
    """
    pre_order_visit(requested_schema, _SchemaCompatibilityVisitor(provided_schema))


class _SchemaCompatibilityVisitor(PreOrderSchemaVisitor[bool]):
    provided_schema: Schema

    def __init__(self, provided_schema: Schema):
        from rich.console import Console
        from rich.table import Table as RichTable

        self.provided_schema = provided_schema
        self.rich_table = RichTable(show_header=True, header_style="bold")
        self.rich_table.add_column("")
        self.rich_table.add_column("Table field")
        self.rich_table.add_column("Dataframe field")
        self.console = Console(record=True)

    def _is_field_compatible(self, lhs: NestedField) -> bool:
        # Validate nullability first.
        # An optional field can be missing in the provided schema
        # But a required field must exist as a required field
        try:
            rhs = self.provided_schema.find_field(lhs.field_id)
        except ValueError:
            if lhs.required:
                self.rich_table.add_row("❌", str(lhs), "Missing")
                return False
            else:
                self.rich_table.add_row("✅", str(lhs), "Missing")
                return True

        if lhs.required and not rhs.required:
            self.rich_table.add_row("❌", str(lhs), str(rhs))
            return False

        # Check type compatibility
        if lhs.field_type == rhs.field_type:
            self.rich_table.add_row("✅", str(lhs), str(rhs))
            return True
        # We only check that the parent node is also of the same type.
        # We check the type of the child nodes when we traverse them later.
        elif any(
            (isinstance(lhs.field_type, container_type) and isinstance(rhs.field_type, container_type))
            for container_type in {StructType, MapType, ListType}
        ):
            self.rich_table.add_row("✅", str(lhs), str(rhs))
            return True
        else:
            try:
                # If type can be promoted to the requested schema
                # it is considered compatible
                promote(rhs.field_type, lhs.field_type)
                self.rich_table.add_row("✅", str(lhs), str(rhs))
                return True
            except ResolveError:
                # UnknownType can only be promoted to Primitive types
                if isinstance(rhs.field_type, UnknownType):
                    if not isinstance(lhs.field_type, PrimitiveType):
                        error_msg = (
                            f"Null type (UnknownType) cannot be promoted to non-primitive type {lhs.field_type}. "
                            "UnknownType can only be promoted to primitive types (string, int, boolean, etc.) "
                            "in V3+ tables."
                        )
                    else:
                        error_msg = (
                            f"Null type (UnknownType) cannot be promoted to {lhs.field_type}. "
                            "This may be due to table format version limitations "
                            "(V1/V2 tables don't support UnknownType promotion)."
                        )
                    self.rich_table.add_row("❌", str(lhs), f"{str(rhs)} - {error_msg}")
                else:
                    self.rich_table.add_row("❌", str(lhs), str(rhs))
                return False

    def schema(self, schema: Schema, struct_result: Callable[[], bool]) -> bool:
        if not (result := struct_result()):
            self.console.print(self.rich_table)
            raise ValueError(f"Mismatch in fields:\n{self.console.export_text()}")
        return result

    def struct(self, struct: StructType, field_results: builtins.list[Callable[[], bool]]) -> bool:
        results = [result() for result in field_results]
        return all(results)

    def field(self, field: NestedField, field_result: Callable[[], bool]) -> bool:
        # Skip child validation for missing optional fields (#2797)
        is_compatible = self._is_field_compatible(field)
        if field.field_id not in self.provided_schema._lazy_id_to_field:
            return is_compatible
        return is_compatible and field_result()

    def list(self, list_type: ListType, element_result: Callable[[], bool]) -> bool:
        return self._is_field_compatible(list_type.element_field) and element_result()

    def map(self, map_type: MapType, key_result: Callable[[], bool], value_result: Callable[[], bool]) -> bool:
        return all(
            [
                self._is_field_compatible(map_type.key_field),
                self._is_field_compatible(map_type.value_field),
                key_result(),
                value_result(),
            ]
        )

    def primitive(self, primitive: PrimitiveType) -> bool:
        return True
