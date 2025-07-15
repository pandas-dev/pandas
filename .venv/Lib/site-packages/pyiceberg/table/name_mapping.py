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
"""
Contains everything around the name mapping.

More information can be found on here:
https://iceberg.apache.org/spec/#name-mapping-serialization
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import ChainMap
from functools import cached_property, singledispatch
from typing import Any, Dict, Generic, Iterator, List, Optional, TypeVar, Union

from pydantic import Field, conlist, field_validator, model_serializer

from pyiceberg.schema import P, PartnerAccessor, Schema, SchemaVisitor, SchemaWithPartnerVisitor, visit, visit_with_partner
from pyiceberg.typedef import IcebergBaseModel, IcebergRootModel
from pyiceberg.types import IcebergType, ListType, MapType, NestedField, PrimitiveType, StructType


class MappedField(IcebergBaseModel):
    field_id: Optional[int] = Field(alias="field-id", default=None)
    names: List[str] = conlist(str)
    fields: List[MappedField] = Field(default_factory=list)

    @field_validator("fields", mode="before")
    @classmethod
    def convert_null_to_empty_List(cls, v: Any) -> Any:
        return v or []

    @model_serializer
    def ser_model(self) -> Dict[str, Any]:
        """Set custom serializer to leave out the field when it is empty."""
        serialized: Dict[str, Any] = {"names": self.names}
        if self.field_id is not None:
            serialized["field-id"] = self.field_id
        if len(self.fields) > 0:
            serialized["fields"] = self.fields
        return serialized

    def __len__(self) -> int:
        """Return the number of fields."""
        return len(self.fields)

    def __str__(self) -> str:
        """Convert the mapped-field into a nicely formatted string."""
        # Otherwise the UTs fail because the order of the set can change
        fields_str = ", ".join([str(e) for e in self.fields]) or ""
        fields_str = " " + fields_str if fields_str else ""
        field_id = "?" if self.field_id is None else (str(self.field_id) or "?")
        return "([" + ", ".join(self.names) + "] -> " + field_id + fields_str + ")"


class NameMapping(IcebergRootModel[List[MappedField]]):
    root: List[MappedField]

    @cached_property
    def _field_by_name(self) -> Dict[str, MappedField]:
        return visit_name_mapping(self, _IndexByName())

    def __len__(self) -> int:
        """Return the number of mappings."""
        return len(self.root)

    def __iter__(self) -> Iterator[MappedField]:
        """Iterate over the mapped fields."""
        return iter(self.root)

    def __str__(self) -> str:
        """Convert the name-mapping into a nicely formatted string."""
        if len(self.root) == 0:
            return "[]"
        else:
            return "[\n  " + "\n  ".join([str(e) for e in self.root]) + "\n]"


S = TypeVar("S")
T = TypeVar("T")


class NameMappingVisitor(Generic[S, T], ABC):
    @abstractmethod
    def mapping(self, nm: NameMapping, field_results: S) -> S:
        """Visit a NameMapping."""

    @abstractmethod
    def fields(self, struct: List[MappedField], field_results: List[T]) -> S:
        """Visit a List[MappedField]."""

    @abstractmethod
    def field(self, field: MappedField, field_result: S) -> T:
        """Visit a MappedField."""


class _IndexByName(NameMappingVisitor[Dict[str, MappedField], Dict[str, MappedField]]):
    def mapping(self, nm: NameMapping, field_results: Dict[str, MappedField]) -> Dict[str, MappedField]:
        return field_results

    def fields(self, struct: List[MappedField], field_results: List[Dict[str, MappedField]]) -> Dict[str, MappedField]:
        return dict(ChainMap(*field_results))

    def field(self, field: MappedField, field_result: Dict[str, MappedField]) -> Dict[str, MappedField]:
        result: Dict[str, MappedField] = {
            f"{field_name}.{key}": result_field for key, result_field in field_result.items() for field_name in field.names
        }

        for name in field.names:
            result[name] = field

        return result


@singledispatch
def visit_name_mapping(obj: Union[NameMapping, List[MappedField], MappedField], visitor: NameMappingVisitor[S, T]) -> S:
    """Traverse the name mapping in post-order traversal."""
    raise NotImplementedError(f"Cannot visit non-type: {obj}")


@visit_name_mapping.register(NameMapping)
def _(obj: NameMapping, visitor: NameMappingVisitor[S, T]) -> S:
    return visitor.mapping(obj, visit_name_mapping(obj.root, visitor))


@visit_name_mapping.register(list)
def _(fields: List[MappedField], visitor: NameMappingVisitor[S, T]) -> S:
    results = [visitor.field(field, visit_name_mapping(field.fields, visitor)) for field in fields]
    return visitor.fields(fields, results)


def parse_mapping_from_json(mapping: str) -> NameMapping:
    return NameMapping.model_validate_json(mapping)


class _CreateMapping(SchemaVisitor[List[MappedField]]):
    def schema(self, schema: Schema, struct_result: List[MappedField]) -> List[MappedField]:
        return struct_result

    def struct(self, struct: StructType, field_results: List[List[MappedField]]) -> List[MappedField]:
        return [
            MappedField(field_id=field.field_id, names=[field.name], fields=result)
            for field, result in zip(struct.fields, field_results)
        ]

    def field(self, field: NestedField, field_result: List[MappedField]) -> List[MappedField]:
        return field_result

    def list(self, list_type: ListType, element_result: List[MappedField]) -> List[MappedField]:
        return [MappedField(field_id=list_type.element_id, names=["element"], fields=element_result)]

    def map(self, map_type: MapType, key_result: List[MappedField], value_result: List[MappedField]) -> List[MappedField]:
        return [
            MappedField(field_id=map_type.key_id, names=["key"], fields=key_result),
            MappedField(field_id=map_type.value_id, names=["value"], fields=value_result),
        ]

    def primitive(self, primitive: PrimitiveType) -> List[MappedField]:
        return []


class _UpdateMapping(NameMappingVisitor[List[MappedField], MappedField]):
    _updates: Dict[int, NestedField]
    _adds: Dict[int, List[NestedField]]

    def __init__(self, updates: Dict[int, NestedField], adds: Dict[int, List[NestedField]]):
        self._updates = updates
        self._adds = adds

    @staticmethod
    def _remove_reassigned_names(field: MappedField, assignments: Dict[str, int]) -> Optional[MappedField]:
        removed_names = set()
        for name in field.names:
            if (assigned_id := assignments.get(name)) and assigned_id != field.field_id:
                removed_names.add(name)

        remaining_names = [f for f in field.names if f not in removed_names]
        if remaining_names:
            return MappedField(field_id=field.field_id, names=remaining_names, fields=field.fields)
        else:
            return None

    def _add_new_fields(self, mapped_fields: List[MappedField], parent_id: int) -> List[MappedField]:
        if fields_to_add := self._adds.get(parent_id):
            fields: List[MappedField] = []
            new_fields: List[MappedField] = []

            for add in fields_to_add:
                new_fields.append(
                    MappedField(field_id=add.field_id, names=[add.name], fields=visit(add.field_type, _CreateMapping()))
                )

            reassignments = {f.name: f.field_id for f in fields_to_add}
            fields = [
                updated_field
                for field in mapped_fields
                if (updated_field := self._remove_reassigned_names(field, reassignments)) is not None
            ] + new_fields
            return fields
        else:
            return mapped_fields

    def mapping(self, nm: NameMapping, field_results: List[MappedField]) -> List[MappedField]:
        return self._add_new_fields(field_results, -1)

    def fields(self, struct: List[MappedField], field_results: List[MappedField]) -> List[MappedField]:
        reassignments: Dict[str, int] = {
            update.name: update.field_id
            for f in field_results
            if f.field_id is not None and (update := self._updates.get(f.field_id))
        }
        return [
            updated_field
            for field in field_results
            if (updated_field := self._remove_reassigned_names(field, reassignments)) is not None
        ]

    def field(self, field: MappedField, field_result: List[MappedField]) -> MappedField:
        if field.field_id is None:
            return field
        field_names = field.names
        if (update := self._updates.get(field.field_id)) is not None and update.name not in field_names:
            field_names.append(update.name)

        return MappedField(field_id=field.field_id, names=field_names, fields=self._add_new_fields(field_result, field.field_id))


def create_mapping_from_schema(schema: Schema) -> NameMapping:
    return NameMapping(visit(schema, _CreateMapping()))


def update_mapping(mapping: NameMapping, updates: Dict[int, NestedField], adds: Dict[int, List[NestedField]]) -> NameMapping:
    return NameMapping(visit_name_mapping(mapping, _UpdateMapping(updates, adds)))


class NameMappingAccessor(PartnerAccessor[MappedField]):
    def schema_partner(self, partner: Optional[MappedField]) -> Optional[MappedField]:
        return partner

    def field_partner(
        self, partner_struct: Optional[Union[List[MappedField], MappedField]], _: int, field_name: str
    ) -> Optional[MappedField]:
        if partner_struct is not None:
            if isinstance(partner_struct, MappedField):
                partner_struct = partner_struct.fields

            for field in partner_struct:
                if field_name in field.names:
                    return field

        return None

    def list_element_partner(self, partner_list: Optional[MappedField]) -> Optional[MappedField]:
        if partner_list is not None:
            for field in partner_list.fields:
                if "element" in field.names:
                    return field
        return None

    def map_key_partner(self, partner_map: Optional[MappedField]) -> Optional[MappedField]:
        if partner_map is not None:
            for field in partner_map.fields:
                if "key" in field.names:
                    return field
        return None

    def map_value_partner(self, partner_map: Optional[MappedField]) -> Optional[MappedField]:
        if partner_map is not None:
            for field in partner_map.fields:
                if "value" in field.names:
                    return field
        return None


class NameMappingProjectionVisitor(SchemaWithPartnerVisitor[MappedField, IcebergType]):
    current_path: List[str]

    def __init__(self) -> None:
        # For keeping track where we are in case when a field cannot be found
        self.current_path = []

    def before_field(self, field: NestedField, field_partner: Optional[P]) -> None:
        self.current_path.append(field.name)

    def after_field(self, field: NestedField, field_partner: Optional[P]) -> None:
        self.current_path.pop()

    def before_list_element(self, element: NestedField, element_partner: Optional[P]) -> None:
        self.current_path.append("element")

    def after_list_element(self, element: NestedField, element_partner: Optional[P]) -> None:
        self.current_path.pop()

    def before_map_key(self, key: NestedField, key_partner: Optional[P]) -> None:
        self.current_path.append("key")

    def after_map_key(self, key: NestedField, key_partner: Optional[P]) -> None:
        self.current_path.pop()

    def before_map_value(self, value: NestedField, value_partner: Optional[P]) -> None:
        self.current_path.append("value")

    def after_map_value(self, value: NestedField, value_partner: Optional[P]) -> None:
        self.current_path.pop()

    def schema(self, schema: Schema, schema_partner: Optional[MappedField], struct_result: StructType) -> IcebergType:
        return Schema(*struct_result.fields, schema_id=schema.schema_id)

    def struct(self, struct: StructType, struct_partner: Optional[MappedField], field_results: List[NestedField]) -> IcebergType:
        return StructType(*field_results)

    def field(self, field: NestedField, field_partner: Optional[MappedField], field_result: IcebergType) -> IcebergType:
        if field_partner is None or field_partner.field_id is None:
            raise ValueError(f"Field or field ID missing from NameMapping: {'.'.join(self.current_path)}")

        return NestedField(
            field_id=field_partner.field_id,
            name=field.name,
            field_type=field_result,
            required=field.required,
            doc=field.doc,
            initial_default=field.initial_default,
            initial_write=field.write_default,
        )

    def list(self, list_type: ListType, list_partner: Optional[MappedField], element_result: IcebergType) -> IcebergType:
        if list_partner is None:
            raise ValueError(f"Could not find field with name: {'.'.join(self.current_path)}")

        element_id = next(field for field in list_partner.fields if "element" in field.names).field_id
        return ListType(element_id=element_id, element=element_result, element_required=list_type.element_required)

    def map(
        self, map_type: MapType, map_partner: Optional[MappedField], key_result: IcebergType, value_result: IcebergType
    ) -> IcebergType:
        if map_partner is None:
            raise ValueError(f"Could not find field with name: {'.'.join(self.current_path)}")

        key_id = next(field for field in map_partner.fields if "key" in field.names).field_id
        value_id = next(field for field in map_partner.fields if "value" in field.names).field_id
        return MapType(
            key_id=key_id,
            key_type=key_result,
            value_id=value_id,
            value_type=value_result,
            value_required=map_type.value_required,
        )

    def primitive(self, primitive: PrimitiveType, primitive_partner: Optional[MappedField]) -> PrimitiveType:
        if primitive_partner is None:
            raise ValueError(f"Could not find field with name: {'.'.join(self.current_path)}")

        return primitive


def apply_name_mapping(schema_without_ids: Schema, name_mapping: NameMapping) -> Schema:
    return visit_with_partner(schema_without_ids, name_mapping, NameMappingProjectionVisitor(), NameMappingAccessor())  # type: ignore
