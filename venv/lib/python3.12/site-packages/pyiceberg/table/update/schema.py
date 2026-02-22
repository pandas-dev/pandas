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

import builtins
import itertools
from copy import copy
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Any

from pyiceberg.exceptions import ResolveError, ValidationError
from pyiceberg.expressions import literal  # type: ignore
from pyiceberg.schema import (
    PartnerAccessor,
    Schema,
    SchemaVisitor,
    SchemaWithPartnerVisitor,
    assign_fresh_schema_ids,
    promote,
    visit,
    visit_with_partner,
)
from pyiceberg.table.name_mapping import (
    NameMapping,
    update_mapping,
)
from pyiceberg.table.update import (
    AddSchemaUpdate,
    AssertCurrentSchemaId,
    SetCurrentSchemaUpdate,
    SetPropertiesUpdate,
    TableRequirement,
    TableUpdate,
    UpdatesAndRequirements,
    UpdateTableMetadata,
)
from pyiceberg.typedef import L, TableVersion
from pyiceberg.types import IcebergType, ListType, MapType, NestedField, PrimitiveType, StructType

if TYPE_CHECKING:
    import pyarrow as pa

    from pyiceberg.table import Transaction

TABLE_ROOT_ID = -1


class _MoveOperation(Enum):
    First = 1
    Before = 2
    After = 3


@dataclass
class _Move:
    field_id: int
    full_name: str
    op: _MoveOperation
    other_field_id: int | None = None


class UpdateSchema(UpdateTableMetadata["UpdateSchema"]):
    _schema: Schema
    _last_column_id: itertools.count[int]
    _identifier_field_names: set[str]

    _adds: dict[int, list[NestedField]] = {}
    _updates: dict[int, NestedField] = {}
    _deletes: set[int] = set()
    _moves: dict[int, list[_Move]] = {}

    _added_name_to_id: dict[str, int] = {}
    # Part of https://github.com/apache/iceberg/pull/8393
    _id_to_parent: dict[int, str] = {}
    _allow_incompatible_changes: bool
    _case_sensitive: bool

    def __init__(
        self,
        transaction: Transaction,
        allow_incompatible_changes: bool = False,
        case_sensitive: bool = True,
        schema: Schema | None = None,
        name_mapping: NameMapping | None = None,
    ) -> None:
        super().__init__(transaction)

        if isinstance(schema, Schema):
            self._schema = schema
            self._last_column_id = itertools.count(1 + schema.highest_field_id)
        else:
            self._schema = self._transaction.table_metadata.schema()
            self._last_column_id = itertools.count(1 + self._transaction.table_metadata.last_column_id)

        self._name_mapping = name_mapping
        self._identifier_field_names = self._schema.identifier_field_names()

        self._adds = {}
        self._updates = {}
        self._deletes = set()
        self._moves = {}

        self._added_name_to_id = {}

        def get_column_name(field_id: int) -> str:
            column_name = self._schema.find_column_name(column_id=field_id)
            if column_name is None:
                raise ValueError(f"Could not find field-id: {field_id}")
            return column_name

        self._id_to_parent = {
            field_id: get_column_name(parent_field_id) for field_id, parent_field_id in self._schema._lazy_id_to_parent.items()
        }

        self._allow_incompatible_changes = allow_incompatible_changes
        self._case_sensitive = case_sensitive
        self._transaction = transaction

    def case_sensitive(self, case_sensitive: bool) -> UpdateSchema:
        """Determine if the case of schema needs to be considered when comparing column names.

        Args:
            case_sensitive: When false case is not considered in column name comparisons.

        Returns:
            This for method chaining
        """
        self._case_sensitive = case_sensitive
        return self

    def union_by_name(
        # TODO: Move TableProperties.DEFAULT_FORMAT_VERSION to separate file and set that as format_version default.
        self,
        new_schema: Schema | pa.Schema,
        format_version: TableVersion = 2,
    ) -> UpdateSchema:
        from pyiceberg.catalog import Catalog

        visit_with_partner(
            Catalog._convert_schema_if_needed(new_schema, format_version=format_version),
            -1,
            _UnionByNameVisitor(update_schema=self, existing_schema=self._schema, case_sensitive=self._case_sensitive),
            # type: ignore
            PartnerIdByNameAccessor(partner_schema=self._schema, case_sensitive=self._case_sensitive),
        )
        return self

    def add_column(
        self,
        path: str | tuple[str, ...],
        field_type: IcebergType,
        doc: str | None = None,
        required: bool = False,
        default_value: L | None = None,
    ) -> UpdateSchema:
        """Add a new column to a nested struct or Add a new top-level column.

        Because "." may be interpreted as a column path separator or may be used in field names, it
        is not allowed to add nested column by passing in a string. To add to nested structures or
        to add fields with names that contain "." use a tuple instead to indicate the path.

        If type is a nested type, its field IDs are reassigned when added to the existing schema.

        Args:
            path: Name for the new column.
            field_type: Type for the new column.
            doc: Documentation string for the new column.
            required: Whether the new column is required.
            default_value: Default value for the new column.

        Returns:
            This for method chaining.
        """
        if isinstance(path, str):
            if "." in path:
                raise ValueError(f"Cannot add column with ambiguous name: {path}, provide a tuple instead")
            path = (path,)

        name = path[-1]
        parent = path[:-1]

        full_name = ".".join(path)
        parent_full_path = ".".join(parent)
        parent_id: int = TABLE_ROOT_ID

        if len(parent) > 0:
            parent_field = self._schema.find_field(parent_full_path, self._case_sensitive)
            parent_type = parent_field.field_type
            if isinstance(parent_type, MapType):
                parent_field = parent_type.value_field
            elif isinstance(parent_type, ListType):
                parent_field = parent_type.element_field

            if not parent_field.field_type.is_struct:
                raise ValueError(f"Cannot add column '{name}' to non-struct type: {parent_full_path}")

            parent_id = parent_field.field_id

        existing_field = None
        try:
            existing_field = self._schema.find_field(full_name, self._case_sensitive)
        except ValueError:
            pass

        if existing_field is not None and existing_field.field_id not in self._deletes:
            raise ValueError(f"Cannot add column, name already exists: {full_name}")

        # assign new IDs in order
        new_id = self.assign_new_column_id()
        new_type = assign_fresh_schema_ids(field_type, self.assign_new_column_id)

        if default_value is not None:
            try:
                # To make sure that the value is valid for the type
                initial_default = literal(default_value).to(new_type).value
            except ValueError as e:
                raise ValueError(f"Invalid default value: {e}") from e
        else:
            initial_default = default_value  # type: ignore

        if (required and initial_default is None) and not self._allow_incompatible_changes:
            # Table format version 1 and 2 cannot add required column because there is no initial value
            raise ValueError(f"Incompatible change: cannot add required column: {'.'.join(path)}")

        # update tracking for moves
        self._added_name_to_id[full_name] = new_id
        self._id_to_parent[new_id] = parent_full_path

        field = NestedField(
            field_id=new_id,
            name=name,
            field_type=new_type,
            required=required,
            doc=doc,
            initial_default=initial_default,
            write_default=initial_default,
        )

        if parent_id in self._adds:
            self._adds[parent_id].append(field)
        else:
            self._adds[parent_id] = [field]

        return self

    def delete_column(self, path: str | tuple[str, ...]) -> UpdateSchema:
        """Delete a column from a table.

        Args:
            path: The path to the column.

        Returns:
            The UpdateSchema with the delete operation staged.
        """
        name = (path,) if isinstance(path, str) else path
        full_name = ".".join(name)

        field = self._schema.find_field(full_name, case_sensitive=self._case_sensitive)

        if field.field_id in self._adds:
            raise ValueError(f"Cannot delete a column that has additions: {full_name}")
        if field.field_id in self._updates:
            raise ValueError(f"Cannot delete a column that has updates: {full_name}")

        self._deletes.add(field.field_id)

        return self

    def set_default_value(self, path: str | tuple[str, ...], default_value: L | None) -> UpdateSchema:
        """Set the default value of a column.

        Args:
            path: The path to the column.

        Returns:
            The UpdateSchema with the delete operation staged.
        """
        self._set_column_default_value(path, default_value)

        return self

    def rename_column(self, path_from: str | tuple[str, ...], new_name: str) -> UpdateSchema:
        """Update the name of a column.

        Args:
            path_from: The path to the column to be renamed.
            new_name: The new path of the column.

        Returns:
            The UpdateSchema with the rename operation staged.
        """
        path_from = ".".join(path_from) if isinstance(path_from, tuple) else path_from
        field_from = self._schema.find_field(path_from, self._case_sensitive)

        if field_from.field_id in self._deletes:
            raise ValueError(f"Cannot rename a column that will be deleted: {path_from}")

        if updated := self._updates.get(field_from.field_id):
            self._updates[field_from.field_id] = NestedField(
                field_id=updated.field_id,
                name=new_name,
                field_type=updated.field_type,
                doc=updated.doc,
                required=updated.required,
                initial_default=updated.initial_default,
                write_default=updated.write_default,
            )
        else:
            self._updates[field_from.field_id] = NestedField(
                field_id=field_from.field_id,
                name=new_name,
                field_type=field_from.field_type,
                doc=field_from.doc,
                required=field_from.required,
                initial_default=field_from.initial_default,
                write_default=field_from.write_default,
            )

        # Lookup the field because of casing
        from_field_correct_casing = self._schema.find_column_name(field_from.field_id)
        if from_field_correct_casing in self._identifier_field_names:
            self._identifier_field_names.remove(from_field_correct_casing)
            new_identifier_path = f"{from_field_correct_casing[: -len(field_from.name)]}{new_name}"
            self._identifier_field_names.add(new_identifier_path)

        return self

    def make_column_optional(self, path: str | tuple[str, ...]) -> UpdateSchema:
        """Make a column optional.

        Args:
            path: The path to the field.

        Returns:
            The UpdateSchema with the requirement change staged.
        """
        self._set_column_requirement(path, required=False)
        return self

    def set_identifier_fields(self, *fields: str) -> None:
        self._identifier_field_names = set(fields)

    def _set_column_requirement(self, path: str | tuple[str, ...], required: bool) -> None:
        path = (path,) if isinstance(path, str) else path
        name = ".".join(path)

        field = self._schema.find_field(name, self._case_sensitive)

        if (field.required and required) or (field.optional and not required):
            # if the change is a noop, allow it even if allowIncompatibleChanges is false
            return

        if not self._allow_incompatible_changes and required:
            raise ValueError(f"Cannot change column nullability: {name}: optional -> required")

        if field.field_id in self._deletes:
            raise ValueError(f"Cannot update a column that will be deleted: {name}")

        if updated := self._updates.get(field.field_id):
            self._updates[field.field_id] = NestedField(
                field_id=updated.field_id,
                name=updated.name,
                field_type=updated.field_type,
                doc=updated.doc,
                required=required,
                initial_default=updated.initial_default,
                write_default=updated.write_default,
            )
        else:
            self._updates[field.field_id] = NestedField(
                field_id=field.field_id,
                name=field.name,
                field_type=field.field_type,
                doc=field.doc,
                required=required,
                initial_default=field.initial_default,
                write_default=field.write_default,
            )

    def _set_column_default_value(self, path: str | tuple[str, ...], default_value: Any) -> None:
        path = (path,) if isinstance(path, str) else path
        name = ".".join(path)

        field = self._schema.find_field(name, self._case_sensitive)

        if default_value is not None:
            try:
                # To make sure that the value is valid for the type
                default_value = literal(default_value).to(field.field_type).value
            except ValueError as e:
                raise ValueError(f"Invalid default value: {e}") from e

        if field.required and default_value == field.write_default:
            # if the change is a noop, allow it even if allowIncompatibleChanges is false
            return

        if not self._allow_incompatible_changes and field.required and default_value is None:
            raise ValueError("Cannot change change default-value of a required column to None")

        if field.field_id in self._deletes:
            raise ValueError(f"Cannot update a column that will be deleted: {name}")

        if updated := self._updates.get(field.field_id):
            self._updates[field.field_id] = NestedField(
                field_id=updated.field_id,
                name=updated.name,
                field_type=updated.field_type,
                doc=updated.doc,
                required=updated.required,
                initial_default=updated.initial_default,
                write_default=default_value,
            )
        else:
            self._updates[field.field_id] = NestedField(
                field_id=field.field_id,
                name=field.name,
                field_type=field.field_type,
                doc=field.doc,
                required=field.required,
                initial_default=field.initial_default,
                write_default=default_value,
            )

    def update_column(
        self,
        path: str | tuple[str, ...],
        field_type: IcebergType | None = None,
        required: bool | None = None,
        doc: str | None = None,
    ) -> UpdateSchema:
        """Update the type of column.

        Args:
            path: The path to the field.
            field_type: The new type
            required: If the field should be required
            doc: Documentation describing the column

        Returns:
            The UpdateSchema with the type update staged.
        """
        path = (path,) if isinstance(path, str) else path
        full_name = ".".join(path)

        if field_type is None and required is None and doc is None:
            return self

        field = self._schema.find_field(full_name, self._case_sensitive)

        if field.field_id in self._deletes:
            raise ValueError(f"Cannot update a column that will be deleted: {full_name}")

        if field_type is not None:
            if not field.field_type.is_primitive:
                raise ValidationError(f"Cannot change column type: {field.field_type} is not a primitive")

            if not self._allow_incompatible_changes and field.field_type != field_type:
                try:
                    promote(field.field_type, field_type)
                except ResolveError as e:
                    raise ValidationError(f"Cannot change column type: {full_name}: {field.field_type} -> {field_type}") from e

        # if other updates for the same field exist in one transaction:
        if updated := self._updates.get(field.field_id):
            self._updates[field.field_id] = NestedField(
                field_id=updated.field_id,
                name=updated.name,
                field_type=field_type or updated.field_type,
                doc=doc if doc is not None else updated.doc,
                required=updated.required,
                initial_default=updated.initial_default,
                write_default=updated.write_default,
            )
        else:
            self._updates[field.field_id] = NestedField(
                field_id=field.field_id,
                name=field.name,
                field_type=field_type or field.field_type,
                doc=doc if doc is not None else field.doc,
                required=field.required,
                initial_default=field.initial_default,
                write_default=field.write_default,
            )

        if required is not None:
            self._set_column_requirement(path, required=required)

        return self

    def _find_for_move(self, name: str) -> int | None:
        try:
            return self._schema.find_field(name, self._case_sensitive).field_id
        except ValueError:
            pass

        return self._added_name_to_id.get(name)

    def _move(self, move: _Move) -> None:
        if parent_name := self._id_to_parent.get(move.field_id):
            parent_field = self._schema.find_field(parent_name, case_sensitive=self._case_sensitive)
            if not parent_field.field_type.is_struct:
                raise ValueError(f"Cannot move fields in non-struct type: {parent_field.field_type}")

            if move.op == _MoveOperation.After or move.op == _MoveOperation.Before:
                if move.other_field_id is None:
                    raise ValueError("Expected other field when performing before/after move")

                if self._id_to_parent.get(move.field_id) != self._id_to_parent.get(move.other_field_id):
                    raise ValueError(f"Cannot move field {move.full_name} to a different struct")

            self._moves[parent_field.field_id] = self._moves.get(parent_field.field_id, []) + [move]
        else:
            # In the top level field
            if move.op == _MoveOperation.After or move.op == _MoveOperation.Before:
                if move.other_field_id is None:
                    raise ValueError("Expected other field when performing before/after move")

                if other_struct := self._id_to_parent.get(move.other_field_id):
                    raise ValueError(f"Cannot move field {move.full_name} to a different struct: {other_struct}")

            self._moves[TABLE_ROOT_ID] = self._moves.get(TABLE_ROOT_ID, []) + [move]

    def move_first(self, path: str | tuple[str, ...]) -> UpdateSchema:
        """Move the field to the first position of the parent struct.

        Args:
            path: The path to the field.

        Returns:
            The UpdateSchema with the move operation staged.
        """
        full_name = ".".join(path) if isinstance(path, tuple) else path

        field_id = self._find_for_move(full_name)

        if field_id is None:
            raise ValueError(f"Cannot move missing column: {full_name}")

        self._move(_Move(field_id=field_id, full_name=full_name, op=_MoveOperation.First))

        return self

    def move_before(self, path: str | tuple[str, ...], before_path: str | tuple[str, ...]) -> UpdateSchema:
        """Move the field to before another field.

        Args:
            path: The path to the field.

        Returns:
            The UpdateSchema with the move operation staged.
        """
        full_name = ".".join(path) if isinstance(path, tuple) else path
        field_id = self._find_for_move(full_name)

        if field_id is None:
            raise ValueError(f"Cannot move missing column: {full_name}")

        before_full_name = (
            ".".join(
                before_path,
            )
            if isinstance(before_path, tuple)
            else before_path
        )
        before_field_id = self._find_for_move(before_full_name)

        if before_field_id is None:
            raise ValueError(f"Cannot move {full_name} before missing column: {before_full_name}")

        if field_id == before_field_id:
            raise ValueError(f"Cannot move {full_name} before itself")

        self._move(_Move(field_id=field_id, full_name=full_name, other_field_id=before_field_id, op=_MoveOperation.Before))

        return self

    def move_after(self, path: str | tuple[str, ...], after_name: str | tuple[str, ...]) -> UpdateSchema:
        """Move the field to after another field.

        Args:
            path: The path to the field.

        Returns:
            The UpdateSchema with the move operation staged.
        """
        full_name = ".".join(path) if isinstance(path, tuple) else path

        field_id = self._find_for_move(full_name)

        if field_id is None:
            raise ValueError(f"Cannot move missing column: {full_name}")

        after_path = ".".join(after_name) if isinstance(after_name, tuple) else after_name
        after_field_id = self._find_for_move(after_path)

        if after_field_id is None:
            raise ValueError(f"Cannot move {full_name} after missing column: {after_path}")

        if field_id == after_field_id:
            raise ValueError(f"Cannot move {full_name} after itself")

        self._move(_Move(field_id=field_id, full_name=full_name, other_field_id=after_field_id, op=_MoveOperation.After))

        return self

    def _commit(self) -> UpdatesAndRequirements:
        """Apply the pending changes and commit."""
        from pyiceberg.table import TableProperties

        new_schema = self._apply()

        existing_schema_id = next(
            (schema.schema_id for schema in self._transaction.table_metadata.schemas if schema == new_schema), None
        )

        requirements: tuple[TableRequirement, ...] = ()
        updates: tuple[TableUpdate, ...] = ()

        # Check if it is different current schema ID
        if existing_schema_id != self._schema.schema_id:
            requirements += (AssertCurrentSchemaId(current_schema_id=self._schema.schema_id),)
            if existing_schema_id is None:
                updates += (
                    AddSchemaUpdate(schema=new_schema),
                    SetCurrentSchemaUpdate(schema_id=-1),
                )
            else:
                updates += (SetCurrentSchemaUpdate(schema_id=existing_schema_id),)

            if name_mapping := self._name_mapping:
                updated_name_mapping = update_mapping(name_mapping, self._updates, self._adds)
                updates += (
                    SetPropertiesUpdate(updates={TableProperties.DEFAULT_NAME_MAPPING: updated_name_mapping.model_dump_json()}),
                )

        return updates, requirements

    def _apply(self) -> Schema:
        """Apply the pending changes to the original schema and returns the result.

        Returns:
            the result Schema when all pending updates are applied
        """
        struct = visit(self._schema, _ApplyChanges(self._adds, self._updates, self._deletes, self._moves))
        if struct is None:
            # Should never happen
            raise ValueError("Could not apply changes")

        # Check the field-ids
        new_schema = Schema(*struct.fields)
        from pyiceberg.partitioning import validate_partition_name

        for spec in self._transaction.table_metadata.partition_specs:
            for partition_field in spec.fields:
                validate_partition_name(
                    partition_field.name, partition_field.transform, partition_field.source_id, new_schema, set()
                )
        field_ids = set()
        for name in self._identifier_field_names:
            try:
                field = new_schema.find_field(name, case_sensitive=self._case_sensitive)
            except ValueError as e:
                raise ValueError(
                    f"Cannot find identifier field {name}. In case of deletion, update the identifier fields first."
                ) from e

            field_ids.add(field.field_id)

        if txn := self._transaction:
            next_schema_id = 1 + (
                max(schema.schema_id for schema in txn.table_metadata.schemas) if txn.table_metadata is not None else 0
            )
        else:
            next_schema_id = 0

        return Schema(*struct.fields, schema_id=next_schema_id, identifier_field_ids=field_ids)

    def assign_new_column_id(self) -> int:
        return next(self._last_column_id)


class _ApplyChanges(SchemaVisitor[IcebergType | None]):
    _adds: dict[int, builtins.list[NestedField]]
    _updates: dict[int, NestedField]
    _deletes: set[int]
    _moves: dict[int, builtins.list[_Move]]

    def __init__(
        self,
        adds: dict[int, builtins.list[NestedField]],
        updates: dict[int, NestedField],
        deletes: set[int],
        moves: dict[int, builtins.list[_Move]],
    ) -> None:
        self._adds = adds
        self._updates = updates
        self._deletes = deletes
        self._moves = moves

    def schema(self, schema: Schema, struct_result: IcebergType | None) -> IcebergType | None:
        added = self._adds.get(TABLE_ROOT_ID)
        moves = self._moves.get(TABLE_ROOT_ID)

        if added is not None or moves is not None:
            if not isinstance(struct_result, StructType):
                raise ValueError(f"Cannot add fields to non-struct: {struct_result}")

            if new_fields := _add_and_move_fields(struct_result.fields, added or [], moves or []):
                return StructType(*new_fields)

        return struct_result

    def struct(self, struct: StructType, field_results: builtins.list[IcebergType | None]) -> IcebergType | None:
        has_changes = False
        new_fields = []

        for idx, result_type in enumerate(field_results):
            result_type = field_results[idx]

            # Has been deleted
            if result_type is None:
                has_changes = True
                continue

            field = struct.fields[idx]

            name = field.name
            doc = field.doc
            required = field.required
            write_default = field.write_default

            # There is an update
            if update := self._updates.get(field.field_id):
                name = update.name
                doc = update.doc
                required = update.required
                write_default = update.write_default

            if (
                field.name == name
                and field.field_type == result_type
                and field.required == required
                and field.doc == doc
                and field.write_default == write_default
            ):
                new_fields.append(field)
            else:
                has_changes = True
                new_fields.append(
                    NestedField(
                        field_id=field.field_id,
                        name=name,
                        field_type=result_type,
                        required=required,
                        doc=doc,
                        initial_default=field.initial_default,
                        write_default=write_default,
                    )
                )

        if has_changes:
            return StructType(*new_fields)

        return struct

    def field(self, field: NestedField, field_result: IcebergType | None) -> IcebergType | None:
        # the API validates deletes, updates, and additions don't conflict handle deletes
        if field.field_id in self._deletes:
            return None

        # handle updates
        if (update := self._updates.get(field.field_id)) and field.field_type != update.field_type:
            return update.field_type

        if isinstance(field_result, StructType):
            # handle add & moves
            added = self._adds.get(field.field_id)
            moves = self._moves.get(field.field_id)
            if added is not None or moves is not None:
                if not isinstance(field.field_type, StructType):
                    raise ValueError(f"Cannot add fields to non-struct: {field}")

                if new_fields := _add_and_move_fields(field_result.fields, added or [], moves or []):
                    return StructType(*new_fields)

        return field_result

    def list(self, list_type: ListType, element_result: IcebergType | None) -> IcebergType | None:
        element_type = self.field(list_type.element_field, element_result)
        if element_type is None:
            raise ValueError(f"Cannot delete element type from list: {element_result}")

        element_required = list_type.element_required
        if update := self._updates.get(list_type.element_id):
            element_required = update.required

        return ListType(element_id=list_type.element_id, element=element_type, element_required=element_required)

    def map(self, map_type: MapType, key_result: IcebergType | None, value_result: IcebergType | None) -> IcebergType | None:
        key_id: int = map_type.key_field.field_id

        if key_id in self._deletes:
            raise ValueError(f"Cannot delete map keys: {map_type}")

        if key_id in self._updates:
            raise ValueError(f"Cannot update map keys: {map_type}")

        if key_id in self._adds:
            raise ValueError(f"Cannot add fields to map keys: {map_type}")

        if map_type.key_type != key_result:
            raise ValueError(f"Cannot alter map keys: {map_type}")

        value_field: NestedField = map_type.value_field
        value_type = self.field(value_field, value_result)
        if value_type is None:
            raise ValueError(f"Cannot delete value type from map: {value_field}")

        value_required = map_type.value_required
        if update := self._updates.get(map_type.value_id):
            value_required = update.required

        return MapType(
            key_id=map_type.key_id,
            key_type=map_type.key_type,
            value_id=map_type.value_id,
            value_type=value_type,
            value_required=value_required,
        )

    def primitive(self, primitive: PrimitiveType) -> IcebergType | None:
        return primitive


class _UnionByNameVisitor(SchemaWithPartnerVisitor[int, bool]):
    update_schema: UpdateSchema
    existing_schema: Schema
    case_sensitive: bool

    def __init__(self, update_schema: UpdateSchema, existing_schema: Schema, case_sensitive: bool) -> None:
        self.update_schema = update_schema
        self.existing_schema = existing_schema
        self.case_sensitive = case_sensitive

    def schema(self, schema: Schema, partner_id: int | None, struct_result: bool) -> bool:
        return struct_result

    def struct(self, struct: StructType, partner_id: int | None, missing_positions: builtins.list[bool]) -> bool:
        if partner_id is None:
            return True

        fields = struct.fields
        partner_struct = self._find_field_type(partner_id)

        if not partner_struct.is_struct:
            raise ValueError(f"Expected a struct, got: {partner_struct}")

        for pos, missing in enumerate(missing_positions):
            if missing:
                self._add_column(partner_id, fields[pos])
            else:
                field = fields[pos]
                if nested_field := partner_struct.field_by_name(field.name, case_sensitive=self.case_sensitive):
                    self._update_column(field, nested_field)

        return False

    def _add_column(self, parent_id: int, field: NestedField) -> None:
        if parent_name := self.existing_schema.find_column_name(parent_id):
            path: tuple[str, ...] = (parent_name, field.name)
        else:
            path = (field.name,)

        self.update_schema.add_column(path=path, field_type=field.field_type, required=field.required, doc=field.doc)

    def _update_column(self, field: NestedField, existing_field: NestedField) -> None:
        full_name = self.existing_schema.find_column_name(existing_field.field_id)

        if full_name is None:
            raise ValueError(f"Could not find field: {existing_field}")

        if field.optional and existing_field.required:
            self.update_schema.make_column_optional(full_name)

        if field.field_type.is_primitive and field.field_type != existing_field.field_type:
            try:
                # If the current type is wider than the new type, then
                # we perform a noop
                _ = promote(field.field_type, existing_field.field_type)
            except ResolveError:
                # If this is not the case, perform the type evolution
                self.update_schema.update_column(full_name, field_type=field.field_type)

        if field.doc is not None and field.doc != existing_field.doc:
            self.update_schema.update_column(full_name, doc=field.doc)

    def _find_field_type(self, field_id: int) -> IcebergType:
        if field_id == -1:
            return self.existing_schema.as_struct()
        else:
            return self.existing_schema.find_field(field_id).field_type

    def field(self, field: NestedField, partner_id: int | None, field_result: bool) -> bool:
        return partner_id is None

    def list(self, list_type: ListType, list_partner_id: int | None, element_missing: bool) -> bool:
        if list_partner_id is None:
            return True

        if element_missing:
            raise ValueError("Error traversing schemas: element is missing, but list is present")

        partner_list_type = self._find_field_type(list_partner_id)
        if not isinstance(partner_list_type, ListType):
            raise ValueError(f"Expected list-type, got: {partner_list_type}")

        self._update_column(list_type.element_field, partner_list_type.element_field)

        return False

    def map(self, map_type: MapType, map_partner_id: int | None, key_missing: bool, value_missing: bool) -> bool:
        if map_partner_id is None:
            return True

        if key_missing:
            raise ValueError("Error traversing schemas: key is missing, but map is present")

        if value_missing:
            raise ValueError("Error traversing schemas: value is missing, but map is present")

        partner_map_type = self._find_field_type(map_partner_id)
        if not isinstance(partner_map_type, MapType):
            raise ValueError(f"Expected map-type, got: {partner_map_type}")

        self._update_column(map_type.key_field, partner_map_type.key_field)
        self._update_column(map_type.value_field, partner_map_type.value_field)

        return False

    def primitive(self, primitive: PrimitiveType, primitive_partner_id: int | None) -> bool:
        return primitive_partner_id is None


class PartnerIdByNameAccessor(PartnerAccessor[int]):
    partner_schema: Schema
    case_sensitive: bool

    def __init__(self, partner_schema: Schema, case_sensitive: bool) -> None:
        self.partner_schema = partner_schema
        self.case_sensitive = case_sensitive

    def schema_partner(self, partner: int | None) -> int | None:
        return -1

    def field_partner(self, partner_field_id: int | None, field_id: int, field_name: str) -> int | None:
        if partner_field_id is not None:
            if partner_field_id == -1:
                struct = self.partner_schema.as_struct()
            else:
                struct = self.partner_schema.find_field(partner_field_id).field_type
                if not struct.is_struct:
                    raise ValueError(f"Expected StructType: {struct}")

            if field := struct.field_by_name(name=field_name, case_sensitive=self.case_sensitive):
                return field.field_id

        return None

    def list_element_partner(self, partner_list_id: int | None) -> int | None:
        if partner_list_id is not None and (field := self.partner_schema.find_field(partner_list_id)):
            if not isinstance(field.field_type, ListType):
                raise ValueError(f"Expected ListType: {field}")
            return field.field_type.element_field.field_id
        else:
            return None

    def map_key_partner(self, partner_map_id: int | None) -> int | None:
        if partner_map_id is not None and (field := self.partner_schema.find_field(partner_map_id)):
            if not isinstance(field.field_type, MapType):
                raise ValueError(f"Expected MapType: {field}")
            return field.field_type.key_field.field_id
        else:
            return None

    def map_value_partner(self, partner_map_id: int | None) -> int | None:
        if partner_map_id is not None and (field := self.partner_schema.find_field(partner_map_id)):
            if not isinstance(field.field_type, MapType):
                raise ValueError(f"Expected MapType: {field}")
            return field.field_type.value_field.field_id
        else:
            return None


def _add_fields(fields: tuple[NestedField, ...], adds: list[NestedField] | None) -> tuple[NestedField, ...]:
    adds = adds or []
    return fields + tuple(adds)


def _move_fields(fields: tuple[NestedField, ...], moves: list[_Move]) -> tuple[NestedField, ...]:
    reordered = list(copy(fields))
    for move in moves:
        # Find the field that we're about to move
        field = next(field for field in reordered if field.field_id == move.field_id)
        # Remove the field that we're about to move from the list
        reordered = [field for field in reordered if field.field_id != move.field_id]

        if move.op == _MoveOperation.First:
            reordered = [field] + reordered
        elif move.op == _MoveOperation.Before or move.op == _MoveOperation.After:
            other_field_id = move.other_field_id
            other_field_pos = next(i for i, field in enumerate(reordered) if field.field_id == other_field_id)
            if move.op == _MoveOperation.Before:
                reordered.insert(other_field_pos, field)
            else:
                reordered.insert(other_field_pos + 1, field)
        else:
            raise ValueError(f"Unknown operation: {move.op}")

    return tuple(reordered)


def _add_and_move_fields(
    fields: tuple[NestedField, ...], adds: list[NestedField], moves: list[_Move]
) -> tuple[NestedField, ...] | None:
    if len(adds) > 0:
        # always apply adds first so that added fields can be moved
        added = _add_fields(fields, adds)
        if len(moves) > 0:
            return _move_fields(added, moves)
        else:
            return added
    elif len(moves) > 0:
        return _move_fields(fields, moves)
    return None if len(adds) == 0 else tuple(*fields, *adds)
