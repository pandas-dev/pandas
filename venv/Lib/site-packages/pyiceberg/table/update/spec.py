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

from typing import TYPE_CHECKING, Any

from pyiceberg.expressions import (
    Reference,
)
from pyiceberg.partitioning import (
    INITIAL_PARTITION_SPEC_ID,
    PARTITION_FIELD_ID_START,
    PartitionField,
    PartitionSpec,
    _PartitionNameGenerator,
    _visit_partition_field,
)
from pyiceberg.schema import Schema
from pyiceberg.table.update import (
    AddPartitionSpecUpdate,
    AssertDefaultSpecId,
    AssertLastAssignedPartitionId,
    SetDefaultSpecUpdate,
    TableRequirement,
    TableUpdate,
    UpdatesAndRequirements,
    UpdateTableMetadata,
)
from pyiceberg.transforms import IdentityTransform, TimeTransform, Transform, VoidTransform, parse_transform

if TYPE_CHECKING:
    from pyiceberg.table import Transaction


class UpdateSpec(UpdateTableMetadata["UpdateSpec"]):
    _transaction: Transaction
    _name_to_field: dict[str, PartitionField] = {}
    _name_to_added_field: dict[str, PartitionField] = {}
    _transform_to_field: dict[tuple[int, str], PartitionField] = {}
    _transform_to_added_field: dict[tuple[int, str], PartitionField] = {}
    _renames: dict[str, str] = {}
    _added_time_fields: dict[int, PartitionField] = {}
    _case_sensitive: bool
    _adds: list[PartitionField]
    _deletes: set[int]
    _last_assigned_partition_id: int

    def __init__(self, transaction: Transaction, case_sensitive: bool = True) -> None:
        super().__init__(transaction)
        self._name_to_field = {field.name: field for field in transaction.table_metadata.spec().fields}
        self._name_to_added_field = {}
        self._transform_to_field = {
            (field.source_id, repr(field.transform)): field for field in transaction.table_metadata.spec().fields
        }
        self._transform_to_added_field = {}
        self._adds = []
        self._deletes = set()
        self._last_assigned_partition_id = transaction.table_metadata.last_partition_id or PARTITION_FIELD_ID_START - 1
        self._renames = {}
        self._transaction = transaction
        self._case_sensitive = case_sensitive
        self._added_time_fields = {}

    def add_field(
        self,
        source_column_name: str,
        transform: str | Transform[Any, Any],
        partition_field_name: str | None = None,
    ) -> UpdateSpec:
        ref = Reference(source_column_name)
        bound_ref = ref.bind(self._transaction.table_metadata.schema(), self._case_sensitive)
        if isinstance(transform, str):
            transform = parse_transform(transform)
        # verify transform can actually bind it
        output_type = bound_ref.field.field_type
        if not transform.can_transform(output_type):
            raise ValueError(f"{transform} cannot transform {output_type} values from {bound_ref.field.name}")

        transform_key = (bound_ref.field.field_id, repr(transform))
        existing_partition_field = self._transform_to_field.get(transform_key)
        if existing_partition_field and self._is_duplicate_partition(transform, existing_partition_field):
            raise ValueError(f"Duplicate partition field for ${ref.name}=${ref}, ${existing_partition_field} already exists")

        added = self._transform_to_added_field.get(transform_key)
        if added:
            raise ValueError(f"Already added partition: {added.name}")

        new_field = self._partition_field((bound_ref.field.field_id, transform), partition_field_name)
        if new_field.name in self._name_to_added_field:
            raise ValueError(f"Already added partition field with name: {new_field.name}")

        if isinstance(new_field.transform, TimeTransform):
            existing_time_field = self._added_time_fields.get(new_field.source_id)
            if existing_time_field:
                raise ValueError(f"Cannot add time partition field: {new_field.name} conflicts with {existing_time_field.name}")
            self._added_time_fields[new_field.source_id] = new_field
        self._transform_to_added_field[transform_key] = new_field

        existing_partition_field = self._name_to_field.get(new_field.name)
        if existing_partition_field and new_field.field_id not in self._deletes:
            if isinstance(existing_partition_field.transform, VoidTransform):
                self.rename_field(
                    existing_partition_field.name, existing_partition_field.name + "_" + str(existing_partition_field.field_id)
                )
            else:
                raise ValueError(f"Cannot add duplicate partition field name: {existing_partition_field.name}")

        self._name_to_added_field[new_field.name] = new_field
        self._adds.append(new_field)
        return self

    def add_identity(self, source_column_name: str) -> UpdateSpec:
        return self.add_field(source_column_name, IdentityTransform(), None)

    def remove_field(self, name: str) -> UpdateSpec:
        added = self._name_to_added_field.get(name)
        if added:
            raise ValueError(f"Cannot delete newly added field {name}")
        renamed = self._renames.get(name)
        if renamed:
            raise ValueError(f"Cannot rename and delete field {name}")
        field = self._name_to_field.get(name)
        if not field:
            raise ValueError(f"No such partition field: {name}")

        self._deletes.add(field.field_id)
        return self

    def rename_field(self, name: str, new_name: str) -> UpdateSpec:
        existing_field = self._name_to_field.get(new_name)
        if existing_field and isinstance(existing_field.transform, VoidTransform):
            return self.rename_field(name, name + "_" + str(existing_field.field_id))
        added = self._name_to_added_field.get(name)
        if added:
            raise ValueError("Cannot rename recently added partitions")
        field = self._name_to_field.get(name)
        if not field:
            raise ValueError(f"Cannot find partition field {name}")
        if field.field_id in self._deletes:
            raise ValueError(f"Cannot delete and rename partition field {name}")
        self._renames[name] = new_name
        return self

    def _commit(self) -> UpdatesAndRequirements:
        new_spec = self._apply()
        updates: tuple[TableUpdate, ...] = ()
        requirements: tuple[TableRequirement, ...] = ()

        if self._transaction.table_metadata.default_spec_id != new_spec.spec_id:
            if new_spec.spec_id not in self._transaction.table_metadata.specs():
                updates = (
                    AddPartitionSpecUpdate(spec=new_spec),
                    SetDefaultSpecUpdate(spec_id=-1),
                )
            else:
                updates = (SetDefaultSpecUpdate(spec_id=new_spec.spec_id),)

            required_last_assigned_partitioned_id = self._transaction.table_metadata.last_partition_id
            default_spec_id = self._transaction.table_metadata.default_spec_id
            requirements = (
                AssertLastAssignedPartitionId(last_assigned_partition_id=required_last_assigned_partitioned_id),
                AssertDefaultSpecId(default_spec_id=default_spec_id),
            )

        return updates, requirements

    def _apply(self) -> PartitionSpec:
        def _check_and_add_partition_name(
            schema: Schema, name: str, source_id: int, transform: Transform[Any, Any], partition_names: set[str]
        ) -> None:
            from pyiceberg.partitioning import validate_partition_name

            validate_partition_name(name, transform, source_id, schema, partition_names)
            partition_names.add(name)

        def _add_new_field(
            schema: Schema, source_id: int, field_id: int, name: str, transform: Transform[Any, Any], partition_names: set[str]
        ) -> PartitionField:
            _check_and_add_partition_name(schema, name, source_id, transform, partition_names)
            return PartitionField(source_id, field_id, transform, name)

        partition_fields = []
        partition_names: set[str] = set()
        for field in self._transaction.table_metadata.spec().fields:
            if field.field_id not in self._deletes:
                renamed = self._renames.get(field.name)
                if renamed:
                    new_field = _add_new_field(
                        self._transaction.table_metadata.schema(),
                        field.source_id,
                        field.field_id,
                        renamed,
                        field.transform,
                        partition_names,
                    )
                else:
                    new_field = _add_new_field(
                        self._transaction.table_metadata.schema(),
                        field.source_id,
                        field.field_id,
                        field.name,
                        field.transform,
                        partition_names,
                    )
                partition_fields.append(new_field)
            elif self._transaction.table_metadata.format_version == 1:
                renamed = self._renames.get(field.name)
                if renamed:
                    new_field = _add_new_field(
                        self._transaction.table_metadata.schema(),
                        field.source_id,
                        field.field_id,
                        renamed,
                        VoidTransform(),
                        partition_names,
                    )
                else:
                    new_field = _add_new_field(
                        self._transaction.table_metadata.schema(),
                        field.source_id,
                        field.field_id,
                        field.name,
                        VoidTransform(),
                        partition_names,
                    )

                partition_fields.append(new_field)

        for added_field in self._adds:
            _check_and_add_partition_name(
                self._transaction.table_metadata.schema(),
                added_field.name,
                added_field.source_id,
                added_field.transform,
                partition_names,
            )
            new_field = PartitionField(
                source_id=added_field.source_id,
                field_id=added_field.field_id,
                transform=added_field.transform,
                name=added_field.name,
            )
            partition_fields.append(new_field)

        # Reuse spec id or create a new one.
        new_spec = PartitionSpec(*partition_fields)
        new_spec_id = INITIAL_PARTITION_SPEC_ID
        for spec in self._transaction.table_metadata.specs().values():
            if new_spec.compatible_with(spec):
                new_spec_id = spec.spec_id
                break
            elif new_spec_id <= spec.spec_id:
                new_spec_id = spec.spec_id + 1
        return PartitionSpec(*partition_fields, spec_id=new_spec_id)

    def _partition_field(self, transform_key: tuple[int, Transform[Any, Any]], name: str | None) -> PartitionField:
        if self._transaction.table_metadata.format_version == 2:
            source_id, transform = transform_key
            historical_fields = []
            for spec in self._transaction.table_metadata.specs().values():
                for field in spec.fields:
                    historical_fields.append(field)

            for field in historical_fields:
                if field.source_id == source_id and repr(field.transform) == repr(transform):
                    if name is None or field.name == name:
                        return PartitionField(source_id, field.field_id, transform, field.name)

        new_field_id = self._new_field_id()
        if name is None:
            tmp_field = PartitionField(transform_key[0], new_field_id, transform_key[1], "unassigned_field_name")
            name = _visit_partition_field(self._transaction.table_metadata.schema(), tmp_field, _PartitionNameGenerator())
        return PartitionField(transform_key[0], new_field_id, transform_key[1], name)

    def _new_field_id(self) -> int:
        self._last_assigned_partition_id += 1
        return self._last_assigned_partition_id

    def _is_duplicate_partition(self, transform: Transform[Any, Any], partition_field: PartitionField) -> bool:
        return partition_field.field_id not in self._deletes and partition_field.transform == transform
