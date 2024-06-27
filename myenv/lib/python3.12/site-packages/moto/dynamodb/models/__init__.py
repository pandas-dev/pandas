import copy
import re
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.exceptions import JsonRESTError
from moto.core.utils import unix_time
from moto.dynamodb.comparisons import get_expected, get_filter_expression
from moto.dynamodb.exceptions import (
    BackupNotFoundException,
    ConditionalCheckFailed,
    DeletionProtectedException,
    ItemSizeTooLarge,
    ItemSizeToUpdateTooLarge,
    MockValidationException,
    MultipleTransactionsException,
    ResourceInUseException,
    ResourceNotFoundException,
    SourceTableNotFoundException,
    StreamAlreadyEnabledException,
    TableAlreadyExistsException,
    TableNotFoundException,
    TooManyTransactionsException,
    TransactionCanceledException,
    TransactWriteSingleOpException,
)
from moto.dynamodb.models.dynamo_type import DynamoType, Item
from moto.dynamodb.models.table import (
    Backup,
    GlobalSecondaryIndex,
    RestoredPITTable,
    RestoredTable,
    Table,
)
from moto.dynamodb.models.table_import import TableImport
from moto.dynamodb.parsing import partiql
from moto.dynamodb.parsing.executors import UpdateExpressionExecutor
from moto.dynamodb.parsing.expressions import UpdateExpressionParser  # type: ignore
from moto.dynamodb.parsing.validators import UpdateExpressionValidator


class DynamoDBBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.tables: Dict[str, Table] = OrderedDict()
        self.backups: Dict[str, Backup] = OrderedDict()
        self.table_imports: Dict[str, TableImport] = {}

    @staticmethod
    def default_vpc_endpoint_service(
        service_region: str, zones: List[str]
    ) -> List[Dict[str, str]]:
        """Default VPC endpoint service."""
        # No 'vpce' in the base endpoint DNS name
        return BaseBackend.default_vpc_endpoint_service_factory(
            service_region,
            zones,
            "dynamodb",
            "Gateway",
            private_dns_names=False,
            base_endpoint_dns_names=[f"dynamodb.{service_region}.amazonaws.com"],
        )

    def create_table(
        self,
        name: str,
        schema: List[Dict[str, str]],
        throughput: Optional[Dict[str, int]],
        attr: List[Dict[str, str]],
        global_indexes: Optional[List[Dict[str, Any]]],
        indexes: Optional[List[Dict[str, Any]]],
        streams: Optional[Dict[str, Any]],
        billing_mode: str,
        sse_specification: Optional[Dict[str, Any]],
        tags: List[Dict[str, str]],
        deletion_protection_enabled: bool,
    ) -> Table:
        if name in self.tables:
            raise ResourceInUseException(f"Table already exists: {name}")
        table = Table(
            name,
            account_id=self.account_id,
            region=self.region_name,
            schema=schema,
            throughput=throughput,
            attr=attr,
            global_indexes=global_indexes,
            indexes=indexes,
            streams=streams,
            billing_mode=billing_mode,
            sse_specification=sse_specification,
            tags=tags,
            deletion_protection_enabled=deletion_protection_enabled,
        )
        self.tables[name] = table
        return table

    def delete_table(self, name: str) -> Table:
        if name not in self.tables:
            raise ResourceNotFoundException
        table_for_deletion = self.tables.get(name)
        if isinstance(table_for_deletion, Table):
            if table_for_deletion.deletion_protection_enabled:
                raise DeletionProtectedException(name)
        return self.tables.pop(name)

    def describe_endpoints(self) -> List[Dict[str, Union[int, str]]]:
        return [
            {
                "Address": f"dynamodb.{self.region_name}.amazonaws.com",
                "CachePeriodInMinutes": 1440,
            }
        ]

    def tag_resource(self, table_arn: str, tags: List[Dict[str, str]]) -> None:
        for table in self.tables:
            if self.tables[table].table_arn == table_arn:
                self.tables[table].tags.extend(tags)

    def untag_resource(self, table_arn: str, tag_keys: List[str]) -> None:
        for table in self.tables:
            if self.tables[table].table_arn == table_arn:
                self.tables[table].tags = [
                    tag for tag in self.tables[table].tags if tag["Key"] not in tag_keys
                ]

    def list_tags_of_resource(self, table_arn: str) -> List[Dict[str, str]]:
        for table in self.tables:
            if self.tables[table].table_arn == table_arn:
                return self.tables[table].tags
        raise ResourceNotFoundException

    def list_tables(
        self, limit: int, exclusive_start_table_name: str
    ) -> Tuple[List[str], Optional[str]]:
        all_tables: List[str] = list(self.tables.keys())

        if exclusive_start_table_name:
            try:
                last_table_index = all_tables.index(exclusive_start_table_name)
            except ValueError:
                start = len(all_tables)
            else:
                start = last_table_index + 1
        else:
            start = 0

        if limit:
            tables = all_tables[start : start + limit]
        else:
            tables = all_tables[start:]

        if limit and len(all_tables) > start + limit:
            return tables, tables[-1]
        return tables, None

    def describe_table(self, name: str) -> Dict[str, Any]:
        # We can't use get_table() here, because the error message is slightly different for this operation
        if name not in self.tables:
            raise ResourceNotFoundException(table_name=name)
        return self.tables[name].describe(base_key="Table")

    def update_table(
        self,
        name: str,
        attr_definitions: List[Dict[str, str]],
        global_index: List[Dict[str, Any]],
        throughput: Dict[str, Any],
        billing_mode: str,
        stream_spec: Dict[str, Any],
        deletion_protection_enabled: bool,
    ) -> Table:
        table = self.get_table(name)
        if attr_definitions:
            table.attr = attr_definitions
        if global_index:
            table = self.update_table_global_indexes(name, global_index)
        if throughput:
            table = self.update_table_throughput(name, throughput)
        if billing_mode:
            table = self.update_table_billing_mode(name, billing_mode)
        if stream_spec:
            table = self.update_table_streams(name, stream_spec)
        if deletion_protection_enabled:
            table = self.update_table_deletion_protection_enabled(
                name, deletion_protection_enabled
            )
        return table

    def update_table_throughput(self, name: str, throughput: Dict[str, int]) -> Table:
        table = self.tables[name]
        table.throughput = throughput
        return table

    def update_table_billing_mode(self, name: str, billing_mode: str) -> Table:
        table = self.tables[name]
        table.billing_mode = billing_mode
        return table

    def update_table_deletion_protection_enabled(
        self, name: str, deletion_protection_enabled: bool
    ) -> Table:
        table = self.tables[name]
        table.deletion_protection_enabled = deletion_protection_enabled
        return table

    def update_table_streams(
        self, name: str, stream_specification: Dict[str, Any]
    ) -> Table:
        table = self.tables[name]
        if (
            stream_specification.get("StreamEnabled")
            or stream_specification.get("StreamViewType")
        ) and table.latest_stream_label:
            raise StreamAlreadyEnabledException
        table.set_stream_specification(stream_specification)
        return table

    def update_table_global_indexes(
        self, name: str, global_index_updates: List[Dict[str, Any]]
    ) -> Table:
        table = self.tables[name]
        gsis_by_name = dict((i.name, i) for i in table.global_indexes)
        for gsi_update in global_index_updates:
            gsi_to_create = gsi_update.get("Create")
            gsi_to_update = gsi_update.get("Update")
            gsi_to_delete = gsi_update.get("Delete")

            if gsi_to_delete:
                index_name = gsi_to_delete["IndexName"]
                if index_name not in gsis_by_name:
                    raise ValueError(
                        "Global Secondary Index does not exist, but tried to delete: %s"
                        % gsi_to_delete["IndexName"]
                    )

                del gsis_by_name[index_name]

            if gsi_to_update:
                index_name = gsi_to_update["IndexName"]
                if index_name not in gsis_by_name:
                    raise ValueError(
                        "Global Secondary Index does not exist, but tried to update: %s"
                        % index_name
                    )
                gsis_by_name[index_name].update(gsi_to_update)

            if gsi_to_create:
                if gsi_to_create["IndexName"] in gsis_by_name:
                    raise ValueError(
                        "Global Secondary Index already exists: %s"
                        % gsi_to_create["IndexName"]
                    )

                gsis_by_name[gsi_to_create["IndexName"]] = GlobalSecondaryIndex.create(
                    gsi_to_create, table.table_key_attrs
                )

        table.global_indexes = list(gsis_by_name.values())
        return table

    def put_item(
        self,
        table_name: str,
        item_attrs: Dict[str, Any],
        expected: Optional[Dict[str, Any]] = None,
        condition_expression: Optional[str] = None,
        expression_attribute_names: Optional[Dict[str, Any]] = None,
        expression_attribute_values: Optional[Dict[str, Any]] = None,
        overwrite: bool = False,
        return_values_on_condition_check_failure: Optional[str] = None,
    ) -> Item:
        table = self.get_table(table_name)
        return table.put_item(
            item_attrs,
            expected,
            condition_expression,
            expression_attribute_names,
            expression_attribute_values,
            overwrite,
            return_values_on_condition_check_failure,
        )

    def get_table_keys_name(
        self, table_name: str, keys: Dict[str, Any]
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Given a set of keys, extracts the key and range key
        """
        table = self.tables.get(table_name)
        if not table:
            return None, None
        else:
            if len(keys) == 1:
                for key in keys:
                    if key in table.hash_key_names:
                        return key, None

            potential_hash, potential_range = None, None
            for key in set(keys):
                if key in table.hash_key_names:
                    potential_hash = key
                elif key in table.range_key_names:
                    potential_range = key
            return potential_hash, potential_range

    def get_keys_value(
        self, table: Table, keys: Dict[str, Any]
    ) -> Tuple[DynamoType, Optional[DynamoType]]:
        if table.hash_key_attr not in keys or (
            table.has_range_key and table.range_key_attr not in keys
        ):
            # "Table has a range key, but no range key was passed into get_item"
            raise MockValidationException("Validation Exception")
        hash_key = DynamoType(keys[table.hash_key_attr])
        range_key = (
            DynamoType(keys[table.range_key_attr]) if table.range_key_attr else None
        )
        return hash_key, range_key

    def get_schema(
        self, table_name: str, index_name: Optional[str]
    ) -> List[Dict[str, Any]]:
        table = self.get_table(table_name)
        if index_name:
            all_indexes = (table.global_indexes or []) + (table.indexes or [])
            indexes_by_name = dict((i.name, i) for i in all_indexes)
            if index_name not in indexes_by_name:
                all_index_names = ", ".join(indexes_by_name.keys())
                raise ResourceNotFoundException(
                    f"Invalid index: {index_name} for table: {table_name}. Available indexes are: {all_index_names}"
                )

            return indexes_by_name[index_name].schema
        else:
            return table.schema

    def get_table(self, table_name: str) -> Table:
        if table_name not in self.tables:
            raise ResourceNotFoundException()
        return self.tables[table_name]

    def get_item(
        self,
        table_name: str,
        keys: Dict[str, Any],
        projection_expressions: Optional[List[List[str]]] = None,
    ) -> Optional[Item]:
        table = self.get_table(table_name)
        hash_key, range_key = self.get_keys_value(table, keys)
        return table.get_item(hash_key, range_key, projection_expressions)

    def query(
        self,
        table_name: str,
        hash_key_dict: Dict[str, Any],
        range_comparison: Optional[str],
        range_value_dicts: List[Dict[str, Any]],
        limit: int,
        exclusive_start_key: Dict[str, Any],
        scan_index_forward: bool,
        projection_expressions: Optional[List[List[str]]],
        index_name: Optional[str] = None,
        consistent_read: bool = False,
        expr_names: Optional[Dict[str, str]] = None,
        expr_values: Optional[Dict[str, Dict[str, str]]] = None,
        filter_expression: Optional[str] = None,
        **filter_kwargs: Any,
    ) -> Tuple[List[Item], int, Optional[Dict[str, Any]]]:
        table = self.get_table(table_name)

        hash_key = DynamoType(hash_key_dict)
        range_values = [DynamoType(range_value) for range_value in range_value_dicts]

        filter_expression_op = get_filter_expression(
            filter_expression, expr_names, expr_values
        )

        return table.query(
            hash_key,
            range_comparison,
            range_values,
            limit,
            exclusive_start_key,
            scan_index_forward,
            projection_expressions,
            index_name,
            consistent_read,
            filter_expression_op,
            **filter_kwargs,
        )

    def scan(
        self,
        table_name: str,
        filters: Dict[str, Any],
        limit: int,
        exclusive_start_key: Dict[str, Any],
        filter_expression: Optional[str],
        expr_names: Dict[str, Any],
        expr_values: Dict[str, Any],
        index_name: str,
        consistent_read: bool,
        projection_expression: Optional[List[List[str]]],
    ) -> Tuple[List[Item], int, Optional[Dict[str, Any]]]:
        table = self.get_table(table_name)

        scan_filters: Dict[str, Any] = {}
        for key, (comparison_operator, comparison_values) in filters.items():
            dynamo_types = [DynamoType(value) for value in comparison_values]
            scan_filters[key] = (comparison_operator, dynamo_types)

        filter_expression_op = get_filter_expression(
            filter_expression, expr_names, expr_values
        )

        return table.scan(
            scan_filters,
            limit,
            exclusive_start_key,
            filter_expression_op,
            index_name,
            consistent_read,
            projection_expression,
        )

    def update_item(
        self,
        table_name: str,
        key: Dict[str, Dict[str, Any]],
        update_expression: str,
        expression_attribute_names: Dict[str, Any],
        expression_attribute_values: Dict[str, Any],
        attribute_updates: Optional[Dict[str, Any]] = None,
        expected: Optional[Dict[str, Any]] = None,
        condition_expression: Optional[str] = None,
        return_values_on_condition_check_failure: Optional[str] = None,
    ) -> Item:
        table = self.get_table(table_name)

        # Support spaces between operators in an update expression
        # E.g. `a = b + c` -> `a=b+c`
        if update_expression:
            # Parse expression to get validation errors
            update_expression_ast = UpdateExpressionParser.make(update_expression)
            update_expression = re.sub(r"\s*([=\+-])\s*", "\\1", update_expression)
            update_expression_ast.validate()

        if all([table.hash_key_attr in key, table.range_key_attr in key]):
            # Covers cases where table has hash and range keys, ``key`` param
            # will be a dict
            hash_value = DynamoType(key[table.hash_key_attr])
            range_value = DynamoType(key[table.range_key_attr])  # type: ignore[index]
        elif table.hash_key_attr in key:
            # Covers tables that have a range key where ``key`` param is a dict
            hash_value = DynamoType(key[table.hash_key_attr])
            range_value = None
        else:
            # Covers other cases
            hash_value = DynamoType(key)
            range_value = None

        item = table.get_item(hash_value, range_value)
        orig_item = copy.deepcopy(item)

        if not expected:
            expected = {}

        if not get_expected(expected).expr(item):
            raise ConditionalCheckFailed
        condition_op = get_filter_expression(
            condition_expression,
            expression_attribute_names,
            expression_attribute_values,
        )
        if not condition_op.expr(item):
            if (
                return_values_on_condition_check_failure == "ALL_OLD"
                and item is not None
            ):
                raise ConditionalCheckFailed(item=item.to_json()["Attributes"])
            else:
                raise ConditionalCheckFailed

        # Update does not fail on new items, so create one
        if item is None:
            if update_expression:
                # Validate AST before creating anything
                item = Item(hash_value, range_value, attrs={})
                UpdateExpressionValidator(
                    update_expression_ast,
                    expression_attribute_names=expression_attribute_names,
                    expression_attribute_values=expression_attribute_values,
                    item=item,
                    table=table,
                ).validate()
            data: Dict[str, Any] = {
                table.hash_key_attr: {hash_value.type: hash_value.value}
            }
            if range_value:
                data.update(
                    {table.range_key_attr: {range_value.type: range_value.value}}  # type: ignore[dict-item]
                )

            table.put_item(data)
            item = table.get_item(hash_value, range_value)

        if attribute_updates:
            item.validate_no_empty_key_values(attribute_updates, table.attribute_keys)  # type: ignore[union-attr]

        if update_expression:
            validator = UpdateExpressionValidator(
                update_expression_ast,
                expression_attribute_names=expression_attribute_names,
                expression_attribute_values=expression_attribute_values,
                item=item,  # type: ignore[arg-type]
                table=table,
            )
            validated_ast = validator.validate()
            validated_ast.normalize()
            try:
                UpdateExpressionExecutor(
                    validated_ast,
                    item,  # type: ignore[arg-type]
                    expression_attribute_names,
                ).execute()
            except ItemSizeTooLarge:
                raise ItemSizeToUpdateTooLarge()
        else:
            item.update_with_attribute_updates(attribute_updates)  # type: ignore
        if table.stream_shard is not None:
            table.stream_shard.add(orig_item, item)
        return item  # type: ignore[return-value]

    def delete_item(
        self,
        table_name: str,
        key: Dict[str, Any],
        expression_attribute_names: Optional[Dict[str, Any]] = None,
        expression_attribute_values: Optional[Dict[str, Any]] = None,
        condition_expression: Optional[str] = None,
    ) -> Optional[Item]:
        table = self.get_table(table_name)

        hash_value, range_value = self.get_keys_value(table, key)
        item = table.get_item(hash_value, range_value)

        condition_op = get_filter_expression(
            condition_expression,
            expression_attribute_names,
            expression_attribute_values,
        )
        if not condition_op.expr(item):
            raise ConditionalCheckFailed

        return table.delete_item(hash_value, range_value)

    def update_time_to_live(self, table_name: str, ttl_spec: Dict[str, Any]) -> None:
        table = self.tables.get(table_name)
        if table is None:
            raise JsonRESTError("ResourceNotFound", "Table not found")

        if "Enabled" not in ttl_spec or "AttributeName" not in ttl_spec:
            raise JsonRESTError(
                "InvalidParameterValue",
                "TimeToLiveSpecification does not contain Enabled and AttributeName",
            )

        if ttl_spec["Enabled"]:
            table.ttl["TimeToLiveStatus"] = "ENABLED"
        else:
            table.ttl["TimeToLiveStatus"] = "DISABLED"
        table.ttl["AttributeName"] = ttl_spec["AttributeName"]

    def describe_time_to_live(self, table_name: str) -> Dict[str, Any]:
        table = self.tables.get(table_name)
        if table is None:
            raise JsonRESTError("ResourceNotFound", "Table not found")

        return table.ttl

    def transact_write_items(self, transact_items: List[Dict[str, Any]]) -> None:
        if len(transact_items) > 100:
            raise TooManyTransactionsException()
        # Create a backup in case any of the transactions fail
        original_table_state = copy.deepcopy(self.tables)
        target_items: Set[Tuple[str, str]] = set()

        def check_unicity(table_name: str, key: Dict[str, Any]) -> None:
            item = (str(table_name), str(key))
            if item in target_items:
                raise MultipleTransactionsException()
            target_items.add(item)

        errors: List[
            Union[Tuple[str, str, Dict[str, Any]], Tuple[None, None, None]]
        ] = []  # [(Code, Message, Item), ..]
        for item in transact_items:
            original_item: Optional[Dict[str, Any]] = None
            # check transact writes are not performing multiple operations
            # in the same item
            if len(list(item.keys())) > 1:
                raise TransactWriteSingleOpException

            try:
                if "ConditionCheck" in item:
                    item = item["ConditionCheck"]
                    key = item["Key"]
                    table_name = item["TableName"]
                    check_unicity(table_name, key)
                    condition_expression = item.get("ConditionExpression", None)
                    expression_attribute_names = item.get(
                        "ExpressionAttributeNames", None
                    )
                    expression_attribute_values = item.get(
                        "ExpressionAttributeValues", None
                    )
                    current = self.get_item(table_name, key)

                    condition_op = get_filter_expression(
                        condition_expression,
                        expression_attribute_names,
                        expression_attribute_values,
                    )
                    if not condition_op.expr(current):
                        raise ConditionalCheckFailed()
                elif "Put" in item:
                    item = item["Put"]
                    attrs = item["Item"]
                    table_name = item["TableName"]
                    table = self.get_table(table_name)
                    key = {table.hash_key_attr: attrs[table.hash_key_attr]}
                    if table.range_key_attr is not None:
                        key[table.range_key_attr] = attrs[table.range_key_attr]
                    check_unicity(table_name, key)
                    condition_expression = item.get("ConditionExpression", None)
                    expression_attribute_names = item.get(
                        "ExpressionAttributeNames", None
                    )
                    expression_attribute_values = item.get(
                        "ExpressionAttributeValues", None
                    )

                    return_values_on_condition_check_failure = item.get(
                        "ReturnValuesOnConditionCheckFailure", None
                    )
                    current = self.get_item(table_name, attrs)
                    if (
                        return_values_on_condition_check_failure == "ALL_OLD"
                        and current
                    ):
                        original_item = current.to_json()["Attributes"]

                    self.put_item(
                        table_name,
                        attrs,
                        condition_expression=condition_expression,
                        expression_attribute_names=expression_attribute_names,
                        expression_attribute_values=expression_attribute_values,
                    )
                elif "Delete" in item:
                    item = item["Delete"]
                    key = item["Key"]
                    table_name = item["TableName"]
                    check_unicity(table_name, key)
                    condition_expression = item.get("ConditionExpression", None)
                    expression_attribute_names = item.get(
                        "ExpressionAttributeNames", None
                    )
                    expression_attribute_values = item.get(
                        "ExpressionAttributeValues", None
                    )
                    self.delete_item(
                        table_name,
                        key,
                        condition_expression=condition_expression,
                        expression_attribute_names=expression_attribute_names,
                        expression_attribute_values=expression_attribute_values,
                    )
                elif "Update" in item:
                    item = item["Update"]
                    key = item["Key"]
                    table_name = item["TableName"]
                    check_unicity(table_name, key)
                    update_expression = item["UpdateExpression"]
                    condition_expression = item.get("ConditionExpression", None)
                    expression_attribute_names = item.get(
                        "ExpressionAttributeNames", None
                    )
                    expression_attribute_values = item.get(
                        "ExpressionAttributeValues", None
                    )
                    self.update_item(
                        table_name,
                        key,
                        update_expression=update_expression,
                        condition_expression=condition_expression,
                        expression_attribute_names=expression_attribute_names,
                        expression_attribute_values=expression_attribute_values,
                    )
                else:
                    raise ValueError
                errors.append((None, None, None))
            except MultipleTransactionsException:
                # Rollback to the original state, and reraise the error
                self.tables = original_table_state
                raise MultipleTransactionsException()
            except Exception as e:  # noqa: E722 Do not use bare except
                errors.append((type(e).__name__, e.message, original_item))  # type: ignore
        if any([code is not None for code, _, _ in errors]):
            # Rollback to the original state, and reraise the errors
            self.tables = original_table_state
            raise TransactionCanceledException(errors)

    def describe_continuous_backups(self, table_name: str) -> Dict[str, Any]:
        try:
            table = self.get_table(table_name)
        except ResourceNotFoundException:
            raise TableNotFoundException(table_name)

        return table.continuous_backups

    def update_continuous_backups(
        self, table_name: str, point_in_time_spec: Dict[str, Any]
    ) -> Dict[str, Any]:
        try:
            table = self.get_table(table_name)
        except ResourceNotFoundException:
            raise TableNotFoundException(table_name)

        if (
            point_in_time_spec["PointInTimeRecoveryEnabled"]
            and table.continuous_backups["PointInTimeRecoveryDescription"][
                "PointInTimeRecoveryStatus"
            ]
            == "DISABLED"
        ):
            table.continuous_backups["PointInTimeRecoveryDescription"] = {
                "PointInTimeRecoveryStatus": "ENABLED",
                "EarliestRestorableDateTime": unix_time(),
                "LatestRestorableDateTime": unix_time(),
            }
        elif not point_in_time_spec["PointInTimeRecoveryEnabled"]:
            table.continuous_backups["PointInTimeRecoveryDescription"] = {
                "PointInTimeRecoveryStatus": "DISABLED"
            }

        return table.continuous_backups

    def get_backup(self, backup_arn: str) -> Backup:
        if backup_arn not in self.backups:
            raise BackupNotFoundException(backup_arn)
        return self.backups[backup_arn]

    def list_backups(self, table_name: str) -> List[Backup]:
        backups = list(self.backups.values())
        if table_name is not None:
            backups = [backup for backup in backups if backup.table.name == table_name]
        return backups

    def create_backup(self, table_name: str, backup_name: str) -> Backup:
        try:
            table = self.get_table(table_name)
        except ResourceNotFoundException:
            raise TableNotFoundException(table_name)
        backup = Backup(self.account_id, self.region_name, backup_name, table)
        self.backups[backup.arn] = backup
        return backup

    def delete_backup(self, backup_arn: str) -> Backup:
        backup = self.get_backup(backup_arn)
        if backup is None:
            raise KeyError()
        backup_deleted = self.backups.pop(backup_arn)
        backup_deleted.status = "DELETED"
        return backup_deleted

    def describe_backup(self, backup_arn: str) -> Backup:
        backup = self.get_backup(backup_arn)
        if backup is None:
            raise KeyError()
        return backup

    def restore_table_from_backup(
        self, target_table_name: str, backup_arn: str
    ) -> RestoredTable:
        backup = self.get_backup(backup_arn)
        if target_table_name in self.tables:
            raise TableAlreadyExistsException(target_table_name)
        new_table = RestoredTable(
            target_table_name,
            account_id=self.account_id,
            region=self.region_name,
            backup=backup,
        )
        self.tables[target_table_name] = new_table
        return new_table

    def restore_table_to_point_in_time(
        self, target_table_name: str, source_table_name: str
    ) -> RestoredPITTable:
        """
        Currently this only accepts the source and target table elements, and will
        copy all items from the source without respect to other arguments.
        """

        try:
            source = self.get_table(source_table_name)
        except ResourceNotFoundException:
            raise SourceTableNotFoundException(source_table_name)
        if target_table_name in self.tables:
            raise TableAlreadyExistsException(target_table_name)
        new_table = RestoredPITTable(
            target_table_name,
            account_id=self.account_id,
            region=self.region_name,
            source=source,
        )
        self.tables[target_table_name] = new_table
        return new_table

    ######################
    # LIST of methods where the logic completely resides in responses.py
    # Duplicated here so that the implementation coverage script is aware
    # TODO: Move logic here
    ######################

    def batch_get_item(self) -> None:
        pass

    def batch_write_item(self) -> None:
        pass

    def transact_get_items(self) -> None:
        pass

    def execute_statement(
        self, statement: str, parameters: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Pagination is not yet implemented.

        Parsing is highly experimental - please raise an issue if you find any bugs.
        """
        # We need to execute a statement - but we don't know which table
        # Just pass all tables to PartiQL
        source_data: Dict[str, str] = dict()
        for table in self.tables.values():
            source_data[table.name] = [  # type: ignore
                item.to_json()["Attributes"] for item in table.all_items()
            ]

        return_data, updates_per_table = partiql.query(
            statement, source_data, parameters
        )

        for table_name, updates in updates_per_table.items():
            table = self.tables[table_name]
            for before, after in updates:
                if after is None and before is not None:
                    # DELETE
                    hash_key = DynamoType(before[table.hash_key_attr])
                    if table.range_key_attr:
                        range_key = DynamoType(before[table.range_key_attr])
                    else:
                        range_key = None
                    table.delete_item(hash_key, range_key)
                elif before is None and after is not None:
                    # CREATE
                    table.put_item(after)
                elif before is not None and after is not None:
                    # UPDATE
                    table.put_item(after)

        return return_data

    def execute_transaction(
        self, statements: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Please see the documentation for `execute_statement` to see the limitations of what is supported.
        """
        responses = []
        for stmt in statements:
            items = self.execute_statement(
                statement=stmt["Statement"], parameters=stmt.get("Parameters", [])
            )
            responses.extend([{"Item": item} for item in items])
        return responses

    def batch_execute_statement(
        self, statements: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Please see the documentation for `execute_statement` to see the limitations of what is supported.
        """
        responses = []
        # Validation
        for stmt in statements:
            metadata = partiql.get_query_metadata(stmt["Statement"])
            table_name = metadata.get_table_names()[0]
            response = {}
            filter_keys = metadata.get_filter_names()
            if table_name not in self.tables:
                response["Error"] = {
                    "Code": "ResourceNotFound",
                    "Message": "Requested resource not found",
                }
            else:
                response["TableName"] = table_name
                if metadata.is_select_query():
                    table = self.tables[table_name]
                    for required_attr in table.table_key_attrs:
                        if required_attr not in filter_keys:
                            response["Error"] = {
                                "Code": "ValidationError",
                                "Message": "Select statements within BatchExecuteStatement must specify the primary key in the where clause.",
                            }
            responses.append(response)

        # Execution
        for idx, stmt in enumerate(statements):
            if "Error" in responses[idx]:
                continue
            try:
                items = self.execute_statement(
                    statement=stmt["Statement"], parameters=stmt.get("Parameters", [])
                )
                # Statements should always contain a HashKey and SortKey
                # An item with those keys may not exist
                if items:
                    # But if it does, it will always only contain one item at most
                    responses[idx]["Item"] = items[0]
            except Exception as e:
                responses[idx] = {"Error": {"Code": e.name, "Message": e.message}}  # type: ignore
        return responses

    def import_table(
        self,
        s3_source: Dict[str, str],
        input_format: Optional[str],
        compression_type: Optional[str],
        table_name: str,
        billing_mode: str,
        throughput: Optional[Dict[str, int]],
        key_schema: List[Dict[str, str]],
        global_indexes: Optional[List[Dict[str, Any]]],
        attrs: List[Dict[str, str]],
    ) -> TableImport:
        """
        Only InputFormat=DYNAMODB_JSON is supported so far.
        InputCompressionType=ZSTD is not supported.
        Other parameters that are not supported: InputFormatOptions, CloudWatchLogGroupArn
        """
        table_import = TableImport(
            account_id=self.account_id,
            s3_source=s3_source,
            region_name=self.region_name,
            table_name=table_name,
            billing_mode=billing_mode,
            throughput=throughput,
            key_schema=key_schema,
            global_indexes=global_indexes,
            attrs=attrs,
            compression_type=compression_type,
        )
        self.table_imports[table_import.arn] = table_import
        table_import.start()
        return table_import

    def describe_import(self, import_arn: str) -> TableImport:
        return self.table_imports[import_arn]


dynamodb_backends = BackendDict(DynamoDBBackend, "dynamodb")
