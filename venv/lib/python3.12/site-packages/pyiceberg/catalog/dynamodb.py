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
import uuid
from time import time
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    Union,
)

import boto3

from pyiceberg.catalog import (
    BOTOCORE_SESSION,
    ICEBERG,
    METADATA_LOCATION,
    PREVIOUS_METADATA_LOCATION,
    TABLE_TYPE,
    MetastoreCatalog,
    PropertiesUpdateSummary,
)
from pyiceberg.exceptions import (
    ConditionalCheckFailedException,
    GenericDynamoDbError,
    NamespaceAlreadyExistsError,
    NamespaceNotEmptyError,
    NoSuchIcebergTableError,
    NoSuchNamespaceError,
    NoSuchPropertyException,
    NoSuchTableError,
    TableAlreadyExistsError,
)
from pyiceberg.io import AWS_ACCESS_KEY_ID, AWS_REGION, AWS_SECRET_ACCESS_KEY, AWS_SESSION_TOKEN, load_file_io
from pyiceberg.partitioning import UNPARTITIONED_PARTITION_SPEC, PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.serializers import FromInputFile
from pyiceberg.table import CommitTableResponse, Table, TableProperties
from pyiceberg.table.locations import load_location_provider
from pyiceberg.table.metadata import new_table_metadata
from pyiceberg.table.sorting import UNSORTED_SORT_ORDER, SortOrder
from pyiceberg.table.update import (
    TableRequirement,
    TableUpdate,
)
from pyiceberg.typedef import EMPTY_DICT, Identifier, Properties
from pyiceberg.utils.properties import get_first_property_value

if TYPE_CHECKING:
    import pyarrow as pa
    from mypy_boto3_dynamodb.client import DynamoDBClient


DYNAMODB_CLIENT = "dynamodb"

DYNAMODB_COL_IDENTIFIER = "identifier"
DYNAMODB_COL_NAMESPACE = "namespace"
DYNAMODB_COL_VERSION = "v"
DYNAMODB_COL_UPDATED_AT = "updated_at"
DYNAMODB_COL_CREATED_AT = "created_at"
DYNAMODB_NAMESPACE = "NAMESPACE"
DYNAMODB_NAMESPACE_GSI = "namespace-identifier"
DYNAMODB_PAY_PER_REQUEST = "PAY_PER_REQUEST"

DYNAMODB_TABLE_NAME = "table-name"
DYNAMODB_TABLE_NAME_DEFAULT = "iceberg"

PROPERTY_KEY_PREFIX = "p."

ACTIVE = "ACTIVE"
ITEM = "Item"

DYNAMODB_PROFILE_NAME = "dynamodb.profile-name"
DYNAMODB_REGION = "dynamodb.region"
DYNAMODB_ACCESS_KEY_ID = "dynamodb.access-key-id"
DYNAMODB_SECRET_ACCESS_KEY = "dynamodb.secret-access-key"
DYNAMODB_SESSION_TOKEN = "dynamodb.session-token"


class DynamoDbCatalog(MetastoreCatalog):
    def __init__(self, name: str, client: Optional["DynamoDBClient"] = None, **properties: str):
        """Dynamodb catalog.

        Args:
            name: Name to identify the catalog.
            client: An optional boto3 dynamodb client.
            properties: Properties for dynamodb client construction and configuration.
        """
        super().__init__(name, **properties)
        if client is not None:
            self.dynamodb = client
        else:
            session = boto3.Session(
                profile_name=properties.get(DYNAMODB_PROFILE_NAME),
                region_name=get_first_property_value(properties, DYNAMODB_REGION, AWS_REGION),
                botocore_session=properties.get(BOTOCORE_SESSION),
                aws_access_key_id=get_first_property_value(properties, DYNAMODB_ACCESS_KEY_ID, AWS_ACCESS_KEY_ID),
                aws_secret_access_key=get_first_property_value(properties, DYNAMODB_SECRET_ACCESS_KEY, AWS_SECRET_ACCESS_KEY),
                aws_session_token=get_first_property_value(properties, DYNAMODB_SESSION_TOKEN, AWS_SESSION_TOKEN),
            )
            self.dynamodb = session.client(DYNAMODB_CLIENT)

        self.dynamodb_table_name = self.properties.get(DYNAMODB_TABLE_NAME, DYNAMODB_TABLE_NAME_DEFAULT)
        self._ensure_catalog_table_exists_or_create()

    def _ensure_catalog_table_exists_or_create(self) -> None:
        if self._dynamodb_table_exists():
            return None

        try:
            self.dynamodb.create_table(
                TableName=self.dynamodb_table_name,
                AttributeDefinitions=CREATE_CATALOG_ATTRIBUTE_DEFINITIONS,
                KeySchema=CREATE_CATALOG_KEY_SCHEMA,
                GlobalSecondaryIndexes=CREATE_CATALOG_GLOBAL_SECONDARY_INDEXES,
                BillingMode=DYNAMODB_PAY_PER_REQUEST,
            )
        except (
            self.dynamodb.exceptions.ResourceInUseException,
            self.dynamodb.exceptions.LimitExceededException,
            self.dynamodb.exceptions.InternalServerError,
        ) as e:
            raise GenericDynamoDbError(e.message) from e

    def _dynamodb_table_exists(self) -> bool:
        try:
            response = self.dynamodb.describe_table(TableName=self.dynamodb_table_name)
        except self.dynamodb.exceptions.ResourceNotFoundException:
            return False
        except self.dynamodb.exceptions.InternalServerError as e:
            raise GenericDynamoDbError(e.message) from e

        if response["Table"]["TableStatus"] != ACTIVE:
            raise GenericDynamoDbError(f"DynamoDB table for catalog {self.dynamodb_table_name} is not {ACTIVE}")
        else:
            return True

    def create_table(
        self,
        identifier: Union[str, Identifier],
        schema: Union[Schema, "pa.Schema"],
        location: Optional[str] = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
    ) -> Table:
        """
        Create an Iceberg table.

        Args:
            identifier: Table identifier.
            schema: Table's schema.
            location: Location for the table. Optional Argument.
            partition_spec: PartitionSpec for the table.
            sort_order: SortOrder for the table.
            properties: Table properties that can be a string based dictionary.

        Returns:
            Table: the created table instance.

        Raises:
            AlreadyExistsError: If a table with the name already exists.
            ValueError: If the identifier is invalid, or no path is given to store metadata.

        """
        schema: Schema = self._convert_schema_if_needed(  # type: ignore
            schema,
            int(properties.get(TableProperties.FORMAT_VERSION, TableProperties.DEFAULT_FORMAT_VERSION)),  # type: ignore
        )

        database_name, table_name = self.identifier_to_database_and_table(identifier)

        location = self._resolve_table_location(location, database_name, table_name)
        provider = load_location_provider(table_location=location, table_properties=properties)
        metadata_location = provider.new_table_metadata_file_location()

        metadata = new_table_metadata(
            location=location, schema=schema, partition_spec=partition_spec, sort_order=sort_order, properties=properties
        )
        io = load_file_io(properties=self.properties, location=metadata_location)
        self._write_metadata(metadata, io, metadata_location)

        self._ensure_namespace_exists(database_name=database_name)

        try:
            self._put_dynamo_item(
                item=_get_create_table_item(
                    database_name=database_name, table_name=table_name, properties=properties, metadata_location=metadata_location
                ),
                condition_expression=f"attribute_not_exists({DYNAMODB_COL_IDENTIFIER})",
            )
        except ConditionalCheckFailedException as e:
            raise TableAlreadyExistsError(f"Table {database_name}.{table_name} already exists") from e

        return self.load_table(identifier=identifier)

    def register_table(self, identifier: Union[str, Identifier], metadata_location: str) -> Table:
        """Register a new table using existing metadata.

        Args:
            identifier Union[str, Identifier]: Table identifier for the table
            metadata_location str: The location to the metadata

        Returns:
            Table: The newly registered table

        Raises:
            TableAlreadyExistsError: If the table already exists
        """
        raise NotImplementedError

    def commit_table(
        self, table: Table, requirements: Tuple[TableRequirement, ...], updates: Tuple[TableUpdate, ...]
    ) -> CommitTableResponse:
        """Commit updates to a table.

        Args:
            table (Table): The table to be updated.
            requirements: (Tuple[TableRequirement, ...]): Table requirements.
            updates: (Tuple[TableUpdate, ...]): Table updates.

        Returns:
            CommitTableResponse: The updated metadata.

        Raises:
            NoSuchTableError: If a table with the given identifier does not exist.
            CommitFailedException: Requirement not met, or a conflict with a concurrent commit.
        """
        raise NotImplementedError

    def load_table(self, identifier: Union[str, Identifier]) -> Table:
        """
        Load the table's metadata and returns the table instance.

        You can also use this method to check for table existence using 'try catalog.table() except TableNotFoundError'.
        Note: This method doesn't scan data stored in the table.

        Args:
            identifier: Table identifier.

        Returns:
            Table: the table instance with its metadata.

        Raises:
            NoSuchTableError: If a table with the name does not exist, or the identifier is invalid.
        """
        database_name, table_name = self.identifier_to_database_and_table(identifier, NoSuchTableError)
        dynamo_table_item = self._get_iceberg_table_item(database_name=database_name, table_name=table_name)
        return self._convert_dynamo_table_item_to_iceberg_table(dynamo_table_item=dynamo_table_item)

    def drop_table(self, identifier: Union[str, Identifier]) -> None:
        """Drop a table.

        Args:
            identifier: Table identifier.

        Raises:
            NoSuchTableError: If a table with the name does not exist, or the identifier is invalid.
        """
        database_name, table_name = self.identifier_to_database_and_table(identifier, NoSuchTableError)

        try:
            self._delete_dynamo_item(
                namespace=database_name,
                identifier=f"{database_name}.{table_name}",
                condition_expression=f"attribute_exists({DYNAMODB_COL_IDENTIFIER})",
            )
        except ConditionalCheckFailedException as e:
            raise NoSuchTableError(f"Table does not exist: {database_name}.{table_name}") from e

    def rename_table(self, from_identifier: Union[str, Identifier], to_identifier: Union[str, Identifier]) -> Table:
        """Rename a fully classified table name.

        This method can only rename Iceberg tables in AWS Glue.

        Args:
            from_identifier: Existing table identifier.
            to_identifier: New table identifier.

        Returns:
            Table: the updated table instance with its metadata.

        Raises:
            ValueError: When from table identifier is invalid.
            NoSuchTableError: When a table with the name does not exist.
            NoSuchIcebergTableError: When from table is not a valid iceberg table.
            NoSuchPropertyException: When from table miss some required properties.
            NoSuchNamespaceError: When the destination namespace doesn't exist.
        """
        from_database_name, from_table_name = self.identifier_to_database_and_table(from_identifier, NoSuchTableError)
        to_database_name, to_table_name = self.identifier_to_database_and_table(to_identifier)

        from_table_item = self._get_iceberg_table_item(database_name=from_database_name, table_name=from_table_name)

        try:
            # Verify that from_identifier is a valid iceberg table
            self._convert_dynamo_table_item_to_iceberg_table(dynamo_table_item=from_table_item)
        except NoSuchPropertyException as e:
            raise NoSuchPropertyException(
                f"Failed to rename table {from_database_name}.{from_table_name} since it is missing required properties"
            ) from e
        except NoSuchIcebergTableError as e:
            raise NoSuchIcebergTableError(
                f"Failed to rename table {from_database_name}.{from_table_name} since it is not a valid iceberg table"
            ) from e

        self._ensure_namespace_exists(database_name=from_database_name)
        self._ensure_namespace_exists(database_name=to_database_name)

        try:
            self._put_dynamo_item(
                item=_get_rename_table_item(
                    from_dynamo_table_item=from_table_item, to_database_name=to_database_name, to_table_name=to_table_name
                ),
                condition_expression=f"attribute_not_exists({DYNAMODB_COL_IDENTIFIER})",
            )
        except ConditionalCheckFailedException as e:
            raise TableAlreadyExistsError(f"Table {to_database_name}.{to_table_name} already exists") from e

        try:
            self.drop_table(from_identifier)
        except (NoSuchTableError, GenericDynamoDbError) as e:
            log_message = f"Failed to drop old table {from_database_name}.{from_table_name}. "

            try:
                self.drop_table(to_identifier)
                log_message += f"Rolled back table creation for {to_database_name}.{to_table_name}."
            except (NoSuchTableError, GenericDynamoDbError):
                log_message += (
                    f"Failed to roll back table creation for {to_database_name}.{to_table_name}. Please clean up manually"
                )

            raise ValueError(log_message) from e

        return self.load_table(to_identifier)

    def create_namespace(self, namespace: Union[str, Identifier], properties: Properties = EMPTY_DICT) -> None:
        """Create a namespace in the catalog.

        Args:
            namespace: Namespace identifier.
            properties: A string dictionary of properties for the given namespace.

        Raises:
            ValueError: If the identifier is invalid.
            AlreadyExistsError: If a namespace with the given name already exists.
        """
        database_name = self.identifier_to_database(namespace)

        try:
            self._put_dynamo_item(
                item=_get_create_database_item(database_name=database_name, properties=properties),
                condition_expression=f"attribute_not_exists({DYNAMODB_COL_NAMESPACE})",
            )
        except ConditionalCheckFailedException as e:
            raise NamespaceAlreadyExistsError(f"Database {database_name} already exists") from e

    def drop_namespace(self, namespace: Union[str, Identifier]) -> None:
        """Drop a namespace.

        A Glue namespace can only be dropped if it is empty.

        Args:
            namespace: Namespace identifier.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist, or the identifier is invalid.
            NamespaceNotEmptyError: If the namespace is not empty.
        """
        database_name = self.identifier_to_database(namespace, NoSuchNamespaceError)
        table_identifiers = self.list_tables(namespace=database_name)

        if len(table_identifiers) > 0:
            raise NamespaceNotEmptyError(f"Database {database_name} is not empty")

        try:
            self._delete_dynamo_item(
                namespace=database_name,
                identifier=DYNAMODB_NAMESPACE,
                condition_expression=f"attribute_exists({DYNAMODB_COL_IDENTIFIER})",
            )
        except ConditionalCheckFailedException as e:
            raise NoSuchNamespaceError(f"Database does not exist: {database_name}") from e

    def list_tables(self, namespace: Union[str, Identifier]) -> List[Identifier]:
        """List Iceberg tables under the given namespace in the catalog.

        Args:
            namespace (str | Identifier): Namespace identifier to search.

        Returns:
            List[Identifier]: list of table identifiers.
        """
        database_name = self.identifier_to_database(namespace, NoSuchNamespaceError)

        paginator = self.dynamodb.get_paginator("query")

        try:
            page_iterator = paginator.paginate(
                TableName=self.dynamodb_table_name,
                IndexName=DYNAMODB_NAMESPACE_GSI,
                KeyConditionExpression=f"{DYNAMODB_COL_NAMESPACE} = :namespace ",
                ExpressionAttributeValues={
                    ":namespace": {
                        "S": database_name,
                    }
                },
            )
        except (
            self.dynamodb.exceptions.ProvisionedThroughputExceededException,
            self.dynamodb.exceptions.RequestLimitExceeded,
            self.dynamodb.exceptions.InternalServerError,
            self.dynamodb.exceptions.ResourceNotFoundException,
        ) as e:
            raise GenericDynamoDbError(e.message) from e

        table_identifiers = []
        for page in page_iterator:
            for item in page["Items"]:
                _dict = _convert_dynamo_item_to_regular_dict(item)
                identifier_col = _dict[DYNAMODB_COL_IDENTIFIER]
                if identifier_col == DYNAMODB_NAMESPACE:
                    continue

                table_identifiers.append(self.identifier_to_tuple(identifier_col))

        return table_identifiers

    def list_namespaces(self, namespace: Union[str, Identifier] = ()) -> List[Identifier]:
        """List top-level namespaces from the catalog.

        We do not support hierarchical namespace.

        Returns:
            List[Identifier]: a List of namespace identifiers.
        """
        # Hierarchical namespace is not supported. Return an empty list
        if namespace:
            return []

        paginator = self.dynamodb.get_paginator("query")

        try:
            page_iterator = paginator.paginate(
                TableName=self.dynamodb_table_name,
                ConsistentRead=True,
                KeyConditionExpression=f"{DYNAMODB_COL_IDENTIFIER} = :identifier",
                ExpressionAttributeValues={
                    ":identifier": {
                        "S": DYNAMODB_NAMESPACE,
                    }
                },
            )
        except (
            self.dynamodb.exceptions.ProvisionedThroughputExceededException,
            self.dynamodb.exceptions.RequestLimitExceeded,
            self.dynamodb.exceptions.InternalServerError,
            self.dynamodb.exceptions.ResourceNotFoundException,
        ) as e:
            raise GenericDynamoDbError(e.message) from e

        database_identifiers = []
        for page in page_iterator:
            for item in page["Items"]:
                _dict = _convert_dynamo_item_to_regular_dict(item)
                namespace_col = _dict[DYNAMODB_COL_NAMESPACE]
                database_identifiers.append(self.identifier_to_tuple(namespace_col))

        return database_identifiers

    def load_namespace_properties(self, namespace: Union[str, Identifier]) -> Properties:
        """
        Get properties for a namespace.

        Args:
            namespace: Namespace identifier.

        Returns:
            Properties: Properties for the given namespace.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist, or identifier is invalid.
        """
        database_name = self.identifier_to_database(namespace, NoSuchNamespaceError)
        namespace_item = self._get_iceberg_namespace_item(database_name=database_name)
        namespace_dict = _convert_dynamo_item_to_regular_dict(namespace_item)
        return _get_namespace_properties(namespace_dict=namespace_dict)

    def update_namespace_properties(
        self, namespace: Union[str, Identifier], removals: Optional[Set[str]] = None, updates: Properties = EMPTY_DICT
    ) -> PropertiesUpdateSummary:
        """
        Remove or update provided property keys for a namespace.

        Args:
            namespace: Namespace identifier
            removals: Set of property keys that need to be removed. Optional Argument.
            updates: Properties to be updated for the given namespace.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist， or identifier is invalid.
            ValueError: If removals and updates have overlapping keys.
        """
        database_name = self.identifier_to_database(namespace, NoSuchNamespaceError)
        namespace_item = self._get_iceberg_namespace_item(database_name=database_name)
        namespace_dict = _convert_dynamo_item_to_regular_dict(namespace_item)
        current_properties = _get_namespace_properties(namespace_dict=namespace_dict)

        properties_update_summary, updated_properties = self._get_updated_props_and_update_summary(
            current_properties=current_properties, removals=removals, updates=updates
        )

        try:
            self._put_dynamo_item(
                item=_get_update_database_item(
                    namespace_item=namespace_item,
                    updated_properties=updated_properties,
                ),
                condition_expression=f"attribute_exists({DYNAMODB_COL_NAMESPACE})",
            )
        except ConditionalCheckFailedException as e:
            raise NoSuchNamespaceError(f"Database {database_name} does not exist") from e

        return properties_update_summary

    def list_views(self, namespace: Union[str, Identifier]) -> List[Identifier]:
        raise NotImplementedError

    def drop_view(self, identifier: Union[str, Identifier]) -> None:
        raise NotImplementedError

    def view_exists(self, identifier: Union[str, Identifier]) -> bool:
        raise NotImplementedError

    def _get_iceberg_table_item(self, database_name: str, table_name: str) -> Dict[str, Any]:
        try:
            return self._get_dynamo_item(identifier=f"{database_name}.{table_name}", namespace=database_name)
        except ValueError as e:
            raise NoSuchTableError(f"Table does not exist: {database_name}.{table_name}") from e

    def _get_iceberg_namespace_item(self, database_name: str) -> Dict[str, Any]:
        try:
            return self._get_dynamo_item(identifier=DYNAMODB_NAMESPACE, namespace=database_name)
        except ValueError as e:
            raise NoSuchNamespaceError(f"Namespace does not exist: {database_name}") from e

    def _ensure_namespace_exists(self, database_name: str) -> Dict[str, Any]:
        return self._get_iceberg_namespace_item(database_name)

    def _get_dynamo_item(self, identifier: str, namespace: str) -> Dict[str, Any]:
        try:
            response = self.dynamodb.get_item(
                TableName=self.dynamodb_table_name,
                ConsistentRead=True,
                Key={
                    DYNAMODB_COL_IDENTIFIER: {
                        "S": identifier,
                    },
                    DYNAMODB_COL_NAMESPACE: {
                        "S": namespace,
                    },
                },
            )
            if ITEM in response:
                return response[ITEM]
            else:
                raise ValueError(f"Item not found. identifier: {identifier} - namespace: {namespace}")
        except self.dynamodb.exceptions.ResourceNotFoundException as e:
            raise ValueError(f"Item not found. identifier: {identifier} - namespace: {namespace}") from e
        except (
            self.dynamodb.exceptions.ProvisionedThroughputExceededException,
            self.dynamodb.exceptions.RequestLimitExceeded,
            self.dynamodb.exceptions.InternalServerError,
        ) as e:
            raise GenericDynamoDbError(e.message) from e

    def _put_dynamo_item(self, item: Dict[str, Any], condition_expression: str) -> None:
        try:
            self.dynamodb.put_item(TableName=self.dynamodb_table_name, Item=item, ConditionExpression=condition_expression)
        except self.dynamodb.exceptions.ConditionalCheckFailedException as e:
            raise ConditionalCheckFailedException(f"Condition expression check failed: {condition_expression} - {item}") from e
        except (
            self.dynamodb.exceptions.ProvisionedThroughputExceededException,
            self.dynamodb.exceptions.RequestLimitExceeded,
            self.dynamodb.exceptions.InternalServerError,
            self.dynamodb.exceptions.ResourceNotFoundException,
            self.dynamodb.exceptions.ItemCollectionSizeLimitExceededException,
            self.dynamodb.exceptions.TransactionConflictException,
        ) as e:
            raise GenericDynamoDbError(e.message) from e

    def _delete_dynamo_item(self, namespace: str, identifier: str, condition_expression: str) -> None:
        try:
            self.dynamodb.delete_item(
                TableName=self.dynamodb_table_name,
                Key={
                    DYNAMODB_COL_IDENTIFIER: {
                        "S": identifier,
                    },
                    DYNAMODB_COL_NAMESPACE: {
                        "S": namespace,
                    },
                },
                ConditionExpression=condition_expression,
            )
        except self.dynamodb.exceptions.ConditionalCheckFailedException as e:
            raise ConditionalCheckFailedException(
                f"Condition expression check failed: {condition_expression} - {identifier}"
            ) from e
        except (
            self.dynamodb.exceptions.ProvisionedThroughputExceededException,
            self.dynamodb.exceptions.RequestLimitExceeded,
            self.dynamodb.exceptions.InternalServerError,
            self.dynamodb.exceptions.ResourceNotFoundException,
            self.dynamodb.exceptions.ItemCollectionSizeLimitExceededException,
            self.dynamodb.exceptions.TransactionConflictException,
        ) as e:
            raise GenericDynamoDbError(e.message) from e

    def _convert_dynamo_table_item_to_iceberg_table(self, dynamo_table_item: Dict[str, Any]) -> Table:
        table_dict = _convert_dynamo_item_to_regular_dict(dynamo_table_item)

        for prop in [_add_property_prefix(prop) for prop in (TABLE_TYPE, METADATA_LOCATION)] + [
            DYNAMODB_COL_IDENTIFIER,
            DYNAMODB_COL_NAMESPACE,
            DYNAMODB_COL_CREATED_AT,
        ]:
            if prop not in table_dict.keys():
                raise NoSuchPropertyException(f"Iceberg required property {prop} is missing: {dynamo_table_item}")

        table_type = table_dict[_add_property_prefix(TABLE_TYPE)]
        identifier = table_dict[DYNAMODB_COL_IDENTIFIER]
        metadata_location = table_dict[_add_property_prefix(METADATA_LOCATION)]
        database_name, table_name = self.identifier_to_database_and_table(identifier, NoSuchTableError)

        if table_type.lower() != ICEBERG:
            raise NoSuchIcebergTableError(
                f"Property table_type is {table_type}, expected {ICEBERG}: {database_name}.{table_name}"
            )

        io = load_file_io(properties=self.properties, location=metadata_location)
        file = io.new_input(metadata_location)
        metadata = FromInputFile.table_metadata(file)
        return Table(
            identifier=(database_name, table_name),
            metadata=metadata,
            metadata_location=metadata_location,
            io=self._load_file_io(metadata.properties, metadata_location),
            catalog=self,
        )

    def _get_default_warehouse_location(self, database_name: str, table_name: str) -> str:
        """Override the default warehouse location to follow Hive-style conventions."""
        return self._get_hive_style_warehouse_location(database_name, table_name)


def _get_create_table_item(database_name: str, table_name: str, properties: Properties, metadata_location: str) -> Dict[str, Any]:
    current_timestamp_ms = str(round(time() * 1000))
    _dict = {
        DYNAMODB_COL_IDENTIFIER: {
            "S": f"{database_name}.{table_name}",
        },
        DYNAMODB_COL_NAMESPACE: {
            "S": database_name,
        },
        DYNAMODB_COL_VERSION: {
            "S": str(uuid.uuid4()),
        },
        DYNAMODB_COL_CREATED_AT: {
            "N": current_timestamp_ms,
        },
        DYNAMODB_COL_UPDATED_AT: {
            "N": current_timestamp_ms,
        },
    }

    for key, val in properties.items():
        _dict[_add_property_prefix(key)] = {"S": val}

    _dict[_add_property_prefix(TABLE_TYPE)] = {"S": ICEBERG.upper()}
    _dict[_add_property_prefix(METADATA_LOCATION)] = {"S": metadata_location}
    _dict[_add_property_prefix(PREVIOUS_METADATA_LOCATION)] = {"S": ""}

    return _dict


def _get_rename_table_item(from_dynamo_table_item: Dict[str, Any], to_database_name: str, to_table_name: str) -> Dict[str, Any]:
    _dict = from_dynamo_table_item
    current_timestamp_ms = str(round(time() * 1000))
    _dict[DYNAMODB_COL_IDENTIFIER]["S"] = f"{to_database_name}.{to_table_name}"
    _dict[DYNAMODB_COL_NAMESPACE]["S"] = to_database_name
    _dict[DYNAMODB_COL_VERSION]["S"] = str(uuid.uuid4())
    _dict[DYNAMODB_COL_UPDATED_AT]["N"] = current_timestamp_ms
    return _dict


def _get_create_database_item(database_name: str, properties: Properties) -> Dict[str, Any]:
    current_timestamp_ms = str(round(time() * 1000))
    _dict = {
        DYNAMODB_COL_IDENTIFIER: {
            "S": DYNAMODB_NAMESPACE,
        },
        DYNAMODB_COL_NAMESPACE: {
            "S": database_name,
        },
        DYNAMODB_COL_VERSION: {
            "S": str(uuid.uuid4()),
        },
        DYNAMODB_COL_CREATED_AT: {
            "N": current_timestamp_ms,
        },
        DYNAMODB_COL_UPDATED_AT: {
            "N": current_timestamp_ms,
        },
    }

    for key, val in properties.items():
        _dict[_add_property_prefix(key)] = {"S": val}

    return _dict


def _get_update_database_item(namespace_item: Dict[str, Any], updated_properties: Properties) -> Dict[str, Any]:
    current_timestamp_ms = str(round(time() * 1000))

    _dict = {
        DYNAMODB_COL_IDENTIFIER: namespace_item[DYNAMODB_COL_IDENTIFIER],
        DYNAMODB_COL_NAMESPACE: namespace_item[DYNAMODB_COL_NAMESPACE],
        DYNAMODB_COL_VERSION: {
            "S": str(uuid.uuid4()),
        },
        DYNAMODB_COL_CREATED_AT: namespace_item[DYNAMODB_COL_CREATED_AT],
        DYNAMODB_COL_UPDATED_AT: {
            "N": current_timestamp_ms,
        },
    }

    for key, val in updated_properties.items():
        _dict[_add_property_prefix(key)] = {"S": val}

    return _dict


CREATE_CATALOG_ATTRIBUTE_DEFINITIONS = [
    {
        "AttributeName": DYNAMODB_COL_IDENTIFIER,
        "AttributeType": "S",
    },
    {
        "AttributeName": DYNAMODB_COL_NAMESPACE,
        "AttributeType": "S",
    },
]

CREATE_CATALOG_KEY_SCHEMA = [
    {
        "AttributeName": DYNAMODB_COL_IDENTIFIER,
        "KeyType": "HASH",
    },
    {
        "AttributeName": DYNAMODB_COL_NAMESPACE,
        "KeyType": "RANGE",
    },
]


CREATE_CATALOG_GLOBAL_SECONDARY_INDEXES = [
    {
        "IndexName": DYNAMODB_NAMESPACE_GSI,
        "KeySchema": [
            {
                "AttributeName": DYNAMODB_COL_NAMESPACE,
                "KeyType": "HASH",
            },
            {
                "AttributeName": DYNAMODB_COL_IDENTIFIER,
                "KeyType": "RANGE",
            },
        ],
        "Projection": {
            "ProjectionType": "KEYS_ONLY",
        },
    }
]


def _get_namespace_properties(namespace_dict: Dict[str, str]) -> Properties:
    return {_remove_property_prefix(key): val for key, val in namespace_dict.items() if key.startswith(PROPERTY_KEY_PREFIX)}


def _convert_dynamo_item_to_regular_dict(dynamo_json: Dict[str, Any]) -> Dict[str, str]:
    """Convert a dynamo json to a regular json.

    Example of a dynamo json:
    {
        "AlbumTitle": {
            "S": "Songs About Life",
        },
        "Artist": {
            "S": "Acme Band",
        },
        "SongTitle": {
            "S": "Happy Day",
        }
    }

    Converted to regular json:
    {
        "AlbumTitle": "Songs About Life",
        "Artist": "Acme Band",
        "SongTitle": "Happy Day"
    }

    Only "S" and "N" data types are supported since those are the only ones that Iceberg is utilizing.
    """
    regular_json = {}
    for column_name, val_dict in dynamo_json.items():
        keys = list(val_dict.keys())

        if len(keys) != 1:
            raise ValueError(f"Expecting only 1 key: {keys}")

        data_type = keys[0]
        if data_type not in ("S", "N"):
            raise ValueError("Only S and N data types are supported.")

        values = list(val_dict.values())
        if len(values) != 1:
            raise ValueError(f"Expecting only 1 value: {values}")

        column_value = values[0]
        regular_json[column_name] = column_value

    return regular_json


def _add_property_prefix(prop: str) -> str:
    return PROPERTY_KEY_PREFIX + prop


def _remove_property_prefix(prop: str) -> str:
    return prop.lstrip(PROPERTY_KEY_PREFIX)
