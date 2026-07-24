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
import json
from typing import TYPE_CHECKING, Any, Union

from google.api_core.exceptions import NotFound
from google.cloud.bigquery import Client, Dataset, DatasetReference, TableReference
from google.cloud.bigquery import Table as BQTable
from google.cloud.bigquery.external_config import ExternalCatalogDatasetOptions, ExternalCatalogTableOptions
from google.cloud.bigquery.schema import SerDeInfo, StorageDescriptor
from google.cloud.exceptions import Conflict
from google.oauth2 import service_account

from pyiceberg.catalog import WAREHOUSE_LOCATION, MetastoreCatalog, PropertiesUpdateSummary
from pyiceberg.exceptions import NamespaceAlreadyExistsError, NoSuchNamespaceError, NoSuchTableError, TableAlreadyExistsError
from pyiceberg.io import load_file_io
from pyiceberg.partitioning import UNPARTITIONED_PARTITION_SPEC, PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.serializers import FromInputFile
from pyiceberg.table import CommitTableResponse, Table
from pyiceberg.table.locations import load_location_provider
from pyiceberg.table.metadata import TableMetadata, new_table_metadata
from pyiceberg.table.snapshots import TOTAL_DATA_FILES, TOTAL_FILE_SIZE, TOTAL_RECORDS
from pyiceberg.table.sorting import UNSORTED_SORT_ORDER, SortOrder
from pyiceberg.table.update import TableRequirement, TableUpdate
from pyiceberg.typedef import EMPTY_DICT, Identifier, Properties
from pyiceberg.utils.config import Config

if TYPE_CHECKING:
    import pyarrow as pa

GCP_PROJECT_ID = "gcp.bigquery.project-id"
GCP_LOCATION = "gcp.bigquery.location"
GCP_CREDENTIALS_FILE = "gcp.bigquery.credential-file"
GCP_CREDENTIALS_INFO = "gcp.bigquery.credentials-info"

METADATA_LOCATION_PROP = "metadata_location"
PREVIOUS_METADATA_LOCATION_PROP = "previous_metadata_location"
TABLE_TYPE_PROP = "table_type"
ICEBERG_TABLE_TYPE_VALUE = "ICEBERG"

HIVE_SERIALIZATION_LIBRARY = "org.apache.iceberg.mr.hive.HiveIcebergSerDe"
HIVE_FILE_INPUT_FORMAT = "org.apache.iceberg.mr.hive.HiveIcebergInputFormat"
HIVE_FILE_OUTPUT_FORMAT = "org.apache.iceberg.mr.hive.HiveIcebergOutputFormat"


class BigQueryMetastoreCatalog(MetastoreCatalog):
    def __init__(self, name: str, **properties: str):
        super().__init__(name, **properties)

        project_id: str | None = self.properties.get(GCP_PROJECT_ID)
        location: str | None = self.properties.get(GCP_LOCATION)
        credentials_file: str | None = self.properties.get(GCP_CREDENTIALS_FILE)
        credentials_info_str: str | None = self.properties.get(GCP_CREDENTIALS_INFO)

        if not project_id:
            raise ValueError(f"Missing property: {GCP_PROJECT_ID}")

        # BigQuery requires current-snapshot-id to be present for tables to be created.
        if not Config().get_bool("legacy-current-snapshot-id"):
            raise ValueError("legacy-current-snapshot-id must be enabled to work with BigQuery.")

        if credentials_file and credentials_info_str:
            raise ValueError("Cannot specify both `gcp.bigquery.credentials-file` and `gcp.bigquery.credentials-info`")

        gcp_credentials = None
        if credentials_file:
            gcp_credentials = service_account.Credentials.from_service_account_file(credentials_file)
        elif credentials_info_str:
            try:
                credentials_info_dict = json.loads(credentials_info_str)
                gcp_credentials = service_account.Credentials.from_service_account_info(credentials_info_dict)
            except json.JSONDecodeError as e:
                raise ValueError(f"Invalid JSON string for {GCP_CREDENTIALS_INFO}: {e}") from e
            except TypeError as e:  # from_service_account_info can raise TypeError for bad structure
                raise ValueError(f"Invalid credentials structure for {GCP_CREDENTIALS_INFO}: {e}") from e

        self.client: Client = Client(
            project=project_id,
            credentials=gcp_credentials,
            location=location,
        )

        self.location = location
        self.project_id = project_id

    def create_table(
        self,
        identifier: str | Identifier,
        schema: Union[Schema, "pa.Schema"],
        location: str | None = None,
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
        schema: Schema = self._convert_schema_if_needed(schema)  # type: ignore

        dataset_name, table_name = self.identifier_to_database_and_table(identifier)

        dataset_ref = DatasetReference(project=self.project_id, dataset_id=dataset_name)
        location = self._resolve_table_location(location, dataset_name, table_name)
        provider = load_location_provider(table_location=location, table_properties=properties)
        metadata_location = provider.new_table_metadata_file_location()

        metadata = new_table_metadata(
            location=location, schema=schema, partition_spec=partition_spec, sort_order=sort_order, properties=properties
        )

        io = load_file_io(properties=self.properties, location=metadata_location)
        self._write_metadata(metadata, io, metadata_location)

        dataset_ref = DatasetReference(project=self.project_id, dataset_id=dataset_name)

        try:
            table = self._make_new_table(
                metadata, metadata_location, TableReference(dataset_ref=dataset_ref, table_id=table_name)
            )
            self.client.create_table(table)
        except Conflict as e:
            raise TableAlreadyExistsError(f"Table {table_name} already exists") from e

        return self.load_table(identifier=identifier)

    def create_namespace(self, namespace: str | Identifier, properties: Properties = EMPTY_DICT) -> None:
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
            dataset_ref = DatasetReference(project=self.project_id, dataset_id=database_name)
            dataset = Dataset(dataset_ref=dataset_ref)
            dataset.external_catalog_dataset_options = self._create_external_catalog_dataset_options(
                self._get_default_warehouse_location_for_dataset(database_name), properties, dataset_ref
            )
            self.client.create_dataset(dataset)
        except Conflict as e:
            raise NamespaceAlreadyExistsError("Namespace {database_name} already exists") from e

    def load_table(self, identifier: str | Identifier) -> Table:
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
        dataset_name, table_name = self.identifier_to_database_and_table(identifier, NoSuchTableError)

        try:
            table_ref = TableReference(
                dataset_ref=DatasetReference(project=self.project_id, dataset_id=dataset_name),
                table_id=table_name,
            )
            table = self.client.get_table(table_ref)
            return self._convert_bigquery_table_to_iceberg_table(identifier, table)
        except NotFound as e:
            raise NoSuchTableError(f"Table does not exist: {dataset_name}.{table_name}") from e

    def drop_table(self, identifier: str | Identifier) -> None:
        """Drop a table.

        Args:
            identifier: Table identifier.

        Raises:
            NoSuchTableError: If a table with the name does not exist, or the identifier is invalid.
        """
        dataset_name, table_name = self.identifier_to_database_and_table(identifier, NoSuchTableError)

        try:
            table_ref = TableReference(
                dataset_ref=DatasetReference(project=self.project_id, dataset_id=dataset_name),
                table_id=table_name,
            )
            self.client.delete_table(table_ref)
        except NoSuchTableError as e:
            raise NoSuchTableError(f"Table does not exist: {dataset_name}.{table_name}") from e

    def commit_table(
        self, table: Table, requirements: tuple[TableRequirement, ...], updates: tuple[TableUpdate, ...]
    ) -> CommitTableResponse:
        raise NotImplementedError

    def rename_table(self, from_identifier: str | Identifier, to_identifier: str | Identifier) -> Table:
        raise NotImplementedError

    def drop_namespace(self, namespace: str | Identifier) -> None:
        database_name = self.identifier_to_database(namespace)

        try:
            dataset_ref = DatasetReference(project=self.project_id, dataset_id=database_name)
            dataset = Dataset(dataset_ref=dataset_ref)
            self.client.delete_dataset(dataset)
        except NotFound as e:
            raise NoSuchNamespaceError(f"Namespace {namespace} does not exist.") from e

    def list_tables(self, namespace: str | Identifier) -> list[Identifier]:
        database_name = self.identifier_to_database(namespace)
        iceberg_tables: list[Identifier] = []
        try:
            dataset_ref = DatasetReference(project=self.project_id, dataset_id=database_name)
            # The list_tables method returns an iterator of TableListItem
            bq_tables_iterator = self.client.list_tables(dataset=dataset_ref)

            for bq_table_list_item in bq_tables_iterator:
                iceberg_tables.append((database_name, bq_table_list_item.table_id))
        except NotFound:
            raise NoSuchNamespaceError(f"Namespace (dataset) '{database_name}' not found.") from None
        return iceberg_tables

    def list_namespaces(self, namespace: str | Identifier = ()) -> list[Identifier]:
        # Since this catalog only supports one-level namespaces, it always returns an empty list unless
        # passed an empty namespace to list all namespaces within the catalog.
        if namespace:
            raise NoSuchNamespaceError(f"Namespace (dataset) '{namespace}' not found.") from None

        # List top-level datasets
        datasets_iterator = self.client.list_datasets()
        return [(dataset.dataset_id,) for dataset in datasets_iterator]

    def register_table(self, identifier: str | Identifier, metadata_location: str) -> Table:
        """Register a new table using existing metadata.

        Args:
            identifier (Union[str, Identifier]): Table identifier for the table
            metadata_location (str): The location to the metadata

        Returns:
            Table: The newly registered table

        Raises:
            TableAlreadyExistsError: If the table already exists
        """
        dataset_name, table_name = self.identifier_to_database_and_table(identifier)

        dataset_ref = DatasetReference(project=self.project_id, dataset_id=dataset_name)

        io = self._load_file_io(location=metadata_location)
        file = io.new_input(metadata_location)
        metadata = FromInputFile.table_metadata(file)

        try:
            table = self._make_new_table(
                metadata, metadata_location, TableReference(dataset_ref=dataset_ref, table_id=table_name)
            )
            self.client.create_table(table)
        except Conflict as e:
            raise TableAlreadyExistsError(f"Table {table_name} already exists") from e

        return self.load_table(identifier=identifier)

    def list_views(self, namespace: str | Identifier) -> list[Identifier]:
        raise NotImplementedError

    def drop_view(self, identifier: str | Identifier) -> None:
        raise NotImplementedError

    def view_exists(self, identifier: str | Identifier) -> bool:
        raise NotImplementedError

    def load_namespace_properties(self, namespace: str | Identifier) -> Properties:
        dataset_name = self.identifier_to_database(namespace)

        try:
            dataset = self.client.get_dataset(DatasetReference(project=self.project_id, dataset_id=dataset_name))

            if dataset and dataset.external_catalog_dataset_options:
                return dataset.external_catalog_dataset_options.to_api_repr()
        except NotFound as e:
            raise NoSuchNamespaceError(f"Namespace {namespace} not found") from e
        return {}

    def update_namespace_properties(
        self, namespace: str | Identifier, removals: set[str] | None = None, updates: Properties = EMPTY_DICT
    ) -> PropertiesUpdateSummary:
        raise NotImplementedError

    def _make_new_table(self, metadata: TableMetadata, metadata_file_location: str, table_ref: TableReference) -> BQTable:
        """To make the table queryable from Hive, the user would likely be setting the HIVE_ENGINE_ENABLED parameter."""
        table = BQTable(table_ref)

        # In Python, you typically set the external data configuration directly.
        # BigQueryMetastoreUtils.create_external_catalog_table_options is mapped to
        # constructing the external_data_configuration for the Table object.
        external_config_options = self._create_external_catalog_table_options(
            metadata.location,
            self._create_table_parameters(metadata_file_location=metadata_file_location, table_metadata=metadata),
        )

        # Apply the external configuration to the Table object.
        # This will depend on the exact structure returned by create_external_catalog_table_options.
        # A common way to set up an external table in BigQuery Python client is:
        table.external_catalog_table_options = external_config_options

        return table

    def _create_external_catalog_table_options(self, location: str, parameters: dict[str, Any]) -> ExternalCatalogTableOptions:
        # This structure directly maps to what BigQuery's ExternalConfig expects for Hive.
        return ExternalCatalogTableOptions(
            storage_descriptor=StorageDescriptor(
                location_uri=location,
                input_format=HIVE_FILE_INPUT_FORMAT,
                output_format=HIVE_FILE_OUTPUT_FORMAT,
                serde_info=SerDeInfo(serialization_library=HIVE_SERIALIZATION_LIBRARY),
            ),
            parameters=parameters,
        )

    def _create_external_catalog_dataset_options(
        self, default_storage_location: str, metadataParameters: dict[str, Any], dataset_ref: DatasetReference
    ) -> ExternalCatalogDatasetOptions:
        return ExternalCatalogDatasetOptions(
            default_storage_location_uri=self._get_default_warehouse_location_for_dataset(dataset_ref.dataset_id),
            parameters=metadataParameters,
        )

    def _convert_bigquery_table_to_iceberg_table(self, identifier: str | Identifier, table: BQTable) -> Table:
        dataset_name, table_name = self.identifier_to_database_and_table(identifier, NoSuchTableError)
        metadata_location = ""
        if table.external_catalog_table_options and table.external_catalog_table_options.parameters:
            metadata_location = table.external_catalog_table_options.parameters[METADATA_LOCATION_PROP]
        io = load_file_io(properties=self.properties, location=metadata_location)
        file = io.new_input(metadata_location)
        metadata = FromInputFile.table_metadata(file)

        return Table(
            identifier=(dataset_name, table_name),
            metadata=metadata,
            metadata_location=metadata_location,
            io=self._load_file_io(metadata.properties, metadata_location),
            catalog=self,
        )

    def _create_table_parameters(self, metadata_file_location: str, table_metadata: TableMetadata) -> dict[str, Any]:
        parameters: dict[str, Any] = table_metadata.properties
        if table_metadata.table_uuid:
            parameters["uuid"] = str(table_metadata.table_uuid)
        parameters[METADATA_LOCATION_PROP] = metadata_file_location
        parameters[TABLE_TYPE_PROP] = ICEBERG_TABLE_TYPE_VALUE
        parameters["EXTERNAL"] = True

        # Add Hive-style basic statistics from snapshot metadata if it exists.
        snapshot = table_metadata.current_snapshot()
        if snapshot:
            summary = snapshot.summary
            if summary:
                if summary.get(TOTAL_DATA_FILES):
                    parameters["numFiles"] = summary.get(TOTAL_DATA_FILES)

                if summary.get(TOTAL_RECORDS):
                    parameters["numRows"] = summary.get(TOTAL_RECORDS)

                if summary.get(TOTAL_FILE_SIZE):
                    parameters["totalSize"] = summary.get(TOTAL_FILE_SIZE)

        return parameters

    def _default_storage_location(self, location: str | None, dataset_ref: DatasetReference) -> str | None:
        if location:
            return location
        dataset = self.client.get_dataset(dataset_ref)
        if dataset and dataset.external_catalog_dataset_options:
            return dataset.external_catalog_dataset_options.default_storage_location_uri

        raise ValueError("Could not find default storage location")

    def _get_default_warehouse_location_for_dataset(self, database_name: str) -> str:
        if warehouse_path := self.properties.get(WAREHOUSE_LOCATION):
            warehouse_path = warehouse_path.rstrip("/")
            return f"{warehouse_path}/{database_name}.db"

        raise ValueError("No default path is set, please specify a location when creating a table")
