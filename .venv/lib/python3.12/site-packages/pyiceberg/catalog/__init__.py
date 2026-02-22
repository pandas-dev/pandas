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

import importlib
import logging
import re
import uuid
from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from typing import (
    TYPE_CHECKING,
    Any,
    cast,
)

from pyiceberg.exceptions import (
    NamespaceAlreadyExistsError,
    NoSuchNamespaceError,
    NoSuchTableError,
    NotInstalledError,
    TableAlreadyExistsError,
)
from pyiceberg.io import FileIO, load_file_io
from pyiceberg.manifest import ManifestFile
from pyiceberg.partitioning import UNPARTITIONED_PARTITION_SPEC, PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.serializers import ToOutputFile
from pyiceberg.table import (
    DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE,
    CommitTableResponse,
    CreateTableTransaction,
    StagedTable,
    Table,
    TableProperties,
)
from pyiceberg.table.locations import load_location_provider
from pyiceberg.table.metadata import TableMetadata, TableMetadataV1, new_table_metadata
from pyiceberg.table.sorting import UNSORTED_SORT_ORDER, SortOrder
from pyiceberg.table.update import (
    TableRequirement,
    TableUpdate,
    update_table_metadata,
)
from pyiceberg.typedef import (
    EMPTY_DICT,
    Identifier,
    Properties,
    RecursiveDict,
    TableVersion,
)
from pyiceberg.utils.config import Config, merge_config
from pyiceberg.utils.properties import property_as_bool

if TYPE_CHECKING:
    import pyarrow as pa

logger = logging.getLogger(__name__)

_ENV_CONFIG = Config()

TOKEN = "token"
TYPE = "type"
PY_CATALOG_IMPL = "py-catalog-impl"
ICEBERG = "iceberg"
TABLE_TYPE = "table_type"
WAREHOUSE_LOCATION = "warehouse"
METADATA_LOCATION = "metadata_location"
PREVIOUS_METADATA_LOCATION = "previous_metadata_location"
MANIFEST = "manifest"
MANIFEST_LIST = "manifest list"
PREVIOUS_METADATA = "previous metadata"
METADATA = "metadata"
URI = "uri"
LOCATION = "location"
EXTERNAL_TABLE = "EXTERNAL_TABLE"
BOTOCORE_SESSION = "botocore_session"

TABLE_METADATA_FILE_NAME_REGEX = re.compile(
    r"""
    (\d+)              # version number
    -                  # separator
    ([\w-]{36})        # UUID (36 characters, including hyphens)
    (?:\.\w+)?         # optional codec name
    \.metadata\.json   # file extension
    """,
    re.X,
)


class CatalogType(Enum):
    REST = "rest"
    HIVE = "hive"
    GLUE = "glue"
    DYNAMODB = "dynamodb"
    SQL = "sql"
    IN_MEMORY = "in-memory"
    BIGQUERY = "bigquery"


def load_rest(name: str, conf: Properties) -> Catalog:
    from pyiceberg.catalog.rest import RestCatalog

    return RestCatalog(name, **conf)


def load_hive(name: str, conf: Properties) -> Catalog:
    try:
        from pyiceberg.catalog.hive import HiveCatalog

        return HiveCatalog(name, **conf)
    except ImportError as exc:
        raise NotInstalledError("Apache Hive support not installed: pip install 'pyiceberg[hive]'") from exc


def load_glue(name: str, conf: Properties) -> Catalog:
    try:
        from pyiceberg.catalog.glue import GlueCatalog

        return GlueCatalog(name, **conf)
    except ImportError as exc:
        raise NotInstalledError("AWS glue support not installed: pip install 'pyiceberg[glue]'") from exc


def load_dynamodb(name: str, conf: Properties) -> Catalog:
    try:
        from pyiceberg.catalog.dynamodb import DynamoDbCatalog

        return DynamoDbCatalog(name, **conf)
    except ImportError as exc:
        raise NotInstalledError("AWS DynamoDB support not installed: pip install 'pyiceberg[dynamodb]'") from exc


def load_sql(name: str, conf: Properties) -> Catalog:
    try:
        from pyiceberg.catalog.sql import SqlCatalog

        return SqlCatalog(name, **conf)
    except ImportError as exc:
        raise NotInstalledError(
            "SQLAlchemy support not installed: pip install 'pyiceberg[sql-postgres]' or pip install 'pyiceberg[sql-sqlite]'"
        ) from exc


def load_in_memory(name: str, conf: Properties) -> Catalog:
    try:
        from pyiceberg.catalog.memory import InMemoryCatalog

        return InMemoryCatalog(name, **conf)
    except ImportError as exc:
        raise NotInstalledError("SQLAlchemy support not installed: pip install 'pyiceberg[sql-sqlite]'") from exc


def load_bigquery(name: str, conf: Properties) -> Catalog:
    try:
        from pyiceberg.catalog.bigquery_metastore import BigQueryMetastoreCatalog

        return BigQueryMetastoreCatalog(name, **conf)
    except ImportError as exc:
        raise NotInstalledError("BigQuery support not installed: pip install 'pyiceberg[bigquery]'") from exc


AVAILABLE_CATALOGS: dict[CatalogType, Callable[[str, Properties], Catalog]] = {
    CatalogType.REST: load_rest,
    CatalogType.HIVE: load_hive,
    CatalogType.GLUE: load_glue,
    CatalogType.DYNAMODB: load_dynamodb,
    CatalogType.SQL: load_sql,
    CatalogType.IN_MEMORY: load_in_memory,
    CatalogType.BIGQUERY: load_bigquery,
}


def infer_catalog_type(name: str, catalog_properties: RecursiveDict) -> CatalogType | None:
    """Try to infer the type based on the dict.

    Args:
        name: Name of the catalog.
        catalog_properties: Catalog properties.

    Returns:
        The inferred type based on the provided properties.

    Raises:
        ValueError: Raises a ValueError in case properties are missing, or the wrong type.
    """
    if uri := catalog_properties.get(URI):
        if isinstance(uri, str):
            if uri.startswith("http"):
                return CatalogType.REST
            elif uri.startswith("thrift"):
                return CatalogType.HIVE
            elif uri.startswith(("sqlite", "postgresql")):
                return CatalogType.SQL
            else:
                raise ValueError(f"Could not infer the catalog type from the uri: {uri}")
        else:
            raise ValueError(f"Expects the URI to be a string, got: {type(uri)}")
    raise ValueError(
        f"URI missing, please provide using --uri, the config or environment variable PYICEBERG_CATALOG__{name.upper()}__URI"
    )


def load_catalog(name: str | None = None, **properties: str | None) -> Catalog:
    """Load the catalog based on the properties.

    Will look up the properties from the config, based on the name.

    Args:
        name: The name of the catalog.
        properties: The properties that are used next to the configuration.

    Returns:
        An initialized Catalog.

    Raises:
        ValueError: Raises a ValueError in case properties are missing or malformed,
            or if it could not determine the catalog based on the properties.
    """
    if name is None:
        name = _ENV_CONFIG.get_default_catalog_name()

    env = _ENV_CONFIG.get_catalog_config(name)
    conf: RecursiveDict = merge_config(env or {}, cast(RecursiveDict, properties))

    catalog_type: CatalogType | None
    provided_catalog_type = conf.get(TYPE)

    if catalog_impl := properties.get(PY_CATALOG_IMPL):
        if provided_catalog_type:
            raise ValueError(
                "Must not set both catalog type and py-catalog-impl configurations, "
                f"but found type {provided_catalog_type} and py-catalog-impl {catalog_impl}"
            )

        if catalog := _import_catalog(name, catalog_impl, properties):
            logger.info("Loaded Catalog: %s", catalog_impl)
            return catalog
        else:
            raise ValueError(f"Could not initialize Catalog: {catalog_impl}")

    catalog_type = None
    if provided_catalog_type and isinstance(provided_catalog_type, str):
        catalog_type = CatalogType(provided_catalog_type.lower())
    elif not provided_catalog_type:
        catalog_type = infer_catalog_type(name, conf)

    if catalog_type:
        return AVAILABLE_CATALOGS[catalog_type](name, cast(dict[str, str], conf))

    raise ValueError(f"Could not initialize catalog with the following properties: {properties}")


def list_catalogs() -> list[str]:
    return _ENV_CONFIG.get_known_catalogs()


def delete_files(io: FileIO, files_to_delete: set[str], file_type: str) -> None:
    """Delete files.

    Log warnings if failing to delete any file.

    Args:
        io: The FileIO used to delete the object.
        files_to_delete: A set of file paths to be deleted.
        file_type: The type of the file.
    """
    for file in files_to_delete:
        try:
            io.delete(file)
        except OSError:
            logger.warning(f"Failed to delete {file_type} file {file}", exc_info=logger.isEnabledFor(logging.DEBUG))


def delete_data_files(io: FileIO, manifests_to_delete: list[ManifestFile]) -> None:
    """Delete data files linked to given manifests.

    Log warnings if failing to delete any file.

    Args:
        io: The FileIO used to delete the object.
        manifests_to_delete: A list of manifest contains paths of data files to be deleted.
    """
    deleted_files: dict[str, bool] = {}
    for manifest_file in manifests_to_delete:
        for entry in manifest_file.fetch_manifest_entry(io, discard_deleted=False):
            path = entry.data_file.file_path
            if not deleted_files.get(path, False):
                try:
                    io.delete(path)
                except OSError:
                    logger.warning(f"Failed to delete data file {path}", exc_info=logger.isEnabledFor(logging.DEBUG))
                deleted_files[path] = True


def _import_catalog(name: str, catalog_impl: str, properties: Properties) -> Catalog | None:
    try:
        path_parts = catalog_impl.split(".")
        if len(path_parts) < 2:
            raise ValueError(f"py-catalog-impl should be full path (module.CustomCatalog), got: {catalog_impl}")
        module_name, class_name = ".".join(path_parts[:-1]), path_parts[-1]
        module = importlib.import_module(module_name)
        class_ = getattr(module, class_name)
        return class_(name, **properties)
    except ModuleNotFoundError:
        logger.warning(f"Could not initialize Catalog: {catalog_impl}", exc_info=logger.isEnabledFor(logging.DEBUG))
        return None


@dataclass
class PropertiesUpdateSummary:
    removed: list[str]
    updated: list[str]
    missing: list[str]


class Catalog(ABC):
    """Base Catalog for table operations like - create, drop, load, list and others.

    The catalog table APIs accept a table identifier, which is fully classified table name. The identifier can be a string or
    tuple of strings. If the identifier is a string, it is split into a tuple on '.'. If it is a tuple, it is used as-is.

    The catalog namespace APIs follow a similar convention wherein they also accept a namespace identifier that can be a string
    or tuple of strings.

    Attributes:
        name (str): Name of the catalog.
        properties (Properties): Catalog properties.
    """

    name: str
    properties: Properties

    def __init__(self, name: str, **properties: str):
        self.name = name
        self.properties = properties

    @abstractmethod
    def create_table(
        self,
        identifier: str | Identifier,
        schema: Schema | pa.Schema,
        location: str | None = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
    ) -> Table:
        """Create a table.

        Args:
            identifier (str | Identifier): Table identifier.
            schema (Schema): Table's schema.
            location (str | None): Location for the table. Optional Argument.
            partition_spec (PartitionSpec): PartitionSpec for the table.
            sort_order (SortOrder): SortOrder for the table.
            properties (Properties): Table properties that can be a string based dictionary.

        Returns:
            Table: the created table instance.

        Raises:
            TableAlreadyExistsError: If a table with the name already exists.
        """

    @abstractmethod
    def create_table_transaction(
        self,
        identifier: str | Identifier,
        schema: Schema | pa.Schema,
        location: str | None = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
    ) -> CreateTableTransaction:
        """Create a CreateTableTransaction.

        Args:
            identifier (str | Identifier): Table identifier.
            schema (Schema): Table's schema.
            location (str | None): Location for the table. Optional Argument.
            partition_spec (PartitionSpec): PartitionSpec for the table.
            sort_order (SortOrder): SortOrder for the table.
            properties (Properties): Table properties that can be a string based dictionary.

        Returns:
            CreateTableTransaction: createTableTransaction instance.
        """

    def create_table_if_not_exists(
        self,
        identifier: str | Identifier,
        schema: Schema | pa.Schema,
        location: str | None = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
    ) -> Table:
        """Create a table if it does not exist.

        Args:
            identifier (str | Identifier): Table identifier.
            schema (Schema): Table's schema.
            location (str | None): Location for the table. Optional Argument.
            partition_spec (PartitionSpec): PartitionSpec for the table.
            sort_order (SortOrder): SortOrder for the table.
            properties (Properties): Table properties that can be a string based dictionary.

        Returns:
            Table: the created table instance if the table does not exist, else the existing
            table instance.
        """
        try:
            return self.create_table(identifier, schema, location, partition_spec, sort_order, properties)
        except TableAlreadyExistsError:
            return self.load_table(identifier)

    @abstractmethod
    def load_table(self, identifier: str | Identifier) -> Table:
        """Load the table's metadata and returns the table instance.

        You can also use this method to check for table existence using 'try catalog.table() except NoSuchTableError'.
        Note: This method doesn't scan data stored in the table.

        Args:
            identifier (str | Identifier): Table identifier.

        Returns:
            Table: the table instance with its metadata.

        Raises:
            NoSuchTableError: If a table with the name does not exist.
        """

    @abstractmethod
    def table_exists(self, identifier: str | Identifier) -> bool:
        """Check if a table exists.

        Args:
            identifier (str | Identifier): Table identifier.

        Returns:
            bool: True if the table exists, False otherwise.
        """

    @abstractmethod
    def view_exists(self, identifier: str | Identifier) -> bool:
        """Check if a view exists.

        Args:
            identifier (str | Identifier): View identifier.

        Returns:
            bool: True if the view exists, False otherwise.
        """

    @abstractmethod
    def namespace_exists(self, namespace: str | Identifier) -> bool:
        """Check if a namespace exists.

        Args:
            namespace (str | Identifier): Namespace identifier.

        Returns:
            bool: True if the namespace exists, False otherwise.
        """

    @abstractmethod
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

    @abstractmethod
    def drop_table(self, identifier: str | Identifier) -> None:
        """Drop a table.

        Args:
            identifier (str | Identifier): Table identifier.

        Raises:
            NoSuchTableError: If a table with the name does not exist.
        """

    @abstractmethod
    def purge_table(self, identifier: str | Identifier) -> None:
        """Drop a table and purge all data and metadata files.

        Note: This method only logs warning rather than raise exception when encountering file deletion failure.

        Args:
            identifier (str | Identifier): Table identifier.

        Raises:
            NoSuchTableError: If a table with the name does not exist, or the identifier is invalid.
        """

    @abstractmethod
    def rename_table(self, from_identifier: str | Identifier, to_identifier: str | Identifier) -> Table:
        """Rename a fully classified table name.

        Args:
            from_identifier (str | Identifier): Existing table identifier.
            to_identifier (str | Identifier): New table identifier.

        Returns:
            Table: the updated table instance with its metadata.

        Raises:
            NoSuchTableError: If a table with the name does not exist.
        """

    @abstractmethod
    def commit_table(
        self, table: Table, requirements: tuple[TableRequirement, ...], updates: tuple[TableUpdate, ...]
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
            CommitStateUnknownException: Failed due to an internal exception on the side of the catalog.
        """

    @abstractmethod
    def create_namespace(self, namespace: str | Identifier, properties: Properties = EMPTY_DICT) -> None:
        """Create a namespace in the catalog.

        Args:
            namespace (str | Identifier): Namespace identifier.
            properties (Properties): A string dictionary of properties for the given namespace.

        Raises:
            NamespaceAlreadyExistsError: If a namespace with the given name already exists.
        """

    def create_namespace_if_not_exists(self, namespace: str | Identifier, properties: Properties = EMPTY_DICT) -> None:
        """Create a namespace if it does not exist.

        Args:
            namespace (str | Identifier): Namespace identifier.
            properties (Properties): A string dictionary of properties for the given namespace.
        """
        try:
            self.create_namespace(namespace, properties)
        except NamespaceAlreadyExistsError:
            pass

    @abstractmethod
    def drop_namespace(self, namespace: str | Identifier) -> None:
        """Drop a namespace.

        Args:
            namespace (str | Identifier): Namespace identifier.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
            NamespaceNotEmptyError: If the namespace is not empty.
        """

    @abstractmethod
    def list_tables(self, namespace: str | Identifier) -> list[Identifier]:
        """List tables under the given namespace in the catalog.

        Args:
            namespace (str | Identifier): Namespace identifier to search.

        Returns:
            List[Identifier]: list of table identifiers.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
        """

    @abstractmethod
    def list_namespaces(self, namespace: str | Identifier = ()) -> list[Identifier]:
        """List namespaces from the given namespace. If not given, list top-level namespaces from the catalog.

        Args:
            namespace (str | Identifier): Namespace identifier to search.

        Returns:
            List[Identifier]: a List of namespace identifiers.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
        """

    @abstractmethod
    def list_views(self, namespace: str | Identifier) -> list[Identifier]:
        """List views under the given namespace in the catalog.

        Args:
            namespace (str | Identifier): Namespace identifier to search.

        Returns:
            List[Identifier]: list of table identifiers.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
        """

    @abstractmethod
    def load_namespace_properties(self, namespace: str | Identifier) -> Properties:
        """Get properties for a namespace.

        Args:
            namespace (str | Identifier): Namespace identifier.

        Returns:
            Properties: Properties for the given namespace.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
        """

    @abstractmethod
    def update_namespace_properties(
        self, namespace: str | Identifier, removals: set[str] | None = None, updates: Properties = EMPTY_DICT
    ) -> PropertiesUpdateSummary:
        """Remove provided property keys and updates properties for a namespace.

        Args:
            namespace (str | Identifier): Namespace identifier.
            removals (Set[str]): Set of property keys that need to be removed. Optional Argument.
            updates (Properties): Properties to be updated for the given namespace.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
            ValueError: If removals and updates have overlapping keys.
        """

    @abstractmethod
    def drop_view(self, identifier: str | Identifier) -> None:
        """Drop a view.

        Args:
            identifier (str | Identifier): View identifier.

        Raises:
            NoSuchViewError: If a view with the given name does not exist.
        """

    @staticmethod
    def identifier_to_tuple(identifier: str | Identifier) -> Identifier:
        """Parse an identifier to a tuple.

        If the identifier is a string, it is split into a tuple on '.'. If it is a tuple, it is used as-is.

        Args:
            identifier (str | Identifier): an identifier, either a string or tuple of strings.

        Returns:
            Identifier: a tuple of strings.
        """
        return identifier if isinstance(identifier, tuple) else tuple(str.split(identifier, "."))

    @staticmethod
    def table_name_from(identifier: str | Identifier) -> str:
        """Extract table name from a table identifier.

        Args:
            identifier (str | Identifier): a table identifier.

        Returns:
            str: Table name.
        """
        return Catalog.identifier_to_tuple(identifier)[-1]

    @staticmethod
    def namespace_from(identifier: str | Identifier) -> Identifier:
        """Extract table namespace from a table identifier.

        Args:
            identifier (Union[str, Identifier]): a table identifier.

        Returns:
            Identifier: Namespace identifier.
        """
        return Catalog.identifier_to_tuple(identifier)[:-1]

    @staticmethod
    def namespace_to_string(identifier: str | Identifier, err: type[ValueError] | type[NoSuchNamespaceError] = ValueError) -> str:
        """Transform a namespace identifier into a string.

        Args:
            identifier (Union[str, Identifier]): a namespace identifier.
            err (Union[Type[ValueError], Type[NoSuchNamespaceError]]): the error type to raise when identifier is empty.

        Returns:
            Identifier: Namespace identifier.
        """
        tuple_identifier = Catalog.identifier_to_tuple(identifier)
        if len(tuple_identifier) < 1:
            raise err("Empty namespace identifier")

        # Check if any segment of the tuple is an empty string
        if any(segment.strip() == "" for segment in tuple_identifier):
            raise err("Namespace identifier contains an empty segment or a segment with only whitespace")

        return ".".join(segment.strip() for segment in tuple_identifier)

    def supports_server_side_planning(self) -> bool:
        """Check if the catalog supports server-side scan planning."""
        return False

    @staticmethod
    def identifier_to_database(
        identifier: str | Identifier, err: type[ValueError] | type[NoSuchNamespaceError] = ValueError
    ) -> str:
        tuple_identifier = Catalog.identifier_to_tuple(identifier)
        if len(tuple_identifier) != 1:
            raise err(f"Invalid database, hierarchical namespaces are not supported: {identifier}")

        return tuple_identifier[0]

    @staticmethod
    def identifier_to_database_and_table(
        identifier: str | Identifier,
        err: type[ValueError] | type[NoSuchTableError] | type[NoSuchNamespaceError] = ValueError,
    ) -> tuple[str, str]:
        tuple_identifier = Catalog.identifier_to_tuple(identifier)
        if len(tuple_identifier) != 2:
            raise err(f"Invalid path, hierarchical namespaces are not supported: {identifier}")

        return tuple_identifier[0], tuple_identifier[1]

    def _load_file_io(self, properties: Properties = EMPTY_DICT, location: str | None = None) -> FileIO:
        return load_file_io({**self.properties, **properties}, location)

    @staticmethod
    def _convert_schema_if_needed(
        schema: Schema | pa.Schema, format_version: TableVersion = TableProperties.DEFAULT_FORMAT_VERSION
    ) -> Schema:
        if isinstance(schema, Schema):
            return schema
        try:
            import pyarrow as pa

            from pyiceberg.io.pyarrow import _ConvertToIcebergWithoutIDs, visit_pyarrow

            downcast_ns_timestamp_to_us = Config().get_bool(DOWNCAST_NS_TIMESTAMP_TO_US_ON_WRITE) or False
            if isinstance(schema, pa.Schema):
                schema: Schema = visit_pyarrow(  # type: ignore
                    schema,
                    _ConvertToIcebergWithoutIDs(
                        downcast_ns_timestamp_to_us=downcast_ns_timestamp_to_us, format_version=format_version
                    ),
                )
                return schema
        except ModuleNotFoundError:
            pass
        raise ValueError(f"{type(schema)=}, but it must be pyiceberg.schema.Schema or pyarrow.Schema")

    @staticmethod
    def _delete_old_metadata(io: FileIO, base: TableMetadata, metadata: TableMetadata) -> None:
        """Delete oldest metadata if config is set to true."""
        delete_after_commit: bool = property_as_bool(
            metadata.properties,
            TableProperties.METADATA_DELETE_AFTER_COMMIT_ENABLED,
            TableProperties.METADATA_DELETE_AFTER_COMMIT_ENABLED_DEFAULT,
        )

        if delete_after_commit:
            removed_previous_metadata_files: set[str] = {log.metadata_file for log in base.metadata_log}
            current_metadata_files: set[str] = {log.metadata_file for log in metadata.metadata_log}
            removed_previous_metadata_files.difference_update(current_metadata_files)
            delete_files(io, removed_previous_metadata_files, METADATA)

    def close(self) -> None:  # noqa: B027
        """Close the catalog and release any resources.

        This method should be called when the catalog is no longer needed to ensure
        proper cleanup of resources like database connections, file handles, etc.

        Default implementation does nothing. Override in subclasses that need cleanup.
        """

    def __enter__(self) -> Catalog:
        """Enter the context manager.

        Returns:
            Catalog: The catalog instance.
        """
        return self

    def __exit__(self, exc_type: type | None, exc_val: BaseException | None, exc_tb: Any | None) -> None:
        """Exit the context manager and close the catalog.

        Args:
            exc_type: Exception type if an exception occurred.
            exc_val: Exception value if an exception occurred.
            exc_tb: Exception traceback if an exception occurred.
        """
        self.close()

    def __repr__(self) -> str:
        """Return the string representation of the Catalog class."""
        return f"{self.name} ({self.__class__})"


class MetastoreCatalog(Catalog, ABC):
    def __init__(self, name: str, **properties: str):
        super().__init__(name, **properties)

    def create_table_transaction(
        self,
        identifier: str | Identifier,
        schema: Schema | pa.Schema,
        location: str | None = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
    ) -> CreateTableTransaction:
        return CreateTableTransaction(
            self._create_staged_table(identifier, schema, location, partition_spec, sort_order, properties)
        )

    def table_exists(self, identifier: str | Identifier) -> bool:
        try:
            self.load_table(identifier)
            return True
        except NoSuchTableError:
            return False

    def namespace_exists(self, namespace: str | Identifier) -> bool:
        """Check if a namespace exists.

        Args:
            namespace (str | Identifier): Namespace identifier.

        Returns:
            bool: True if the namespace exists, False otherwise.
        """
        try:
            self.load_namespace_properties(namespace)
            return True
        except NoSuchNamespaceError:
            return False

    def purge_table(self, identifier: str | Identifier) -> None:
        table = self.load_table(identifier)
        self.drop_table(identifier)
        io = load_file_io(self.properties, table.metadata_location)
        metadata = table.metadata
        manifest_lists_to_delete = set()
        manifests_to_delete: list[ManifestFile] = []
        for snapshot in metadata.snapshots:
            manifests_to_delete += snapshot.manifests(io)
            manifest_lists_to_delete.add(snapshot.manifest_list)

        manifest_paths_to_delete = {manifest.manifest_path for manifest in manifests_to_delete}
        prev_metadata_files = {log.metadata_file for log in metadata.metadata_log}

        delete_data_files(io, manifests_to_delete)
        delete_files(io, manifest_paths_to_delete, MANIFEST)
        delete_files(io, manifest_lists_to_delete, MANIFEST_LIST)
        delete_files(io, prev_metadata_files, PREVIOUS_METADATA)
        delete_files(io, {table.metadata_location}, METADATA)

    def _create_staged_table(
        self,
        identifier: str | Identifier,
        schema: Schema | pa.Schema,
        location: str | None = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
    ) -> StagedTable:
        """Create a table and return the table instance without committing the changes.

        Args:
            identifier (str | Identifier): Table identifier.
            schema (Schema): Table's schema.
            location (str | None): Location for the table. Optional Argument.
            partition_spec (PartitionSpec): PartitionSpec for the table.
            sort_order (SortOrder): SortOrder for the table.
            properties (Properties): Table properties that can be a string based dictionary.

        Returns:
            StagedTable: the created staged table instance.
        """
        schema: Schema = self._convert_schema_if_needed(  # type: ignore
            schema,
            int(properties.get(TableProperties.FORMAT_VERSION, TableProperties.DEFAULT_FORMAT_VERSION)),  # type: ignore
        )

        database_name, table_name = self.identifier_to_database_and_table(identifier)

        location = self._resolve_table_location(location, database_name, table_name)
        provider = load_location_provider(location, properties)
        metadata_location = provider.new_table_metadata_file_location()
        metadata = new_table_metadata(
            location=location, schema=schema, partition_spec=partition_spec, sort_order=sort_order, properties=properties
        )
        io = self._load_file_io(properties=properties, location=metadata_location)
        return StagedTable(
            identifier=(database_name, table_name),
            metadata=metadata,
            metadata_location=metadata_location,
            io=io,
            catalog=self,
        )

    def _update_and_stage_table(
        self,
        current_table: Table | None,
        table_identifier: Identifier,
        requirements: tuple[TableRequirement, ...],
        updates: tuple[TableUpdate, ...],
    ) -> StagedTable:
        for requirement in requirements:
            requirement.validate(current_table.metadata if current_table else None)

        updated_metadata = update_table_metadata(
            base_metadata=current_table.metadata if current_table else self._empty_table_metadata(),
            updates=updates,
            enforce_validation=current_table is None,
            metadata_location=current_table.metadata_location if current_table else None,
        )

        new_metadata_version = self._parse_metadata_version(current_table.metadata_location) + 1 if current_table else 0
        provider = load_location_provider(updated_metadata.location, updated_metadata.properties)
        new_metadata_location = provider.new_table_metadata_file_location(new_metadata_version)

        return StagedTable(
            identifier=table_identifier,
            metadata=updated_metadata,
            metadata_location=new_metadata_location,
            io=self._load_file_io(properties=updated_metadata.properties, location=new_metadata_location),
            catalog=self,
        )

    def _get_updated_props_and_update_summary(
        self, current_properties: Properties, removals: set[str] | None, updates: Properties
    ) -> tuple[PropertiesUpdateSummary, Properties]:
        self._check_for_overlap(updates=updates, removals=removals)
        updated_properties = dict(current_properties)

        removed: set[str] = set()
        updated: set[str] = set()

        if removals:
            for key in removals:
                if key in updated_properties:
                    updated_properties.pop(key)
                    removed.add(key)
        if updates:
            for key, value in updates.items():
                updated_properties[key] = value
                updated.add(key)

        expected_to_change = (removals or set()).difference(removed)
        properties_update_summary = PropertiesUpdateSummary(
            removed=list(removed or []), updated=list(updated or []), missing=list(expected_to_change)
        )

        return properties_update_summary, updated_properties

    def _resolve_table_location(self, location: str | None, database_name: str, table_name: str) -> str:
        if not location:
            return self._get_default_warehouse_location(database_name, table_name)
        return location.rstrip("/")

    def _get_default_warehouse_location(self, database_name: str, table_name: str) -> str:
        """Return the default warehouse location using the convention of `warehousePath/databaseName/tableName`."""
        database_properties = self.load_namespace_properties(database_name)
        if database_location := database_properties.get(LOCATION):
            database_location = database_location.rstrip("/")
            return f"{database_location}/{table_name}"

        if warehouse_path := self.properties.get(WAREHOUSE_LOCATION):
            warehouse_path = warehouse_path.rstrip("/")
            return f"{warehouse_path}/{database_name}/{table_name}"

        raise ValueError("No default path is set, please specify a location when creating a table")

    def _get_hive_style_warehouse_location(self, database_name: str, table_name: str) -> str:
        """Return the default warehouse location following the Hive convention of `warehousePath/databaseName.db/tableName`."""
        database_properties = self.load_namespace_properties(database_name)
        if database_location := database_properties.get(LOCATION):
            database_location = database_location.rstrip("/")
            return f"{database_location}/{table_name}"

        if warehouse_path := self.properties.get(WAREHOUSE_LOCATION):
            warehouse_path = warehouse_path.rstrip("/")
            return f"{warehouse_path}/{database_name}.db/{table_name}"

        raise ValueError("No default path is set, please specify a location when creating a table")

    @staticmethod
    def _write_metadata(metadata: TableMetadata, io: FileIO, metadata_path: str) -> None:
        ToOutputFile.table_metadata(metadata, io.new_output(metadata_path))

    @staticmethod
    def _parse_metadata_version(metadata_location: str) -> int:
        """Parse the version from the metadata location.

        The version is the first part of the file name, before the first dash.
        For example, the version of the metadata file
        `s3://bucket/db/tb/metadata/00001-6c97e413-d51b-4538-ac70-12fe2a85cb83.metadata.json`
        is 1.
        If the path does not comply with the pattern, the version is defaulted to be -1, ensuring
        that the next metadata file is treated as having version 0.

        Args:
            metadata_location (str): The location of the metadata file.

        Returns:
            int: The version of the metadata file. -1 if the file name does not have valid version string
        """
        file_name = metadata_location.split("/")[-1]
        if file_name_match := TABLE_METADATA_FILE_NAME_REGEX.fullmatch(file_name):
            try:
                uuid.UUID(file_name_match.group(2))
            except ValueError:
                return -1
            return int(file_name_match.group(1))
        else:
            return -1

    @staticmethod
    def _check_for_overlap(removals: set[str] | None, updates: Properties) -> None:
        if updates and removals:
            overlap = set(removals) & set(updates.keys())
            if overlap:
                raise ValueError(f"Updates and deletes have an overlap: {overlap}")

    @staticmethod
    def _empty_table_metadata() -> TableMetadata:
        """Return an empty TableMetadata instance.

        It is used to build a TableMetadata from a sequence of initial TableUpdates.
        It is a V1 TableMetadata because there will be a UpgradeFormatVersionUpdate in
        initial changes to bump the metadata to the target version.

        Returns:
            TableMetadata: An empty TableMetadata instance.
        """
        return TableMetadataV1.model_construct(last_column_id=-1, schema=Schema())
