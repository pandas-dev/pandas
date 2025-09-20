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

from typing import (
    TYPE_CHECKING,
    List,
    Optional,
    Set,
    Tuple,
    Union,
)

from sqlalchemy import (
    String,
    create_engine,
    delete,
    insert,
    select,
    union,
    update,
)
from sqlalchemy.exc import IntegrityError, NoResultFound, OperationalError, ProgrammingError
from sqlalchemy.orm import (
    DeclarativeBase,
    Mapped,
    MappedAsDataclass,
    Session,
    mapped_column,
)

from pyiceberg.catalog import (
    METADATA_LOCATION,
    URI,
    Catalog,
    MetastoreCatalog,
    PropertiesUpdateSummary,
)
from pyiceberg.exceptions import (
    CommitFailedException,
    NamespaceAlreadyExistsError,
    NamespaceNotEmptyError,
    NoSuchNamespaceError,
    NoSuchPropertyException,
    NoSuchTableError,
    TableAlreadyExistsError,
)
from pyiceberg.io import load_file_io
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
from pyiceberg.types import strtobool

if TYPE_CHECKING:
    import pyarrow as pa

DEFAULT_ECHO_VALUE = "false"
DEFAULT_POOL_PRE_PING_VALUE = "false"
DEFAULT_INIT_CATALOG_TABLES = "true"


class SqlCatalogBaseTable(MappedAsDataclass, DeclarativeBase):
    pass


class IcebergTables(SqlCatalogBaseTable):
    __tablename__ = "iceberg_tables"

    catalog_name: Mapped[str] = mapped_column(String(255), nullable=False, primary_key=True)
    table_namespace: Mapped[str] = mapped_column(String(255), nullable=False, primary_key=True)
    table_name: Mapped[str] = mapped_column(String(255), nullable=False, primary_key=True)
    metadata_location: Mapped[Optional[str]] = mapped_column(String(1000), nullable=True)
    previous_metadata_location: Mapped[Optional[str]] = mapped_column(String(1000), nullable=True)


class IcebergNamespaceProperties(SqlCatalogBaseTable):
    __tablename__ = "iceberg_namespace_properties"
    # Catalog minimum Namespace Properties
    NAMESPACE_MINIMAL_PROPERTIES = {"exists": "true"}

    catalog_name: Mapped[str] = mapped_column(String(255), nullable=False, primary_key=True)
    namespace: Mapped[str] = mapped_column(String(255), nullable=False, primary_key=True)
    property_key: Mapped[str] = mapped_column(String(255), nullable=False, primary_key=True)
    property_value: Mapped[str] = mapped_column(String(1000), nullable=False)


class SqlCatalog(MetastoreCatalog):
    """Implementation of a SQL based catalog.

    In the `JDBCCatalog` implementation, a `Namespace` is composed of a list of strings separated by dots: `'ns1.ns2.ns3'`.
    And you can have as many levels as you want, but you need at least one.  The `SqlCatalog` honors the same convention.

    In the `JDBCCatalog` implementation, a `TableIdentifier` is composed of an optional `Namespace` and a table name.
    When a `Namespace` is present, the full name will be `'ns1.ns2.ns3.table'`.  A valid `TableIdentifier` could be `'name'` (no namespace).
    The `SqlCatalog` has a different convention where a `TableIdentifier` requires a `Namespace`.
    """

    def __init__(self, name: str, **properties: str):
        super().__init__(name, **properties)

        if not (uri_prop := self.properties.get(URI)):
            raise NoSuchPropertyException("SQL connection URI is required")

        echo_str = str(self.properties.get("echo", DEFAULT_ECHO_VALUE)).lower()
        echo = strtobool(echo_str) if echo_str != "debug" else "debug"
        pool_pre_ping = strtobool(self.properties.get("pool_pre_ping", DEFAULT_POOL_PRE_PING_VALUE))
        init_catalog_tables = strtobool(self.properties.get("init_catalog_tables", DEFAULT_INIT_CATALOG_TABLES))

        self.engine = create_engine(uri_prop, echo=echo, pool_pre_ping=pool_pre_ping)

        if init_catalog_tables:
            self._ensure_tables_exist()

    def _ensure_tables_exist(self) -> None:
        with Session(self.engine) as session:
            for table in [IcebergTables, IcebergNamespaceProperties]:
                stmt = select(1).select_from(table)
                try:
                    session.scalar(stmt)
                except (
                    OperationalError,
                    ProgrammingError,
                ):  # sqlalchemy returns OperationalError in case of sqlite and ProgrammingError with postgres.
                    self.create_tables()
                    return

    def create_tables(self) -> None:
        SqlCatalogBaseTable.metadata.create_all(self.engine)

    def destroy_tables(self) -> None:
        SqlCatalogBaseTable.metadata.drop_all(self.engine)

    def _convert_orm_to_iceberg(self, orm_table: IcebergTables) -> Table:
        # Check for expected properties.
        if not (metadata_location := orm_table.metadata_location):
            raise NoSuchTableError(f"Table property {METADATA_LOCATION} is missing")
        if not (table_namespace := orm_table.table_namespace):
            raise NoSuchTableError(f"Table property {IcebergTables.table_namespace} is missing")
        if not (table_name := orm_table.table_name):
            raise NoSuchTableError(f"Table property {IcebergTables.table_name} is missing")

        io = load_file_io(properties=self.properties, location=metadata_location)
        file = io.new_input(metadata_location)
        metadata = FromInputFile.table_metadata(file)
        return Table(
            identifier=Catalog.identifier_to_tuple(table_namespace) + (table_name,),
            metadata=metadata,
            metadata_location=metadata_location,
            io=self._load_file_io(metadata.properties, metadata_location),
            catalog=self,
        )

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

        namespace_identifier = Catalog.namespace_from(identifier)
        table_name = Catalog.table_name_from(identifier)
        if not self._namespace_exists(namespace_identifier):
            raise NoSuchNamespaceError(f"Namespace does not exist: {namespace_identifier}")

        namespace = Catalog.namespace_to_string(namespace_identifier)
        location = self._resolve_table_location(location, namespace, table_name)
        location_provider = load_location_provider(table_location=location, table_properties=properties)
        metadata_location = location_provider.new_table_metadata_file_location()
        metadata = new_table_metadata(
            location=location, schema=schema, partition_spec=partition_spec, sort_order=sort_order, properties=properties
        )
        io = load_file_io(properties=self.properties, location=metadata_location)
        self._write_metadata(metadata, io, metadata_location)

        with Session(self.engine) as session:
            try:
                session.add(
                    IcebergTables(
                        catalog_name=self.name,
                        table_namespace=namespace,
                        table_name=table_name,
                        metadata_location=metadata_location,
                        previous_metadata_location=None,
                    )
                )
                session.commit()
            except IntegrityError as e:
                raise TableAlreadyExistsError(f"Table {namespace}.{table_name} already exists") from e

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
            NoSuchNamespaceError: If namespace does not exist
        """
        namespace_tuple = Catalog.namespace_from(identifier)
        namespace = Catalog.namespace_to_string(namespace_tuple)
        table_name = Catalog.table_name_from(identifier)
        if not self._namespace_exists(namespace):
            raise NoSuchNamespaceError(f"Namespace does not exist: {namespace}")

        with Session(self.engine) as session:
            try:
                session.add(
                    IcebergTables(
                        catalog_name=self.name,
                        table_namespace=namespace,
                        table_name=table_name,
                        metadata_location=metadata_location,
                        previous_metadata_location=None,
                    )
                )
                session.commit()
            except IntegrityError as e:
                raise TableAlreadyExistsError(f"Table {namespace}.{table_name} already exists") from e

        return self.load_table(identifier=identifier)

    def load_table(self, identifier: Union[str, Identifier]) -> Table:
        """Load the table's metadata and return the table instance.

        You can also use this method to check for table existence using 'try catalog.table() except NoSuchTableError'.
        Note: This method doesn't scan data stored in the table.

        Args:
            identifier (str | Identifier): Table identifier.

        Returns:
            Table: the table instance with its metadata.

        Raises:
            NoSuchTableError: If a table with the name does not exist.
        """
        namespace_tuple = Catalog.namespace_from(identifier)
        namespace = Catalog.namespace_to_string(namespace_tuple)
        table_name = Catalog.table_name_from(identifier)
        with Session(self.engine) as session:
            stmt = select(IcebergTables).where(
                IcebergTables.catalog_name == self.name,
                IcebergTables.table_namespace == namespace,
                IcebergTables.table_name == table_name,
            )
            result = session.scalar(stmt)
        if result:
            return self._convert_orm_to_iceberg(result)
        raise NoSuchTableError(f"Table does not exist: {namespace}.{table_name}")

    def drop_table(self, identifier: Union[str, Identifier]) -> None:
        """Drop a table.

        Args:
            identifier (str | Identifier): Table identifier.

        Raises:
            NoSuchTableError: If a table with the name does not exist.
        """
        namespace_tuple = Catalog.namespace_from(identifier)
        namespace = Catalog.namespace_to_string(namespace_tuple)
        table_name = Catalog.table_name_from(identifier)
        with Session(self.engine) as session:
            if self.engine.dialect.supports_sane_rowcount:
                res = session.execute(
                    delete(IcebergTables).where(
                        IcebergTables.catalog_name == self.name,
                        IcebergTables.table_namespace == namespace,
                        IcebergTables.table_name == table_name,
                    )
                )
                if res.rowcount < 1:
                    raise NoSuchTableError(f"Table does not exist: {namespace}.{table_name}")
            else:
                try:
                    tbl = (
                        session.query(IcebergTables)
                        .with_for_update(of=IcebergTables)
                        .filter(
                            IcebergTables.catalog_name == self.name,
                            IcebergTables.table_namespace == namespace,
                            IcebergTables.table_name == table_name,
                        )
                        .one()
                    )
                    session.delete(tbl)
                except NoResultFound as e:
                    raise NoSuchTableError(f"Table does not exist: {namespace}.{table_name}") from e
            session.commit()

    def rename_table(self, from_identifier: Union[str, Identifier], to_identifier: Union[str, Identifier]) -> Table:
        """Rename a fully classified table name.

        Args:
            from_identifier (str | Identifier): Existing table identifier.
            to_identifier (str | Identifier): New table identifier.

        Returns:
            Table: the updated table instance with its metadata.

        Raises:
            NoSuchTableError: If a table with the name does not exist.
            TableAlreadyExistsError: If a table with the new name already exist.
            NoSuchNamespaceError: If the target namespace does not exist.
        """
        from_namespace_tuple = Catalog.namespace_from(from_identifier)
        from_namespace = Catalog.namespace_to_string(from_namespace_tuple)
        from_table_name = Catalog.table_name_from(from_identifier)
        to_namespace_tuple = Catalog.namespace_from(to_identifier)
        to_namespace = Catalog.namespace_to_string(to_namespace_tuple)
        to_table_name = Catalog.table_name_from(to_identifier)
        if not self._namespace_exists(to_namespace):
            raise NoSuchNamespaceError(f"Namespace does not exist: {to_namespace}")
        with Session(self.engine) as session:
            try:
                if self.engine.dialect.supports_sane_rowcount:
                    stmt = (
                        update(IcebergTables)
                        .where(
                            IcebergTables.catalog_name == self.name,
                            IcebergTables.table_namespace == from_namespace,
                            IcebergTables.table_name == from_table_name,
                        )
                        .values(table_namespace=to_namespace, table_name=to_table_name)
                    )
                    result = session.execute(stmt)
                    if result.rowcount < 1:
                        raise NoSuchTableError(f"Table does not exist: {from_table_name}")
                else:
                    try:
                        tbl = (
                            session.query(IcebergTables)
                            .with_for_update(of=IcebergTables)
                            .filter(
                                IcebergTables.catalog_name == self.name,
                                IcebergTables.table_namespace == from_namespace,
                                IcebergTables.table_name == from_table_name,
                            )
                            .one()
                        )
                        tbl.table_namespace = to_namespace
                        tbl.table_name = to_table_name
                    except NoResultFound as e:
                        raise NoSuchTableError(f"Table does not exist: {from_table_name}") from e
                session.commit()
            except IntegrityError as e:
                raise TableAlreadyExistsError(f"Table {to_namespace}.{to_table_name} already exists") from e
        return self.load_table(to_identifier)

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
        table_identifier = table.name()
        namespace_tuple = Catalog.namespace_from(table_identifier)
        namespace = Catalog.namespace_to_string(namespace_tuple)
        table_name = Catalog.table_name_from(table_identifier)

        current_table: Optional[Table]
        try:
            current_table = self.load_table(table_identifier)
        except NoSuchTableError:
            current_table = None

        updated_staged_table = self._update_and_stage_table(current_table, table.name(), requirements, updates)
        if current_table and updated_staged_table.metadata == current_table.metadata:
            # no changes, do nothing
            return CommitTableResponse(metadata=current_table.metadata, metadata_location=current_table.metadata_location)
        self._write_metadata(
            metadata=updated_staged_table.metadata,
            io=updated_staged_table.io,
            metadata_path=updated_staged_table.metadata_location,
        )

        with Session(self.engine) as session:
            if current_table:
                # table exists, update it
                if self.engine.dialect.supports_sane_rowcount:
                    stmt = (
                        update(IcebergTables)
                        .where(
                            IcebergTables.catalog_name == self.name,
                            IcebergTables.table_namespace == namespace,
                            IcebergTables.table_name == table_name,
                            IcebergTables.metadata_location == current_table.metadata_location,
                        )
                        .values(
                            metadata_location=updated_staged_table.metadata_location,
                            previous_metadata_location=current_table.metadata_location,
                        )
                    )
                    result = session.execute(stmt)
                    if result.rowcount < 1:
                        raise CommitFailedException(f"Table has been updated by another process: {namespace}.{table_name}")
                else:
                    try:
                        tbl = (
                            session.query(IcebergTables)
                            .with_for_update(of=IcebergTables)
                            .filter(
                                IcebergTables.catalog_name == self.name,
                                IcebergTables.table_namespace == namespace,
                                IcebergTables.table_name == table_name,
                                IcebergTables.metadata_location == current_table.metadata_location,
                            )
                            .one()
                        )
                        tbl.metadata_location = updated_staged_table.metadata_location
                        tbl.previous_metadata_location = current_table.metadata_location
                    except NoResultFound as e:
                        raise CommitFailedException(f"Table has been updated by another process: {namespace}.{table_name}") from e
                session.commit()
            else:
                # table does not exist, create it
                try:
                    session.add(
                        IcebergTables(
                            catalog_name=self.name,
                            table_namespace=namespace,
                            table_name=table_name,
                            metadata_location=updated_staged_table.metadata_location,
                            previous_metadata_location=None,
                        )
                    )
                    session.commit()
                except IntegrityError as e:
                    raise TableAlreadyExistsError(f"Table {namespace}.{table_name} already exists") from e

        return CommitTableResponse(
            metadata=updated_staged_table.metadata, metadata_location=updated_staged_table.metadata_location
        )

    def _namespace_exists(self, identifier: Union[str, Identifier]) -> bool:
        namespace_tuple = Catalog.identifier_to_tuple(identifier)
        namespace = Catalog.namespace_to_string(namespace_tuple, NoSuchNamespaceError)
        namespace_starts_with = namespace.replace("!", "!!").replace("_", "!_").replace("%", "!%") + ".%"

        with Session(self.engine) as session:
            stmt = (
                select(IcebergTables)
                .where(
                    IcebergTables.catalog_name == self.name,
                    (IcebergTables.table_namespace == namespace)
                    | (IcebergTables.table_namespace.like(namespace_starts_with, escape="!")),
                )
                .limit(1)
            )
            result = session.execute(stmt).all()
            if result:
                return True
            stmt = (
                select(IcebergNamespaceProperties)
                .where(
                    IcebergNamespaceProperties.catalog_name == self.name,
                    (IcebergNamespaceProperties.namespace == namespace)
                    | (IcebergNamespaceProperties.namespace.like(namespace_starts_with, escape="!")),
                )
                .limit(1)
            )
            result = session.execute(stmt).all()
            if result:
                return True
        return False

    def create_namespace(self, namespace: Union[str, Identifier], properties: Properties = EMPTY_DICT) -> None:
        """Create a namespace in the catalog.

        Args:
            namespace (str | Identifier): Namespace identifier.
            properties (Properties): A string dictionary of properties for the given namespace.

        Raises:
            NamespaceAlreadyExistsError: If a namespace with the given name already exists.
        """
        if self._namespace_exists(namespace):
            raise NamespaceAlreadyExistsError(f"Namespace {namespace} already exists")

        if not properties:
            properties = IcebergNamespaceProperties.NAMESPACE_MINIMAL_PROPERTIES
        create_properties = properties if properties else IcebergNamespaceProperties.NAMESPACE_MINIMAL_PROPERTIES
        with Session(self.engine) as session:
            for key, value in create_properties.items():
                session.add(
                    IcebergNamespaceProperties(
                        catalog_name=self.name,
                        namespace=Catalog.namespace_to_string(namespace, NoSuchNamespaceError),
                        property_key=key,
                        property_value=value,
                    )
                )
            session.commit()

    def drop_namespace(self, namespace: Union[str, Identifier]) -> None:
        """Drop a namespace.

        Args:
            namespace (str | Identifier): Namespace identifier.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
            NamespaceNotEmptyError: If the namespace is not empty.
        """
        if not self._namespace_exists(namespace):
            raise NoSuchNamespaceError(f"Namespace does not exist: {namespace}")

        namespace_str = Catalog.namespace_to_string(namespace)
        if tables := self.list_tables(namespace):
            raise NamespaceNotEmptyError(f"Namespace {namespace_str} is not empty. {len(tables)} tables exist.")

        with Session(self.engine) as session:
            session.execute(
                delete(IcebergNamespaceProperties).where(
                    IcebergNamespaceProperties.catalog_name == self.name,
                    IcebergNamespaceProperties.namespace == namespace_str,
                )
            )
            session.commit()

    def list_tables(self, namespace: Union[str, Identifier]) -> List[Identifier]:
        """List tables under the given namespace in the catalog.

        Args:
            namespace (str | Identifier): Namespace identifier to search.

        Returns:
            List[Identifier]: list of table identifiers.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
        """
        if namespace and not self._namespace_exists(namespace):
            raise NoSuchNamespaceError(f"Namespace does not exist: {namespace}")

        namespace = Catalog.namespace_to_string(namespace)
        stmt = select(IcebergTables).where(IcebergTables.catalog_name == self.name, IcebergTables.table_namespace == namespace)
        with Session(self.engine) as session:
            result = session.scalars(stmt)
            return [(Catalog.identifier_to_tuple(table.table_namespace) + (table.table_name,)) for table in result]

    def list_namespaces(self, namespace: Union[str, Identifier] = ()) -> List[Identifier]:
        """List namespaces from the given namespace. If not given, list top-level namespaces from the catalog.

        Args:
            namespace (str | Identifier): Namespace identifier to search.

        Returns:
            List[Identifier]: a List of namespace identifiers.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
        """
        if namespace and not self._namespace_exists(namespace):
            raise NoSuchNamespaceError(f"Namespace does not exist: {namespace}")

        table_stmt = select(IcebergTables.table_namespace).where(IcebergTables.catalog_name == self.name)
        namespace_stmt = select(IcebergNamespaceProperties.namespace).where(IcebergNamespaceProperties.catalog_name == self.name)
        if namespace:
            namespace_like = Catalog.namespace_to_string(namespace, NoSuchNamespaceError) + "%"
            table_stmt = table_stmt.where(IcebergTables.table_namespace.like(namespace_like))
            namespace_stmt = namespace_stmt.where(IcebergNamespaceProperties.namespace.like(namespace_like))
        stmt = union(
            table_stmt,
            namespace_stmt,
        )
        with Session(self.engine) as session:
            namespace_tuple = Catalog.identifier_to_tuple(namespace)
            sub_namespaces_level_length = len(namespace_tuple) + 1

            namespaces = list(
                {  # only get distinct namespaces
                    ns[:sub_namespaces_level_length]  # truncate to the required level
                    for ns in {Catalog.identifier_to_tuple(ns) for ns in session.execute(stmt).scalars()}
                    if len(ns) >= sub_namespaces_level_length  # only get sub namespaces/children
                    and ns[: sub_namespaces_level_length - 1] == namespace_tuple
                    # exclude fuzzy matches when `namespace` contains `%` or `_`
                }
            )

            return namespaces

    def load_namespace_properties(self, namespace: Union[str, Identifier]) -> Properties:
        """Get properties for a namespace.

        Args:
            namespace (str | Identifier): Namespace identifier.

        Returns:
            Properties: Properties for the given namespace.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
        """
        namespace_str = Catalog.namespace_to_string(namespace)
        if not self._namespace_exists(namespace):
            raise NoSuchNamespaceError(f"Namespace {namespace_str} does not exists")

        stmt = select(IcebergNamespaceProperties).where(
            IcebergNamespaceProperties.catalog_name == self.name, IcebergNamespaceProperties.namespace == namespace_str
        )
        with Session(self.engine) as session:
            result = session.scalars(stmt)
            return {props.property_key: props.property_value for props in result}

    def update_namespace_properties(
        self, namespace: Union[str, Identifier], removals: Optional[Set[str]] = None, updates: Properties = EMPTY_DICT
    ) -> PropertiesUpdateSummary:
        """Remove provided property keys and update properties for a namespace.

        Args:
            namespace (str | Identifier): Namespace identifier.
            removals (Set[str]): Set of property keys that need to be removed. Optional Argument.
            updates (Properties): Properties to be updated for the given namespace.

        Raises:
            NoSuchNamespaceError: If a namespace with the given name does not exist.
            ValueError: If removals and updates have overlapping keys.
        """
        namespace_str = Catalog.namespace_to_string(namespace)
        if not self._namespace_exists(namespace):
            raise NoSuchNamespaceError(f"Namespace {namespace_str} does not exists")

        current_properties = self.load_namespace_properties(namespace=namespace)
        properties_update_summary = self._get_updated_props_and_update_summary(
            current_properties=current_properties, removals=removals, updates=updates
        )[0]

        with Session(self.engine) as session:
            if removals:
                delete_stmt = delete(IcebergNamespaceProperties).where(
                    IcebergNamespaceProperties.catalog_name == self.name,
                    IcebergNamespaceProperties.namespace == namespace_str,
                    IcebergNamespaceProperties.property_key.in_(removals),
                )
                session.execute(delete_stmt)

            if updates:
                # SQLAlchemy does not (yet) support engine agnostic UPSERT
                # https://docs.sqlalchemy.org/en/20/orm/queryguide/dml.html#orm-upsert-statements
                # This is not a problem since it runs in a single transaction
                delete_stmt = delete(IcebergNamespaceProperties).where(
                    IcebergNamespaceProperties.catalog_name == self.name,
                    IcebergNamespaceProperties.namespace == namespace_str,
                    IcebergNamespaceProperties.property_key.in_(set(updates.keys())),
                )
                session.execute(delete_stmt)
                insert_stmt_values = [
                    {
                        IcebergNamespaceProperties.catalog_name: self.name,
                        IcebergNamespaceProperties.namespace: namespace_str,
                        IcebergNamespaceProperties.property_key: property_key,
                        IcebergNamespaceProperties.property_value: property_value,
                    }
                    for property_key, property_value in updates.items()
                ]
                insert_stmt = insert(IcebergNamespaceProperties).values(insert_stmt_values)
                session.execute(insert_stmt)
            session.commit()
        return properties_update_summary

    def list_views(self, namespace: Union[str, Identifier]) -> List[Identifier]:
        raise NotImplementedError

    def view_exists(self, identifier: Union[str, Identifier]) -> bool:
        raise NotImplementedError

    def drop_view(self, identifier: Union[str, Identifier]) -> None:
        raise NotImplementedError

    def close(self) -> None:
        """Close the catalog and release database connections.

        This method closes the SQLAlchemy engine and disposes of all connection pools.
        This ensures that any cached connections are properly closed, which is especially
        important for blobfuse scenarios where file handles need to be closed for
        data to be flushed to persistent storage.
        """
        if hasattr(self, "engine"):
            self.engine.dispose()
