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
from enum import Enum
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

from pydantic import Field, field_validator
from requests import HTTPError, Session
from tenacity import RetryCallState, retry, retry_if_exception_type, stop_after_attempt

from pyiceberg import __version__
from pyiceberg.catalog import (
    BOTOCORE_SESSION,
    TOKEN,
    URI,
    WAREHOUSE_LOCATION,
    Catalog,
    PropertiesUpdateSummary,
)
from pyiceberg.catalog.rest.auth import AuthManager, AuthManagerAdapter, AuthManagerFactory, LegacyOAuth2AuthManager
from pyiceberg.catalog.rest.response import _handle_non_200_response
from pyiceberg.exceptions import (
    AuthorizationExpiredError,
    CommitFailedException,
    CommitStateUnknownException,
    NamespaceAlreadyExistsError,
    NamespaceNotEmptyError,
    NoSuchIdentifierError,
    NoSuchNamespaceError,
    NoSuchTableError,
    NoSuchViewError,
    TableAlreadyExistsError,
    UnauthorizedError,
)
from pyiceberg.io import AWS_ACCESS_KEY_ID, AWS_REGION, AWS_SECRET_ACCESS_KEY, AWS_SESSION_TOKEN
from pyiceberg.partitioning import UNPARTITIONED_PARTITION_SPEC, PartitionSpec, assign_fresh_partition_spec_ids
from pyiceberg.schema import Schema, assign_fresh_schema_ids
from pyiceberg.table import (
    CommitTableRequest,
    CommitTableResponse,
    CreateTableTransaction,
    StagedTable,
    Table,
    TableIdentifier,
    TableProperties,
)
from pyiceberg.table.metadata import TableMetadata
from pyiceberg.table.sorting import UNSORTED_SORT_ORDER, SortOrder, assign_fresh_sort_order_ids
from pyiceberg.table.update import (
    TableRequirement,
    TableUpdate,
)
from pyiceberg.typedef import EMPTY_DICT, UTF8, IcebergBaseModel, Identifier, Properties
from pyiceberg.types import transform_dict_value_to_str
from pyiceberg.utils.deprecated import deprecation_message
from pyiceberg.utils.properties import get_first_property_value, get_header_properties, property_as_bool

if TYPE_CHECKING:
    import pyarrow as pa


class Endpoints:
    get_config: str = "config"
    list_namespaces: str = "namespaces"
    create_namespace: str = "namespaces"
    load_namespace_metadata: str = "namespaces/{namespace}"
    drop_namespace: str = "namespaces/{namespace}"
    update_namespace_properties: str = "namespaces/{namespace}/properties"
    namespace_exists: str = "namespaces/{namespace}"
    list_tables: str = "namespaces/{namespace}/tables"
    create_table: str = "namespaces/{namespace}/tables"
    register_table = "namespaces/{namespace}/register"
    load_table: str = "namespaces/{namespace}/tables/{table}"
    update_table: str = "namespaces/{namespace}/tables/{table}"
    drop_table: str = "namespaces/{namespace}/tables/{table}"
    table_exists: str = "namespaces/{namespace}/tables/{table}"
    get_token: str = "oauth/tokens"
    rename_table: str = "tables/rename"
    list_views: str = "namespaces/{namespace}/views"
    drop_view: str = "namespaces/{namespace}/views/{view}"
    view_exists: str = "namespaces/{namespace}/views/{view}"


class IdentifierKind(Enum):
    TABLE = "table"
    VIEW = "view"


ACCESS_DELEGATION_DEFAULT = "vended-credentials"
AUTHORIZATION_HEADER = "Authorization"
BEARER_PREFIX = "Bearer"
CATALOG_SCOPE = "catalog"
CLIENT_ID = "client_id"
PREFIX = "prefix"
CLIENT_SECRET = "client_secret"
CLIENT_CREDENTIALS = "client_credentials"
CREDENTIAL = "credential"
GRANT_TYPE = "grant_type"
SCOPE = "scope"
AUDIENCE = "audience"
RESOURCE = "resource"
TOKEN_EXCHANGE = "urn:ietf:params:oauth:grant-type:token-exchange"
SEMICOLON = ":"
KEY = "key"
CERT = "cert"
CLIENT = "client"
CA_BUNDLE = "cabundle"
SSL = "ssl"
SIGV4 = "rest.sigv4-enabled"
SIGV4_REGION = "rest.signing-region"
SIGV4_SERVICE = "rest.signing-name"
OAUTH2_SERVER_URI = "oauth2-server-uri"
SNAPSHOT_LOADING_MODE = "snapshot-loading-mode"
AUTH = "auth"
CUSTOM = "custom"

NAMESPACE_SEPARATOR = b"\x1f".decode(UTF8)


def _retry_hook(retry_state: RetryCallState) -> None:
    rest_catalog: RestCatalog = retry_state.args[0]
    rest_catalog._refresh_token()  # pylint: disable=protected-access


_RETRY_ARGS = {
    "retry": retry_if_exception_type((AuthorizationExpiredError, UnauthorizedError)),
    "stop": stop_after_attempt(2),
    "before_sleep": _retry_hook,
    "reraise": True,
}


class TableResponse(IcebergBaseModel):
    metadata_location: Optional[str] = Field(alias="metadata-location", default=None)
    metadata: TableMetadata
    config: Properties = Field(default_factory=dict)


class CreateTableRequest(IcebergBaseModel):
    name: str = Field()
    location: Optional[str] = Field()
    table_schema: Schema = Field(alias="schema")
    partition_spec: Optional[PartitionSpec] = Field(alias="partition-spec")
    write_order: Optional[SortOrder] = Field(alias="write-order")
    stage_create: bool = Field(alias="stage-create", default=False)
    properties: Dict[str, str] = Field(default_factory=dict)

    # validators
    @field_validator("properties", mode="before")
    def transform_properties_dict_value_to_str(cls, properties: Properties) -> Dict[str, str]:
        return transform_dict_value_to_str(properties)


class RegisterTableRequest(IcebergBaseModel):
    name: str
    metadata_location: str = Field(..., alias="metadata-location")


class ConfigResponse(IcebergBaseModel):
    defaults: Optional[Properties] = Field(default_factory=dict)
    overrides: Optional[Properties] = Field(default_factory=dict)


class ListNamespaceResponse(IcebergBaseModel):
    namespaces: List[Identifier] = Field()


class NamespaceResponse(IcebergBaseModel):
    namespace: Identifier = Field()
    properties: Properties = Field()


class UpdateNamespacePropertiesResponse(IcebergBaseModel):
    removed: List[str] = Field()
    updated: List[str] = Field()
    missing: List[str] = Field()


class ListTableResponseEntry(IcebergBaseModel):
    name: str = Field()
    namespace: Identifier = Field()


class ListViewResponseEntry(IcebergBaseModel):
    name: str = Field()
    namespace: Identifier = Field()


class ListTablesResponse(IcebergBaseModel):
    identifiers: List[ListTableResponseEntry] = Field()


class ListViewsResponse(IcebergBaseModel):
    identifiers: List[ListViewResponseEntry] = Field()


class RestCatalog(Catalog):
    uri: str
    _session: Session

    def __init__(self, name: str, **properties: str):
        """Rest Catalog.

        You either need to provide a client_id and client_secret, or an already valid token.

        Args:
            name: Name to identify the catalog.
            properties: Properties that are passed along to the configuration.
        """
        super().__init__(name, **properties)
        self.uri = properties[URI]
        self._fetch_config()
        self._session = self._create_session()

    def _create_session(self) -> Session:
        """Create a request session with provided catalog configuration."""
        session = Session()

        # Sets the client side and server side SSL cert verification, if provided as properties.
        if ssl_config := self.properties.get(SSL):
            if ssl_ca_bundle := ssl_config.get(CA_BUNDLE):
                session.verify = ssl_ca_bundle
            if ssl_client := ssl_config.get(CLIENT):
                if all(k in ssl_client for k in (CERT, KEY)):
                    session.cert = (ssl_client[CERT], ssl_client[KEY])
                elif ssl_client_cert := ssl_client.get(CERT):
                    session.cert = ssl_client_cert

        if auth_config := self.properties.get(AUTH):
            auth_type = auth_config.get("type")
            if auth_type is None:
                raise ValueError("auth.type must be defined")
            auth_type_config = auth_config.get(auth_type, {})
            auth_impl = auth_config.get("impl")

            if auth_type == CUSTOM and not auth_impl:
                raise ValueError("auth.impl must be specified when using custom auth.type")

            if auth_type != CUSTOM and auth_impl:
                raise ValueError("auth.impl can only be specified when using custom auth.type")

            session.auth = AuthManagerAdapter(AuthManagerFactory.create(auth_impl or auth_type, auth_type_config))
        else:
            session.auth = AuthManagerAdapter(self._create_legacy_oauth2_auth_manager(session))

        # Set HTTP headers
        self._config_headers(session)

        # Configure SigV4 Request Signing
        if property_as_bool(self.properties, SIGV4, False):
            self._init_sigv4(session)

        return session

    def _create_legacy_oauth2_auth_manager(self, session: Session) -> AuthManager:
        """Create the LegacyOAuth2AuthManager by fetching required properties.

        This will be removed in PyIceberg 1.0
        """
        client_credentials = self.properties.get(CREDENTIAL)
        # We want to call `self.auth_url` only when we are using CREDENTIAL
        # with the legacy OAUTH2 flow as it will raise a DeprecationWarning
        auth_url = self.auth_url if client_credentials is not None else None

        auth_config = {
            "session": session,
            "auth_url": auth_url,
            "credential": client_credentials,
            "initial_token": self.properties.get(TOKEN),
            "optional_oauth_params": self._extract_optional_oauth_params(),
        }

        return AuthManagerFactory.create("legacyoauth2", auth_config)

    def _check_valid_namespace_identifier(self, identifier: Union[str, Identifier]) -> Identifier:
        """Check if the identifier has at least one element."""
        identifier_tuple = Catalog.identifier_to_tuple(identifier)
        if len(identifier_tuple) < 1:
            raise NoSuchNamespaceError(f"Empty namespace identifier: {identifier}")
        return identifier_tuple

    def url(self, endpoint: str, prefixed: bool = True, **kwargs: Any) -> str:
        """Construct the endpoint.

        Args:
            endpoint: Resource identifier that points to the REST catalog.
            prefixed: If the prefix return by the config needs to be appended.

        Returns:
            The base url of the rest catalog.
        """
        url = self.uri
        url = url + "v1/" if url.endswith("/") else url + "/v1/"

        if prefixed:
            url += self.properties.get(PREFIX, "")
            url = url if url.endswith("/") else url + "/"

        return url + endpoint.format(**kwargs)

    @property
    def auth_url(self) -> str:
        self._warn_oauth_tokens_deprecation()

        if url := self.properties.get(OAUTH2_SERVER_URI):
            return url
        else:
            return self.url(Endpoints.get_token, prefixed=False)

    def _warn_oauth_tokens_deprecation(self) -> None:
        has_oauth_server_uri = OAUTH2_SERVER_URI in self.properties
        has_credential = CREDENTIAL in self.properties
        has_init_token = TOKEN in self.properties
        has_sigv4_enabled = property_as_bool(self.properties, SIGV4, False)

        if not has_oauth_server_uri and (has_init_token or has_credential) and not has_sigv4_enabled:
            deprecation_message(
                deprecated_in="0.8.0",
                removed_in="1.0.0",
                help_message="Iceberg REST client is missing the OAuth2 server URI "
                f"configuration and defaults to {self.uri}{Endpoints.get_token}. "
                "This automatic fallback will be removed in a future Iceberg release."
                f"It is recommended to configure the OAuth2 endpoint using the '{OAUTH2_SERVER_URI}'"
                "property to be prepared. This warning will disappear if the OAuth2"
                "endpoint is explicitly configured. See https://github.com/apache/iceberg/issues/10537",
            )

    def _extract_optional_oauth_params(self) -> Dict[str, str]:
        optional_oauth_param = {SCOPE: self.properties.get(SCOPE) or CATALOG_SCOPE}
        set_of_optional_params = {AUDIENCE, RESOURCE}
        for param in set_of_optional_params:
            if param_value := self.properties.get(param):
                optional_oauth_param[param] = param_value

        return optional_oauth_param

    def _fetch_config(self) -> None:
        params = {}
        if warehouse_location := self.properties.get(WAREHOUSE_LOCATION):
            params[WAREHOUSE_LOCATION] = warehouse_location

        with self._create_session() as session:
            response = session.get(self.url(Endpoints.get_config, prefixed=False), params=params)
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {})
        config_response = ConfigResponse.model_validate_json(response.text)

        config = config_response.defaults
        config.update(self.properties)
        config.update(config_response.overrides)
        self.properties = config

        # Update URI based on overrides
        self.uri = config[URI]

    def _identifier_to_validated_tuple(self, identifier: Union[str, Identifier]) -> Identifier:
        identifier_tuple = self.identifier_to_tuple(identifier)
        if len(identifier_tuple) <= 1:
            raise NoSuchIdentifierError(f"Missing namespace or invalid identifier: {'.'.join(identifier_tuple)}")
        return identifier_tuple

    def _split_identifier_for_path(
        self, identifier: Union[str, Identifier, TableIdentifier], kind: IdentifierKind = IdentifierKind.TABLE
    ) -> Properties:
        if isinstance(identifier, TableIdentifier):
            return {"namespace": NAMESPACE_SEPARATOR.join(identifier.namespace.root), kind.value: identifier.name}
        identifier_tuple = self._identifier_to_validated_tuple(identifier)

        return {"namespace": NAMESPACE_SEPARATOR.join(identifier_tuple[:-1]), kind.value: identifier_tuple[-1]}

    def _split_identifier_for_json(self, identifier: Union[str, Identifier]) -> Dict[str, Union[Identifier, str]]:
        identifier_tuple = self._identifier_to_validated_tuple(identifier)
        return {"namespace": identifier_tuple[:-1], "name": identifier_tuple[-1]}

    def _init_sigv4(self, session: Session) -> None:
        from urllib import parse

        import boto3
        from botocore.auth import SigV4Auth
        from botocore.awsrequest import AWSRequest
        from requests import PreparedRequest
        from requests.adapters import HTTPAdapter

        class SigV4Adapter(HTTPAdapter):
            def __init__(self, **properties: str):
                super().__init__()
                self._properties = properties
                self._boto_session = boto3.Session(
                    region_name=get_first_property_value(self._properties, AWS_REGION),
                    botocore_session=self._properties.get(BOTOCORE_SESSION),
                    aws_access_key_id=get_first_property_value(self._properties, AWS_ACCESS_KEY_ID),
                    aws_secret_access_key=get_first_property_value(self._properties, AWS_SECRET_ACCESS_KEY),
                    aws_session_token=get_first_property_value(self._properties, AWS_SESSION_TOKEN),
                )

            def add_headers(self, request: PreparedRequest, **kwargs: Any) -> None:  # pylint: disable=W0613
                credentials = self._boto_session.get_credentials().get_frozen_credentials()
                region = self._properties.get(SIGV4_REGION, self._boto_session.region_name)
                service = self._properties.get(SIGV4_SERVICE, "execute-api")

                url = str(request.url).split("?")[0]
                query = str(parse.urlsplit(request.url).query)
                params = dict(parse.parse_qsl(query))

                # remove the connection header as it will be updated after signing
                del request.headers["connection"]

                aws_request = AWSRequest(
                    method=request.method, url=url, params=params, data=request.body, headers=dict(request.headers)
                )

                SigV4Auth(credentials, service, region).add_auth(aws_request)
                original_header = request.headers
                signed_headers = aws_request.headers
                relocated_headers = {}

                # relocate headers if there is a conflict with signed headers
                for header, value in original_header.items():
                    if header in signed_headers and signed_headers[header] != value:
                        relocated_headers[f"Original-{header}"] = value

                request.headers.update(relocated_headers)
                request.headers.update(signed_headers)

        session.mount(self.uri, SigV4Adapter(**self.properties))

    def _response_to_table(self, identifier_tuple: Tuple[str, ...], table_response: TableResponse) -> Table:
        return Table(
            identifier=identifier_tuple,
            metadata_location=table_response.metadata_location,  # type: ignore
            metadata=table_response.metadata,
            io=self._load_file_io(
                {**table_response.metadata.properties, **table_response.config}, table_response.metadata_location
            ),
            catalog=self,
            config=table_response.config,
        )

    def _response_to_staged_table(self, identifier_tuple: Tuple[str, ...], table_response: TableResponse) -> StagedTable:
        return StagedTable(
            identifier=identifier_tuple,
            metadata_location=table_response.metadata_location,  # type: ignore
            metadata=table_response.metadata,
            io=self._load_file_io(
                {**table_response.metadata.properties, **table_response.config}, table_response.metadata_location
            ),
            catalog=self,
        )

    def _refresh_token(self) -> None:
        # Reactive token refresh is atypical - we should proactively refresh tokens in a separate thread
        # instead of retrying on Auth Exceptions. Keeping refresh behavior for the LegacyOAuth2AuthManager
        # for backward compatibility
        auth_manager = self._session.auth.auth_manager  # type: ignore[union-attr]
        if isinstance(auth_manager, LegacyOAuth2AuthManager):
            auth_manager._refresh_token()

    def _config_headers(self, session: Session) -> None:
        header_properties = get_header_properties(self.properties)
        session.headers.update(header_properties)
        session.headers["Content-type"] = "application/json"
        session.headers["User-Agent"] = f"PyIceberg/{__version__}"
        session.headers.setdefault("X-Iceberg-Access-Delegation", ACCESS_DELEGATION_DEFAULT)

    def _create_table(
        self,
        identifier: Union[str, Identifier],
        schema: Union[Schema, "pa.Schema"],
        location: Optional[str] = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
        stage_create: bool = False,
    ) -> TableResponse:
        iceberg_schema = self._convert_schema_if_needed(
            schema,
            int(properties.get(TableProperties.FORMAT_VERSION, TableProperties.DEFAULT_FORMAT_VERSION)),  # type: ignore
        )
        fresh_schema = assign_fresh_schema_ids(iceberg_schema)
        fresh_partition_spec = assign_fresh_partition_spec_ids(partition_spec, iceberg_schema, fresh_schema)
        fresh_sort_order = assign_fresh_sort_order_ids(sort_order, iceberg_schema, fresh_schema)

        namespace_and_table = self._split_identifier_for_path(identifier)
        if location:
            location = location.rstrip("/")
        request = CreateTableRequest(
            name=namespace_and_table["table"],
            location=location,
            table_schema=fresh_schema,
            partition_spec=fresh_partition_spec,
            write_order=fresh_sort_order,
            stage_create=stage_create,
            properties=properties,
        )
        serialized_json = request.model_dump_json().encode(UTF8)
        response = self._session.post(
            self.url(Endpoints.create_table, namespace=namespace_and_table["namespace"]),
            data=serialized_json,
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {409: TableAlreadyExistsError, 404: NoSuchNamespaceError})
        return TableResponse.model_validate_json(response.text)

    @retry(**_RETRY_ARGS)
    def create_table(
        self,
        identifier: Union[str, Identifier],
        schema: Union[Schema, "pa.Schema"],
        location: Optional[str] = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
    ) -> Table:
        table_response = self._create_table(
            identifier=identifier,
            schema=schema,
            location=location,
            partition_spec=partition_spec,
            sort_order=sort_order,
            properties=properties,
            stage_create=False,
        )
        return self._response_to_table(self.identifier_to_tuple(identifier), table_response)

    @retry(**_RETRY_ARGS)
    def create_table_transaction(
        self,
        identifier: Union[str, Identifier],
        schema: Union[Schema, "pa.Schema"],
        location: Optional[str] = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
    ) -> CreateTableTransaction:
        table_response = self._create_table(
            identifier=identifier,
            schema=schema,
            location=location,
            partition_spec=partition_spec,
            sort_order=sort_order,
            properties=properties,
            stage_create=True,
        )
        staged_table = self._response_to_staged_table(self.identifier_to_tuple(identifier), table_response)
        return CreateTableTransaction(staged_table)

    @retry(**_RETRY_ARGS)
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
        namespace_and_table = self._split_identifier_for_path(identifier)
        request = RegisterTableRequest(
            name=namespace_and_table["table"],
            metadata_location=metadata_location,
        )
        serialized_json = request.model_dump_json().encode(UTF8)
        response = self._session.post(
            self.url(Endpoints.register_table, namespace=namespace_and_table["namespace"]),
            data=serialized_json,
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {409: TableAlreadyExistsError})

        table_response = TableResponse.model_validate_json(response.text)
        return self._response_to_table(self.identifier_to_tuple(identifier), table_response)

    @retry(**_RETRY_ARGS)
    def list_tables(self, namespace: Union[str, Identifier]) -> List[Identifier]:
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace_concat = NAMESPACE_SEPARATOR.join(namespace_tuple)
        response = self._session.get(self.url(Endpoints.list_tables, namespace=namespace_concat))
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError})
        return [(*table.namespace, table.name) for table in ListTablesResponse.model_validate_json(response.text).identifiers]

    @retry(**_RETRY_ARGS)
    def load_table(self, identifier: Union[str, Identifier]) -> Table:
        params = {}
        if mode := self.properties.get(SNAPSHOT_LOADING_MODE):
            if mode in {"all", "refs"}:
                params["snapshots"] = mode
            else:
                raise ValueError("Invalid snapshot-loading-mode: {}")

        response = self._session.get(
            self.url(Endpoints.load_table, prefixed=True, **self._split_identifier_for_path(identifier)), params=params
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchTableError})

        table_response = TableResponse.model_validate_json(response.text)
        return self._response_to_table(self.identifier_to_tuple(identifier), table_response)

    @retry(**_RETRY_ARGS)
    def drop_table(self, identifier: Union[str, Identifier], purge_requested: bool = False) -> None:
        response = self._session.delete(
            self.url(Endpoints.drop_table, prefixed=True, **self._split_identifier_for_path(identifier)),
            params={"purgeRequested": purge_requested},
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchTableError})

    @retry(**_RETRY_ARGS)
    def purge_table(self, identifier: Union[str, Identifier]) -> None:
        self.drop_table(identifier=identifier, purge_requested=True)

    @retry(**_RETRY_ARGS)
    def rename_table(self, from_identifier: Union[str, Identifier], to_identifier: Union[str, Identifier]) -> Table:
        payload = {
            "source": self._split_identifier_for_json(from_identifier),
            "destination": self._split_identifier_for_json(to_identifier),
        }
        response = self._session.post(self.url(Endpoints.rename_table), json=payload)
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchTableError, 409: TableAlreadyExistsError})

        return self.load_table(to_identifier)

    def _remove_catalog_name_from_table_request_identifier(self, table_request: CommitTableRequest) -> CommitTableRequest:
        if table_request.identifier.namespace.root[0] == self.name:
            return table_request.model_copy(
                update={
                    "identifier": TableIdentifier(
                        namespace=table_request.identifier.namespace.root[1:], name=table_request.identifier.name
                    )
                }
            )
        return table_request

    @retry(**_RETRY_ARGS)
    def list_views(self, namespace: Union[str, Identifier]) -> List[Identifier]:
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace_concat = NAMESPACE_SEPARATOR.join(namespace_tuple)
        response = self._session.get(self.url(Endpoints.list_views, namespace=namespace_concat))
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError})
        return [(*view.namespace, view.name) for view in ListViewsResponse.model_validate_json(response.text).identifiers]

    @retry(**_RETRY_ARGS)
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
            CommitStateUnknownException: Failed due to an internal exception on the side of the catalog.
        """
        identifier = table.name()
        table_identifier = TableIdentifier(namespace=identifier[:-1], name=identifier[-1])
        table_request = CommitTableRequest(identifier=table_identifier, requirements=requirements, updates=updates)

        headers = self._session.headers
        if table_token := table.config.get(TOKEN):
            headers[AUTHORIZATION_HEADER] = f"{BEARER_PREFIX} {table_token}"

        response = self._session.post(
            self.url(Endpoints.update_table, prefixed=True, **self._split_identifier_for_path(table_request.identifier)),
            data=table_request.model_dump_json().encode(UTF8),
            headers=headers,
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(
                exc,
                {
                    409: CommitFailedException,
                    500: CommitStateUnknownException,
                    502: CommitStateUnknownException,
                    504: CommitStateUnknownException,
                },
            )
        return CommitTableResponse.model_validate_json(response.text)

    @retry(**_RETRY_ARGS)
    def create_namespace(self, namespace: Union[str, Identifier], properties: Properties = EMPTY_DICT) -> None:
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        payload = {"namespace": namespace_tuple, "properties": properties}
        response = self._session.post(self.url(Endpoints.create_namespace), json=payload)
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {409: NamespaceAlreadyExistsError})

    @retry(**_RETRY_ARGS)
    def drop_namespace(self, namespace: Union[str, Identifier]) -> None:
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace = NAMESPACE_SEPARATOR.join(namespace_tuple)
        response = self._session.delete(self.url(Endpoints.drop_namespace, namespace=namespace))
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError, 409: NamespaceNotEmptyError})

    @retry(**_RETRY_ARGS)
    def list_namespaces(self, namespace: Union[str, Identifier] = ()) -> List[Identifier]:
        namespace_tuple = self.identifier_to_tuple(namespace)
        response = self._session.get(
            self.url(
                f"{Endpoints.list_namespaces}?parent={NAMESPACE_SEPARATOR.join(namespace_tuple)}"
                if namespace_tuple
                else Endpoints.list_namespaces
            ),
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError})

        return ListNamespaceResponse.model_validate_json(response.text).namespaces

    @retry(**_RETRY_ARGS)
    def load_namespace_properties(self, namespace: Union[str, Identifier]) -> Properties:
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace = NAMESPACE_SEPARATOR.join(namespace_tuple)
        response = self._session.get(self.url(Endpoints.load_namespace_metadata, namespace=namespace))
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError})

        return NamespaceResponse.model_validate_json(response.text).properties

    @retry(**_RETRY_ARGS)
    def update_namespace_properties(
        self, namespace: Union[str, Identifier], removals: Optional[Set[str]] = None, updates: Properties = EMPTY_DICT
    ) -> PropertiesUpdateSummary:
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace = NAMESPACE_SEPARATOR.join(namespace_tuple)
        payload = {"removals": list(removals or []), "updates": updates}
        response = self._session.post(self.url(Endpoints.update_namespace_properties, namespace=namespace), json=payload)
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError})
        parsed_response = UpdateNamespacePropertiesResponse.model_validate_json(response.text)
        return PropertiesUpdateSummary(
            removed=parsed_response.removed,
            updated=parsed_response.updated,
            missing=parsed_response.missing,
        )

    @retry(**_RETRY_ARGS)
    def namespace_exists(self, namespace: Union[str, Identifier]) -> bool:
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace = NAMESPACE_SEPARATOR.join(namespace_tuple)
        response = self._session.head(self.url(Endpoints.namespace_exists, namespace=namespace))

        if response.status_code == 404:
            return False
        elif response.status_code in (200, 204):
            return True

        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {})

        return False

    @retry(**_RETRY_ARGS)
    def table_exists(self, identifier: Union[str, Identifier]) -> bool:
        """Check if a table exists.

        Args:
            identifier (str | Identifier): Table identifier.

        Returns:
            bool: True if the table exists, False otherwise.
        """
        response = self._session.head(
            self.url(Endpoints.load_table, prefixed=True, **self._split_identifier_for_path(identifier))
        )

        if response.status_code == 404:
            return False
        elif response.status_code in (200, 204):
            return True

        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {})

        return False

    @retry(**_RETRY_ARGS)
    def view_exists(self, identifier: Union[str, Identifier]) -> bool:
        """Check if a view exists.

        Args:
            identifier (str | Identifier): View identifier.

        Returns:
            bool: True if the view exists, False otherwise.
        """
        response = self._session.head(
            self.url(Endpoints.view_exists, prefixed=True, **self._split_identifier_for_path(identifier, IdentifierKind.VIEW)),
        )
        if response.status_code == 404:
            return False
        elif response.status_code in [200, 204]:
            return True

        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {})

        return False

    @retry(**_RETRY_ARGS)
    def drop_view(self, identifier: Union[str]) -> None:
        response = self._session.delete(
            self.url(Endpoints.drop_view, prefixed=True, **self._split_identifier_for_path(identifier, IdentifierKind.VIEW)),
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchViewError})

    def close(self) -> None:
        """Close the catalog and release Session connection adapters.

        This method closes mounted HttpAdapters' pooled connections and any active Proxy pooled connections.
        """
        self._session.close()
