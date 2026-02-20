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
from collections import deque
from enum import Enum
from typing import (
    TYPE_CHECKING,
    Any,
    Union,
)
from urllib.parse import quote, unquote

from pydantic import ConfigDict, Field, TypeAdapter, field_validator
from requests import HTTPError, Session
from tenacity import RetryCallState, retry, retry_if_exception_type, stop_after_attempt

from pyiceberg import __version__
from pyiceberg.catalog import BOTOCORE_SESSION, TOKEN, URI, WAREHOUSE_LOCATION, Catalog, PropertiesUpdateSummary
from pyiceberg.catalog.rest.auth import AUTH_MANAGER, AuthManager, AuthManagerAdapter, AuthManagerFactory, LegacyOAuth2AuthManager
from pyiceberg.catalog.rest.response import _handle_non_200_response
from pyiceberg.catalog.rest.scan_planning import (
    FetchScanTasksRequest,
    PlanCancelled,
    PlanCompleted,
    PlanFailed,
    PlanningResponse,
    PlanSubmitted,
    PlanTableScanRequest,
    ScanTasks,
)
from pyiceberg.exceptions import (
    AuthorizationExpiredError,
    CommitFailedException,
    CommitStateUnknownException,
    NamespaceAlreadyExistsError,
    NamespaceNotEmptyError,
    NoSuchIdentifierError,
    NoSuchNamespaceError,
    NoSuchPlanTaskError,
    NoSuchTableError,
    NoSuchViewError,
    TableAlreadyExistsError,
    UnauthorizedError,
)
from pyiceberg.io import AWS_ACCESS_KEY_ID, AWS_REGION, AWS_SECRET_ACCESS_KEY, AWS_SESSION_TOKEN, FileIO, load_file_io
from pyiceberg.partitioning import UNPARTITIONED_PARTITION_SPEC, PartitionSpec, assign_fresh_partition_spec_ids
from pyiceberg.schema import Schema, assign_fresh_schema_ids
from pyiceberg.table import (
    CommitTableRequest,
    CommitTableResponse,
    CreateTableTransaction,
    FileScanTask,
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


class HttpMethod(str, Enum):
    GET = "GET"
    HEAD = "HEAD"
    POST = "POST"
    DELETE = "DELETE"


class Endpoint(IcebergBaseModel):
    model_config = ConfigDict(frozen=True)

    http_method: HttpMethod = Field()
    path: str = Field()

    @field_validator("path", mode="before")
    @classmethod
    def _validate_path(cls, raw_path: str) -> str:
        raw_path = raw_path.strip()
        if not raw_path:
            raise ValueError("Invalid path: empty")
        return raw_path

    def __str__(self) -> str:
        """Return the string representation of the Endpoint class."""
        return f"{self.http_method.value} {self.path}"

    @classmethod
    def from_string(cls, endpoint: str) -> "Endpoint":
        elements = endpoint.strip().split(None, 1)
        if len(elements) != 2:
            raise ValueError(f"Invalid endpoint (must consist of two elements separated by a single space): {endpoint}")
        return cls(http_method=HttpMethod(elements[0].upper()), path=elements[1])


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
    register_table: str = "namespaces/{namespace}/register"
    load_table: str = "namespaces/{namespace}/tables/{table}"
    update_table: str = "namespaces/{namespace}/tables/{table}"
    drop_table: str = "namespaces/{namespace}/tables/{table}"
    table_exists: str = "namespaces/{namespace}/tables/{table}"
    get_token: str = "oauth/tokens"
    rename_table: str = "tables/rename"
    list_views: str = "namespaces/{namespace}/views"
    drop_view: str = "namespaces/{namespace}/views/{view}"
    view_exists: str = "namespaces/{namespace}/views/{view}"
    plan_table_scan: str = "namespaces/{namespace}/tables/{table}/plan"
    fetch_scan_tasks: str = "namespaces/{namespace}/tables/{table}/tasks"


API_PREFIX = "/v1/{prefix}"


class Capability:
    V1_LIST_NAMESPACES = Endpoint(http_method=HttpMethod.GET, path=f"{API_PREFIX}/{Endpoints.list_namespaces}")
    V1_LOAD_NAMESPACE = Endpoint(http_method=HttpMethod.GET, path=f"{API_PREFIX}/{Endpoints.load_namespace_metadata}")
    V1_NAMESPACE_EXISTS = Endpoint(http_method=HttpMethod.HEAD, path=f"{API_PREFIX}/{Endpoints.namespace_exists}")
    V1_UPDATE_NAMESPACE = Endpoint(http_method=HttpMethod.POST, path=f"{API_PREFIX}/{Endpoints.update_namespace_properties}")
    V1_CREATE_NAMESPACE = Endpoint(http_method=HttpMethod.POST, path=f"{API_PREFIX}/{Endpoints.create_namespace}")
    V1_DELETE_NAMESPACE = Endpoint(http_method=HttpMethod.DELETE, path=f"{API_PREFIX}/{Endpoints.drop_namespace}")

    V1_LIST_TABLES = Endpoint(http_method=HttpMethod.GET, path=f"{API_PREFIX}/{Endpoints.list_tables}")
    V1_LOAD_TABLE = Endpoint(http_method=HttpMethod.GET, path=f"{API_PREFIX}/{Endpoints.load_table}")
    V1_TABLE_EXISTS = Endpoint(http_method=HttpMethod.HEAD, path=f"{API_PREFIX}/{Endpoints.table_exists}")
    V1_CREATE_TABLE = Endpoint(http_method=HttpMethod.POST, path=f"{API_PREFIX}/{Endpoints.create_table}")
    V1_UPDATE_TABLE = Endpoint(http_method=HttpMethod.POST, path=f"{API_PREFIX}/{Endpoints.update_table}")
    V1_DELETE_TABLE = Endpoint(http_method=HttpMethod.DELETE, path=f"{API_PREFIX}/{Endpoints.drop_table}")
    V1_RENAME_TABLE = Endpoint(http_method=HttpMethod.POST, path=f"{API_PREFIX}/{Endpoints.rename_table}")
    V1_REGISTER_TABLE = Endpoint(http_method=HttpMethod.POST, path=f"{API_PREFIX}/{Endpoints.register_table}")

    V1_LIST_VIEWS = Endpoint(http_method=HttpMethod.GET, path=f"{API_PREFIX}/{Endpoints.list_views}")
    V1_VIEW_EXISTS = Endpoint(http_method=HttpMethod.HEAD, path=f"{API_PREFIX}/{Endpoints.view_exists}")
    V1_DELETE_VIEW = Endpoint(http_method=HttpMethod.DELETE, path=f"{API_PREFIX}/{Endpoints.drop_view}")
    V1_SUBMIT_TABLE_SCAN_PLAN = Endpoint(http_method=HttpMethod.POST, path=f"{API_PREFIX}/{Endpoints.plan_table_scan}")
    V1_TABLE_SCAN_PLAN_TASKS = Endpoint(http_method=HttpMethod.POST, path=f"{API_PREFIX}/{Endpoints.fetch_scan_tasks}")


# Default endpoints for backwards compatibility with legacy servers that don't return endpoints
# in ConfigResponse. Only includes namespace and table endpoints.
DEFAULT_ENDPOINTS: frozenset[Endpoint] = frozenset(
    (
        Capability.V1_LIST_NAMESPACES,
        Capability.V1_LOAD_NAMESPACE,
        Capability.V1_CREATE_NAMESPACE,
        Capability.V1_UPDATE_NAMESPACE,
        Capability.V1_DELETE_NAMESPACE,
        Capability.V1_LIST_TABLES,
        Capability.V1_LOAD_TABLE,
        Capability.V1_CREATE_TABLE,
        Capability.V1_UPDATE_TABLE,
        Capability.V1_DELETE_TABLE,
        Capability.V1_RENAME_TABLE,
        Capability.V1_REGISTER_TABLE,
    )
)

# View endpoints conditionally added based on VIEW_ENDPOINTS_SUPPORTED property.
VIEW_ENDPOINTS: frozenset[Endpoint] = frozenset(
    (
        Capability.V1_LIST_VIEWS,
        Capability.V1_DELETE_VIEW,
    )
)


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
EMPTY_BODY_SHA256: str = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
OAUTH2_SERVER_URI = "oauth2-server-uri"
SNAPSHOT_LOADING_MODE = "snapshot-loading-mode"
AUTH = "auth"
CUSTOM = "custom"
REST_SCAN_PLANNING_ENABLED = "rest-scan-planning-enabled"
REST_SCAN_PLANNING_ENABLED_DEFAULT = False
# for backwards compatibility with older REST servers where it can be assumed that a particular
# server supports view endpoints but doesn't send the "endpoints" field in the ConfigResponse
VIEW_ENDPOINTS_SUPPORTED = "view-endpoints-supported"
VIEW_ENDPOINTS_SUPPORTED_DEFAULT = False

NAMESPACE_SEPARATOR_PROPERTY = "namespace-separator"
DEFAULT_NAMESPACE_SEPARATOR = b"\x1f".decode(UTF8)


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
    metadata_location: str | None = Field(alias="metadata-location", default=None)
    metadata: TableMetadata
    config: Properties = Field(default_factory=dict)


class CreateTableRequest(IcebergBaseModel):
    name: str = Field()
    location: str | None = Field()
    table_schema: Schema = Field(alias="schema")
    partition_spec: PartitionSpec | None = Field(alias="partition-spec")
    write_order: SortOrder | None = Field(alias="write-order")
    stage_create: bool = Field(alias="stage-create", default=False)
    properties: dict[str, str] = Field(default_factory=dict)

    # validators
    @field_validator("properties", mode="before")
    def transform_properties_dict_value_to_str(cls, properties: Properties) -> dict[str, str]:
        return transform_dict_value_to_str(properties)


class RegisterTableRequest(IcebergBaseModel):
    name: str
    metadata_location: str = Field(..., alias="metadata-location")


class ConfigResponse(IcebergBaseModel):
    defaults: Properties | None = Field(default_factory=dict)
    overrides: Properties | None = Field(default_factory=dict)
    endpoints: set[Endpoint] | None = Field(default=None)

    @field_validator("endpoints", mode="before")
    @classmethod
    def _parse_endpoints(cls, v: list[str] | None) -> set[Endpoint] | None:
        if v is None:
            return None
        return {Endpoint.from_string(s) for s in v}


class ListNamespaceResponse(IcebergBaseModel):
    namespaces: list[Identifier] = Field()


class NamespaceResponse(IcebergBaseModel):
    namespace: Identifier = Field()
    properties: Properties = Field()


class UpdateNamespacePropertiesResponse(IcebergBaseModel):
    removed: list[str] = Field()
    updated: list[str] = Field()
    missing: list[str] = Field()


class ListTableResponseEntry(IcebergBaseModel):
    name: str = Field()
    namespace: Identifier = Field()


class ListViewResponseEntry(IcebergBaseModel):
    name: str = Field()
    namespace: Identifier = Field()


class ListTablesResponse(IcebergBaseModel):
    identifiers: list[ListTableResponseEntry] = Field()


class ListViewsResponse(IcebergBaseModel):
    identifiers: list[ListViewResponseEntry] = Field()


_PLANNING_RESPONSE_ADAPTER = TypeAdapter(PlanningResponse)


class RestCatalog(Catalog):
    uri: str
    _session: Session
    _auth_manager: AuthManager | None
    _supported_endpoints: set[Endpoint]
    _namespace_separator: str

    def __init__(self, name: str, **properties: str):
        """Rest Catalog.

        You either need to provide a client_id and client_secret, or an already valid token.

        Args:
            name: Name to identify the catalog.
            properties: Properties that are passed along to the configuration.
        """
        super().__init__(name, **properties)
        self._auth_manager: AuthManager | None = None
        self.uri = properties[URI]
        self._fetch_config()
        self._session = self._create_session()

    def _create_session(self) -> Session:
        """Create a request session with provided catalog configuration."""
        session = Session()

        # Set HTTP headers
        self._config_headers(session)

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

            self._auth_manager = AuthManagerFactory.create(auth_impl or auth_type, auth_type_config)
            session.auth = AuthManagerAdapter(self._auth_manager)
        else:
            self._auth_manager = self._create_legacy_oauth2_auth_manager(session)
            session.auth = AuthManagerAdapter(self._auth_manager)

        # Configure SigV4 Request Signing
        if property_as_bool(self.properties, SIGV4, False):
            self._init_sigv4(session)

        return session

    def _load_file_io(self, properties: Properties = EMPTY_DICT, location: str | None = None) -> FileIO:
        merged_properties = {**self.properties, **properties}
        if self._auth_manager:
            merged_properties[AUTH_MANAGER] = self._auth_manager
        return load_file_io(merged_properties, location)

    def supports_server_side_planning(self) -> bool:
        """Check if the catalog supports server-side scan planning."""
        return Capability.V1_SUBMIT_TABLE_SCAN_PLAN in self._supported_endpoints and property_as_bool(
            self.properties, REST_SCAN_PLANNING_ENABLED, REST_SCAN_PLANNING_ENABLED_DEFAULT
        )

    @retry(**_RETRY_ARGS)
    def _plan_table_scan(self, identifier: str | Identifier, request: PlanTableScanRequest) -> PlanningResponse:
        """Submit a scan plan request to the REST server.

        Args:
            identifier: Table identifier.
            request: The scan plan request parameters.

        Returns:
            PlanningResponse the result of the scan plan request representing the status

        Raises:
            NoSuchTableError: If a table with the given identifier does not exist.
        """
        self._check_endpoint(Capability.V1_SUBMIT_TABLE_SCAN_PLAN)
        response = self._session.post(
            self.url(Endpoints.plan_table_scan, prefixed=True, **self._split_identifier_for_path(identifier)),
            data=request.model_dump_json(by_alias=True, exclude_none=True).encode(UTF8),
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchTableError})

        return _PLANNING_RESPONSE_ADAPTER.validate_json(response.text)

    @retry(**_RETRY_ARGS)
    def _fetch_scan_tasks(self, identifier: str | Identifier, plan_task: str) -> ScanTasks:
        """Fetch additional scan tasks using a plan task token.

        Args:
            identifier: Table identifier.
            plan_task: The plan task token from a previous response.

        Returns:
            ScanTasks containing file scan tasks and possibly more plan-task tokens.

        Raises:
            NoSuchPlanTaskError: If a plan task with the given identifier or task does not exist.
        """
        self._check_endpoint(Capability.V1_TABLE_SCAN_PLAN_TASKS)
        request = FetchScanTasksRequest(plan_task=plan_task)
        response = self._session.post(
            self.url(Endpoints.fetch_scan_tasks, prefixed=True, **self._split_identifier_for_path(identifier)),
            data=request.model_dump_json(by_alias=True).encode(UTF8),
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchPlanTaskError})

        return ScanTasks.model_validate_json(response.text)

    def plan_scan(self, identifier: str | Identifier, request: PlanTableScanRequest) -> list[FileScanTask]:
        """Plan a table scan and return FileScanTasks.

        Handles the full scan planning lifecycle including pagination.

        Args:
            identifier: Table identifier.
            request: The scan plan request parameters.

        Returns:
            List of FileScanTask objects ready for execution.

        Raises:
            RuntimeError: If planning fails, is cancelled, or returns unexpected response.
            NotImplementedError: If async planning is required but not yet supported.
        """
        response = self._plan_table_scan(identifier, request)

        if isinstance(response, PlanFailed):
            error_msg = response.error.message if response.error else "unknown error"
            raise RuntimeError(f"Received status: failed: {error_msg}")

        if isinstance(response, PlanCancelled):
            raise RuntimeError("Received status: cancelled")

        if isinstance(response, PlanSubmitted):
            # TODO: implement polling for async planning
            raise NotImplementedError(f"Async scan planning not yet supported for planId: {response.plan_id}")

        if not isinstance(response, PlanCompleted):
            raise RuntimeError(f"Invalid planStatus for response: {type(response).__name__}")

        tasks: list[FileScanTask] = []

        # Collect tasks from initial response
        for task in response.file_scan_tasks:
            tasks.append(FileScanTask.from_rest_response(task, response.delete_files))

        # Fetch and collect from additional batches
        pending_tasks = deque(response.plan_tasks)
        while pending_tasks:
            plan_task = pending_tasks.popleft()
            batch = self._fetch_scan_tasks(identifier, plan_task)
            for task in batch.file_scan_tasks:
                tasks.append(FileScanTask.from_rest_response(task, batch.delete_files))
            pending_tasks.extend(batch.plan_tasks)

        return tasks

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

    def _check_valid_namespace_identifier(self, identifier: str | Identifier) -> Identifier:
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

    def _check_endpoint(self, endpoint: Endpoint) -> None:
        """Check if an endpoint is supported by the server.

        Args:
            endpoint: The endpoint to check against the set of supported endpoints

        Raises:
            NotImplementedError: If the endpoint is not supported.
        """
        if endpoint not in self._supported_endpoints:
            raise NotImplementedError(f"Server does not support endpoint: {endpoint}")

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

    def _extract_optional_oauth_params(self) -> dict[str, str]:
        optional_oauth_param = {SCOPE: self.properties.get(SCOPE) or CATALOG_SCOPE}
        set_of_optional_params = {AUDIENCE, RESOURCE}
        for param in set_of_optional_params:
            if param_value := self.properties.get(param):
                optional_oauth_param[param] = param_value

        return optional_oauth_param

    def _encode_namespace_path(self, namespace: Identifier) -> str:
        """
        Encode a namespace for use as a path parameter in a URL.

        Each part of the namespace is URL-encoded using `urllib.parse.quote`
        (ensuring characters like '/' are encoded) and then joined by the
        configured namespace separator.
        """
        return self._namespace_separator.join(quote(part, safe="") for part in namespace)

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

        # Determine supported endpoints
        endpoints = config_response.endpoints
        if endpoints:
            self._supported_endpoints = set(endpoints)
        else:
            # Use default endpoints for legacy servers that don't return endpoints
            self._supported_endpoints = set(DEFAULT_ENDPOINTS)
            # Conditionally add view endpoints based on config
            if property_as_bool(self.properties, VIEW_ENDPOINTS_SUPPORTED, VIEW_ENDPOINTS_SUPPORTED_DEFAULT):
                self._supported_endpoints.update(VIEW_ENDPOINTS)

        separator_from_properties = self.properties.get(NAMESPACE_SEPARATOR_PROPERTY, DEFAULT_NAMESPACE_SEPARATOR)
        if not separator_from_properties:
            raise ValueError("Namespace separator cannot be an empty string")
        self._namespace_separator = unquote(separator_from_properties)

    def _identifier_to_validated_tuple(self, identifier: str | Identifier) -> Identifier:
        identifier_tuple = self.identifier_to_tuple(identifier)
        if len(identifier_tuple) <= 1:
            raise NoSuchIdentifierError(f"Missing namespace or invalid identifier: {'.'.join(identifier_tuple)}")
        return identifier_tuple

    def _split_identifier_for_path(
        self, identifier: str | Identifier | TableIdentifier, kind: IdentifierKind = IdentifierKind.TABLE
    ) -> Properties:
        if isinstance(identifier, TableIdentifier):
            return {
                "namespace": self._encode_namespace_path(tuple(identifier.namespace.root)),
                kind.value: quote(identifier.name, safe=""),
            }
        identifier_tuple = self._identifier_to_validated_tuple(identifier)

        # Use quote to ensure that '/' aren't treated as path separators.
        return {
            "namespace": self._encode_namespace_path(identifier_tuple[:-1]),
            kind.value: quote(identifier_tuple[-1], safe=""),
        }

    def _split_identifier_for_json(self, identifier: str | Identifier) -> dict[str, Identifier | str]:
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
                if "connection" in request.headers:
                    del request.headers["connection"]
                # For empty bodies, explicitly set the content hash header to the SHA256 of an empty string
                if not request.body:
                    request.headers["x-amz-content-sha256"] = EMPTY_BODY_SHA256

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

    def _response_to_table(self, identifier_tuple: tuple[str, ...], table_response: TableResponse) -> Table:
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

    def _response_to_staged_table(self, identifier_tuple: tuple[str, ...], table_response: TableResponse) -> StagedTable:
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
        session.headers["X-Client-Version"] = f"PyIceberg {__version__}"
        session.headers.setdefault("X-Iceberg-Access-Delegation", ACCESS_DELEGATION_DEFAULT)

    def _create_table(
        self,
        identifier: str | Identifier,
        schema: Union[Schema, "pa.Schema"],
        location: str | None = None,
        partition_spec: PartitionSpec = UNPARTITIONED_PARTITION_SPEC,
        sort_order: SortOrder = UNSORTED_SORT_ORDER,
        properties: Properties = EMPTY_DICT,
        stage_create: bool = False,
    ) -> TableResponse:
        self._check_endpoint(Capability.V1_CREATE_TABLE)
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
            name=self._identifier_to_validated_tuple(identifier)[-1],
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
        identifier: str | Identifier,
        schema: Union[Schema, "pa.Schema"],
        location: str | None = None,
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
        identifier: str | Identifier,
        schema: Union[Schema, "pa.Schema"],
        location: str | None = None,
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
        self._check_endpoint(Capability.V1_REGISTER_TABLE)
        namespace_and_table = self._split_identifier_for_path(identifier)
        request = RegisterTableRequest(
            name=self._identifier_to_validated_tuple(identifier)[-1],
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
    def list_tables(self, namespace: str | Identifier) -> list[Identifier]:
        self._check_endpoint(Capability.V1_LIST_TABLES)
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace_concat = self._encode_namespace_path(namespace_tuple)
        response = self._session.get(self.url(Endpoints.list_tables, namespace=namespace_concat))
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError})
        return [(*table.namespace, table.name) for table in ListTablesResponse.model_validate_json(response.text).identifiers]

    @retry(**_RETRY_ARGS)
    def load_table(self, identifier: str | Identifier) -> Table:
        self._check_endpoint(Capability.V1_LOAD_TABLE)
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
    def drop_table(self, identifier: str | Identifier, purge_requested: bool = False) -> None:
        self._check_endpoint(Capability.V1_DELETE_TABLE)
        response = self._session.delete(
            self.url(Endpoints.drop_table, prefixed=True, **self._split_identifier_for_path(identifier)),
            params={"purgeRequested": purge_requested},
        )
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchTableError})

    @retry(**_RETRY_ARGS)
    def purge_table(self, identifier: str | Identifier) -> None:
        self.drop_table(identifier=identifier, purge_requested=True)

    @retry(**_RETRY_ARGS)
    def rename_table(self, from_identifier: str | Identifier, to_identifier: str | Identifier) -> Table:
        self._check_endpoint(Capability.V1_RENAME_TABLE)
        payload = {
            "source": self._split_identifier_for_json(from_identifier),
            "destination": self._split_identifier_for_json(to_identifier),
        }

        # Ensure that namespaces exist on source and destination.
        source_namespace = self._split_identifier_for_json(from_identifier)["namespace"]
        if not self.namespace_exists(source_namespace):
            raise NoSuchNamespaceError(f"Source namespace does not exist: {source_namespace}")

        destination_namespace = self._split_identifier_for_json(to_identifier)["namespace"]
        if not self.namespace_exists(destination_namespace):
            raise NoSuchNamespaceError(f"Destination namespace does not exist: {destination_namespace}")

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
    def list_views(self, namespace: str | Identifier) -> list[Identifier]:
        if Capability.V1_LIST_VIEWS not in self._supported_endpoints:
            return []
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace_concat = self._encode_namespace_path(namespace_tuple)
        response = self._session.get(self.url(Endpoints.list_views, namespace=namespace_concat))
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError})
        return [(*view.namespace, view.name) for view in ListViewsResponse.model_validate_json(response.text).identifiers]

    @retry(**_RETRY_ARGS)
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
        self._check_endpoint(Capability.V1_UPDATE_TABLE)
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
    def create_namespace(self, namespace: str | Identifier, properties: Properties = EMPTY_DICT) -> None:
        self._check_endpoint(Capability.V1_CREATE_NAMESPACE)
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        payload = {"namespace": namespace_tuple, "properties": properties}
        response = self._session.post(self.url(Endpoints.create_namespace), json=payload)
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {409: NamespaceAlreadyExistsError})

    @retry(**_RETRY_ARGS)
    def drop_namespace(self, namespace: str | Identifier) -> None:
        self._check_endpoint(Capability.V1_DELETE_NAMESPACE)
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace = self._encode_namespace_path(namespace_tuple)
        response = self._session.delete(self.url(Endpoints.drop_namespace, namespace=namespace))
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError, 409: NamespaceNotEmptyError})

    @retry(**_RETRY_ARGS)
    def list_namespaces(self, namespace: str | Identifier = ()) -> list[Identifier]:
        self._check_endpoint(Capability.V1_LIST_NAMESPACES)
        namespace_tuple = self.identifier_to_tuple(namespace)
        response = self._session.get(
            self.url(
                f"{Endpoints.list_namespaces}?parent={self._encode_namespace_path(namespace_tuple)}"
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
    def load_namespace_properties(self, namespace: str | Identifier) -> Properties:
        self._check_endpoint(Capability.V1_LOAD_NAMESPACE)
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace = self._encode_namespace_path(namespace_tuple)
        response = self._session.get(self.url(Endpoints.load_namespace_metadata, namespace=namespace))
        try:
            response.raise_for_status()
        except HTTPError as exc:
            _handle_non_200_response(exc, {404: NoSuchNamespaceError})

        return NamespaceResponse.model_validate_json(response.text).properties

    @retry(**_RETRY_ARGS)
    def update_namespace_properties(
        self, namespace: str | Identifier, removals: set[str] | None = None, updates: Properties = EMPTY_DICT
    ) -> PropertiesUpdateSummary:
        self._check_endpoint(Capability.V1_UPDATE_NAMESPACE)
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace = self._encode_namespace_path(namespace_tuple)
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
    def namespace_exists(self, namespace: str | Identifier) -> bool:
        namespace_tuple = self._check_valid_namespace_identifier(namespace)
        namespace = self._encode_namespace_path(namespace_tuple)

        # fallback in order to work with older rest catalog implementations
        if Capability.V1_NAMESPACE_EXISTS not in self._supported_endpoints:
            try:
                self.load_namespace_properties(namespace_tuple)
                return True
            except NoSuchNamespaceError:
                return False

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
    def table_exists(self, identifier: str | Identifier) -> bool:
        """Check if a table exists.

        Args:
            identifier (str | Identifier): Table identifier.

        Returns:
            bool: True if the table exists, False otherwise.
        """
        # fallback in order to work with older rest catalog implementations
        if Capability.V1_TABLE_EXISTS not in self._supported_endpoints:
            try:
                self.load_table(identifier)
                return True
            except NoSuchTableError:
                return False

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
    def view_exists(self, identifier: str | Identifier) -> bool:
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
    def drop_view(self, identifier: str) -> None:
        self._check_endpoint(Capability.V1_DELETE_VIEW)
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
