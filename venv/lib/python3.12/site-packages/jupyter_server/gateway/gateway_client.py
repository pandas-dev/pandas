"""A kernel gateway client."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import asyncio
import json
import logging
import os
import typing as ty
from abc import ABC, ABCMeta, abstractmethod
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from http.cookies import SimpleCookie
from socket import gaierror

from jupyter_events import EventLogger
from tornado import web
from tornado.httpclient import AsyncHTTPClient, HTTPClientError, HTTPResponse
from traitlets import (
    Bool,
    Float,
    Instance,
    Int,
    TraitError,
    Type,
    Unicode,
    default,
    observe,
    validate,
)
from traitlets.config import LoggingConfigurable, SingletonConfigurable

from jupyter_server import DEFAULT_EVENTS_SCHEMA_PATH, JUPYTER_SERVER_EVENTS_URI

ERROR_STATUS = "error"
SUCCESS_STATUS = "success"
STATUS_KEY = "status"
STATUS_CODE_KEY = "status_code"
MESSAGE_KEY = "msg"

if ty.TYPE_CHECKING:
    from http.cookies import Morsel


class GatewayTokenRenewerMeta(ABCMeta, type(LoggingConfigurable)):  # type: ignore[misc]
    """The metaclass necessary for proper ABC behavior in a Configurable."""


class GatewayTokenRenewerBase(  # type:ignore[misc]
    ABC, LoggingConfigurable, metaclass=GatewayTokenRenewerMeta
):
    """
    Abstract base class for refreshing tokens used between this server and a Gateway
    server.  Implementations requiring additional configuration can extend their class
    with appropriate configuration values or convey those values via appropriate
    environment variables relative to the implementation.
    """

    @abstractmethod
    def get_token(
        self,
        auth_header_key: str,
        auth_scheme: ty.Union[str, None],
        auth_token: str,
        **kwargs: ty.Any,
    ) -> str:
        """
        Given the current authorization header key, scheme, and token, this method returns
        a (potentially renewed) token for use against the Gateway server.
        """


class NoOpTokenRenewer(GatewayTokenRenewerBase):  # type:ignore[misc]
    """NoOpTokenRenewer is the default value to the GatewayClient trait
    `gateway_token_renewer` and merely returns the provided token.
    """

    def get_token(
        self,
        auth_header_key: str,
        auth_scheme: ty.Union[str, None],
        auth_token: str,
        **kwargs: ty.Any,
    ) -> str:
        """This implementation simply returns the current authorization token."""
        return auth_token


class GatewayClient(SingletonConfigurable):
    """This class manages the configuration.  It's its own singleton class so
    that we can share these values across all objects.  It also contains some
    options.
    helper methods to build request arguments out of the various config
    """

    event_schema_id = JUPYTER_SERVER_EVENTS_URI + "/gateway_client/v1"
    event_logger = Instance(EventLogger).tag(config=True)

    @default("event_logger")
    def _default_event_logger(self):
        if self.parent and hasattr(self.parent, "event_logger"):
            # Event logger is attached from serverapp.
            return self.parent.event_logger
        else:
            # If parent does not have an event logger, create one.
            logger = EventLogger()
            schema_path = DEFAULT_EVENTS_SCHEMA_PATH / "gateway_client" / "v1.yaml"
            logger.register_event_schema(schema_path)
            self.log.info("Event is registered in GatewayClient.")
            return logger

    def emit(self, data):
        """Emit event using the core event schema from Jupyter Server's Gateway Client."""
        self.event_logger.emit(schema_id=self.event_schema_id, data=data)

    url = Unicode(
        default_value=None,
        allow_none=True,
        config=True,
        help="""The url of the Kernel or Enterprise Gateway server where
kernel specifications are defined and kernel management takes place.
If defined, this Notebook server acts as a proxy for all kernel
management and kernel specification retrieval.  (JUPYTER_GATEWAY_URL env var)
        """,
    )

    url_env = "JUPYTER_GATEWAY_URL"

    @default("url")
    def _url_default(self):
        return os.environ.get(self.url_env)

    @validate("url")
    def _url_validate(self, proposal):
        value = proposal["value"]
        # Ensure value, if present, starts with 'http'
        if value is not None and len(value) > 0 and not str(value).lower().startswith("http"):
            message = "GatewayClient url must start with 'http': '%r'" % value
            self.emit(data={STATUS_KEY: ERROR_STATUS, STATUS_CODE_KEY: 400, MESSAGE_KEY: message})
            raise TraitError(message)
        return value

    ws_url = Unicode(
        default_value=None,
        allow_none=True,
        config=True,
        help="""The websocket url of the Kernel or Enterprise Gateway server.  If not provided, this value
will correspond to the value of the Gateway url with 'ws' in place of 'http'.  (JUPYTER_GATEWAY_WS_URL env var)
        """,
    )

    ws_url_env = "JUPYTER_GATEWAY_WS_URL"

    @default("ws_url")
    def _ws_url_default(self):
        default_value = os.environ.get(self.ws_url_env)
        if self.url is not None and default_value is None and self.gateway_enabled:
            default_value = self.url.lower().replace("http", "ws")
        return default_value

    @validate("ws_url")
    def _ws_url_validate(self, proposal):
        value = proposal["value"]
        # Ensure value, if present, starts with 'ws'
        if value is not None and len(value) > 0 and not str(value).lower().startswith("ws"):
            message = "GatewayClient ws_url must start with 'ws': '%r'" % value
            self.emit(data={STATUS_KEY: ERROR_STATUS, STATUS_CODE_KEY: 400, MESSAGE_KEY: message})
            raise TraitError(message)
        return value

    kernels_endpoint_default_value = "/api/kernels"
    kernels_endpoint_env = "JUPYTER_GATEWAY_KERNELS_ENDPOINT"
    kernels_endpoint = Unicode(
        default_value=kernels_endpoint_default_value,
        config=True,
        help="""The gateway API endpoint for accessing kernel resources (JUPYTER_GATEWAY_KERNELS_ENDPOINT env var)""",
    )

    @default("kernels_endpoint")
    def _kernels_endpoint_default(self):
        return os.environ.get(self.kernels_endpoint_env, self.kernels_endpoint_default_value)

    kernelspecs_endpoint_default_value = "/api/kernelspecs"
    kernelspecs_endpoint_env = "JUPYTER_GATEWAY_KERNELSPECS_ENDPOINT"
    kernelspecs_endpoint = Unicode(
        default_value=kernelspecs_endpoint_default_value,
        config=True,
        help="""The gateway API endpoint for accessing kernelspecs (JUPYTER_GATEWAY_KERNELSPECS_ENDPOINT env var)""",
    )

    @default("kernelspecs_endpoint")
    def _kernelspecs_endpoint_default(self):
        return os.environ.get(
            self.kernelspecs_endpoint_env, self.kernelspecs_endpoint_default_value
        )

    kernelspecs_resource_endpoint_default_value = "/kernelspecs"
    kernelspecs_resource_endpoint_env = "JUPYTER_GATEWAY_KERNELSPECS_RESOURCE_ENDPOINT"
    kernelspecs_resource_endpoint = Unicode(
        default_value=kernelspecs_resource_endpoint_default_value,
        config=True,
        help="""The gateway endpoint for accessing kernelspecs resources
(JUPYTER_GATEWAY_KERNELSPECS_RESOURCE_ENDPOINT env var)""",
    )

    @default("kernelspecs_resource_endpoint")
    def _kernelspecs_resource_endpoint_default(self):
        return os.environ.get(
            self.kernelspecs_resource_endpoint_env,
            self.kernelspecs_resource_endpoint_default_value,
        )

    connect_timeout_default_value = 40.0
    connect_timeout_env = "JUPYTER_GATEWAY_CONNECT_TIMEOUT"
    connect_timeout = Float(
        default_value=connect_timeout_default_value,
        config=True,
        help="""The time allowed for HTTP connection establishment with the Gateway server.
(JUPYTER_GATEWAY_CONNECT_TIMEOUT env var)""",
    )

    @default("connect_timeout")
    def _connect_timeout_default(self):
        return float(os.environ.get(self.connect_timeout_env, self.connect_timeout_default_value))

    request_timeout_default_value = 42.0
    request_timeout_env = "JUPYTER_GATEWAY_REQUEST_TIMEOUT"
    request_timeout = Float(
        default_value=request_timeout_default_value,
        config=True,
        help="""The time allowed for HTTP request completion. (JUPYTER_GATEWAY_REQUEST_TIMEOUT env var)""",
    )

    @default("request_timeout")
    def _request_timeout_default(self):
        return float(os.environ.get(self.request_timeout_env, self.request_timeout_default_value))

    client_key = Unicode(
        default_value=None,
        allow_none=True,
        config=True,
        help="""The filename for client SSL key, if any.  (JUPYTER_GATEWAY_CLIENT_KEY env var)
        """,
    )
    client_key_env = "JUPYTER_GATEWAY_CLIENT_KEY"

    @default("client_key")
    def _client_key_default(self):
        return os.environ.get(self.client_key_env)

    client_cert = Unicode(
        default_value=None,
        allow_none=True,
        config=True,
        help="""The filename for client SSL certificate, if any.  (JUPYTER_GATEWAY_CLIENT_CERT env var)
        """,
    )
    client_cert_env = "JUPYTER_GATEWAY_CLIENT_CERT"

    @default("client_cert")
    def _client_cert_default(self):
        return os.environ.get(self.client_cert_env)

    ca_certs = Unicode(
        default_value=None,
        allow_none=True,
        config=True,
        help="""The filename of CA certificates or None to use defaults.  (JUPYTER_GATEWAY_CA_CERTS env var)
        """,
    )
    ca_certs_env = "JUPYTER_GATEWAY_CA_CERTS"

    @default("ca_certs")
    def _ca_certs_default(self):
        return os.environ.get(self.ca_certs_env)

    http_user = Unicode(
        default_value=None,
        allow_none=True,
        config=True,
        help="""The username for HTTP authentication. (JUPYTER_GATEWAY_HTTP_USER env var)
        """,
    )
    http_user_env = "JUPYTER_GATEWAY_HTTP_USER"

    @default("http_user")
    def _http_user_default(self):
        return os.environ.get(self.http_user_env)

    http_pwd = Unicode(
        default_value=None,
        allow_none=True,
        config=True,
        help="""The password for HTTP authentication.  (JUPYTER_GATEWAY_HTTP_PWD env var)
        """,
    )
    http_pwd_env = "JUPYTER_GATEWAY_HTTP_PWD"  # noqa: S105

    @default("http_pwd")
    def _http_pwd_default(self):
        return os.environ.get(self.http_pwd_env)

    headers_default_value = "{}"
    headers_env = "JUPYTER_GATEWAY_HEADERS"
    headers = Unicode(
        default_value=headers_default_value,
        allow_none=True,
        config=True,
        help="""Additional HTTP headers to pass on the request.  This value will be converted to a dict.
          (JUPYTER_GATEWAY_HEADERS env var)
        """,
    )

    @default("headers")
    def _headers_default(self):
        return os.environ.get(self.headers_env, self.headers_default_value)

    auth_header_key_default_value = "Authorization"
    auth_header_key = Unicode(
        config=True,
        help="""The authorization header's key name (typically 'Authorization') used in the HTTP headers. The
header will be formatted as::

{'{auth_header_key}': '{auth_scheme} {auth_token}'}

If the authorization header key takes a single value, `auth_scheme` should be set to None and
'auth_token' should be configured to use the appropriate value.

(JUPYTER_GATEWAY_AUTH_HEADER_KEY env var)""",
    )
    auth_header_key_env = "JUPYTER_GATEWAY_AUTH_HEADER_KEY"

    @default("auth_header_key")
    def _auth_header_key_default(self):
        return os.environ.get(self.auth_header_key_env, self.auth_header_key_default_value)

    auth_token_default_value = ""
    auth_token = Unicode(
        default_value=None,
        allow_none=True,
        config=True,
        help="""The authorization token used in the HTTP headers. The header will be formatted as::

{'{auth_header_key}': '{auth_scheme} {auth_token}'}

(JUPYTER_GATEWAY_AUTH_TOKEN env var)""",
    )
    auth_token_env = "JUPYTER_GATEWAY_AUTH_TOKEN"  # noqa: S105

    @default("auth_token")
    def _auth_token_default(self):
        return os.environ.get(self.auth_token_env, self.auth_token_default_value)

    auth_scheme_default_value = "token"  # This value is purely for backwards compatibility
    auth_scheme = Unicode(
        allow_none=True,
        config=True,
        help="""The auth scheme, added as a prefix to the authorization token used in the HTTP headers.
(JUPYTER_GATEWAY_AUTH_SCHEME env var)""",
    )
    auth_scheme_env = "JUPYTER_GATEWAY_AUTH_SCHEME"

    @default("auth_scheme")
    def _auth_scheme_default(self):
        return os.environ.get(self.auth_scheme_env, self.auth_scheme_default_value)

    validate_cert_default_value = True
    validate_cert_env = "JUPYTER_GATEWAY_VALIDATE_CERT"
    validate_cert = Bool(
        default_value=validate_cert_default_value,
        config=True,
        help="""For HTTPS requests, determines if server's certificate should be validated or not.
(JUPYTER_GATEWAY_VALIDATE_CERT env var)""",
    )

    @default("validate_cert")
    def _validate_cert_default(self):
        return bool(
            os.environ.get(self.validate_cert_env, str(self.validate_cert_default_value))
            not in ["no", "false"]
        )

    allowed_envs_default_value = ""
    allowed_envs_env = "JUPYTER_GATEWAY_ALLOWED_ENVS"
    allowed_envs = Unicode(
        default_value=allowed_envs_default_value,
        config=True,
        help="""A comma-separated list of environment variable names that will be included, along with
their values, in the kernel startup request.  The corresponding `client_envs` configuration
value must also be set on the Gateway server - since that configuration value indicates which
environmental values to make available to the kernel. (JUPYTER_GATEWAY_ALLOWED_ENVS env var)""",
    )

    @default("allowed_envs")
    def _allowed_envs_default(self):
        return os.environ.get(
            self.allowed_envs_env,
            os.environ.get("JUPYTER_GATEWAY_ENV_WHITELIST", self.allowed_envs_default_value),
        )

    env_whitelist = Unicode(
        default_value=allowed_envs_default_value,
        config=True,
        help="""Deprecated, use `GatewayClient.allowed_envs`""",
    )

    gateway_retry_interval_default_value = 1.0
    gateway_retry_interval_env = "JUPYTER_GATEWAY_RETRY_INTERVAL"
    gateway_retry_interval = Float(
        default_value=gateway_retry_interval_default_value,
        config=True,
        help="""The time allowed for HTTP reconnection with the Gateway server for the first time.
Next will be JUPYTER_GATEWAY_RETRY_INTERVAL multiplied by two in factor of numbers of retries
but less than JUPYTER_GATEWAY_RETRY_INTERVAL_MAX.
(JUPYTER_GATEWAY_RETRY_INTERVAL env var)""",
    )

    @default("gateway_retry_interval")
    def _gateway_retry_interval_default(self):
        return float(
            os.environ.get(
                self.gateway_retry_interval_env,
                self.gateway_retry_interval_default_value,
            )
        )

    gateway_retry_interval_max_default_value = 30.0
    gateway_retry_interval_max_env = "JUPYTER_GATEWAY_RETRY_INTERVAL_MAX"
    gateway_retry_interval_max = Float(
        default_value=gateway_retry_interval_max_default_value,
        config=True,
        help="""The maximum time allowed for HTTP reconnection retry with the Gateway server.
(JUPYTER_GATEWAY_RETRY_INTERVAL_MAX env var)""",
    )

    @default("gateway_retry_interval_max")
    def _gateway_retry_interval_max_default(self):
        return float(
            os.environ.get(
                self.gateway_retry_interval_max_env,
                self.gateway_retry_interval_max_default_value,
            )
        )

    gateway_retry_max_default_value = 5
    gateway_retry_max_env = "JUPYTER_GATEWAY_RETRY_MAX"
    gateway_retry_max = Int(
        default_value=gateway_retry_max_default_value,
        config=True,
        help="""The maximum retries allowed for HTTP reconnection with the Gateway server.
(JUPYTER_GATEWAY_RETRY_MAX env var)""",
    )

    @default("gateway_retry_max")
    def _gateway_retry_max_default(self):
        return int(os.environ.get(self.gateway_retry_max_env, self.gateway_retry_max_default_value))

    gateway_token_renewer_class_default_value = (
        "jupyter_server.gateway.gateway_client.NoOpTokenRenewer"  # noqa: S105
    )
    gateway_token_renewer_class_env = "JUPYTER_GATEWAY_TOKEN_RENEWER_CLASS"  # noqa: S105
    gateway_token_renewer_class = Type(
        klass=GatewayTokenRenewerBase,
        config=True,
        help="""The class to use for Gateway token renewal. (JUPYTER_GATEWAY_TOKEN_RENEWER_CLASS env var)""",
    )

    @default("gateway_token_renewer_class")
    def _gateway_token_renewer_class_default(self):
        return os.environ.get(
            self.gateway_token_renewer_class_env, self.gateway_token_renewer_class_default_value
        )

    launch_timeout_pad_default_value = 2.0
    launch_timeout_pad_env = "JUPYTER_GATEWAY_LAUNCH_TIMEOUT_PAD"
    launch_timeout_pad = Float(
        default_value=launch_timeout_pad_default_value,
        config=True,
        help="""Timeout pad to be ensured between KERNEL_LAUNCH_TIMEOUT and request_timeout
such that request_timeout >= KERNEL_LAUNCH_TIMEOUT + launch_timeout_pad.
(JUPYTER_GATEWAY_LAUNCH_TIMEOUT_PAD env var)""",
    )

    @default("launch_timeout_pad")
    def _launch_timeout_pad_default(self):
        return float(
            os.environ.get(
                self.launch_timeout_pad_env,
                self.launch_timeout_pad_default_value,
            )
        )

    accept_cookies_value = False
    accept_cookies_env = "JUPYTER_GATEWAY_ACCEPT_COOKIES"
    accept_cookies = Bool(
        default_value=accept_cookies_value,
        config=True,
        help="""Accept and manage cookies sent by the service side. This is often useful
        for load balancers to decide which backend node to use.
        (JUPYTER_GATEWAY_ACCEPT_COOKIES env var)""",
    )

    @default("accept_cookies")
    def _accept_cookies_default(self):
        return bool(
            os.environ.get(self.accept_cookies_env, str(self.accept_cookies_value).lower())
            not in ["no", "false"]
        )

    _deprecated_traits = {
        "env_whitelist": ("allowed_envs", "2.0"),
    }

    # Method copied from
    # https://github.com/jupyterhub/jupyterhub/blob/d1a85e53dccfc7b1dd81b0c1985d158cc6b61820/jupyterhub/auth.py#L143-L161
    @observe(*list(_deprecated_traits))
    def _deprecated_trait(self, change):
        """observer for deprecated traits"""
        old_attr = change.name
        new_attr, version = self._deprecated_traits[old_attr]
        new_value = getattr(self, new_attr)
        if new_value != change.new:
            # only warn if different
            # protects backward-compatible config from warnings
            # if they set the same value under both names
            self.log.warning(
                f"{self.__class__.__name__}.{old_attr} is deprecated in jupyter_server "
                f"{version}, use {self.__class__.__name__}.{new_attr} instead"
            )
            setattr(self, new_attr, change.new)

    @property
    def gateway_enabled(self):
        return bool(self.url is not None and len(self.url) > 0)

    # Ensure KERNEL_LAUNCH_TIMEOUT has a default value.
    KERNEL_LAUNCH_TIMEOUT = int(os.environ.get("KERNEL_LAUNCH_TIMEOUT", 40))

    _connection_args: dict[str, ty.Any]  # initialized on first use

    gateway_token_renewer: GatewayTokenRenewerBase

    def __init__(self, **kwargs):
        """Initialize a gateway client."""
        super().__init__(**kwargs)
        self._connection_args = {}  # initialized on first use
        self.gateway_token_renewer = self.gateway_token_renewer_class(parent=self, log=self.log)  # type:ignore[abstract]

        # store of cookies with store time
        self._cookies: dict[str, tuple[Morsel[ty.Any], datetime]] = {}

    def init_connection_args(self):
        """Initialize arguments used on every request.  Since these are primarily static values,
        we'll perform this operation once.
        """
        # Ensure that request timeout and KERNEL_LAUNCH_TIMEOUT are in sync, taking the
        #  greater value of the two and taking into account the following relation:
        #  request_timeout = KERNEL_LAUNCH_TIME + padding
        minimum_request_timeout = (
            float(GatewayClient.KERNEL_LAUNCH_TIMEOUT) + self.launch_timeout_pad
        )
        if self.request_timeout < minimum_request_timeout:
            self.request_timeout = minimum_request_timeout
        elif self.request_timeout > minimum_request_timeout:
            GatewayClient.KERNEL_LAUNCH_TIMEOUT = int(
                self.request_timeout - self.launch_timeout_pad
            )
        # Ensure any adjustments are reflected in env.
        os.environ["KERNEL_LAUNCH_TIMEOUT"] = str(GatewayClient.KERNEL_LAUNCH_TIMEOUT)

        if self.headers:
            self._connection_args["headers"] = json.loads(self.headers)
            if self.auth_header_key not in self._connection_args["headers"]:
                self._connection_args["headers"].update(
                    {f"{self.auth_header_key}": f"{self.auth_scheme} {self.auth_token}"}
                )
        self._connection_args["connect_timeout"] = self.connect_timeout
        self._connection_args["request_timeout"] = self.request_timeout
        self._connection_args["validate_cert"] = self.validate_cert
        if self.client_cert:
            self._connection_args["client_cert"] = self.client_cert
            self._connection_args["client_key"] = self.client_key
            if self.ca_certs:
                self._connection_args["ca_certs"] = self.ca_certs
        if self.http_user:
            self._connection_args["auth_username"] = self.http_user
        if self.http_pwd:
            self._connection_args["auth_password"] = self.http_pwd

    def load_connection_args(self, **kwargs):
        """Merges the static args relative to the connection, with the given keyword arguments.  If static
        args have yet to be initialized, we'll do that here.

        """
        if len(self._connection_args) == 0:
            self.init_connection_args()

        # Give token renewal a shot at renewing the token
        prev_auth_token = self.auth_token
        if self.auth_token is not None:
            try:
                self.auth_token = self.gateway_token_renewer.get_token(
                    self.auth_header_key, self.auth_scheme, self.auth_token
                )
            except Exception as ex:
                self.log.error(
                    f"An exception occurred attempting to renew the "
                    f"Gateway authorization token using an instance of class "
                    f"'{self.gateway_token_renewer_class}'.  The request will "
                    f"proceed using the current token value.  Exception was: {ex}"
                )
                self.auth_token = prev_auth_token

        for arg, value in self._connection_args.items():
            if arg == "headers":
                given_value = kwargs.setdefault(arg, {})
                if isinstance(given_value, dict):
                    given_value.update(value)
                    # Ensure the auth header is current
                    given_value.update(
                        {f"{self.auth_header_key}": f"{self.auth_scheme} {self.auth_token}"}
                    )
            else:
                kwargs[arg] = value

        if self.accept_cookies:
            self._update_cookie_header(kwargs)

        return kwargs

    def update_cookies(self, cookie: SimpleCookie) -> None:
        """Update cookies from existing requests for load balancers"""
        if not self.accept_cookies:
            return

        store_time = datetime.now(tz=timezone.utc)
        for key, item in cookie.items():
            # Convert "expires" arg into "max-age" to facilitate expiration management.
            # As "max-age" has precedence, ignore "expires" when "max-age" exists.
            if item.get("expires") and not item.get("max-age"):
                expire_timedelta = parsedate_to_datetime(item["expires"]) - store_time
                item["max-age"] = str(expire_timedelta.total_seconds())

            self._cookies[key] = (item, store_time)

    def _clear_expired_cookies(self) -> None:
        """Clear expired cookies."""
        check_time = datetime.now(tz=timezone.utc)
        expired_keys = []

        for key, (morsel, store_time) in self._cookies.items():
            cookie_max_age = morsel.get("max-age")
            if not cookie_max_age:
                continue
            expired_timedelta = check_time - store_time
            if expired_timedelta.total_seconds() > float(cookie_max_age):
                expired_keys.append(key)

        for key in expired_keys:
            self._cookies.pop(key)

    def _update_cookie_header(self, connection_args: dict[str, ty.Any]) -> None:
        """Update a cookie header."""
        self._clear_expired_cookies()

        gateway_cookie_values = "; ".join(
            f"{name}={morsel.coded_value}" for name, (morsel, _time) in self._cookies.items()
        )
        if gateway_cookie_values:
            headers = connection_args.get("headers", {})

            # As headers are case-insensitive, we get existing name of cookie header,
            #  or use "Cookie" by default.
            cookie_header_name = next(
                (header_key for header_key in headers if header_key.lower() == "cookie"),
                "Cookie",
            )
            existing_cookie = headers.get(cookie_header_name)

            # merge gateway-managed cookies with cookies already in arguments
            if existing_cookie:
                gateway_cookie_values = existing_cookie + "; " + gateway_cookie_values
            headers[cookie_header_name] = gateway_cookie_values

            connection_args["headers"] = headers


class RetryableHTTPClient:
    """
    Inspired by urllib.util.Retry (https://urllib3.readthedocs.io/en/stable/reference/urllib3.util.html),
    this class is initialized with desired retry characteristics, uses a recursive method `fetch()` against an instance
    of `AsyncHTTPClient` which tracks the current retry count across applicable request retries.
    """

    MAX_RETRIES_DEFAULT = 2
    MAX_RETRIES_CAP = 10  # The upper limit to max_retries value.
    max_retries: int = int(os.getenv("JUPYTER_GATEWAY_MAX_REQUEST_RETRIES", MAX_RETRIES_DEFAULT))
    max_retries = max(0, min(max_retries, MAX_RETRIES_CAP))  # Enforce boundaries
    retried_methods: set[str] = {"GET", "DELETE"}
    retried_errors: set[int] = {502, 503, 504, 599}
    retried_exceptions: set[type] = {ConnectionError}
    backoff_factor: float = 0.1

    def __init__(self):
        """Initialize the retryable http client."""
        self.retry_count: int = 0
        self.client: AsyncHTTPClient = AsyncHTTPClient()

    async def fetch(self, endpoint: str, **kwargs: ty.Any) -> HTTPResponse:
        """
        Retryable AsyncHTTPClient.fetch() method.  When the request fails, this method will
        recurse up to max_retries times if the condition deserves a retry.
        """
        self.retry_count = 0
        return await self._fetch(endpoint, **kwargs)

    async def _fetch(self, endpoint: str, **kwargs: ty.Any) -> HTTPResponse:
        """
        Performs the fetch against the contained AsyncHTTPClient instance and determines
        if retry is necessary on any exceptions.  If so, retry is performed recursively.
        """
        try:
            response: HTTPResponse = await self.client.fetch(endpoint, **kwargs)
        except Exception as e:
            is_retryable: bool = await self._is_retryable(kwargs["method"], e)
            if not is_retryable:
                raise e
            logging.getLogger("ServerApp").info(
                f"Attempting retry ({self.retry_count}) against "
                f"endpoint '{endpoint}'.  Retried error: '{e!r}'"
            )
            response = await self._fetch(endpoint, **kwargs)
        return response

    async def _is_retryable(self, method: str, exception: Exception) -> bool:
        """Determines if the given exception is retryable based on object's configuration."""

        if method not in self.retried_methods:
            return False
        if self.retry_count == self.max_retries:
            return False

        # Determine if error is retryable...
        if isinstance(exception, HTTPClientError):
            hce: HTTPClientError = exception
            if hce.code not in self.retried_errors:
                return False
        elif not any(isinstance(exception, error) for error in self.retried_exceptions):
            return False

        # Is retryable, wait for backoff, then increment count
        await asyncio.sleep(self.backoff_factor * (2**self.retry_count))
        self.retry_count += 1
        return True


async def gateway_request(endpoint: str, **kwargs: ty.Any) -> HTTPResponse:
    """Make an async request to kernel gateway endpoint, returns a response"""
    gateway_client = GatewayClient.instance()
    kwargs = gateway_client.load_connection_args(**kwargs)
    rhc = RetryableHTTPClient()
    try:
        response = await rhc.fetch(endpoint, **kwargs)
        gateway_client.emit(
            data={STATUS_KEY: SUCCESS_STATUS, STATUS_CODE_KEY: 200, MESSAGE_KEY: "success"}
        )
    # Trap a set of common exceptions so that we can inform the user that their Gateway url is incorrect
    # or the server is not running.
    # NOTE: We do this here since this handler is called during the server's startup and subsequent refreshes
    # of the tree view.
    except HTTPClientError as e:
        gateway_client.emit(
            data={STATUS_KEY: ERROR_STATUS, STATUS_CODE_KEY: e.code, MESSAGE_KEY: str(e.message)}
        )
        error_reason = (
            f"Exception while attempting to connect to Gateway server url '{gateway_client.url}'"
        )
        error_message = e.message
        if e.response:
            try:
                error_payload = json.loads(e.response.body)
                error_reason = error_payload.get("reason") or error_reason
                error_message = error_payload.get("message") or error_message
            except json.decoder.JSONDecodeError:
                error_reason = e.response.body.decode()

        raise web.HTTPError(
            e.code,
            f"Error from Gateway: [{error_message}] {error_reason}. "
            "Ensure gateway url is valid and the Gateway instance is running.",
        ) from e
    except ConnectionError as e:
        gateway_client.emit(
            data={STATUS_KEY: ERROR_STATUS, STATUS_CODE_KEY: 503, MESSAGE_KEY: str(e)}
        )
        raise web.HTTPError(
            503,
            f"ConnectionError was received from Gateway server url '{gateway_client.url}'.  "
            "Check to be sure the Gateway instance is running.",
        ) from e
    except gaierror as e:
        gateway_client.emit(
            data={STATUS_KEY: ERROR_STATUS, STATUS_CODE_KEY: 404, MESSAGE_KEY: str(e)}
        )
        raise web.HTTPError(
            404,
            f"The Gateway server specified in the gateway_url '{gateway_client.url}' doesn't "
            f"appear to be valid.  Ensure gateway url is valid and the Gateway instance is running.",
        ) from e
    except Exception as e:
        gateway_client.emit(
            data={STATUS_KEY: ERROR_STATUS, STATUS_CODE_KEY: 505, MESSAGE_KEY: str(e)}
        )
        logging.getLogger("ServerApp").error(
            "Exception while trying to launch kernel via Gateway URL %s: %s",
            gateway_client.url,
            e,
        )
        raise e

    if gateway_client.accept_cookies:
        # Update cookies on GatewayClient from server if configured.
        cookie_values = response.headers.get("Set-Cookie")
        if cookie_values:
            cookie: SimpleCookie = SimpleCookie()
            cookie.load(cookie_values)
            gateway_client.update_cookies(cookie)
    return response
