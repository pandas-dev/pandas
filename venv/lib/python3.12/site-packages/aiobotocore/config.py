import collections.abc
import copy
import ssl
import sys
from typing import Optional, TypedDict, Union

import botocore.client
from aiohttp import SocketFactoryType
from aiohttp.abc import AbstractResolver
from botocore.exceptions import ParamValidationError

from ._constants import DEFAULT_KEEPALIVE_TIMEOUT
from .endpoint import DEFAULT_HTTP_SESSION_CLS
from .httpsession import AIOHTTPSession
from .httpxsession import HttpxSession

if sys.version_info >= (3, 11):
    from typing import NotRequired
else:
    from typing_extensions import NotRequired

TIMEOUT_ARGS = frozenset(
    ('keepalive_timeout', 'write_timeout', 'pool_timeout')
)


class _ConnectorArgs(TypedDict):
    use_dns_cache: NotRequired[bool]
    ttl_dns_cache: NotRequired[Optional[int]]
    keepalive_timeout: NotRequired[Optional[float]]
    write_timeout: NotRequired[Optional[float]]
    pool_timeout: NotRequired[Optional[float]]
    force_close: NotRequired[bool]
    ssl_context: NotRequired[ssl.SSLContext]
    resolver: NotRequired[AbstractResolver]
    socket_factory: NotRequired[Optional[SocketFactoryType]]


_HttpSessionType = Union[AIOHTTPSession, HttpxSession]


class AioConfig(botocore.client.Config):
    def __init__(
        self,
        connector_args: Optional[_ConnectorArgs] = None,
        http_session_cls: type[_HttpSessionType] = DEFAULT_HTTP_SESSION_CLS,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self.connector_args: _ConnectorArgs = (
            copy.copy(connector_args) if connector_args else {}
        )
        self.http_session_cls: type[_HttpSessionType] = http_session_cls
        self._validate_connector_args(
            self.connector_args, self.http_session_cls
        )

        if 'keepalive_timeout' not in self.connector_args:
            self.connector_args['keepalive_timeout'] = (
                DEFAULT_KEEPALIVE_TIMEOUT
            )

    def merge(self, other_config):
        # Adapted from parent class
        config_options = copy.copy(self._user_provided_options)
        config_options.update(other_config._user_provided_options)
        return AioConfig(self.connector_args, **config_options)

    @staticmethod
    def _validate_connector_args(
        connector_args: _ConnectorArgs,
        http_session_cls: type[_HttpSessionType],
    ) -> None:
        for k, v in connector_args.items():
            # verify_ssl is handled by verify parameter to create_client
            if k == 'use_dns_cache':
                if http_session_cls is HttpxSession:
                    raise ParamValidationError(
                        report='Httpx does not support dns caching. https://github.com/encode/httpx/discussions/2211'
                    )
                if not isinstance(v, bool):
                    raise ParamValidationError(
                        report=f'{k} value must be a boolean'
                    )
            elif k == 'ttl_dns_cache':
                if v is not None and not isinstance(v, int):
                    raise ParamValidationError(
                        report=f'{k} value must be an int or None'
                    )
            elif k in TIMEOUT_ARGS:
                if v is not None and not isinstance(v, (float, int)):
                    raise ParamValidationError(
                        report=f'{k} value must be a float/int or None'
                    )
            elif k == 'force_close':
                if http_session_cls is HttpxSession:
                    raise ParamValidationError(
                        report=f'Httpx backend does not currently support {k}.'
                    )
                if not isinstance(v, bool):
                    raise ParamValidationError(
                        report=f'{k} value must be a boolean'
                    )
            # limit is handled by max_pool_connections
            elif k == 'ssl_context':
                if not isinstance(v, ssl.SSLContext):
                    raise ParamValidationError(
                        report=f'{k} must be an SSLContext instance'
                    )
            elif k == "resolver":
                if http_session_cls is HttpxSession:
                    raise ParamValidationError(
                        report=f'Httpx backend does not support {k}.'
                    )
                if not isinstance(v, AbstractResolver):
                    raise ParamValidationError(
                        report=f'{k} must be an instance of a AbstractResolver'
                    )
            elif k == "socket_factory":
                if http_session_cls is HttpxSession:
                    raise ParamValidationError(
                        report=f'Httpx backend does not support {k}.'
                    )
                if v is not None and not isinstance(
                    v, collections.abc.Callable
                ):
                    raise ParamValidationError(
                        report=f'{k} must be a callable'
                    )
            else:
                raise ParamValidationError(report=f'invalid connector_arg:{k}')
