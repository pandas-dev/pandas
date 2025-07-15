import copy

import botocore.client
from botocore.exceptions import ParamValidationError

from aiobotocore.endpoint import DEFAULT_HTTP_SESSION_CLS
from aiobotocore.httpxsession import HttpxSession

# AWS has a 20 second idle timeout:
#   https://web.archive.org/web/20150926192339/https://forums.aws.amazon.com/message.jspa?messageID=215367
# and aiohttp default timeout is 30s so we set it to something
# reasonable here
DEFAULT_KEEPALIVE_TIMEOUT = 12

TIMEOUT_ARGS = frozenset(
    ('keepalive_timeout', 'write_timeout', 'pool_timeout')
)


class AioConfig(botocore.client.Config):
    def __init__(
        self,
        connector_args=None,
        http_session_cls=DEFAULT_HTTP_SESSION_CLS,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self._validate_connector_args(connector_args, http_session_cls)
        self.connector_args = copy.copy(connector_args)
        self.http_session_cls = http_session_cls
        if not self.connector_args:
            self.connector_args = dict()

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
    def _validate_connector_args(connector_args, http_session_cls):
        if connector_args is None:
            return

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
                import ssl

                if not isinstance(v, ssl.SSLContext):
                    raise ParamValidationError(
                        report=f'{k} must be an SSLContext instance'
                    )
            elif k == "resolver":
                from aiohttp.abc import AbstractResolver

                if http_session_cls is HttpxSession:
                    raise ParamValidationError(
                        report=f'Httpx backend does not support {k}.'
                    )
                if not isinstance(v, AbstractResolver):
                    raise ParamValidationError(
                        report=f'{k} must be an instance of a AbstractResolver'
                    )
            else:
                raise ParamValidationError(report=f'invalid connector_arg:{k}')
