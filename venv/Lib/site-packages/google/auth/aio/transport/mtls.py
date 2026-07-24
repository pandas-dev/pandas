# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Helper functions for mTLS in async for discovery of certs.
"""

import asyncio
import inspect
import logging
import ssl
from typing import Optional

from google.auth import exceptions
from google.auth.transport._mtls_helper import secure_cert_key_paths
import google.auth.transport.mtls

_LOGGER = logging.getLogger(__name__)


def make_client_cert_ssl_context(
    cert_bytes: bytes, key_bytes: bytes, passphrase: Optional[bytes] = None
) -> ssl.SSLContext:
    """Creates an SSLContext with the given client certificate and key.
    This function writes the certificate and key to temporary files so that
    ssl.create_default_context can load them, as the ssl module requires
    file paths for client certificates. These temporary files are deleted
    immediately after the SSL context is created.
    Args:
        cert_bytes (bytes): The client certificate content in PEM format.
        key_bytes (bytes): The client private key content in PEM format.
        passphrase (Optional[bytes]): The passphrase for the private key, if any.
    Returns:
        ssl.SSLContext: The configured SSL context with client certificate.

    Raises:
        google.auth.exceptions.TransportError: If there is an error loading the certificate.
    """
    try:
        with secure_cert_key_paths(cert_bytes, key_bytes, passphrase=passphrase) as (
            cert_path,
            key_path,
            passphrase_val,
        ):
            context = ssl.create_default_context(ssl.Purpose.SERVER_AUTH)
            if cert_path:
                password = passphrase_val
                context.load_cert_chain(
                    certfile=cert_path,
                    keyfile=key_path,
                    password=password,
                )
            return context
    except (ssl.SSLError, OSError, IOError, ValueError, RuntimeError, TypeError) as exc:
        raise exceptions.TransportError(
            "Failed to load client certificate and key for mTLS."
        ) from exc


async def _run_in_executor(func, *args):
    """Run a blocking function in an executor to avoid blocking the event loop.

    This implements the non-blocking execution strategy for disk I/O operations.
    """
    try:
        # For python versions 3.9 and newer versions
        return await asyncio.to_thread(func, *args)
    except AttributeError:
        # Fallback for older Python versions
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, func, *args)


def default_client_cert_source():
    """Get a callback which returns the default client SSL credentials.

    Returns:
        Awaitable[Callable[[], Tuple[bytes, bytes]]]: A callback which returns the default
            client certificate bytes and private key bytes, both in PEM format.

    Raises:
        google.auth.exceptions.DefaultClientCertSourceError: If the default
            client SSL credentials don't exist or are malformed.
    """
    if not google.auth.transport.mtls.has_default_client_cert_source(
        include_context_aware=False
    ):
        raise exceptions.MutualTLSChannelError(
            "Default client cert source doesn't exist"
        )

    async def callback():
        try:
            _, cert_bytes, key_bytes = await get_client_cert_and_key()
        except (OSError, RuntimeError, ValueError) as caught_exc:
            new_exc = exceptions.MutualTLSChannelError(caught_exc)
            raise new_exc from caught_exc

        return cert_bytes, key_bytes

    return callback


async def get_client_ssl_credentials(
    certificate_config_path=None,
):
    """Returns the client side certificate, private key and passphrase.

    We look for certificates and keys with the following order of priority:
        1. Certificate and key specified by certificate_config.json.
               Currently, only X.509 workload certificates are supported.

    Args:
        certificate_config_path (str): The certificate_config.json file path.

    Returns:
        Tuple[bool, bytes, bytes, bytes]:
            A boolean indicating if cert, key and passphrase are obtained, the
            cert bytes and key bytes both in PEM format, and passphrase bytes.

    Raises:
        google.auth.exceptions.ClientCertError: if problems occurs when getting
            the cert, key and passphrase.
    """

    # Attempt to retrieve X.509 Workload cert and key.
    cert, key = await _run_in_executor(
        google.auth.transport._mtls_helper._get_workload_cert_and_key,
        certificate_config_path,
        False,
    )

    if cert and key:
        return True, cert, key, None

    return False, None, None, None


async def get_client_cert_and_key(client_cert_callback=None):
    """Returns the client side certificate and private key. The function first
    tries to get certificate and key from client_cert_callback; if the callback
    is None or doesn't provide certificate and key, the function tries application
    default SSL credentials.

    Args:
        client_cert_callback (Optional[Callable[[], (bytes, bytes)]]): An
            optional callback which returns client certificate bytes and private
            key bytes both in PEM format.

    Returns:
        Tuple[bool, bytes, bytes]:
            A boolean indicating if cert and key are obtained, the cert bytes
            and key bytes both in PEM format.

    Raises:
        google.auth.exceptions.ClientCertError: if problems occurs when getting
            the cert and key.
    """
    if client_cert_callback:
        result = client_cert_callback()
        if inspect.isawaitable(result):
            cert, key = await result
        else:
            cert, key = result
        return True, cert, key

    has_cert, cert, key, _ = await get_client_ssl_credentials()
    return has_cert, cert, key
