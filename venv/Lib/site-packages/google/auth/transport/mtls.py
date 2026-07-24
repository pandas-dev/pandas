# Copyright 2020 Google LLC
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

"""Utilites for mutual TLS."""

import enum
import logging
from os import getenv
import ssl
from typing import Optional

from google.auth import environment_vars
from google.auth import exceptions
from google.auth.transport import _mtls_helper

_LOGGER = logging.getLogger(__name__)


class UseMtlsEndpointMode(enum.Enum):
    ALWAYS = "always"
    NEVER = "never"
    AUTO = "auto"


def has_default_client_cert_source(include_context_aware=True):
    """Check if default client SSL credentials exists on the device.

    Args:
       include_context_aware (bool): include_context_aware indicates if context_aware
       path location will be checked or should it be skipped.

    Returns:
        bool: indicating if the default client cert source exists.
    """
    cert_path = _mtls_helper._get_cert_config_path(
        include_context_aware=include_context_aware
    )
    if cert_path is not None:
        return True
    if (
        include_context_aware
        and _mtls_helper._check_config_path(_mtls_helper.CONTEXT_AWARE_METADATA_PATH)
        is not None
    ):
        return True

    return False


def default_client_cert_source():
    """Get a callback which returns the default client SSL credentials.

    Returns:
        Callable[[], [bytes, bytes]]: A callback which returns the default
            client certificate bytes and private key bytes, both in PEM format.

    Raises:
        google.auth.exceptions.MutualTLSChannelError: If the default
            client SSL credentials don't exist or are malformed.
    """
    if not has_default_client_cert_source(include_context_aware=True):
        raise exceptions.MutualTLSChannelError(
            "Default client cert source doesn't exist"
        )

    def callback():
        try:
            _, cert_bytes, key_bytes = _mtls_helper.get_client_cert_and_key()
        except (OSError, RuntimeError, ValueError) as caught_exc:
            new_exc = exceptions.MutualTLSChannelError(caught_exc)
            raise new_exc from caught_exc

        return cert_bytes, key_bytes

    return callback


def default_client_encrypted_cert_source(cert_path, key_path):
    """Get a callback which returns the default encrpyted client SSL credentials.

    Args:
        cert_path (str): The cert file path. The default client certificate will
            be written to this file when the returned callback is called.
        key_path (str): The key file path. The default encrypted client key will
            be written to this file when the returned callback is called.

    Returns:
        Callable[[], [str, str, bytes]]: A callback which generates the default
            client certificate, encrpyted private key and passphrase. It writes
            the certificate and private key into the cert_path and key_path, and
            returns the cert_path, key_path and passphrase bytes.

    Raises:
        google.auth.exceptions.MutualTLSChannelError: If any problem
            occurs when loading or saving the client certificate and key.
    """
    if not has_default_client_cert_source(include_context_aware=True):
        raise exceptions.MutualTLSChannelError(
            "Default client encrypted cert source doesn't exist"
        )

    def callback():
        try:
            (
                _,
                cert_bytes,
                key_bytes,
                passphrase_bytes,
            ) = _mtls_helper.get_client_ssl_credentials(generate_encrypted_key=True)
            with open(cert_path, "wb") as cert_file:
                cert_file.write(cert_bytes)
            with open(key_path, "wb") as key_file:
                key_file.write(key_bytes)
        except (exceptions.ClientCertError, OSError) as caught_exc:
            new_exc = exceptions.MutualTLSChannelError(caught_exc)
            raise new_exc from caught_exc

        return cert_path, key_path, passphrase_bytes

    return callback


def should_use_client_cert():
    """Returns boolean for whether the client certificate should be used for mTLS.

    This is a wrapper around _mtls_helper.check_use_client_cert().
    If GOOGLE_API_USE_CLIENT_CERTIFICATE is set to true or false, a corresponding
    bool value will be returned
    If GOOGLE_API_USE_CLIENT_CERTIFICATE is unset, the value will be inferred by
    reading a file pointed at by GOOGLE_API_CERTIFICATE_CONFIG or
    CLOUDSDK_CONTEXT_AWARE_CERTIFICATE_CONFIG_FILE_PATH, or the default path
    like ~/.config/gcloud/certificate_config.json, and verifying it
    contains a "workload" section. If so, the function will return True,
    otherwise False.

    Returns:
       bool: indicating whether the client certificate should be used for mTLS.
    """
    return _mtls_helper.check_use_client_cert()


def _load_client_cert_into_context(
    ctx: ssl.SSLContext,
    cert_bytes: bytes,
    key_bytes: bytes,
    passphrase: Optional[bytes] = None,
) -> None:
    """Load a client certificate and key into an SSL context.

    Args:
        ctx (ssl.SSLContext): The SSL context to load the certificate and key into.
        cert_bytes (bytes): The client certificate bytes in PEM format.
        key_bytes (bytes): The client private key bytes in PEM format.
        passphrase (Optional[bytes]): The passphrase for the client private key.

    Raises:
        google.auth.exceptions.MutualTLSChannelError: If the SSL context is invalid,
            or if loading the certificate and key fails.
    """
    if not isinstance(ctx, ssl.SSLContext):
        raise exceptions.MutualTLSChannelError(
            "Failed to load client certificate and key for mTLS. The provided context "
            "object is invalid or does not support loading certificate chains."
        )

    try:
        with _mtls_helper.secure_cert_key_paths(
            cert_bytes, key_bytes, passphrase=passphrase
        ) as (
            cert_path,
            key_path,
            passphrase_val,
        ):
            if cert_path is None or key_path is None:
                raise exceptions.MutualTLSChannelError(
                    "Failed to generate temporary file paths for the client certificate and key."
                )
            ctx.load_cert_chain(
                certfile=cert_path, keyfile=key_path, password=passphrase_val
            )
    except (
        ssl.SSLError,
        OSError,
        ValueError,
        RuntimeError,
        TypeError,
    ) as caught_exc:
        new_exc = exceptions.MutualTLSChannelError(caught_exc)
        raise new_exc from caught_exc


def load_default_client_cert(ctx: ssl.SSLContext) -> bool:
    """Load the default client certificate and key into an SSL context if configured.

    If client certificates are enabled and a default client certificate source is
    found, the certificate and key are loaded into the SSL context.

    Args:
        ctx (ssl.SSLContext): The SSL context to load the default client certificate
            and key into.

    Returns:
        bool: True if client certificates are enabled and the default client
            certificate was successfully loaded. False if client certificates
            are disabled or if no default certificate source is configured.

    Raises:
        google.auth.exceptions.MutualTLSChannelError: If the default client certificate
            or key is malformed.
    """
    if not should_use_client_cert() or not has_default_client_cert_source():
        return False
    try:
        (
            has_cert,
            cert_bytes,
            key_bytes,
            passphrase,
        ) = _mtls_helper.get_client_ssl_credentials()
    except (
        exceptions.ClientCertError,
        OSError,
        RuntimeError,
        ValueError,
    ) as caught_exc:
        new_exc = exceptions.MutualTLSChannelError(caught_exc)
        raise new_exc from caught_exc
    else:
        if not has_cert:
            return False
        _load_client_cert_into_context(ctx, cert_bytes, key_bytes, passphrase)
        return True


def get_default_ssl_context() -> Optional[ssl.SSLContext]:
    """Get a default SSL context loaded with the default client certificate.

    Returns:
        ssl.SSLContext: An SSL context loaded with the default client
            certificate, or None if client certificates are not configured
            or available.

    Raises:
        google.auth.exceptions.MutualTLSChannelError: If the default client certificate
            or key is malformed.
    """
    if not should_use_client_cert() or not has_default_client_cert_source():
        return None

    ctx = ssl.create_default_context(ssl.Purpose.SERVER_AUTH)
    return ctx if load_default_client_cert(ctx) else None


def should_use_mtls_endpoint(
    client_cert_available: Optional[bool] = None,
) -> bool:
    """Determine whether to use an mTLS endpoint.

    This relies on the GOOGLE_API_USE_MTLS_ENDPOINT environment variable. If set to
    "always", returns True. If set to "never", returns False. If set to "auto"
    or unset, returns whether a client certificate is available.

    Args:
        client_cert_available (Optional[bool]): indicating if a client certificate
            is available. If None, this is determined by checking if client
            certificates are enabled using :func:`should_use_client_cert`.

    Returns:
        bool: indicating if an mTLS endpoint should be used.
    """
    if client_cert_available is None:
        client_cert_available = should_use_client_cert()

    use_mtls_endpoint = getenv(environment_vars.GOOGLE_API_USE_MTLS_ENDPOINT)
    use_mtls_endpoint = (use_mtls_endpoint or "auto").strip().lower()
    try:
        mode = UseMtlsEndpointMode(use_mtls_endpoint)
    except ValueError:
        raise exceptions.MutualTLSChannelError(
            f"Unsupported {environment_vars.GOOGLE_API_USE_MTLS_ENDPOINT} value "
            f"'{use_mtls_endpoint}'. Accepted values: never, auto, always."
        )

    if mode == UseMtlsEndpointMode.ALWAYS:
        return True
    if mode == UseMtlsEndpointMode.NEVER:
        return False
    if mode == UseMtlsEndpointMode.AUTO:
        return client_cert_available
