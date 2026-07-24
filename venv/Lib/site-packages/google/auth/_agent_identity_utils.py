# Copyright 2025 Google LLC
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

"""Helpers for Agent Identity credentials."""

import base64
import hashlib
import os
import re
import stat
import time
from urllib.parse import quote, urlparse
import warnings

from google.auth import environment_vars, exceptions

CRYPTOGRAPHY_NOT_FOUND_ERROR = (
    "The cryptography library is required for certificate-based authentication."
    "Please install it with `pip install google-auth[cryptography]`."
)

# SPIFFE trust domain patterns for Agent Identities.
_AGENT_IDENTITY_SPIFFE_TRUST_DOMAIN_PATTERNS = [
    r"^agents\.global\.org-\d+\.system\.id\.goog$",
    r"^agents\.global\.proj-\d+\.system\.id\.goog$",
    r"^agents-nonprod\.global\.org-\d+\.system\.id\.goog$",
    r"^agents-nonprod\.global\.proj-\d+\.system\.id\.goog$",
]

_WELL_KNOWN_CERT_PATH = "/var/run/secrets/workload-spiffe-credentials/certificates.pem"

# Constants for polling the certificate file.
_FAST_POLL_CYCLES = 50
_FAST_POLL_INTERVAL = 0.1  # 100ms
_SLOW_POLL_INTERVAL = 0.5  # 500ms
_TOTAL_TIMEOUT = 30  # seconds

# Calculate the number of slow poll cycles based on the total timeout.
_SLOW_POLL_CYCLES = int(
    (_TOTAL_TIMEOUT - (_FAST_POLL_CYCLES * _FAST_POLL_INTERVAL)) / _SLOW_POLL_INTERVAL
)

_POLLING_INTERVALS = ([_FAST_POLL_INTERVAL] * _FAST_POLL_CYCLES) + (
    [_SLOW_POLL_INTERVAL] * _SLOW_POLL_CYCLES
)


def _is_certificate_file_ready(path):
    """Checks if a file exists, is a regular file, and is not empty."""
    if not path:
        return False
    try:
        # Check if the path points to a regular file and is not empty.
        # stat.S_ISREG is used instead of os.path.isfile to avoid swallowing
        # PermissionError exceptions, which the caller needs to propagate.
        st = os.stat(path)
        return stat.S_ISREG(st.st_mode) and st.st_size > 0
    except PermissionError:
        # Propagate PermissionError to let caller handle it (e.g., return early or fallback)
        raise
    except OSError:
        return False


def get_agent_identity_certificate_path():
    """Gets the agent certificate path from the certificate config file.

    The path to the certificate config file is read from the
    GOOGLE_API_CERTIFICATE_CONFIG environment variable. This function
    can optionally trigger polling to handle cases where the environment
    variable is set before the files are available on the filesystem.

    Returns:
        Optional[str]: The path to the agent's certificate file, or None if unavailable.

    Raises:
        google.auth.exceptions.RefreshError: If the certificate config file
            or the certificate file cannot be found after retries.
    """
    cert_config_path = os.environ.get(environment_vars.GOOGLE_API_CERTIFICATE_CONFIG)

    if not cert_config_path:
        return None

    # We trigger polling only if the config path points to the well-known directory.
    # Cloud Run dynamically generates these files in this directory, and both the
    # config file and the certificate file may experience a brief startup latency.
    # For all other paths, we return early to avoid introducing unnecessary startup
    # delays.
    well_known_dir = os.path.dirname(_WELL_KNOWN_CERT_PATH)
    try:
        abs_cert_path = os.path.abspath(cert_config_path)
        abs_well_known_dir = os.path.abspath(well_known_dir)
        should_poll = (
            os.path.commonpath([abs_well_known_dir, abs_cert_path])
            == abs_well_known_dir
        )
    except ValueError:
        should_poll = False

    return _get_cert_path_with_optional_polling(cert_config_path, should_poll)


def _get_cert_path_with_optional_polling(cert_config_path, should_poll):
    """Gets the certificate path, optionally polling until it is ready.

    Args:
        cert_config_path (str): The path to the certificate configuration file.
        should_poll (bool): If True, the function will poll for the file and
            certificate to be ready. If False, it will check only once and
            return early if they are not immediately available.

    Returns:
        str: The path to the certificate file, or None if unavailable.

    Raises:
        google.auth.exceptions.RefreshError: If the certificate config file
            or the certificate file cannot be found after retries.
    """
    has_logged_config_warning = False
    has_logged_cert_warning = False

    for interval in _POLLING_INTERVALS:
        try:
            cert_path = _parse_cert_path_from_config(cert_config_path)

            if cert_path is None:
                return None

            if _is_certificate_file_ready(cert_path):
                return cert_path

            # The config was parsed, but the cert file is not ready yet
            if not should_poll:
                # If polling is disabled, return early.
                return None

            if not has_logged_cert_warning:
                warnings.warn(
                    f"Certificate file not ready at {cert_path}. Retrying until startup timeout (up to {_TOTAL_TIMEOUT} seconds total)..."
                )
                has_logged_cert_warning = True

        except PermissionError as e:
            warnings.warn(
                f"Permission denied when accessing certificate config or certificate file: {e}. "
                "Token binding protection cannot be enabled. Falling back to unbound tokens."
            )
            return None
        except (IOError, ValueError, KeyError) as e:
            if os.path.exists(cert_config_path):
                # If the file exists but has invalid JSON or is unreadable,
                # we assume it is in its final format and return early (returning None).
                return None

            if not should_poll:
                # If polling is disabled, return early if the file doesn't exist.
                return None

            if not has_logged_config_warning:
                warnings.warn(
                    f"Certificate config file not found or incomplete: {e} (from "
                    f"{environment_vars.GOOGLE_API_CERTIFICATE_CONFIG} environment variable). "
                    f"Retrying until startup timeout (up to {_TOTAL_TIMEOUT} seconds total)..."
                )
                has_logged_config_warning = True

        # Sleep before the next polling attempt.
        time.sleep(interval)

    raise exceptions.RefreshError(
        "Certificate config or certificate file not found after multiple retries. "
        f"Token binding protection is failing. You can turn off this protection by setting "
        f"{environment_vars.GOOGLE_API_PREVENT_AGENT_TOKEN_SHARING_FOR_GCP_SERVICES} to false "
        "to fall back to unbound tokens."
    )


def _parse_cert_path_from_config(cert_config_path):
    """Reads the cert config file and returns the cert_path.

    Args:
        cert_config_path (str): The path to the certificate configuration file.

    Returns:
        Optional[str]: The path to the certificate file, or None if not found
            in the config.

    Raises:
        IOError: If the certificate config file cannot be read.
        ValueError: If the certificate config file contains invalid JSON.
        KeyError: If the certificate config file does not contain the
            expected structure.
    """
    import json

    with open(cert_config_path, "r", encoding="utf-8") as f:
        cert_config = json.load(f)

    cert_configs = (
        cert_config.get("cert_configs") if isinstance(cert_config, dict) else None
    )
    workload_config = (
        cert_configs.get("workload") if isinstance(cert_configs, dict) else None
    )

    if not isinstance(workload_config, dict) or "cert_path" not in workload_config:
        return None

    return workload_config["cert_path"]


def get_and_parse_agent_identity_certificate():
    """Gets and parses the agent identity certificate if not opted out.

    Checks if the user has opted out of certificate-bound tokens. If not,
    it gets the certificate path, reads the file, and parses it.

    Returns:
        The parsed certificate object if found and not opted out, otherwise None.
    """
    # If the user has opted out of cert bound tokens, there is no need to
    # look up the certificate.
    is_opted_out = (
        os.environ.get(
            environment_vars.GOOGLE_API_PREVENT_AGENT_TOKEN_SHARING_FOR_GCP_SERVICES,
            "true",
        ).lower()
        == "false"
    )
    if is_opted_out:
        return None

    # Respect explicit opt-out of mTLS / client certs
    from google.auth.transport import _mtls_helper

    env_override = _mtls_helper._check_use_client_cert_env()
    if env_override is False:
        return None

    cert_path = get_agent_identity_certificate_path()
    if not cert_path:
        return None

    try:
        with open(cert_path, "rb") as cert_file:
            cert_bytes = cert_file.read()
    except PermissionError as e:
        warnings.warn(
            f"Failed to read agent identity certificate file at {cert_path}: {e}. "
            "Token binding protection cannot be enabled. Falling back to unbound tokens."
        )
        return None

    return parse_certificate(cert_bytes)


def parse_certificate(cert_bytes):
    """Parses a PEM-encoded certificate.

    Args:
        cert_bytes (bytes): The PEM-encoded certificate bytes.

    Returns:
        cryptography.x509.Certificate: The parsed certificate object.
    """
    try:
        from cryptography import x509

        return x509.load_pem_x509_certificate(cert_bytes)
    except ImportError as e:
        raise ImportError(CRYPTOGRAPHY_NOT_FOUND_ERROR) from e


def _is_agent_identity_certificate(cert):
    """Checks if a certificate is an Agent Identity certificate.

    This is determined by checking the Subject Alternative Name (SAN) for a
    SPIFFE ID with a trust domain matching Agent Identity patterns.

    Args:
        cert (cryptography.x509.Certificate): The parsed certificate object.

    Returns:
        bool: True if the certificate is an Agent Identity certificate,
            False otherwise.
    """
    try:
        from cryptography import x509
        from cryptography.x509.oid import ExtensionOID

        try:
            ext = cert.extensions.get_extension_for_oid(
                ExtensionOID.SUBJECT_ALTERNATIVE_NAME
            )
        except x509.ExtensionNotFound:
            return False
        uris = ext.value.get_values_for_type(x509.UniformResourceIdentifier)

        for uri in uris:
            parsed_uri = urlparse(uri)
            if parsed_uri.scheme == "spiffe":
                trust_domain = parsed_uri.netloc
                for pattern in _AGENT_IDENTITY_SPIFFE_TRUST_DOMAIN_PATTERNS:
                    if re.match(pattern, trust_domain):
                        return True
        return False
    except ImportError as e:
        raise ImportError(CRYPTOGRAPHY_NOT_FOUND_ERROR) from e


def calculate_certificate_fingerprint(cert):
    """Calculates the URL-encoded, unpadded, base64-encoded SHA256 hash of a
    DER-encoded certificate.

    Args:
        cert (cryptography.x509.Certificate): The parsed certificate object.

    Returns:
        str: The URL-encoded, unpadded, base64-encoded SHA256 fingerprint.
    """
    try:
        from cryptography.hazmat.primitives import serialization

        der_cert = cert.public_bytes(serialization.Encoding.DER)
        fingerprint = hashlib.sha256(der_cert).digest()
        # The certificate fingerprint is generated in two steps to align with GFE's
        # expectations and ensure proper URL transmission:
        # 1. Standard base64 encoding is applied, and padding ('=') is removed.
        # 2. The resulting string is then URL-encoded to handle special characters
        #    ('+', '/') that would otherwise be misinterpreted in URL parameters.
        base64_fingerprint = base64.b64encode(fingerprint).decode("utf-8")
        unpadded_base64_fingerprint = base64_fingerprint.rstrip("=")
        return quote(unpadded_base64_fingerprint)
    except ImportError as e:
        raise ImportError(CRYPTOGRAPHY_NOT_FOUND_ERROR) from e


def should_request_bound_token(cert):
    """Determines if a bound token should be requested.

    This is based on the GOOGLE_API_PREVENT_AGENT_TOKEN_SHARING_FOR_GCP_SERVICES
    environment variable and whether the certificate is an agent identity cert.

    Args:
        cert (cryptography.x509.Certificate): The parsed certificate object.

    Returns:
        bool: True if a bound token should be requested, False otherwise.
    """
    is_agent_cert = _is_agent_identity_certificate(cert)
    is_opted_in = (
        os.environ.get(
            environment_vars.GOOGLE_API_PREVENT_AGENT_TOKEN_SHARING_FOR_GCP_SERVICES,
            "true",
        ).lower()
        == "true"
    )
    if not (is_agent_cert and is_opted_in):
        return False

    # Respect explicit opt-out of mTLS / client certs
    from google.auth.transport import _mtls_helper

    env_override = _mtls_helper._check_use_client_cert_env()
    if env_override is False:
        return False

    return True


def get_cached_cert_fingerprint(cached_cert):
    """Returns the fingerprint of the cached certificate."""
    if cached_cert:
        cert_obj = parse_certificate(cached_cert)
        cached_cert_fingerprint = calculate_certificate_fingerprint(cert_obj)
    else:
        raise ValueError("mTLS connection is not configured.")
    return cached_cert_fingerprint
