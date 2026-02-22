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

"""Helper functions for getting mTLS cert and key."""

import json
import logging
from os import environ, getenv, path
import re
import subprocess

from google.auth import _agent_identity_utils
from google.auth import environment_vars
from google.auth import exceptions

CONTEXT_AWARE_METADATA_PATH = "~/.secureConnect/context_aware_metadata.json"
CERTIFICATE_CONFIGURATION_DEFAULT_PATH = "~/.config/gcloud/certificate_config.json"
_CERT_PROVIDER_COMMAND = "cert_provider_command"
_CERT_REGEX = re.compile(
    b"-----BEGIN CERTIFICATE-----.+-----END CERTIFICATE-----\r?\n?", re.DOTALL
)

# support various format of key files, e.g.
# "-----BEGIN PRIVATE KEY-----...",
# "-----BEGIN EC PRIVATE KEY-----...",
# "-----BEGIN RSA PRIVATE KEY-----..."
# "-----BEGIN ENCRYPTED PRIVATE KEY-----"
_KEY_REGEX = re.compile(
    b"-----BEGIN [A-Z ]*PRIVATE KEY-----.+-----END [A-Z ]*PRIVATE KEY-----\r?\n?",
    re.DOTALL,
)

_LOGGER = logging.getLogger(__name__)


_PASSPHRASE_REGEX = re.compile(
    b"-----BEGIN PASSPHRASE-----(.+)-----END PASSPHRASE-----", re.DOTALL
)

# Temporary patch to accomodate incorrect cert config in Cloud Run prod environment.
_WELL_KNOWN_CLOUD_RUN_CERT_PATH = (
    "/var/run/secrets/workload-spiffe-credentials/certificates.pem"
)
_WELL_KNOWN_CLOUD_RUN_KEY_PATH = (
    "/var/run/secrets/workload-spiffe-credentials/private_key.pem"
)
_INCORRECT_CLOUD_RUN_CERT_PATH = (
    "/var/lib/volumes/certificate/workload-certificates/certificates.pem"
)
_INCORRECT_CLOUD_RUN_KEY_PATH = (
    "/var/lib/volumes/certificate/workload-certificates/private_key.pem"
)


def _check_config_path(config_path):
    """Checks for config file path. If it exists, returns the absolute path with user expansion;
    otherwise returns None.

    Args:
        config_path (str): The config file path for either context_aware_metadata.json or certificate_config.json for example

    Returns:
        str: absolute path if exists and None otherwise.
    """
    config_path = path.expanduser(config_path)
    if not path.exists(config_path):
        _LOGGER.debug("%s is not found.", config_path)
        return None
    return config_path


def _load_json_file(path):
    """Reads and loads JSON from the given path. Used to read both X509 workload certificate and
    secure connect configurations.

    Args:
        path (str): the path to read from.

    Returns:
        Dict[str, str]: The JSON stored at the file.

    Raises:
        google.auth.exceptions.ClientCertError: If failed to parse the file as JSON.
    """
    try:
        with open(path) as f:
            json_data = json.load(f)
    except ValueError as caught_exc:
        new_exc = exceptions.ClientCertError(caught_exc)
        raise new_exc from caught_exc

    return json_data


def _get_workload_cert_and_key(certificate_config_path=None):
    """Read the workload identity cert and key files specified in the certificate config provided.
    If no config path is provided, check the environment variable: "GOOGLE_API_CERTIFICATE_CONFIG"
    first, then the well known gcloud location: "~/.config/gcloud/certificate_config.json".

    Args:
        certificate_config_path (string): The certificate config path. If no path is provided,
        the environment variable will be checked first, then the well known gcloud location.

    Returns:
        Tuple[Optional[bytes], Optional[bytes]]: client certificate bytes in PEM format and key
            bytes in PEM format.

    Raises:
        google.auth.exceptions.ClientCertError: if problems occurs when retrieving
        the certificate or key information.
    """

    cert_path, key_path = _get_workload_cert_and_key_paths(certificate_config_path)

    if cert_path is None and key_path is None:
        return None, None

    return _read_cert_and_key_files(cert_path, key_path)


def _get_cert_config_path(certificate_config_path=None):
    """Get the certificate configuration path based on the following order:

    1: Explicit override, if set
    2: Environment variable, if set
    3: Well-known location

    Returns "None" if the selected config file does not exist.

    Args:
        certificate_config_path (string): The certificate config path. If provided, the well known
        location and environment variable will be ignored.

    Returns:
        The absolute path of the certificate config file, and None if the file does not exist.
    """

    if certificate_config_path is None:
        env_path = environ.get(environment_vars.GOOGLE_API_CERTIFICATE_CONFIG, None)
        if env_path is not None and env_path != "":
            certificate_config_path = env_path
        else:
            certificate_config_path = CERTIFICATE_CONFIGURATION_DEFAULT_PATH

    certificate_config_path = path.expanduser(certificate_config_path)
    if not path.exists(certificate_config_path):
        return None
    return certificate_config_path


def _get_workload_cert_and_key_paths(config_path):
    absolute_path = _get_cert_config_path(config_path)
    if absolute_path is None:
        return None, None

    data = _load_json_file(absolute_path)

    if "cert_configs" not in data:
        raise exceptions.ClientCertError(
            'Certificate config file {} is in an invalid format, a "cert configs" object is expected'.format(
                absolute_path
            )
        )
    cert_configs = data["cert_configs"]

    if "workload" not in cert_configs:
        raise exceptions.ClientCertError(
            'Certificate config file {} is in an invalid format, a "workload" cert config is expected'.format(
                absolute_path
            )
        )
    workload = cert_configs["workload"]

    if "cert_path" not in workload:
        raise exceptions.ClientCertError(
            'Certificate config file {} is in an invalid format, a "cert_path" is expected in the workload cert config'.format(
                absolute_path
            )
        )
    cert_path = workload["cert_path"]

    if "key_path" not in workload:
        raise exceptions.ClientCertError(
            'Certificate config file {} is in an invalid format, a "key_path" is expected in the workload cert config'.format(
                absolute_path
            )
        )
    key_path = workload["key_path"]

    # == BEGIN Temporary Cloud Run PATCH ==
    # See https://github.com/googleapis/google-auth-library-python/issues/1881
    if (cert_path == _INCORRECT_CLOUD_RUN_CERT_PATH) and (
        key_path == _INCORRECT_CLOUD_RUN_KEY_PATH
    ):
        if not path.exists(cert_path) and not path.exists(key_path):
            _LOGGER.debug(
                "Applying Cloud Run certificate path patch. "
                "Configured paths not found: %s, %s. "
                "Using well-known paths: %s, %s",
                cert_path,
                key_path,
                _WELL_KNOWN_CLOUD_RUN_CERT_PATH,
                _WELL_KNOWN_CLOUD_RUN_KEY_PATH,
            )
            cert_path = _WELL_KNOWN_CLOUD_RUN_CERT_PATH
            key_path = _WELL_KNOWN_CLOUD_RUN_KEY_PATH
    # == END Temporary Cloud Run PATCH ==

    return cert_path, key_path


def _read_cert_and_key_files(cert_path, key_path):
    cert_data = _read_cert_file(cert_path)
    key_data = _read_key_file(key_path)

    return cert_data, key_data


def _read_cert_file(cert_path):
    with open(cert_path, "rb") as cert_file:
        cert_data = cert_file.read()

    cert_match = re.findall(_CERT_REGEX, cert_data)
    if len(cert_match) != 1:
        raise exceptions.ClientCertError(
            "Certificate file {} is in an invalid format, a single PEM formatted certificate is expected".format(
                cert_path
            )
        )
    return cert_match[0]


def _read_key_file(key_path):
    with open(key_path, "rb") as key_file:
        key_data = key_file.read()

    key_match = re.findall(_KEY_REGEX, key_data)
    if len(key_match) != 1:
        raise exceptions.ClientCertError(
            "Private key file {} is in an invalid format, a single PEM formatted private key is expected".format(
                key_path
            )
        )

    return key_match[0]


def _run_cert_provider_command(command, expect_encrypted_key=False):
    """Run the provided command, and return client side mTLS cert, key and
    passphrase.

    Args:
        command (List[str]): cert provider command.
        expect_encrypted_key (bool): If encrypted private key is expected.

    Returns:
        Tuple[bytes, bytes, bytes]: client certificate bytes in PEM format, key
            bytes in PEM format and passphrase bytes.

    Raises:
        google.auth.exceptions.ClientCertError: if problems occurs when running
            the cert provider command or generating cert, key and passphrase.
    """
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
    except OSError as caught_exc:
        new_exc = exceptions.ClientCertError(caught_exc)
        raise new_exc from caught_exc

    # Check cert provider command execution error.
    if process.returncode != 0:
        raise exceptions.ClientCertError(
            "Cert provider command returns non-zero status code %s" % process.returncode
        )

    # Extract certificate (chain), key and passphrase.
    cert_match = re.findall(_CERT_REGEX, stdout)
    if len(cert_match) != 1:
        raise exceptions.ClientCertError("Client SSL certificate is missing or invalid")
    key_match = re.findall(_KEY_REGEX, stdout)
    if len(key_match) != 1:
        raise exceptions.ClientCertError("Client SSL key is missing or invalid")
    passphrase_match = re.findall(_PASSPHRASE_REGEX, stdout)

    if expect_encrypted_key:
        if len(passphrase_match) != 1:
            raise exceptions.ClientCertError("Passphrase is missing or invalid")
        if b"ENCRYPTED" not in key_match[0]:
            raise exceptions.ClientCertError("Encrypted private key is expected")
        return cert_match[0], key_match[0], passphrase_match[0].strip()

    if b"ENCRYPTED" in key_match[0]:
        raise exceptions.ClientCertError("Encrypted private key is not expected")
    if len(passphrase_match) > 0:
        raise exceptions.ClientCertError("Passphrase is not expected")
    return cert_match[0], key_match[0], None


def get_client_ssl_credentials(
    generate_encrypted_key=False,
    context_aware_metadata_path=CONTEXT_AWARE_METADATA_PATH,
    certificate_config_path=None,
):
    """Returns the client side certificate, private key and passphrase.

    We look for certificates and keys with the following order of priority:
        1. Certificate and key specified by certificate_config.json.
               Currently, only X.509 workload certificates are supported.
        2. Certificate and key specified by context aware metadata (i.e. SecureConnect).

    Args:
        generate_encrypted_key (bool): If set to True, encrypted private key
            and passphrase will be generated; otherwise, unencrypted private key
            will be generated and passphrase will be None. This option only
            affects keys obtained via context_aware_metadata.json.
        context_aware_metadata_path (str): The context_aware_metadata.json file path.
        certificate_config_path (str): The certificate_config.json file path.

    Returns:
        Tuple[bool, bytes, bytes, bytes]:
            A boolean indicating if cert, key and passphrase are obtained, the
            cert bytes and key bytes both in PEM format, and passphrase bytes.

    Raises:
        google.auth.exceptions.ClientCertError: if problems occurs when getting
            the cert, key and passphrase.
    """

    # 1.  Attempt to retrieve X.509 Workload cert and key.
    cert, key = _get_workload_cert_and_key(certificate_config_path)
    if cert and key:
        return True, cert, key, None

    # 2. Check for context aware metadata json
    metadata_path = _check_config_path(context_aware_metadata_path)

    if metadata_path:
        metadata_json = _load_json_file(metadata_path)

        if _CERT_PROVIDER_COMMAND not in metadata_json:
            raise exceptions.ClientCertError("Cert provider command is not found")

        command = metadata_json[_CERT_PROVIDER_COMMAND]

        if generate_encrypted_key and "--with_passphrase" not in command:
            command.append("--with_passphrase")

        # Execute the command.
        cert, key, passphrase = _run_cert_provider_command(
            command, expect_encrypted_key=generate_encrypted_key
        )
        return True, cert, key, passphrase

    return False, None, None, None


def get_client_cert_and_key(client_cert_callback=None):
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
        cert, key = client_cert_callback()
        return True, cert, key

    has_cert, cert, key, _ = get_client_ssl_credentials(generate_encrypted_key=False)
    return has_cert, cert, key


def decrypt_private_key(key, passphrase):
    """A helper function to decrypt the private key with the given passphrase.
    google-auth library doesn't support passphrase protected private key for
    mutual TLS channel. This helper function can be used to decrypt the
    passphrase protected private key in order to estalish mutual TLS channel.

    For example, if you have a function which produces client cert, passphrase
    protected private key and passphrase, you can convert it to a client cert
    callback function accepted by google-auth::

        from google.auth.transport import _mtls_helper

        def your_client_cert_function():
            return cert, encrypted_key, passphrase

        # callback accepted by google-auth for mutual TLS channel.
        def client_cert_callback():
            cert, encrypted_key, passphrase = your_client_cert_function()
            decrypted_key = _mtls_helper.decrypt_private_key(encrypted_key,
                passphrase)
            return cert, decrypted_key

    Args:
        key (bytes): The private key bytes in PEM format.
        passphrase (bytes): The passphrase bytes.

    Returns:
        bytes: The decrypted private key in PEM format.

    Raises:
        ImportError: If pyOpenSSL is not installed.
        OpenSSL.crypto.Error: If there is any problem decrypting the private key.
    """
    from OpenSSL import crypto

    # First convert encrypted_key_bytes to PKey object
    pkey = crypto.load_privatekey(crypto.FILETYPE_PEM, key, passphrase=passphrase)

    # Then dump the decrypted key bytes
    return crypto.dump_privatekey(crypto.FILETYPE_PEM, pkey)


def check_use_client_cert():
    """Returns boolean for whether the client certificate should be used for mTLS.

    If GOOGLE_API_USE_CLIENT_CERTIFICATE is set to true or false, a corresponding
    bool value will be returned. If the value is set to an unexpected string, it
    will default to False.
    If GOOGLE_API_USE_CLIENT_CERTIFICATE is unset, the value will be inferred
    by reading a file pointed at by GOOGLE_API_CERTIFICATE_CONFIG, and verifying
    it contains a "workload" section. If so, the function will return True,
    otherwise False.

    Returns:
        bool: Whether the client certificate should be used for mTLS connection.
    """
    use_client_cert = getenv("GOOGLE_API_USE_CLIENT_CERTIFICATE")
    # Check if the value of GOOGLE_API_USE_CLIENT_CERTIFICATE is set.
    if use_client_cert:
        return use_client_cert.lower() == "true"
    else:
        # Check if the value of GOOGLE_API_CERTIFICATE_CONFIG is set.
        cert_path = getenv("GOOGLE_API_CERTIFICATE_CONFIG")
        if cert_path:
            try:
                with open(cert_path, "r") as f:
                    content = json.load(f)
                    # verify json has workload key
                    content["cert_configs"]["workload"]
                    return True
            except (
                FileNotFoundError,
                OSError,
                KeyError,
                TypeError,
                json.JSONDecodeError,
            ) as e:
                _LOGGER.debug("error decoding certificate: %s", e)
        return False


def check_parameters_for_unauthorized_response(cached_cert):
    """Returns the cached and current cert fingerprint for reconfiguring mTLS.

    Args:
        cached_cert(bytes): The cached client certificate.

    Returns:
        bytes: The client callback cert bytes.
        bytes: The client callback key bytes.
        str: The base64-encoded SHA256 cached fingerprint.
        str: The base64-encoded SHA256 current cert fingerprint.
    """
    call_cert_bytes, call_key_bytes = call_client_cert_callback()
    cert_obj = _agent_identity_utils.parse_certificate(call_cert_bytes)
    current_cert_fingerprint = _agent_identity_utils.calculate_certificate_fingerprint(
        cert_obj
    )
    if cached_cert:
        cached_fingerprint = _agent_identity_utils.get_cached_cert_fingerprint(
            cached_cert
        )
    else:
        cached_fingerprint = current_cert_fingerprint
    return call_cert_bytes, call_key_bytes, cached_fingerprint, current_cert_fingerprint


def call_client_cert_callback():
    """Calls the client cert callback and returns the certificate and key."""
    _, cert_bytes, key_bytes, passphrase = get_client_ssl_credentials(
        generate_encrypted_key=True
    )
    return cert_bytes, key_bytes
