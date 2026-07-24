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

import contextlib
import json
import logging
import os
from os import environ, getenv, path
import re
import subprocess
import sys
import tempfile
from typing import cast, Generator, List, Optional, Tuple, Union

from google.auth import _agent_identity_utils
from google.auth import _cloud_sdk
from google.auth import environment_vars
from google.auth import exceptions

CONTEXT_AWARE_METADATA_PATH = "~/.secureConnect/context_aware_metadata.json"

# Default gcloud config path, to be used with path.expanduser for cross-platform compatibility.
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


class _MemfdCreationError(OSError):
    """Raised when Linux in-memory virtual file creation (memfd) fails."""

    pass


def _can_read(path: Optional[str]) -> bool:
    if path is None:
        return True
    try:
        with open(path, "rb"):
            pass
        return True
    except OSError:
        return False


@contextlib.contextmanager
def secure_cert_key_paths(
    cert: Union[bytes, str, None],
    key: Union[bytes, str, None],
    passphrase: Optional[bytes] = None,
) -> Generator[Tuple[Optional[str], Optional[str], Optional[bytes]], None, None]:
    """Provides secure file paths for certificate and key.

    This function is implemented as a context manager generator to ensure that
    any temporary resources (such as in-memory virtual files or encrypted physical
    temp files) are automatically cleaned up and securely wiped when the context exits.

    It supports mixed inputs (e.g. passing one as a string path and the other as bytes).
    If a parameter is already a string path or None, it is passed through as-is, and
    only raw bytes are written to temporary storage.

    Args:
        cert (Union[str, bytes, None]): Certificate path, raw PEM content bytes, or None.
        key (Union[str, bytes, None]): Private key path, raw PEM content bytes, or None.
        passphrase (Optional[bytes]): Optional passphrase for the private key.

    Yields:
        Tuple[str, str, Optional[bytes]]: The certificate path, key path, and
            the passphrase needed to load the key (either the user's original,
            or the newly generated one if Tier 3 had to encrypt the key).

    Raises:
        OSError: If temporary file creation or writing fails during the Tier 3 fallback.
    """
    # Normalize PEM strings to bytes so they are written to temporary storage.
    # We check for "-----BEGIN " to distinguish between file paths and PEM payloads.
    if isinstance(cert, str) and "-----BEGIN " in cert:
        cert = cert.encode("utf-8")
    if isinstance(key, str) and "-----BEGIN " in key:
        key = key.encode("utf-8")

    # Tier 1: Pass-through (No-op). If the caller already provided file paths,
    # we yield them directly to avoid any unnecessary file creation.
    if isinstance(cert, str) and isinstance(key, str):
        yield cert, key, passphrase
        return

    # If a value is a string path, it is passed through. If bytes, we will write
    # it to temporary storage. None values are also passed through as-is.
    cert_bytes = cert if isinstance(cert, bytes) else None
    key_bytes = key if isinstance(key, bytes) else None

    # Tier 2: Linux RAM-backed virtual files. If supported by the OS, we write
    # the bytes to anonymous in-memory files using memfd_create. This yields
    # /proc/self/fd/... paths, keeping the private key entirely in memory.
    if sys.platform == "linux" and hasattr(os, "memfd_create"):
        try:
            with _memfd_cert_key_paths(cert_bytes, key_bytes) as (cert_path, key_path):
                # Handle cases where path exists but might be restricted.
                if (cert_path is None or os.path.exists(cert_path)) and (
                    key_path is None or os.path.exists(key_path)
                ):
                    if _can_read(cert_path) and _can_read(key_path):
                        yield cast(str, cert_path or cert), cast(
                            str, key_path or key
                        ), passphrase
                        return
        except _MemfdCreationError:
            pass  # Fallback to Tier 3 on failure.

    # Tier 3: Fallback Encrypted Temp Files. If in-memory files are not supported
    # (macOS/Windows), we write to disk. To protect the key, we encrypt plaintext
    # keys on-the-fly and securely wipe the files with null bytes during cleanup.
    with _tempfile_cert_key_paths(cert_bytes, key_bytes, passphrase) as (
        cert_path,
        key_path,
        new_passphrase,
    ):
        yield cast(str, cert_path or cert), cast(str, key_path or key), new_passphrase


def _encrypt_key_if_plaintext(
    key_bytes: bytes, passphrase: Optional[bytes]
) -> Tuple[bytes, Optional[bytes]]:
    """Encrypts a plaintext PEM key if necessary, returning the bytes and passphrase.

    If the key is already encrypted, or if parsing/encryption fails, the key is
    returned as-is (plaintext) as a fallback. This allows the caller (underlying SSL
    context) to attempt loading the key directly and handle any failures.
    """
    import cryptography
    from cryptography.hazmat.primitives import serialization
    import secrets

    try:
        pkey = serialization.load_pem_private_key(key_bytes, password=None)
        # It's plaintext, encrypt it.
        target_passphrase = passphrase
        if target_passphrase is None:
            target_passphrase = secrets.token_hex(32).encode("utf-8")
        elif isinstance(target_passphrase, str):
            target_passphrase = target_passphrase.encode("utf-8")

        encrypted_content = pkey.private_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PrivateFormat.PKCS8,
            encryption_algorithm=serialization.BestAvailableEncryption(
                target_passphrase
            ),
        )
        return encrypted_content, target_passphrase
    except (ValueError, TypeError, cryptography.exceptions.UnsupportedAlgorithm):
        # Likely already encrypted, invalid, or unsupported algorithm, return as-is.
        return key_bytes, passphrase


def _secure_wipe_and_remove(file_path: str):
    """Overwrites a file with null bytes before deleting it.

    This is an extra security measure to make file recovery harder. However, on modern
    solid-state drives (SSDs), the hardware optimizes where data is written, meaning
    the original private key bytes might still physically remain on the storage chips
    until the drive cleans them up.
    """
    if not os.path.exists(file_path):
        return
    try:
        size = os.path.getsize(file_path)
        with open(file_path, "r+b") as f:
            f.write(b"\0" * size)
            f.flush()
            os.fsync(f.fileno())
    except OSError:
        pass  # Ignore permission/lock errors during cleanup.
    finally:
        try:
            os.remove(file_path)
        except OSError:
            pass


@contextlib.contextmanager
def _memfd_cert_key_paths(
    cert_bytes: Optional[bytes], key_bytes: Optional[bytes]
) -> Generator[Tuple[Optional[str], Optional[str]], None, None]:
    """Creates secure, in-memory virtual files on Linux using memfd_create.

    Yields:
        Tuple[Optional[str], Optional[str]]: In-memory file paths pointing to
            the active descriptors (e.g., '/proc/self/fd/3').
    """
    cleanup_fds = []
    paths: List[Optional[str]] = []

    try:
        try:
            for data, name in [(cert_bytes, "mtls_cert"), (key_bytes, "mtls_key")]:
                if data is not None:
                    # MFD_CLOEXEC prevents FD leaks to spawned subprocesses.
                    fd = os.memfd_create(name, os.MFD_CLOEXEC)  # type: ignore[attr-defined]
                    cleanup_fds.append(fd)
                    with os.fdopen(fd, "wb", closefd=False) as f:
                        f.write(data)
                    paths.append(f"/proc/self/fd/{fd}")
                else:
                    paths.append(None)
        except (OSError, AttributeError) as exc:
            raise _MemfdCreationError(
                "Failed to create in-memory virtual files"
            ) from exc

        cert_path, key_path = paths
        yield cert_path, key_path
    finally:
        # Closing the descriptors automatically frees the RAM allocation.
        for fd in cleanup_fds:
            try:
                os.close(fd)
            except OSError:
                pass


def _write_secure_tempfile(fd: int, data: bytes) -> None:
    """Writes data to a file descriptor, securely flushes to disk, and closes it."""
    try:
        f = os.fdopen(fd, "wb")
    except BaseException:
        try:
            os.close(fd)
        except OSError:
            pass
        raise

    with f:
        f.write(data)
        f.flush()
        try:
            os.fsync(f.fileno())
        except OSError:
            pass


@contextlib.contextmanager
def _tempfile_cert_key_paths(
    cert_bytes: Optional[bytes],
    key_bytes: Optional[bytes],
    passphrase: Optional[bytes],
) -> Generator[Tuple[Optional[str], Optional[str], Optional[bytes]], None, None]:
    """Creates secure temporary file paths on disk, encrypting private keys.

    Yields:
        Tuple[Optional[str], Optional[str], Optional[bytes]]: The temporary file
            paths and the passphrase needed to load the key.
    """
    # Prioritize RAM-backed /dev/shm to avoid writing secrets to physical storage.
    tmp_dir = (
        "/dev/shm"
        if os.path.isdir("/dev/shm") and os.access("/dev/shm", os.W_OK)
        else None
    )
    cleanup_files: List[Optional[str]] = [None, None]
    new_passphrase = passphrase
    cert_data = cert_bytes
    key_data = None
    if key_bytes is not None:
        key_data, new_passphrase = _encrypt_key_if_plaintext(key_bytes, passphrase)

    try:
        for i, data in enumerate([cert_data, key_data]):
            if data is not None:
                try:
                    fd, path = tempfile.mkstemp(dir=tmp_dir)
                except OSError:
                    fd, path = tempfile.mkstemp(dir=None)

                cleanup_files[i] = path
                _write_secure_tempfile(fd, data)

        yield cleanup_files[0], cleanup_files[1], new_passphrase
    finally:
        cert_cleanup_path = cleanup_files[0]
        key_cleanup_path = cleanup_files[1]

        try:
            if key_cleanup_path:
                _secure_wipe_and_remove(key_cleanup_path)
        except Exception:
            pass
        finally:
            if cert_cleanup_path:
                try:
                    if os.path.exists(cert_cleanup_path):
                        os.remove(cert_cleanup_path)
                except OSError:
                    pass


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


def _get_workload_cert_and_key(
    certificate_config_path=None, include_context_aware=True
):
    """Read the workload identity cert and key files specified in the certificate config provided.
    If no config path is provided, check the environment variable: "GOOGLE_API_CERTIFICATE_CONFIG"
    first, then the well known gcloud location: "~/.config/gcloud/certificate_config.json".

    Args:
        certificate_config_path (string): The certificate config path. If no path is provided,
        the environment variable will be checked first, then the well known gcloud location.
        include_context_aware (bool): If context aware metadata path should be checked for the
        SecureConnect mTLS configuration.

    Returns:
        Tuple[Optional[bytes], Optional[bytes]]: client certificate bytes in PEM format and key
            bytes in PEM format.

    Raises:
        google.auth.exceptions.ClientCertError: if problems occurs when retrieving
        the certificate or key information.
    """

    cert_path, key_path = _get_workload_cert_and_key_paths(
        certificate_config_path, include_context_aware
    )

    if cert_path is None and key_path is None:
        return None, None

    return _read_cert_and_key_files(cert_path, key_path)


def _get_cert_config_path(certificate_config_path=None, include_context_aware=True):
    """Get the certificate configuration path based on the following order:

    1: Explicit override, if set
    2: Environment variable, if set
    3: Well-known location

    Returns "None" if the selected config file does not exist.

    Args:
        certificate_config_path (string): The certificate config path. If provided, the well known
        location and environment variable will be ignored.
        include_context_aware (bool): If context aware metadata path should be checked for the
        SecureConnect mTLS configuration.

    Returns:
        The absolute path of the certificate config file, and None if the file does not exist.
    """

    source = "function argument"
    is_explicit = True
    if certificate_config_path is None:
        env_path = environ.get(environment_vars.GOOGLE_API_CERTIFICATE_CONFIG, None)
        if env_path is not None and env_path != "":
            certificate_config_path = env_path
            source = (
                f"environment variable {environment_vars.GOOGLE_API_CERTIFICATE_CONFIG}"
            )
        else:
            env_path = environ.get(
                environment_vars.CLOUDSDK_CONTEXT_AWARE_CERTIFICATE_CONFIG_FILE_PATH,
                None,
            )
            if include_context_aware and env_path is not None and env_path != "":
                certificate_config_path = env_path
                source = f"environment variable {environment_vars.CLOUDSDK_CONTEXT_AWARE_CERTIFICATE_CONFIG_FILE_PATH}"
            else:
                certificate_config_path = os.path.join(
                    _cloud_sdk.get_config_path(), "certificate_config.json"
                )
                is_explicit = False

    certificate_config_path = path.expanduser(certificate_config_path)
    if not path.exists(certificate_config_path):
        if is_explicit:
            _LOGGER.debug(
                "Certificate configuration file explicitly specified via %s at %s does not exist",
                source,
                certificate_config_path,
            )
        return None
    return certificate_config_path


def _get_workload_cert_and_key_paths(config_path, include_context_aware=True):
    absolute_path = _get_cert_config_path(config_path, include_context_aware)
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

    # We return None, None if the expected workload fields are not present.
    # The certificate config might be present for other types of connections (e.g. gECC),
    # and we want to gracefully fallback to testing other mTLS configurations
    # like SecureConnect instead of throwing an exception.

    if "workload" not in cert_configs:
        return None, None
    workload = cert_configs["workload"]

    if "cert_path" not in workload or "key_path" not in workload:
        raise exceptions.ClientCertError(
            'Workload certificate configuration is missing "cert_path" or "key_path" in {}'.format(
                absolute_path
            )
        )
    cert_path = workload["cert_path"]
    key_path = workload["key_path"]

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
        ValueError: If there is any problem decrypting the private key.
    """
    if isinstance(key, str):
        key = key.encode("utf-8")
    if isinstance(passphrase, str):
        passphrase = passphrase.encode("utf-8")

    from cryptography.hazmat.primitives import serialization

    # First convert encrypted_key_bytes to PKey object
    pkey = serialization.load_pem_private_key(key, password=passphrase)

    # Then dump the decrypted key bytes
    return pkey.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.PKCS8,
        encryption_algorithm=serialization.NoEncryption(),
    )


def _check_use_client_cert_env():
    use_client_cert = getenv(
        environment_vars.GOOGLE_API_USE_CLIENT_CERTIFICATE
    ) or getenv(environment_vars.CLOUDSDK_CONTEXT_AWARE_USE_CLIENT_CERTIFICATE)

    if use_client_cert:
        return use_client_cert.lower() == "true"
    return None


def check_use_client_cert():
    """Returns boolean for whether the client certificate should be used for mTLS.

    If GOOGLE_API_USE_CLIENT_CERTIFICATE is set to true or false, a corresponding
    bool value will be returned. If the value is set to an unexpected string, it
    will default to False.
    If GOOGLE_API_USE_CLIENT_CERTIFICATE is unset, the value will be inferred
    as True (auto-enabled) if a workload config file exists (pointed at by
    GOOGLE_API_CERTIFICATE_CONFIG or CLOUDSDK_CONTEXT_AWARE_CERTIFICATE_CONFIG_FILE_PATH,
    or the default path like ~/.config/gcloud/certificate_config.json)
    containing a "workload" section.
    Otherwise, it returns False.

    Returns:
        bool: Whether the client certificate should be used for mTLS connection.
    """
    env_override = _check_use_client_cert_env()
    if env_override is not None:
        return env_override

    # Auto-enablement checks (when GOOGLE_API_USE_CLIENT_CERTIFICATE is not set)

    # Check if a workload config file exists.
    cert_path = _get_cert_config_path(include_context_aware=True)

    if cert_path:
        try:
            with open(cert_path, "r") as f:
                content = json.load(f)
        except (FileNotFoundError, OSError, json.JSONDecodeError) as e:
            _LOGGER.debug(
                "mTLS auto-enablement failed: Could not read/parse certificate file at %s. Error: %s",
                cert_path,
                e,
            )
            return False

        # Structural validation
        if isinstance(content, dict):
            cert_configs = content.get("cert_configs")
            if isinstance(cert_configs, dict) and "workload" in cert_configs:
                return True

        # If we got here, the file exists but the expected structure is missing
        _LOGGER.debug(
            "mTLS auto-enablement failed: Certificate configuration file at %s is missing the required ['cert_configs']['workload'] section.",
            cert_path,
        )
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
