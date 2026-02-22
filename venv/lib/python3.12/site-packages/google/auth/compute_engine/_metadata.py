# Copyright 2016 Google LLC
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

"""Provides helper methods for talking to the Compute Engine metadata server.

See https://cloud.google.com/compute/docs/metadata for more details.
"""

import datetime
import http.client as http_client
import json
import logging
import os
from urllib.parse import urljoin

import requests

from google.auth import _helpers
from google.auth import environment_vars
from google.auth import exceptions
from google.auth import metrics
from google.auth import transport
from google.auth._exponential_backoff import ExponentialBackoff
from google.auth.compute_engine import _mtls


_LOGGER = logging.getLogger(__name__)

_GCE_DEFAULT_MDS_IP = "169.254.169.254"
_GCE_DEFAULT_HOST = "metadata.google.internal"
_GCE_DEFAULT_MDS_HOSTS = [_GCE_DEFAULT_HOST, _GCE_DEFAULT_MDS_IP]

# Environment variable GCE_METADATA_HOST is originally named
# GCE_METADATA_ROOT. For compatibility reasons, here it checks
# the new variable first; if not set, the system falls back
# to the old variable.
_GCE_METADATA_HOST = os.getenv(environment_vars.GCE_METADATA_HOST, None)
if not _GCE_METADATA_HOST:
    _GCE_METADATA_HOST = os.getenv(
        environment_vars.GCE_METADATA_ROOT, _GCE_DEFAULT_HOST
    )


def _validate_gce_mds_configured_environment():
    """Validates the GCE metadata server environment configuration for mTLS.

    mTLS is only supported when connecting to the default metadata server hosts.
    If we are in strict mode (which requires mTLS), ensure that the metadata host
    has not been overridden to a custom value (which means mTLS will fail).

    Raises:
        google.auth.exceptions.MutualTLSChannelError: if the environment
            configuration is invalid for mTLS.
    """
    mode = _mtls._parse_mds_mode()
    if mode == _mtls.MdsMtlsMode.STRICT:
        # mTLS is only supported when connecting to the default metadata host.
        # Raise an exception if we are in strict mode (which requires mTLS)
        # but the metadata host has been overridden to a custom MDS. (which means mTLS will fail)
        if _GCE_METADATA_HOST not in _GCE_DEFAULT_MDS_HOSTS:
            raise exceptions.MutualTLSChannelError(
                "Mutual TLS is required, but the metadata host has been overridden. "
                "mTLS is only supported when connecting to the default metadata host."
            )


def _get_metadata_root(use_mtls: bool):
    """Returns the metadata server root URL."""

    scheme = "https" if use_mtls else "http"
    return "{}://{}/computeMetadata/v1/".format(scheme, _GCE_METADATA_HOST)


def _get_metadata_ip_root(use_mtls: bool):
    """Returns the metadata server IP root URL."""
    scheme = "https" if use_mtls else "http"
    return "{}://{}".format(
        scheme, os.getenv(environment_vars.GCE_METADATA_IP, _GCE_DEFAULT_MDS_IP)
    )


_METADATA_FLAVOR_HEADER = "metadata-flavor"
_METADATA_FLAVOR_VALUE = "Google"
_METADATA_HEADERS = {_METADATA_FLAVOR_HEADER: _METADATA_FLAVOR_VALUE}

# Timeout in seconds to wait for the GCE metadata server when detecting the
# GCE environment.
try:
    _METADATA_DEFAULT_TIMEOUT = int(os.getenv(environment_vars.GCE_METADATA_TIMEOUT, 3))
except ValueError:  # pragma: NO COVER
    _METADATA_DEFAULT_TIMEOUT = 3

# The number of tries to perform when waiting for the GCE metadata server
# when detecting the GCE environment.
try:
    _METADATA_DETECT_RETRIES = int(
        os.getenv(environment_vars.GCE_METADATA_DETECT_RETRIES, 3)
    )
except ValueError:  # pragma: NO COVER
    _METADATA_DETECT_RETRIES = 3

# This is used to disable checking for the GCE metadata server and directly
# assuming it's not available.
_NO_GCE_CHECK = os.getenv(environment_vars.NO_GCE_CHECK) == "true"

# Detect GCE Residency
_GOOGLE = "Google"
_GCE_PRODUCT_NAME_FILE = "/sys/class/dmi/id/product_name"


def is_on_gce(request):
    """Checks to see if the code runs on Google Compute Engine

    Args:
        request (google.auth.transport.Request): A callable used to make
            HTTP requests.

    Returns:
        bool: True if the code runs on Google Compute Engine, False otherwise.
    """
    if _NO_GCE_CHECK:
        return False

    if ping(request):
        return True

    if os.name == "nt":
        # TODO: implement GCE residency detection on Windows
        return False

    # Detect GCE residency on Linux
    return detect_gce_residency_linux()


def detect_gce_residency_linux():
    """Detect Google Compute Engine residency by smbios check on Linux

    Returns:
        bool: True if the GCE product name file is detected, False otherwise.
    """
    try:
        with open(_GCE_PRODUCT_NAME_FILE, "r") as file_obj:
            content = file_obj.read().strip()

    except Exception:
        return False

    return content.startswith(_GOOGLE)


def _prepare_request_for_mds(request, use_mtls=False) -> None:
    """Prepares a request for the metadata server.

    This will check if mTLS should be used and mount the mTLS adapter if needed.

    Args:
        request (google.auth.transport.Request): A callable used to make
            HTTP requests.
        use_mtls (bool): Whether to use mTLS for the request.

    Returns:
        google.auth.transport.Request: A request object to use.
            If mTLS is enabled, the request will have the mTLS adapter mounted.
            Otherwise, the original request will be returned unchanged.
    """
    # Only modify the request if mTLS is enabled.
    if use_mtls:
        # Ensure the request has a session to mount the adapter to.
        if not request.session:
            request.session = requests.Session()

        adapter = _mtls.MdsMtlsAdapter()
        # Mount the adapter for all default GCE metadata hosts.
        for host in _GCE_DEFAULT_MDS_HOSTS:
            request.session.mount(f"https://{host}/", adapter)


def ping(
    request, timeout=_METADATA_DEFAULT_TIMEOUT, retry_count=_METADATA_DETECT_RETRIES
):
    """Checks to see if the metadata server is available.

    Args:
        request (google.auth.transport.Request): A callable used to make
            HTTP requests.
        timeout (int): How long to wait for the metadata server to respond.
        retry_count (int): How many times to attempt connecting to metadata
            server using above timeout.

    Returns:
        bool: True if the metadata server is reachable, False otherwise.
    """
    use_mtls = _mtls.should_use_mds_mtls()
    _prepare_request_for_mds(request, use_mtls=use_mtls)
    # NOTE: The explicit ``timeout`` is a workaround. The underlying
    #       issue is that resolving an unknown host on some networks will take
    #       20-30 seconds; making this timeout short fixes the issue, but
    #       could lead to false negatives in the event that we are on GCE, but
    #       the metadata resolution was particularly slow. The latter case is
    #       "unlikely".
    headers = _METADATA_HEADERS.copy()
    headers[metrics.API_CLIENT_HEADER] = metrics.mds_ping()

    backoff = ExponentialBackoff(total_attempts=retry_count)

    for attempt in backoff:
        try:
            response = request(
                url=_get_metadata_ip_root(use_mtls),
                method="GET",
                headers=headers,
                timeout=timeout,
            )

            metadata_flavor = response.headers.get(_METADATA_FLAVOR_HEADER)
            return (
                response.status == http_client.OK
                and metadata_flavor == _METADATA_FLAVOR_VALUE
            )

        except exceptions.TransportError as e:
            _LOGGER.warning(
                "Compute Engine Metadata server unavailable on "
                "attempt %s of %s. Reason: %s",
                attempt,
                retry_count,
                e,
            )

    return False


def get(
    request,
    path,
    root=None,
    params=None,
    recursive=False,
    retry_count=5,
    headers=None,
    return_none_for_not_found_error=False,
    timeout=_METADATA_DEFAULT_TIMEOUT,
):
    """Fetch a resource from the metadata server.

    Args:
        request (google.auth.transport.Request): A callable used to make
            HTTP requests.
        path (str): The resource to retrieve. For example,
            ``'instance/service-accounts/default'``.
        root (Optional[str]): The full path to the metadata server root. If not
            provided, the default root will be used.
        params (Optional[Mapping[str, str]]): A mapping of query parameter
            keys to values.
        recursive (bool): Whether to do a recursive query of metadata. See
            https://cloud.google.com/compute/docs/metadata#aggcontents for more
            details.
        retry_count (int): How many times to attempt connecting to metadata
            server using above timeout.
        headers (Optional[Mapping[str, str]]): Headers for the request.
        return_none_for_not_found_error (Optional[bool]): If True, returns None
            for 404 error instead of throwing an exception.
        timeout (int): How long to wait, in seconds for the metadata server to respond.

    Returns:
        Union[Mapping, str]: If the metadata server returns JSON, a mapping of
            the decoded JSON is returned. Otherwise, the response content is
            returned as a string.

    Raises:
        google.auth.exceptions.TransportError: if an error occurred while
            retrieving metadata.
        google.auth.exceptions.MutualTLSChannelError: if using mtls and the environment
            configuration is invalid for mTLS (for example, the metadata host
            has been overridden in strict mTLS mode).

    """
    use_mtls = _mtls.should_use_mds_mtls()
    # Prepare the request object for mTLS if needed.
    # This will create a new request object with the mTLS session.
    _prepare_request_for_mds(request, use_mtls=use_mtls)

    if root is None:
        root = _get_metadata_root(use_mtls)

    # mTLS is only supported when connecting to the default metadata host.
    # If we are in strict mode (which requires mTLS), ensure that the metadata host
    # has not been overridden to a non-default host value (which means mTLS will fail).
    _validate_gce_mds_configured_environment()

    base_url = urljoin(root, path)
    query_params = {} if params is None else params

    headers_to_use = _METADATA_HEADERS.copy()
    if headers:
        headers_to_use.update(headers)

    if recursive:
        query_params["recursive"] = "true"

    url = _helpers.update_query(base_url, query_params)

    backoff = ExponentialBackoff(total_attempts=retry_count)
    last_exception = None
    for attempt in backoff:
        try:
            response = request(
                url=url, method="GET", headers=headers_to_use, timeout=timeout
            )
            if response.status in transport.DEFAULT_RETRYABLE_STATUS_CODES:
                _LOGGER.warning(
                    "Compute Engine Metadata server unavailable on "
                    "attempt %s of %s. Response status: %s",
                    attempt,
                    retry_count,
                    response.status,
                )
                last_exception = None
                continue
            else:
                last_exception = None
                break

        except exceptions.TransportError as e:
            _LOGGER.warning(
                "Compute Engine Metadata server unavailable on "
                "attempt %s of %s. Reason: %s",
                attempt,
                retry_count,
                e,
            )
            last_exception = e
    else:
        if last_exception:
            raise exceptions.TransportError(
                "Failed to retrieve {} from the Google Compute Engine "
                "metadata service. Compute Engine Metadata server unavailable. "
                "Last exception: {}".format(url, last_exception)
            ) from last_exception
        else:
            error_details = (
                response.data.decode("utf-8")
                if hasattr(response.data, "decode")
                else response.data
            )
            raise exceptions.TransportError(
                "Failed to retrieve {} from the Google Compute Engine "
                "metadata service. Compute Engine Metadata server unavailable. "
                "Response status: {}\nResponse details:\n{}".format(
                    url, response.status, error_details
                )
            )

    content = _helpers.from_bytes(response.data)

    if response.status == http_client.NOT_FOUND and return_none_for_not_found_error:
        return None

    if response.status == http_client.OK:
        if (
            _helpers.parse_content_type(response.headers["content-type"])
            == "application/json"
        ):
            try:
                return json.loads(content)
            except ValueError as caught_exc:
                new_exc = exceptions.TransportError(
                    "Received invalid JSON from the Google Compute Engine "
                    "metadata service: {:.20}".format(content)
                )
                raise new_exc from caught_exc
        else:
            return content

    raise exceptions.TransportError(
        "Failed to retrieve {} from the Google Compute Engine "
        "metadata service. Status: {} Response:\n{}".format(
            url, response.status, response.data
        ),
        response,
    )


def get_project_id(request):
    """Get the Google Cloud Project ID from the metadata server.

    Args:
        request (google.auth.transport.Request): A callable used to make
            HTTP requests.

    Returns:
        str: The project ID

    Raises:
        google.auth.exceptions.TransportError: if an error occurred while
            retrieving metadata.
    """
    return get(request, "project/project-id")


def get_universe_domain(request):
    """Get the universe domain value from the metadata server.

    Args:
        request (google.auth.transport.Request): A callable used to make
            HTTP requests.

    Returns:
        str: The universe domain value. If the universe domain endpoint is not
        not found, return the default value, which is googleapis.com

    Raises:
        google.auth.exceptions.TransportError: if an error other than
            404 occurs while retrieving metadata.
    """
    universe_domain = get(
        request, "universe/universe-domain", return_none_for_not_found_error=True
    )
    if not universe_domain:
        return "googleapis.com"
    return universe_domain


def get_service_account_info(request, service_account="default"):
    """Get information about a service account from the metadata server.

    Args:
        request (google.auth.transport.Request): A callable used to make
            HTTP requests.
        service_account (str): The string 'default' or a service account email
            address. The determines which service account for which to acquire
            information.

    Returns:
        Mapping: The service account's information, for example::

            {
                'email': '...',
                'scopes': ['scope', ...],
                'aliases': ['default', '...']
            }

    Raises:
        google.auth.exceptions.TransportError: if an error occurred while
            retrieving metadata.
    """
    path = "instance/service-accounts/{0}/".format(service_account)
    # See https://cloud.google.com/compute/docs/metadata#aggcontents
    # for more on the use of 'recursive'.
    return get(request, path, params={"recursive": "true"})


def get_service_account_token(request, service_account="default", scopes=None):
    """Get the OAuth 2.0 access token for a service account.

    Args:
        request (google.auth.transport.Request): A callable used to make
            HTTP requests.
        service_account (str): The string 'default' or a service account email
            address. The determines which service account for which to acquire
            an access token.
        scopes (Optional[Union[str, List[str]]]): Optional string or list of
            strings with auth scopes.
    Returns:
        Tuple[str, datetime]: The access token and its expiration.

    Raises:
        google.auth.exceptions.TransportError: if an error occurred while
            retrieving metadata.
    """
    from google.auth import _agent_identity_utils

    params = {}
    if scopes:
        if not isinstance(scopes, str):
            scopes = ",".join(scopes)
        params["scopes"] = scopes

    cert = _agent_identity_utils.get_and_parse_agent_identity_certificate()
    if cert:
        if _agent_identity_utils.should_request_bound_token(cert):
            fingerprint = _agent_identity_utils.calculate_certificate_fingerprint(cert)
            params["bindCertificateFingerprint"] = fingerprint

    metrics_header = {
        metrics.API_CLIENT_HEADER: metrics.token_request_access_token_mds()
    }

    path = "instance/service-accounts/{0}/token".format(service_account)
    token_json = get(request, path, params=params, headers=metrics_header)
    token_expiry = _helpers.utcnow() + datetime.timedelta(
        seconds=token_json["expires_in"]
    )
    return token_json["access_token"], token_expiry
