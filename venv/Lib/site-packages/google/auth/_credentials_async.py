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


"""Interfaces for credentials."""

import abc
import inspect

from google.auth import _regional_access_boundary_utils
from google.auth import credentials


class Credentials(credentials.Credentials, metaclass=abc.ABCMeta):
    """Async inherited credentials class from google.auth.credentials.
    The added functionality is the before_request call which requires
    async/await syntax.
    All credentials have a :attr:`token` that is used for authentication and
    may also optionally set an :attr:`expiry` to indicate when the token will
    no longer be valid.

    Most credentials will be :attr:`invalid` until :meth:`refresh` is called.
    Credentials can do this automatically before the first HTTP request in
    :meth:`before_request`.

    Although the token and expiration will change as the credentials are
    :meth:`refreshed <refresh>` and used, credentials should be considered
    immutable. Various credentials will accept configuration such as private
    keys, scopes, and other options. These options are not changeable after
    construction. Some classes will provide mechanisms to copy the credentials
    with modifications such as :meth:`ScopedCredentials.with_scopes`.
    """

    async def before_request(self, request, method, url, headers):
        """Performs credential-specific before request logic.

        Refreshes the credentials if necessary, then calls :meth:`apply` to
        apply the token to the authentication header.

        Args:
            request (google.auth.transport.Request): The object used to make
                HTTP requests.
            method (str): The request's HTTP method or the RPC method being
                invoked.
            url (str): The request's URI or the RPC service's URI.
            headers (Mapping): The request's headers.
        """
        # pylint: disable=unused-argument
        # (Subclasses may use these arguments to ascertain information about
        # the http request.)

        if not self.valid:
            if inspect.iscoroutinefunction(self.refresh):
                await self.refresh(request)
            else:
                self.refresh(request)

        if inspect.iscoroutinefunction(self._after_refresh):
            await self._after_refresh(request, method, url, headers)
        else:
            self._after_refresh(request, method, url, headers)

        self.apply(headers)

    def _after_refresh(self, request, method, url, headers):
        """Hook for subclasses to perform actions after refresh but before
        applying credentials to headers.

        Args:
            request (google.auth.transport.Request): The object used to make
                HTTP requests.
            method (str): The request's HTTP method or the RPC method being
                invoked.
            url (str): The request's URI or the RPC service's URI.
            headers (Mapping[str, str]): The request's headers.
        """
        pass


class CredentialsWithQuotaProject(credentials.CredentialsWithQuotaProject):
    """Abstract base for credentials supporting ``with_quota_project`` factory"""


class AnonymousCredentials(credentials.AnonymousCredentials, Credentials):
    """Credentials that do not provide any authentication information.

    These are useful in the case of services that support anonymous access or
    local service emulators that do not use credentials. This class inherits
    from the sync anonymous credentials file, but is kept if async credentials
    is initialized and we would like anonymous credentials.
    """


class ReadOnlyScoped(credentials.ReadOnlyScoped, metaclass=abc.ABCMeta):
    """Interface for credentials whose scopes can be queried.

    OAuth 2.0-based credentials allow limiting access using scopes as described
    in `RFC6749 Section 3.3`_.
    If a credential class implements this interface then the credentials either
    use scopes in their implementation.

    Some credentials require scopes in order to obtain a token. You can check
    if scoping is necessary with :attr:`requires_scopes`::

        if credentials.requires_scopes:
            # Scoping is required.
            credentials = _credentials_async.with_scopes(scopes=['one', 'two'])

    Credentials that require scopes must either be constructed with scopes::

        credentials = SomeScopedCredentials(scopes=['one', 'two'])

    Or must copy an existing instance using :meth:`with_scopes`::

        scoped_credentials = _credentials_async.with_scopes(scopes=['one', 'two'])

    Some credentials have scopes but do not allow or require scopes to be set,
    these credentials can be used as-is.

    .. _RFC6749 Section 3.3: https://tools.ietf.org/html/rfc6749#section-3.3
    """


class Scoped(credentials.Scoped):
    """Interface for credentials whose scopes can be replaced while copying.

    OAuth 2.0-based credentials allow limiting access using scopes as described
    in `RFC6749 Section 3.3`_.
    If a credential class implements this interface then the credentials either
    use scopes in their implementation.

    Some credentials require scopes in order to obtain a token. You can check
    if scoping is necessary with :attr:`requires_scopes`::

        if credentials.requires_scopes:
            # Scoping is required.
            credentials = _credentials_async.create_scoped(['one', 'two'])

    Credentials that require scopes must either be constructed with scopes::

        credentials = SomeScopedCredentials(scopes=['one', 'two'])

    Or must copy an existing instance using :meth:`with_scopes`::

        scoped_credentials = credentials.with_scopes(scopes=['one', 'two'])

    Some credentials have scopes but do not allow or require scopes to be set,
    these credentials can be used as-is.

    .. _RFC6749 Section 3.3: https://tools.ietf.org/html/rfc6749#section-3.3
    """


def with_scopes_if_required(credentials, scopes):
    """Creates a copy of the credentials with scopes if scoping is required.

    This helper function is useful when you do not know (or care to know) the
    specific type of credentials you are using (such as when you use
    :func:`google.auth.default`). This function will call
    :meth:`Scoped.with_scopes` if the credentials are scoped credentials and if
    the credentials require scoping. Otherwise, it will return the credentials
    as-is.

    Args:
        credentials (google.auth.credentials.Credentials): The credentials to
            scope if necessary.
        scopes (Sequence[str]): The list of scopes to use.

    Returns:
        google.auth._credentials_async.Credentials: Either a new set of scoped
            credentials, or the passed in credentials instance if no scoping
            was required.
    """
    if isinstance(credentials, Scoped) and credentials.requires_scopes:
        return credentials.with_scopes(scopes)
    else:
        return credentials


class Signing(credentials.Signing, metaclass=abc.ABCMeta):
    """Interface for credentials that can cryptographically sign messages."""


class CredentialsWithRegionalAccessBoundary(
    Credentials, credentials.CredentialsWithRegionalAccessBoundary
):
    """Async base for credentials supporting regional access boundary configuration."""

    def __init__(self):
        super().__init__()
        self._rab_manager.refresh_manager = (
            _regional_access_boundary_utils._AsyncRegionalAccessBoundaryRefreshManager()
        )

    def __setstate__(self, state):
        super().__setstate__(state)
        self._rab_manager.refresh_manager = (
            _regional_access_boundary_utils._AsyncRegionalAccessBoundaryRefreshManager()
        )

    async def _after_refresh(self, request, method, url, headers):
        """Triggers the Regional Access Boundary lookup asynchronously if necessary."""
        await self._maybe_start_regional_access_boundary_refresh_async(request, url)

    async def _maybe_start_regional_access_boundary_refresh_async(self, request, url):
        """Starts a background refresh or performs a blocking refresh asynchronously.

        Args:
            request (google.auth.aio.transport.Request): The object used to make
                HTTP requests.
            url (str): The URL of the request.
        """
        # Do not perform a lookup if the request is for a regional endpoint.
        if self._is_regional_endpoint(url):
            return

        # A refresh is only needed if the feature is enabled.
        if not self._is_regional_access_boundary_lookup_required():
            return

        # Trigger background or blocking refresh if needed.
        await self._rab_manager.maybe_start_refresh_async(self, request)

    async def _lookup_regional_access_boundary(self, request, fail_fast=False):
        """Calls the Regional Access Boundary lookup API asynchronously.

        Args:
            request (google.auth.aio.transport.Request): The object used to make
                HTTP requests.
            fail_fast (bool): Whether the lookup should fail fast (short timeout, no retries).

        Returns:
            Optional[Dict[str, str]]: The Regional Access Boundary information
                returned by the lookup API, or None if the lookup failed.
        """
        url_builder = self._build_regional_access_boundary_lookup_url
        if inspect.iscoroutinefunction(url_builder):
            url = await url_builder(request=request)
        else:
            url = url_builder(request=request)

        if not url:
            return None

        headers = {}
        self._apply(headers)

        from google.oauth2 import _client_async

        return await _client_async._lookup_regional_access_boundary(
            request, url, headers=headers, fail_fast=fail_fast
        )
