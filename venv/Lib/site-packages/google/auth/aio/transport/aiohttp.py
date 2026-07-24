# Copyright 2024 Google LLC
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

"""Transport adapter for Asynchronous HTTP Requests based on aiohttp."""

import asyncio
import logging
from typing import AsyncGenerator, Mapping, Optional, TYPE_CHECKING, Union

try:
    import aiohttp  # type: ignore
except ImportError as caught_exc:  # pragma: NO COVER
    raise ImportError(
        "The aiohttp library is not installed from please install the aiohttp package to use the aiohttp transport."
    ) from caught_exc

from google.auth import _helpers
from google.auth import exceptions
from google.auth.aio import _helpers as _helpers_async
from google.auth.aio import transport

if TYPE_CHECKING:  # pragma: NO COVER
    from aiohttp import ClientTimeout  # type: ignore

else:
    try:
        from aiohttp import ClientTimeout
    except (ImportError, AttributeError):  # pragma: NO COVER
        ClientTimeout = None

_LOGGER = logging.getLogger(__name__)


class Response(transport.Response):
    """
    Represents an HTTP response and its data. It is returned by ``google.auth.aio.transport.sessions.AsyncAuthorizedSession``.

    Args:
        response (aiohttp.ClientResponse): An instance of aiohttp.ClientResponse.

    Attributes:
        status_code (int): The HTTP status code of the response.
        headers (Mapping[str, str]): The HTTP headers of the response.
    """

    def __init__(self, response: aiohttp.ClientResponse):
        self._response = response

    @property
    @_helpers.copy_docstring(transport.Response)
    def status_code(self) -> int:
        return self._response.status

    @property
    @_helpers.copy_docstring(transport.Response)
    def headers(self) -> Mapping[str, str]:
        return {key: value for key, value in self._response.headers.items()}

    @_helpers.copy_docstring(transport.Response)
    async def content(self, chunk_size: int = 1024) -> AsyncGenerator[bytes, None]:
        try:
            async for chunk in self._response.content.iter_chunked(
                chunk_size
            ):  # pragma: no branch
                yield chunk
        except aiohttp.ClientPayloadError as exc:
            raise exceptions.ResponseError(
                "Failed to read from the payload stream."
            ) from exc

    @_helpers.copy_docstring(transport.Response)
    async def read(self) -> bytes:
        try:
            return await self._response.read()
        except aiohttp.ClientResponseError as exc:
            raise exceptions.ResponseError("Failed to read the response body.") from exc

    @_helpers.copy_docstring(transport.Response)
    async def close(self):
        self._response.close()


class Request(transport.Request):
    """Asynchronous Requests request adapter.

    This class is used internally for making requests using aiohttp
    in a consistent way. If you use :class:`google.auth.aio.transport.sessions.AsyncAuthorizedSession`
    you do not need to construct or use this class directly.

    This class can be useful if you want to configure a Request callable
    with a custom ``aiohttp.ClientSession`` in :class:`AuthorizedSession` or if
    you want to manually refresh a :class:`~google.auth.aio.credentials.Credentials` instance::

        import aiohttp
        import google.auth.aio.transport.aiohttp

        # Default example:
        request = google.auth.aio.transport.aiohttp.Request()
        await credentials.refresh(request)

        # Custom aiohttp Session Example:
        session = session=aiohttp.ClientSession(auto_decompress=False)
        request = google.auth.aio.transport.aiohttp.Request(session=session)
        auth_session = google.auth.aio.transport.sessions.AsyncAuthorizedSession(auth_request=request)

    Args:
        session (aiohttp.ClientSession): An instance :class:`aiohttp.ClientSession` used
            to make HTTP requests. If not specified, a session will be created.

    .. automethod:: __call__
    """

    def __init__(self, session: Optional[aiohttp.ClientSession] = None):
        self._session = session
        self._closed = False

    async def __call__(
        self,
        url: str,
        method: str = "GET",
        body: Optional[bytes] = None,
        headers: Optional[Mapping[str, str]] = None,
        timeout: Union[float, ClientTimeout] = transport._DEFAULT_TIMEOUT_SECONDS,
        **kwargs,
    ) -> transport.Response:
        """
        Make an HTTP request using aiohttp.

        Args:
            url (str): The URL to be requested.
            method (Optional[str]):
                The HTTP method to use for the request. Defaults to 'GET'.
            body (Optional[bytes]):
                The payload or body in HTTP request.
            headers (Optional[Mapping[str, str]]):
                Request headers.
            timeout (float): The number of seconds to wait for a
                response from the server. If not specified or if None, the
                requests default timeout will be used.
            kwargs: Additional arguments passed through to the underlying
                aiohttp :meth:`aiohttp.Session.request` method.

        Returns:
            google.auth.aio.transport.Response: The HTTP response.

        Raises:
            - google.auth.exceptions.TransportError: If the request fails or if the session is closed.
            - google.auth.exceptions.TimeoutError: If the request times out.
        """

        try:
            if self._closed:
                raise exceptions.TransportError("session is closed.")

            if not self._session:
                self._session = aiohttp.ClientSession()

            if isinstance(timeout, aiohttp.ClientTimeout):
                client_timeout = timeout
            else:
                client_timeout = aiohttp.ClientTimeout(total=timeout)
            _helpers.request_log(_LOGGER, method, url, body, headers)
            response = await self._session.request(
                method,
                url,
                data=body,
                headers=headers,
                timeout=client_timeout,
                **kwargs,
            )
            await _helpers_async.response_log_async(_LOGGER, response)
            return Response(response)

        except aiohttp.ClientError as caught_exc:
            client_exc = exceptions.TransportError(f"Failed to send request to {url}.")
            raise client_exc from caught_exc

        except asyncio.TimeoutError as caught_exc:
            if isinstance(timeout, aiohttp.ClientTimeout):
                timeout_seconds = timeout.total
            else:
                timeout_seconds = timeout
            timeout_exc = exceptions.TimeoutError(
                f"Request timed out after {timeout_seconds} seconds."
            )
            raise timeout_exc from caught_exc

    async def close(self) -> None:
        """
        Close the underlying aiohttp session to release the acquired resources.
        """
        if not self._closed and self._session:
            await self._session.close()
        self._closed = True

    def _clone(self) -> "Request":
        """Creates an independent copy of this request adapter.

        Clones the connection settings, trace configurations, and session defaults
        (headers, cookies, basic auth, and timeouts).

        Only standard `aiohttp.TCPConnector` and `aiohttp.UnixConnector` connectors
        are supported. The DNS resolver is not copied to avoid closing shared resolver
        resources.

        Returns:
            google.auth.aio.transport.aiohttp.Request: A new request adapter.

        Raises:
            google.auth.exceptions.TransportError: If the transport is closed, or if the
                session uses an unsupported connector.
        """
        if self._closed:
            raise exceptions.TransportError("Cannot clone a closed transport.")

        if not self._session:
            new_session = aiohttp.ClientSession(
                auto_decompress=False,
                trust_env=True,
            )
            return Request(session=new_session)

        session_kwargs: dict = {
            "auto_decompress": False,
            "trust_env": getattr(self._session, "_trust_env", True),
        }

        # Copy underlying connection pool settings (SSL context, IP bindings, limits).
        orig_connector = getattr(self._session, "_connector", None)
        if orig_connector and not orig_connector.closed:
            if isinstance(orig_connector, aiohttp.TCPConnector):
                # We explicitly do not copy the resolver. The connector
                # owns the resolver, and closing the cloned session would
                # close the shared resolver, breaking the original session.
                session_kwargs["connector"] = aiohttp.TCPConnector(
                    ssl=getattr(orig_connector, "_ssl", None),  # type: ignore
                    limit=getattr(orig_connector, "_limit", 100),
                    limit_per_host=getattr(orig_connector, "_limit_per_host", 0),
                    force_close=getattr(orig_connector, "_force_close", False),
                    local_addr=_helpers_async._get_local_addr(orig_connector),
                )
            elif getattr(aiohttp, "UnixConnector", None) and isinstance(
                orig_connector, getattr(aiohttp, "UnixConnector")
            ):
                path = getattr(orig_connector, "_path", None)
                if path:
                    session_kwargs["connector"] = aiohttp.UnixConnector(
                        path=path,
                        limit=getattr(orig_connector, "_limit", 100),
                        force_close=getattr(orig_connector, "_force_close", False),
                    )
            else:
                raise exceptions.TransportError(
                    f"Unsupported connector type for cloning: {type(orig_connector)}"
                )

        # Preserve distributed tracing configurations.
        trace_configs = getattr(self._session, "_trace_configs", None)
        if trace_configs:
            session_kwargs["trace_configs"] = list(trace_configs)

        # Copy session-level defaults (headers, cookies, auth, timeout).
        for attr_name, kwarg_name in [
            ("_default_headers", "headers"),
            ("_cookie_jar", "cookie_jar"),
            ("_default_auth", "auth"),
            ("_timeout", "timeout"),
            ("_json_serialize", "json_serialize"),
        ]:
            val = getattr(self._session, attr_name, None)
            if val is not None:
                session_kwargs[kwarg_name] = val

        return Request(session=aiohttp.ClientSession(**session_kwargs))  # type: ignore
