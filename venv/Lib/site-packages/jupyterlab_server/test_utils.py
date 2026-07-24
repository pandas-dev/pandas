# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Testing utils."""
from __future__ import annotations

import json
import os
import sys
from http.cookies import SimpleCookie
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import tornado.httpclient
import tornado.web
from openapi_core import V30RequestValidator, V30ResponseValidator
from openapi_core.spec.paths import Spec
from openapi_core.validation.request.datatypes import RequestParameters
from tornado.httpclient import HTTPRequest, HTTPResponse
from werkzeug.datastructures import Headers, ImmutableMultiDict

from jupyterlab_server.spec import get_openapi_spec

HERE = Path(os.path.dirname(__file__)).resolve()

with open(HERE / "test_data" / "app-settings" / "overrides.json", encoding="utf-8") as fid:
    big_unicode_string = json.load(fid)["@jupyterlab/unicode-extension:plugin"]["comment"]


class TornadoOpenAPIRequest:
    """
    Converts a torando request to an OpenAPI one
    """

    def __init__(self, request: HTTPRequest, spec: Spec):
        """Initialize the request."""
        self.request = request
        self.spec = spec
        if request.url is None:
            msg = "Request URL is missing"  # type:ignore[unreachable]
            raise RuntimeError(msg)
        self._url_parsed = urlparse(request.url)

        cookie: SimpleCookie = SimpleCookie()
        cookie.load(request.headers.get("Set-Cookie", ""))
        cookies = {}
        for key, morsel in cookie.items():
            cookies[key] = morsel.value

        # extract the path
        o = urlparse(request.url)

        # gets deduced by path finder against spec
        path: dict = {}

        self.parameters = RequestParameters(
            query=ImmutableMultiDict(parse_qs(o.query)),
            header=dict(request.headers),
            cookie=ImmutableMultiDict(cookies),
            path=path,
        )

    @property
    def content_type(self) -> str:
        return "application/json"

    @property
    def host_url(self) -> str:
        url = self.request.url
        return url[: url.index("/lab")]

    @property
    def path(self) -> str:
        # extract the best matching url
        # work around lack of support for path parameters which can contain slashes
        # https://github.com/OAI/OpenAPI-Specification/issues/892
        url = None
        o = urlparse(self.request.url)
        for path_ in self.spec["paths"]:
            if url:
                continue  # type:ignore[unreachable]
            has_arg = "{" in path_
            path = path_[: path_.index("{")] if has_arg else path_
            if path in o.path:
                u = o.path[o.path.index(path) :]
                if not has_arg and len(u) == len(path):
                    url = u
                if has_arg and not u.endswith("/"):
                    url = u[: len(path)] + r"foo"

        if url is None:
            msg = f"Could not find matching pattern for {o.path}"
            raise ValueError(msg)
        return url

    @property
    def method(self) -> str:
        method = self.request.method
        return method and method.lower() or ""

    @property
    def body(self) -> bytes | None:
        if self.request.body is None:
            return None  # type:ignore[unreachable]
        if not isinstance(self.request.body, bytes):
            msg = "Request body is invalid"  # type:ignore[unreachable]
            raise AssertionError(msg)
        return self.request.body

    @property
    def mimetype(self) -> str:
        # Order matters because all tornado requests
        # include Accept */* which does not necessarily match the content type
        request = self.request
        return (
            request.headers.get("Content-Type")
            or request.headers.get("Accept")
            or "application/json"
        )


class TornadoOpenAPIResponse:
    """A tornado open API response."""

    def __init__(self, response: HTTPResponse):
        """Initialize the response."""
        self.response = response

    @property
    def data(self) -> bytes | None:
        if not isinstance(self.response.body, bytes):
            msg = "Response body is invalid"  # type:ignore[unreachable]
            raise AssertionError(msg)
        return self.response.body

    @property
    def status_code(self) -> int:
        return int(self.response.code)

    @property
    def content_type(self) -> str:
        return "application/json"

    @property
    def mimetype(self) -> str:
        return str(self.response.headers.get("Content-Type", "application/json"))

    @property
    def headers(self) -> Headers:
        return Headers(dict(self.response.headers))


def validate_request(response: HTTPResponse) -> None:
    """Validate an API request"""
    openapi_spec = get_openapi_spec()

    request = TornadoOpenAPIRequest(response.request, openapi_spec)
    V30RequestValidator(openapi_spec).validate(request)

    torn_response = TornadoOpenAPIResponse(response)
    V30ResponseValidator(openapi_spec).validate(request, torn_response)


def maybe_patch_ioloop() -> None:
    """a windows 3.8+ patch for the asyncio loop"""
    if (
        sys.platform.startswith("win")
        and tornado.version_info < (6, 1)
        and sys.version_info >= (3, 8)
    ):
        try:
            from asyncio import WindowsProactorEventLoopPolicy, WindowsSelectorEventLoopPolicy
        except ImportError:
            pass
            # not affected
        else:
            from asyncio import get_event_loop_policy, set_event_loop_policy

            if type(get_event_loop_policy()) is WindowsProactorEventLoopPolicy:
                # WindowsProactorEventLoopPolicy is not compatible with tornado 6
                # fallback to the pre-3.8 default of Selector
                set_event_loop_policy(WindowsSelectorEventLoopPolicy())


def expected_http_error(
    error: Exception, expected_code: int, expected_message: str | None = None
) -> bool:
    """Check that the error matches the expected output error."""
    e = error.value  # type:ignore[attr-defined]
    if isinstance(e, tornado.web.HTTPError):
        if expected_code != e.status_code:
            return False
        if expected_message is not None and expected_message != str(e):
            return False
        return True
    if any(
        [
            isinstance(e, tornado.httpclient.HTTPClientError),
            isinstance(e, tornado.httpclient.HTTPError),
        ]
    ):
        if expected_code != e.code:
            return False
        if expected_message:
            message = json.loads(e.response.body.decode())["message"]
            if expected_message != message:
                return False
        return True

    return False
