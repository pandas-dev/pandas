import gzip
import json as json_module
import re
from json.decoder import JSONDecodeError
from typing import Any
from typing import Callable
from typing import List
from typing import Mapping
from typing import MutableMapping
from typing import Optional
from typing import Pattern
from typing import Tuple
from typing import Union
from urllib.parse import parse_qsl
from urllib.parse import urlparse

from requests import PreparedRequest
from urllib3.util.url import parse_url


def _filter_dict_recursively(
    dict1: Mapping[Any, Any], dict2: Mapping[Any, Any]
) -> Mapping[Any, Any]:
    filtered_dict = {}
    for k, val in dict1.items():
        if k in dict2:
            if isinstance(val, dict):
                val = _filter_dict_recursively(val, dict2[k])
            filtered_dict[k] = val

    return filtered_dict


def body_matcher(params: str, *, allow_blank: bool = False) -> Callable[..., Any]:
    def match(request: PreparedRequest) -> Tuple[bool, str]:
        reason = ""
        if isinstance(request.body, bytes):
            request_body = request.body.decode("utf-8")
        else:
            request_body = str(request.body)
        valid = True if request_body == params else False
        if not valid:
            reason = f"request.body doesn't match {params} doesn't match {request_body}"
        return valid, reason

    return match


def urlencoded_params_matcher(
    params: Optional[Mapping[str, str]], *, allow_blank: bool = False
) -> Callable[..., Any]:
    """
    Matches URL encoded data

    :param params: (dict) data provided to 'data' arg of request
    :return: (func) matcher
    """

    def match(request: PreparedRequest) -> Tuple[bool, str]:
        reason = ""
        request_body = request.body
        qsl_body = (
            dict(parse_qsl(request_body, keep_blank_values=allow_blank))  # type: ignore[type-var]
            if request_body
            else {}
        )
        params_dict = params or {}
        valid = params is None if request_body is None else params_dict == qsl_body
        if not valid:
            reason = (
                f"request.body doesn't match: {qsl_body} doesn't match {params_dict}"
            )

        return valid, reason

    return match


def json_params_matcher(
    params: Optional[Union[Mapping[str, Any], List[Any]]], *, strict_match: bool = True
) -> Callable[..., Any]:
    """Matches JSON encoded data of request body.

    Parameters
    ----------
    params : dict or list
        JSON object provided to 'json' arg of request or a part of it if used in
        conjunction with ``strict_match=False``.
    strict_match : bool, default=True
        Applied only when JSON object is a dictionary.
        If set to ``True``, validates that all keys of JSON object match.
        If set to ``False``, original request may contain additional keys.


    Returns
    -------
    Callable
        Matcher function.

    """

    def match(request: PreparedRequest) -> Tuple[bool, str]:
        reason = ""
        request_body = request.body
        json_params = (params or {}) if not isinstance(params, list) else params
        try:
            if isinstance(request.body, bytes):
                try:
                    request_body = request.body.decode("utf-8")
                except UnicodeDecodeError:
                    request_body = gzip.decompress(request.body).decode("utf-8")
            json_body = json_module.loads(request_body) if request_body else {}

            if (
                not strict_match
                and isinstance(json_body, dict)
                and isinstance(json_params, dict)
            ):
                # filter down to just the params specified in the matcher
                json_body = _filter_dict_recursively(json_body, json_params)

            valid = params is None if request_body is None else json_params == json_body

            if not valid:
                reason = f"request.body doesn't match: {json_body} doesn't match {json_params}"
                if not strict_match:
                    reason += (
                        "\nNote: You use non-strict parameters check, "
                        "to change it use `strict_match=True`."
                    )

        except JSONDecodeError:
            valid = False
            reason = (
                "request.body doesn't match: JSONDecodeError: Cannot parse request.body"
            )

        return valid, reason

    return match


def fragment_identifier_matcher(identifier: Optional[str]) -> Callable[..., Any]:
    def match(request: PreparedRequest) -> Tuple[bool, str]:
        reason = ""
        url_fragment = urlparse(request.url).fragment
        if identifier:
            url_fragment_qsl = sorted(parse_qsl(url_fragment))  # type: ignore[type-var]
            identifier_qsl = sorted(parse_qsl(identifier))
            valid = identifier_qsl == url_fragment_qsl
        else:
            valid = not url_fragment

        if not valid:
            reason = (
                "URL fragment identifier is different: "  # type: ignore[str-bytes-safe]
                f"{identifier} doesn't match {url_fragment}"
            )

        return valid, reason

    return match


def query_param_matcher(
    params: Optional[MutableMapping[str, Any]], *, strict_match: bool = True
) -> Callable[..., Any]:
    """Matcher to match 'params' argument in request.

    Parameters
    ----------
    params : dict
        The same as provided to request or a part of it if used in
        conjunction with ``strict_match=False``.
    strict_match : bool, default=True
        If set to ``True``, validates that all parameters match.
        If set to ``False``, original request may contain additional parameters.


    Returns
    -------
    Callable
        Matcher function.

    """

    params_dict = params or {}

    for k, v in params_dict.items():
        if isinstance(v, (int, float)):
            params_dict[k] = str(v)

    def match(request: PreparedRequest) -> Tuple[bool, str]:
        reason = ""
        request_params = request.params  # type: ignore[attr-defined]
        request_params_dict = request_params or {}

        if not strict_match:
            # filter down to just the params specified in the matcher
            request_params_dict = {
                k: v for k, v in request_params_dict.items() if k in params_dict
            }

        valid = sorted(params_dict.items()) == sorted(request_params_dict.items())

        if not valid:
            reason = f"Parameters do not match. {request_params_dict} doesn't match {params_dict}"
            if not strict_match:
                reason += (
                    "\nYou can use `strict_match=True` to do a strict parameters check."
                )

        return valid, reason

    return match


def query_string_matcher(query: Optional[str]) -> Callable[..., Any]:
    """
    Matcher to match query string part of request

    :param query: (str), same as constructed by request
    :return: (func) matcher
    """

    def match(request: PreparedRequest) -> Tuple[bool, str]:
        reason = ""
        data = parse_url(request.url or "")
        request_query = data.query

        request_qsl = sorted(parse_qsl(request_query)) if request_query else {}
        matcher_qsl = sorted(parse_qsl(query)) if query else {}

        valid = not query if request_query is None else request_qsl == matcher_qsl

        if not valid:
            reason = (
                "Query string doesn't match. "
                f"{dict(request_qsl)} doesn't match {dict(matcher_qsl)}"
            )

        return valid, reason

    return match


def request_kwargs_matcher(kwargs: Optional[Mapping[str, Any]]) -> Callable[..., Any]:
    """
    Matcher to match keyword arguments provided to request

    :param kwargs: (dict), keyword arguments, same as provided to request
    :return: (func) matcher
    """

    def match(request: PreparedRequest) -> Tuple[bool, str]:
        reason = ""
        kwargs_dict = kwargs or {}
        # validate only kwargs that were requested for comparison, skip defaults
        req_kwargs = request.req_kwargs  # type: ignore[attr-defined]
        request_kwargs = {k: v for k, v in req_kwargs.items() if k in kwargs_dict}

        valid = (
            not kwargs_dict
            if not request_kwargs
            else sorted(kwargs_dict.items()) == sorted(request_kwargs.items())
        )

        if not valid:
            reason = (
                f"Arguments don't match: {request_kwargs} doesn't match {kwargs_dict}"
            )

        return valid, reason

    return match


def multipart_matcher(
    files: Mapping[str, Any], data: Optional[Mapping[str, str]] = None
) -> Callable[..., Any]:
    """
    Matcher to match 'multipart/form-data' content-type.
    This function constructs request body and headers from provided 'data' and 'files'
    arguments and compares to actual request

    :param files: (dict), same as provided to request
    :param data: (dict), same as provided to request
    :return: (func) matcher
    """
    if not files:
        raise TypeError("files argument cannot be empty")

    prepared = PreparedRequest()
    prepared.headers = {"Content-Type": ""}  # type: ignore[assignment]
    prepared.prepare_body(data=data, files=files)

    def get_boundary(content_type: str) -> str:
        """
        Parse 'boundary' value from header.

        :param content_type: (str) headers["Content-Type"] value
        :return: (str) boundary value
        """
        if "boundary=" not in content_type:
            return ""

        return content_type.split("boundary=")[1]

    def match(request: PreparedRequest) -> Tuple[bool, str]:
        reason = "multipart/form-data doesn't match. "
        if "Content-Type" not in request.headers:
            return False, reason + "Request is missing the 'Content-Type' header"

        request_boundary = get_boundary(request.headers["Content-Type"])
        prepared_boundary = get_boundary(prepared.headers["Content-Type"])

        # replace boundary value in header and in body, since by default
        # urllib3.filepost.encode_multipart_formdata dynamically calculates
        # random boundary alphanumeric value
        request_content_type = request.headers["Content-Type"]
        prepared_content_type = prepared.headers["Content-Type"].replace(
            prepared_boundary, request_boundary
        )

        request_body = request.body
        prepared_body = prepared.body or ""

        if isinstance(prepared_body, bytes):
            # since headers always come as str, need to convert to bytes
            prepared_boundary = prepared_boundary.encode("utf-8")  # type: ignore[assignment]
            request_boundary = request_boundary.encode("utf-8")  # type: ignore[assignment]

        prepared_body = prepared_body.replace(
            prepared_boundary, request_boundary  # type: ignore[arg-type]
        )

        headers_valid = prepared_content_type == request_content_type
        if not headers_valid:
            return (
                False,
                reason
                + "Request headers['Content-Type'] is different. {} isn't equal to {}".format(
                    request_content_type, prepared_content_type
                ),
            )

        body_valid = prepared_body == request_body
        if not body_valid:
            return (
                False,
                reason
                + "Request body differs. {} aren't equal {}".format(  # type: ignore[str-bytes-safe]
                    request_body, prepared_body
                ),
            )

        return True, ""

    return match


def header_matcher(
    headers: Mapping[str, Union[str, Pattern[str]]], strict_match: bool = False
) -> Callable[..., Any]:
    """
    Matcher to match 'headers' argument in request using the responses library.

    Because ``requests`` will send several standard headers in addition to what
    was specified by your code, request headers that are additional to the ones
    passed to the matcher are ignored by default. You can change this behaviour
    by passing ``strict_match=True``.

    :param headers: (dict), same as provided to request
    :param strict_match: (bool), whether headers in addition to those specified
                         in the matcher should cause the match to fail.
    :return: (func) matcher
    """

    def _compare_with_regex(request_headers: Union[Mapping[Any, Any], Any]) -> bool:
        if strict_match and len(request_headers) != len(headers):
            return False

        for k, v in headers.items():
            if request_headers.get(k) is not None:
                if isinstance(v, re.Pattern):
                    if re.match(v, request_headers[k]) is None:
                        return False
                else:
                    if not v == request_headers[k]:
                        return False
            else:
                return False

        return True

    def match(request: PreparedRequest) -> Tuple[bool, str]:
        request_headers: Union[Mapping[Any, Any], Any] = request.headers or {}

        if not strict_match:
            # filter down to just the headers specified in the matcher
            request_headers = {k: v for k, v in request_headers.items() if k in headers}

        valid = _compare_with_regex(request_headers)

        if not valid:
            return (
                False,
                f"Headers do not match: {request_headers} doesn't match {headers}",
            )

        return valid, ""

    return match
