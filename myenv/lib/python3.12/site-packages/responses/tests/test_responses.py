import inspect
import os
import re
import warnings
from io import BufferedReader
from io import BytesIO
from typing import Any
from typing import List
from typing import Optional
from unittest.mock import Mock
from unittest.mock import patch

import pytest
import requests
import urllib3
from requests.exceptions import ChunkedEncodingError
from requests.exceptions import ConnectionError
from requests.exceptions import HTTPError
from requests.exceptions import RetryError
from urllib3.util.retry import Retry

import responses
from responses import BaseResponse
from responses import Call
from responses import CallbackResponse
from responses import PassthroughResponse
from responses import Response
from responses import matchers
from responses import registries


def assert_reset():
    assert len(responses.mock.registered()) == 0
    assert len(responses.calls) == 0


def assert_response(
    resp: Any, body: Optional[Any] = None, content_type: "Optional[str]" = "text/plain"
) -> None:
    assert resp.status_code == 200
    assert resp.reason == "OK"
    if content_type is not None:
        assert resp.headers["Content-Type"] == content_type
    else:
        assert "Content-Type" not in resp.headers
    assert resp.text == body


def assert_params(resp, expected):
    assert hasattr(resp, "request"), "Missing request"
    assert hasattr(
        resp.request, "params"
    ), "Missing params on request that responses should add"
    assert getattr(resp.request, "params") == expected, "Incorrect parameters"


def test_response():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com", body=b"test")
        resp = requests.get("http://example.com")
        assert_response(resp, "test")
        assert len(responses.calls) == 1
        assert responses.calls[0].request.url == "http://example.com/"
        assert responses.calls[0].response.content == b"test"

        resp = requests.get("http://example.com?foo=bar")
        assert_response(resp, "test")
        assert len(responses.calls) == 2
        assert responses.calls[1].request.url == "http://example.com/?foo=bar"
        assert responses.calls[1].response.content == b"test"

    run()
    assert_reset()


def test_response_encoded():
    @responses.activate
    def run():
        # Path contains urlencoded =/()[]
        url = "http://example.org/foo.bar%3D%2F%28%29%5B%5D"
        responses.add(responses.GET, url, body="it works", status=200)
        resp = requests.get(url)
        assert_response(resp, "it works")

    run()
    assert_reset()


def test_response_with_instance():
    @responses.activate
    def run():
        responses.add(
            responses.Response(method=responses.GET, url="http://example.com")
        )
        resp = requests.get("http://example.com")
        assert_response(resp, "")
        assert len(responses.calls) == 1
        assert responses.calls[0].request.url == "http://example.com/"

        resp = requests.get("http://example.com?foo=bar")
        assert_response(resp, "")
        assert len(responses.calls) == 2
        assert responses.calls[1].request.url == "http://example.com/?foo=bar"

    run()
    assert_reset()


@pytest.mark.parametrize(
    "original,replacement",
    [
        ("http://example.com/two", "http://example.com/two"),
        (
            Response(method=responses.GET, url="http://example.com/two"),
            Response(
                method=responses.GET, url="http://example.com/two", body="testtwo"
            ),
        ),
        (
            re.compile(r"http://example\.com/two"),
            re.compile(r"http://example\.com/two"),
        ),
    ],
)
def test_replace(original, replacement):  # type: ignore[misc]
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com/one", body="test1")

        if isinstance(original, BaseResponse):
            responses.add(original)
        else:
            responses.add(responses.GET, original, body="test2")

        responses.add(responses.GET, "http://example.com/three", body="test3")
        responses.add(
            responses.GET, re.compile(r"http://example\.com/four"), body="test3"
        )

        if isinstance(replacement, BaseResponse):
            responses.replace(replacement)
        else:
            responses.replace(responses.GET, replacement, body="testtwo")

        resp = requests.get("http://example.com/two")
        assert_response(resp, "testtwo")

    run()
    assert_reset()


@pytest.mark.parametrize(
    "original,replacement",
    [
        ("http://example.com/one", re.compile(r"http://example\.com/one")),
        (re.compile(r"http://example\.com/one"), "http://example.com/one"),
    ],
)
def test_replace_error(original, replacement):  # type: ignore[misc]
    @responses.activate
    def run():
        responses.add(responses.GET, original)
        with pytest.raises(ValueError) as excinfo:
            responses.replace(responses.GET, replacement)
        assert "Response is not registered for URL %s" % replacement in str(
            excinfo.value
        )

    run()
    assert_reset()


def test_replace_response_object_error():
    @responses.activate
    def run():
        responses.add(Response(method=responses.GET, url="http://example.com/one"))
        with pytest.raises(ValueError) as excinfo:
            responses.replace(
                Response(method=responses.GET, url="http://example.com/two")
            )
        assert "Response is not registered for URL http://example.com/two" in str(
            excinfo.value
        )

    run()
    assert_reset()


@pytest.mark.parametrize(
    "original,replacement",
    [
        ("http://example.com/two", "http://example.com/two"),
        (
            Response(method=responses.GET, url="http://example.com/two"),
            Response(
                method=responses.GET, url="http://example.com/two", body="testtwo"
            ),
        ),
        (
            re.compile(r"http://example\.com/two"),
            re.compile(r"http://example\.com/two"),
        ),
    ],
)
def test_upsert_replace(original, replacement):  # type: ignore[misc]
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com/one", body="test1")

        if isinstance(original, BaseResponse):
            responses.add(original)
        else:
            responses.add(responses.GET, original, body="test2")

        if isinstance(replacement, BaseResponse):
            responses.upsert(replacement)
        else:
            responses.upsert(responses.GET, replacement, body="testtwo")

        resp = requests.get("http://example.com/two")
        assert_response(resp, "testtwo")

    run()
    assert_reset()


@pytest.mark.parametrize(
    "original,replacement",
    [
        ("http://example.com/two", "http://example.com/two"),
        (
            Response(method=responses.GET, url="http://example.com/two"),
            Response(
                method=responses.GET, url="http://example.com/two", body="testtwo"
            ),
        ),
        (
            re.compile(r"http://example\.com/two"),
            re.compile(r"http://example\.com/two"),
        ),
    ],
)
def test_upsert_add(original, replacement):  # type: ignore[misc]
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com/one", body="test1")

        if isinstance(replacement, BaseResponse):
            responses.upsert(replacement)
        else:
            responses.upsert(responses.GET, replacement, body="testtwo")

        resp = requests.get("http://example.com/two")
        assert_response(resp, "testtwo")

    run()
    assert_reset()


def test_remove():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com/zero")
        responses.add(responses.GET, "http://example.com/one")
        responses.add(responses.GET, "http://example.com/two")
        responses.add(responses.GET, re.compile(r"http://example\.com/three"))
        responses.add(responses.GET, re.compile(r"http://example\.com/four"))
        re.purge()
        responses.remove(responses.GET, "http://example.com/two")
        responses.remove(Response(method=responses.GET, url="http://example.com/zero"))
        responses.remove(responses.GET, re.compile(r"http://example\.com/four"))

        with pytest.raises(ConnectionError):
            requests.get("http://example.com/zero")
        requests.get("http://example.com/one")
        with pytest.raises(ConnectionError):
            requests.get("http://example.com/two")
        requests.get("http://example.com/three")
        with pytest.raises(ConnectionError):
            requests.get("http://example.com/four")

    run()
    assert_reset()


@pytest.mark.parametrize(
    "args1,kwargs1,args2,kwargs2,expected",
    [
        ((responses.GET, "a"), {}, (responses.GET, "a"), {}, True),
        ((responses.GET, "a"), {}, (responses.GET, "b"), {}, False),
        ((responses.GET, "a"), {}, (responses.POST, "a"), {}, False),
        (
            (responses.GET, "a"),
            {"match_querystring": True},
            (responses.GET, "a"),
            {},
            True,
        ),
    ],
)
def test_response_equality(args1, kwargs1, args2, kwargs2, expected):  # type: ignore[misc]
    o1 = BaseResponse(*args1, **kwargs1)
    o2 = BaseResponse(*args2, **kwargs2)
    assert (o1 == o2) is expected
    assert (o1 != o2) is not expected


def test_response_equality_different_objects():
    o1 = BaseResponse(method=responses.GET, url="a")
    o2 = "str"
    assert (o1 == o2) is False
    assert (o1 != o2) is True


def test_connection_error():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com")

        with pytest.raises(ConnectionError):
            requests.get("http://example.com/foo")

        assert len(responses.calls) == 1
        assert responses.calls[0].request.url == "http://example.com/foo"
        assert type(responses.calls[0].response) is ConnectionError
        assert responses.calls[0].response.request

    run()
    assert_reset()


def test_match_querystring():
    @responses.activate
    def run():
        url = "http://example.com?test=1&foo=bar"
        responses.add(responses.GET, url, match_querystring=True, body=b"test")
        resp = requests.get("http://example.com?test=1&foo=bar")
        assert_response(resp, "test")
        resp = requests.get("http://example.com?foo=bar&test=1")
        assert_response(resp, "test")
        resp = requests.get("http://example.com/?foo=bar&test=1")
        assert_response(resp, "test")

    run()
    assert_reset()


def test_match_querystring_empty():
    @responses.activate
    def run():
        responses.add(
            responses.GET, "http://example.com", body=b"test", match_querystring=True
        )
        resp = requests.get("http://example.com")
        assert_response(resp, "test")
        resp = requests.get("http://example.com/")
        assert_response(resp, "test")
        with pytest.raises(ConnectionError):
            requests.get("http://example.com?query=foo")

    run()
    assert_reset()


def test_match_querystring_error():
    @responses.activate
    def run():
        responses.add(
            responses.GET, "http://example.com/?test=1", match_querystring=True
        )

        with pytest.raises(ConnectionError):
            requests.get("http://example.com/foo/?test=2")

    run()
    assert_reset()


def test_match_querystring_regex():
    @responses.activate
    def run():
        """Note that `match_querystring` value shouldn't matter when passing a
        regular expression"""

        responses.add(
            responses.GET,
            re.compile(r"http://example\.com/foo/\?test=1"),
            body="test1",
            match_querystring=True,
        )

        resp = requests.get("http://example.com/foo/?test=1")
        assert_response(resp, "test1")

        responses.add(
            responses.GET,
            re.compile(r"http://example\.com/foo/\?test=2"),
            body="test2",
            match_querystring=False,
        )

        resp = requests.get("http://example.com/foo/?test=2")
        assert_response(resp, "test2")

    run()
    assert_reset()


def test_match_querystring_error_regex():
    @responses.activate
    def run():
        """Note that `match_querystring` value shouldn't matter when passing a
        regular expression"""

        responses.add(
            responses.GET,
            re.compile(r"http://example\.com/foo/\?test=1"),
            match_querystring=True,
        )

        with pytest.raises(ConnectionError):
            requests.get("http://example.com/foo/?test=3")

        responses.add(
            responses.GET,
            re.compile(r"http://example\.com/foo/\?test=2"),
            match_querystring=False,
        )

        with pytest.raises(ConnectionError):
            requests.get("http://example.com/foo/?test=4")

    run()
    assert_reset()


def test_match_querystring_auto_activates():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com?test=1", body=b"test")
        resp = requests.get("http://example.com?test=1")
        assert_response(resp, "test")
        with pytest.raises(ConnectionError):
            requests.get("http://example.com/?test=2")

    run()
    assert_reset()


def test_match_querystring_missing_key():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com?foo=1&bar=2", body=b"test")
        with pytest.raises(ConnectionError):
            requests.get("http://example.com/?foo=1&baz=2")

        with pytest.raises(ConnectionError):
            requests.get("http://example.com/?bar=2&fez=1")

    run()
    assert_reset()


def test_accept_string_body():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(responses.GET, url, body="test")
        resp = requests.get(url)
        assert_response(resp, "test")

    run()
    assert_reset()


def test_accept_json_body():
    @responses.activate
    def run():
        content_type = "application/json"

        url = "http://example.com/"
        responses.add(responses.GET, url, json={"message": "success"})
        resp = requests.get(url)
        assert_response(resp, '{"message": "success"}', content_type)

        url = "http://example.com/1/"
        responses.add(responses.GET, url, json=[])
        resp = requests.get(url)
        assert_response(resp, "[]", content_type)

    run()
    assert_reset()


def test_no_content_type():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(responses.GET, url, body="test", content_type=None)
        resp = requests.get(url)
        assert_response(resp, "test", content_type=None)

    run()
    assert_reset()


def test_arbitrary_status_code():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(responses.GET, url, body="test", status=419)
        resp = requests.get(url)
        assert resp.status_code == 419
        assert resp.reason is None

    run()
    assert_reset()


def test_throw_connection_error_explicit():
    @responses.activate
    def run():
        url = "http://example.com"
        exception = HTTPError("HTTP Error")
        responses.add(responses.GET, url, exception)

        with pytest.raises(HTTPError) as HE:
            requests.get(url)

        assert str(HE.value) == "HTTP Error"

    run()
    assert_reset()


def test_callback():
    body = b"test callback"
    status = 400
    reason = "Bad Request"
    headers = {
        "foo": "bar",
        "Content-Type": "application/json",
        "Content-Length": "13",
    }
    url = "http://example.com/"

    def request_callback(_request):
        return status, headers, body

    @responses.activate
    def run():
        responses.add_callback(responses.GET, url, request_callback)
        resp = requests.get(url)
        assert resp.text == "test callback"
        assert resp.status_code == status
        assert resp.reason == reason
        assert "bar" == resp.headers.get("foo")
        assert "application/json" == resp.headers.get("Content-Type")
        assert "13" == resp.headers.get("Content-Length")

    run()
    assert_reset()


def test_deprecated_package_attributes():
    """Validates that deprecation warning is raised when package attributes are called."""
    # keep separate context manager to avoid leakage
    with pytest.deprecated_call():
        responses.assert_all_requests_are_fired

    with pytest.deprecated_call():
        responses.passthru_prefixes

    with pytest.deprecated_call():
        responses.target


def test_callback_deprecated_stream_argument():
    with pytest.deprecated_call():
        CallbackResponse(responses.GET, "url", lambda x: x, stream=False)


def test_callback_deprecated_match_querystring_argument():
    with pytest.deprecated_call():
        CallbackResponse(responses.GET, "url", lambda x: x, match_querystring=False)


def test_callback_match_querystring_default_false():
    """
    Test to ensure that by default 'match_querystring' in 'add_callback' is set to False
    and does not raise deprecation
    see: https://github.com/getsentry/responses/issues/464 and related PR
    """
    body = b"test callback"
    status = 200
    params = {"hello": "world", "I am": "a big test"}
    headers = {"foo": "bar"}
    url = "http://example.com/"

    def request_callback(_request):
        return status, headers, body

    @responses.activate
    def run():
        responses.add_callback(responses.GET, url, request_callback, content_type=None)
        resp = requests.get(url, params=params)
        assert resp.text == "test callback"
        assert resp.status_code == status
        assert "foo" in resp.headers

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        run()

    assert_reset()


def test_callback_exception_result():
    result = Exception()
    url = "http://example.com/"

    def request_callback(_request):
        return result

    @responses.activate
    def run():
        responses.add_callback(responses.GET, url, request_callback)

        with pytest.raises(Exception) as e:
            requests.get(url)

        assert e.value is result

    run()
    assert_reset()


def test_callback_exception_body():
    body = Exception()
    url = "http://example.com/"

    def request_callback(_request):
        return 200, {}, body

    @responses.activate
    def run():
        responses.add_callback(responses.GET, url, request_callback)

        with pytest.raises(Exception) as e:
            requests.get(url)

        assert e.value is body

    run()
    assert_reset()


def test_callback_no_content_type():
    body = b"test callback"
    status = 400
    reason = "Bad Request"
    headers = {"foo": "bar"}
    url = "http://example.com/"

    def request_callback(_request):
        return status, headers, body

    @responses.activate
    def run():
        responses.add_callback(responses.GET, url, request_callback, content_type=None)
        resp = requests.get(url)
        assert resp.text == "test callback"
        assert resp.status_code == status
        assert resp.reason == reason
        assert "foo" in resp.headers
        assert "Content-Type" not in resp.headers

    run()
    assert_reset()


def test_callback_content_type_dict():
    def request_callback(_request):
        return (
            200,
            {"Content-Type": "application/json"},
            b"foo",
        )

    @responses.activate
    def run():
        responses.add_callback("GET", "http://mockhost/.foo", callback=request_callback)
        resp = requests.get("http://mockhost/.foo")
        assert resp.text == "foo"
        assert resp.headers["content-type"] == "application/json"

    run()
    assert_reset()


def test_callback_matchers():
    def request_callback(_request):
        return (
            200,
            {"Content-Type": "application/json"},
            b"foo",
        )

    @responses.activate
    def run():
        req_data = {"some": "other", "data": "fields"}
        req_files = {"file_name": b"Old World!"}

        responses.add_callback(
            responses.POST,
            url="http://httpbin.org/post",
            match=[matchers.multipart_matcher(req_files, data=req_data)],
            callback=request_callback,
        )
        resp = requests.post("http://httpbin.org/post", data=req_data, files=req_files)
        assert resp.text == "foo"
        assert resp.headers["content-type"] == "application/json"

    run()
    assert_reset()


def test_callback_matchers_fail():
    @responses.activate
    def run():
        req_data = {"some": "other", "data": "fields"}
        req_files = {"file_name": b"Old World!"}

        responses.add_callback(
            responses.POST,
            url="http://httpbin.org/post",
            match=[matchers.multipart_matcher(req_files, data=req_data)],
            callback=lambda x: (
                0,
                {"a": ""},
                "",
            ),
        )
        with pytest.raises(ConnectionError) as exc:
            requests.post(
                "http://httpbin.org/post",
                data={"some": "other", "data": "wrong"},
                files=req_files,
            )

        assert "multipart/form-data doesn't match." in str(exc.value)

    run()
    assert_reset()


def test_callback_content_type_tuple():
    def request_callback(_request):
        return (
            200,
            [("Content-Type", "application/json")],
            b"foo",
        )

    @responses.activate
    def run():
        responses.add_callback("GET", "http://mockhost/.foo", callback=request_callback)
        resp = requests.get("http://mockhost/.foo")
        assert resp.text == "foo"
        assert resp.headers["content-type"] == "application/json"

    run()
    assert_reset()


def test_regular_expression_url():
    @responses.activate
    def run():
        url = re.compile(r"https?://(.*\.)?example.com")
        responses.add(responses.GET, url, body=b"test")

        resp = requests.get("http://example.com")
        assert_response(resp, "test")

        resp = requests.get("https://example.com")
        assert_response(resp, "test")

        resp = requests.get("https://uk.example.com")
        assert_response(resp, "test")

        with pytest.raises(ConnectionError):
            requests.get("https://uk.exaaample.com")

    run()
    assert_reset()


def test_base_response_get_response():
    resp = BaseResponse("GET", ".com")
    with pytest.raises(NotImplementedError):
        resp.get_response(requests.PreparedRequest())


class TestAdapters:
    class CustomAdapter(requests.adapters.HTTPAdapter):
        """Classic custom adapter."""

        def send(self, *a, **k):
            return super().send(*a, **k)

    class PositionalArgsAdapter(requests.adapters.HTTPAdapter):
        """Custom adapter that sends only positional args.
        See https://github.com/getsentry/responses/issues/642 for more into.
        """

        def send(
            self,
            request,
            stream=False,
            timeout=None,
            verify=True,
            cert=None,
            proxies=None,
        ):
            return super().send(request, stream, timeout, verify, cert, proxies)

    class PositionalArgsIncompleteAdapter(requests.adapters.HTTPAdapter):
        """Custom adapter that sends only positional args.
        Not all arguments are forwarded to the send method.
                    See https://github.com/getsentry/responses/issues/642 for more into.
        """

        def send(
            self,
            request,
            stream=False,
            timeout=None,
            verify=True,
            # following args are intentionally not forwarded
            cert=None,
            proxies=None,
        ):
            return super().send(request, stream, timeout, verify)

    @pytest.mark.parametrize(
        "adapter_class",
        (CustomAdapter, PositionalArgsAdapter, PositionalArgsIncompleteAdapter),
    )
    def test_custom_adapter(self, adapter_class):  # type: ignore[misc]
        """Test basic adapter implementation and that responses can patch them properly."""

        @responses.activate
        def run():
            url = "http://example.com"
            responses.add(responses.GET, url, body=b"test adapter")

            # Test that the adapter is actually used
            session = requests.Session()
            adapter = adapter_class()
            session.mount("http://", adapter)
            with patch.object(adapter, "send", side_effect=adapter.send) as mock_send:
                resp = session.get(url, allow_redirects=False)

            assert mock_send.call_count == 1
            assert_response(resp, "test adapter")

        run()


def test_responses_as_context_manager():
    def run():
        with responses.mock:
            responses.add(responses.GET, "http://example.com", body=b"test")
            resp = requests.get("http://example.com")
            assert_response(resp, "test")
            assert len(responses.calls) == 1
            assert responses.calls[0].request.url == "http://example.com/"
            assert responses.calls[0].response.content == b"test"

            resp = requests.get("http://example.com?foo=bar")
            assert_response(resp, "test")
            assert len(responses.calls) == 2
            assert responses.calls[1].request.url == "http://example.com/?foo=bar"
            assert responses.calls[1].response.content == b"test"

    run()
    assert_reset()


def test_activate_doesnt_change_signature():
    def test_function(a, b=None):
        return a, b

    decorated_test_function = responses.activate(test_function)
    assert inspect.signature(test_function) == inspect.signature(
        decorated_test_function
    )

    assert decorated_test_function(1, 2) == test_function(1, 2)
    assert decorated_test_function(3) == test_function(3)


@pytest.fixture
def my_fruit():  # type: ignore[misc]
    return "apple"


@pytest.fixture
def fruit_basket(my_fruit):  # type: ignore[misc]
    return ["banana", my_fruit]


@pytest.mark.usefixtures("my_fruit", "fruit_basket")
class TestFixtures:
    """
    Test that pytest fixtures work well with 'activate' decorator
    """

    def test_function(self, my_fruit, fruit_basket):
        assert my_fruit in fruit_basket
        assert my_fruit == "apple"

    test_function_decorated = responses.activate(test_function)


def test_activate_mock_interaction():
    @patch("sys.stdout")
    def test_function(mock_stdout):  # type: ignore[misc]
        return mock_stdout

    decorated_test_function = responses.activate(test_function)
    assert inspect.signature(test_function) == inspect.signature(
        decorated_test_function
    )

    value = test_function()
    assert isinstance(value, Mock)

    value = decorated_test_function()
    assert isinstance(value, Mock)


def test_activate_doesnt_change_signature_with_return_type():
    def test_function(a, b=None):
        return a, b

    # Add type annotations as they are syntax errors in py2.
    # Use a class to test for import errors in evaled code.
    test_function.__annotations__["return"] = Mock
    test_function.__annotations__["a"] = Mock

    decorated_test_function = responses.activate(test_function)
    assert inspect.signature(test_function) == inspect.signature(
        decorated_test_function
    )

    assert decorated_test_function(1, 2) == test_function(1, 2)
    assert decorated_test_function(3) == test_function(3)


def test_activate_doesnt_change_signature_for_method():
    class TestCase:
        def test_function(self, a, b=None):
            return self, a, b

        decorated_test_function = responses.activate(test_function)

    test_case = TestCase()
    assert test_case.decorated_test_function(1, 2) == test_case.test_function(1, 2)
    assert test_case.decorated_test_function(3) == test_case.test_function(3)


def test_response_cookies():
    body = b"test callback"
    status = 200
    headers = {"set-cookie": "session_id=12345; a=b; c=d"}
    url = "http://example.com/"

    def request_callback(_request):
        return status, headers, body

    @responses.activate
    def run():
        responses.add_callback(responses.GET, url, request_callback)
        resp = requests.get(url)
        assert resp.text == "test callback"
        assert resp.status_code == status
        assert "session_id" in resp.cookies
        assert resp.cookies["session_id"] == "12345"
        assert set(resp.cookies.keys()) == {"session_id"}

    run()
    assert_reset()


def test_response_cookies_secure():
    body = b"test callback"
    status = 200
    headers = {"set-cookie": "session_id=12345; a=b; c=d; secure"}
    url = "http://example.com/"

    def request_callback(_request):
        return status, headers, body

    @responses.activate
    def run():
        responses.add_callback(responses.GET, url, request_callback)
        resp = requests.get(url)
        assert resp.text == "test callback"
        assert resp.status_code == status
        assert "session_id" in resp.cookies
        assert resp.cookies["session_id"] == "12345"
        assert set(resp.cookies.keys()) == {"session_id"}

    run()
    assert_reset()


def test_response_cookies_multiple():
    body = b"test callback"
    status = 200
    headers = [
        ("set-cookie", "1P_JAR=2019-12-31-23; path=/; domain=.example.com; HttpOnly"),
        ("set-cookie", "NID=some=value; path=/; domain=.example.com; secure"),
    ]
    url = "http://example.com/"

    def request_callback(_request):
        return status, headers, body

    @responses.activate
    def run():
        responses.add_callback(responses.GET, url, request_callback)
        resp = requests.get(url)
        assert resp.text == "test callback"
        assert resp.status_code == status
        assert set(resp.cookies.keys()) == {"1P_JAR", "NID"}
        assert resp.cookies["1P_JAR"] == "2019-12-31-23"
        assert resp.cookies["NID"] == "some=value"

    run()
    assert_reset()


@pytest.mark.parametrize("request_stream", (True, False, None))
@pytest.mark.parametrize("responses_stream", (True, False, None))
def test_response_cookies_session(request_stream, responses_stream):  # type: ignore[misc]
    @responses.activate
    def run():
        url = "https://example.com/path"
        responses.add(
            responses.GET,
            url,
            headers=[
                ("Set-cookie", "mycookie=cookieval; path=/; secure"),
            ],
            body="ok",
            stream=responses_stream,
        )
        session = requests.session()
        resp = session.get(url, stream=request_stream)
        assert resp.text == "ok"
        assert resp.status_code == 200

        assert "mycookie" in resp.cookies
        assert resp.cookies["mycookie"] == "cookieval"
        assert set(resp.cookies.keys()) == {"mycookie"}

        assert "mycookie" in session.cookies
        assert session.cookies["mycookie"] == "cookieval"
        assert set(session.cookies.keys()) == {"mycookie"}

    run()
    assert_reset()


def test_response_callback():
    """adds a callback to decorate the response, then checks it"""

    def run():
        def response_callback(response):
            response._is_mocked = True
            return response

        with responses.RequestsMock(response_callback=response_callback) as m:
            m.add(responses.GET, "http://example.com", body=b"test")
            resp = requests.get("http://example.com")
            assert resp.text == "test"
            assert hasattr(resp, "_is_mocked")
            assert getattr(resp, "_is_mocked") is True

    run()
    assert_reset()


def test_response_filebody():
    """Adds the possibility to use actual (binary) files as responses"""

    def run():
        current_file = os.path.abspath(__file__)
        with responses.RequestsMock() as m:
            with open(current_file, encoding="utf-8") as out:
                m.add(responses.GET, "http://example.com", body=out.read(), stream=True)
                resp = requests.get("http://example.com", stream=True)
            with open(current_file, encoding="utf-8") as out:
                assert resp.text == out.read()

    run()
    assert_reset()


def test_use_stream_twice_to_double_raw_io():
    @responses.activate
    def run():
        url = "http://example.com"
        responses.add(responses.GET, url, body=b"42", stream=True)
        resp = requests.get(url, stream=True)
        assert resp.raw.read() == b"42"

    run()
    assert_reset()


def test_assert_all_requests_are_fired():
    def request_callback(_request):
        raise BaseException()

    def run():
        with pytest.raises(AssertionError) as excinfo:
            with responses.RequestsMock(assert_all_requests_are_fired=True) as m:
                m.add(responses.GET, "http://example.com", body=b"test")
        assert "http://example.com" in str(excinfo.value)
        assert responses.GET in str(excinfo.value)

        # check that assert_all_requests_are_fired default to True
        with pytest.raises(AssertionError):
            with responses.RequestsMock() as m:
                m.add(responses.GET, "http://example.com", body=b"test")

        # check that assert_all_requests_are_fired doesn't swallow exceptions
        with pytest.raises(ValueError):
            with responses.RequestsMock() as m:
                m.add(responses.GET, "http://example.com", body=b"test")
                raise ValueError()

        # check that assert_all_requests_are_fired=True doesn't remove urls
        with responses.RequestsMock(assert_all_requests_are_fired=True) as m:
            m.add(responses.GET, "http://example.com", body=b"test")
            assert len(m.registered()) == 1
            requests.get("http://example.com")
            assert len(m.registered()) == 1

        # check that assert_all_requests_are_fired=True counts mocked errors
        with responses.RequestsMock(assert_all_requests_are_fired=True) as m:
            m.add(responses.GET, "http://example.com", body=Exception())
            assert len(m.registered()) == 1
            with pytest.raises(Exception):
                requests.get("http://example.com")
            assert len(m.registered()) == 1

        with responses.RequestsMock(assert_all_requests_are_fired=True) as m:
            m.add_callback(responses.GET, "http://example.com", request_callback)
            assert len(m.registered()) == 1
            with pytest.raises(BaseException):
                requests.get("http://example.com")
            assert len(m.registered()) == 1

    run()
    assert_reset()


def test_assert_all_requests_fired_multiple():
    @responses.activate(assert_all_requests_are_fired=True)
    def test_some_function():
        # Not all mocks are called so we'll get an AssertionError
        responses.add(responses.GET, "http://other_url", json={})
        responses.add(responses.GET, "http://some_api", json={})
        requests.get("http://some_api")

    @responses.activate(assert_all_requests_are_fired=True)
    def test_some_second_function():
        # This should pass as mocks should be reset.
        responses.add(responses.GET, "http://some_api", json={})
        requests.get("http://some_api")

    with pytest.raises(AssertionError):
        test_some_function()
    assert_reset()

    test_some_second_function()
    assert_reset()


def test_allow_redirects_samehost():
    redirecting_url = "http://example.com"
    final_url_path = "/1"
    final_url = f"{redirecting_url}{final_url_path}"
    url_re = re.compile(r"^http://example.com(/)?(\d+)?$")

    def request_callback(request):
        # endpoint of chained redirect
        if request.url.endswith(final_url_path):
            return 200, (), b"test"

        # otherwise redirect to an integer path
        else:
            if request.url.endswith("/0"):
                n = 1
            else:
                n = 0
            redirect_headers = {"location": f"/{n!s}"}
            return 301, redirect_headers, None

    def run():
        # setup redirect
        with responses.mock:
            responses.add_callback(responses.GET, url_re, request_callback)
            resp_no_redirects = requests.get(redirecting_url, allow_redirects=False)
            assert resp_no_redirects.status_code == 301
            assert len(responses.calls) == 1  # 1x300
            assert responses.calls[0][1].status_code == 301
        assert_reset()

        with responses.mock:
            responses.add_callback(responses.GET, url_re, request_callback)
            resp_yes_redirects = requests.get(redirecting_url, allow_redirects=True)
            assert len(responses.calls) == 3  # 2x300 + 1x200
            assert len(resp_yes_redirects.history) == 2
            assert resp_yes_redirects.status_code == 200
            assert final_url == resp_yes_redirects.url
            status_codes = [call[1].status_code for call in responses.calls]
            assert status_codes == [301, 301, 200]
        assert_reset()

    run()
    assert_reset()


def test_path_segments():
    """Test that path segment after ``;`` is preserved.

    Validate compliance with RFC 3986.
    The path is terminated by the first question mark ("?") or
    number sign ("#") character, or by the end of the URI.
    See more about how path should be treated under:
    https://datatracker.ietf.org/doc/html/rfc3986.html#section-3.3
    """

    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com/here/we", status=669)
        responses.add(responses.GET, "http://example.com/here/we;go", status=777)

        resp = requests.get("http://example.com/here/we;go")
        assert resp.status_code == 777

    run()
    assert_reset()


def test_handles_unicode_querystring():
    url = "http://example.com/test?type=2&ie=utf8&query=汉字"

    @responses.activate
    def run():
        responses.add(responses.GET, url, body="test", match_querystring=True)

        resp = requests.get(url)

        assert_response(resp, "test")

    run()
    assert_reset()


def test_handles_unicode_url():
    url = "http://www.संजाल.भारत/hi/वेबसाइट-डिजाइन"

    @responses.activate
    def run():
        responses.add(responses.GET, url, body="test")

        resp = requests.get(url)

        assert_response(resp, "test")

    run()
    assert_reset()


def test_handles_unicode_body():
    url = "http://example.com/test"

    @responses.activate
    def run():
        responses.add(responses.GET, url, body="михољско лето")

        resp = requests.get(url)

        assert_response(resp, "михољско лето", content_type="text/plain; charset=utf-8")

    run()
    assert_reset()


def test_handles_buffered_reader_body():
    url = "http://example.com/test"

    @responses.activate
    def run():
        responses.add(responses.GET, url, body=BufferedReader(BytesIO(b"test")))  # type: ignore

        resp = requests.get(url)

        assert_response(resp, "test")

    run()
    assert_reset()


def test_headers():
    @responses.activate
    def run():
        responses.add(
            responses.GET, "http://example.com", body="", headers={"X-Test": "foo"}
        )
        resp = requests.get("http://example.com")
        assert resp.headers["X-Test"] == "foo"

    run()
    assert_reset()


def test_headers_deduplicated_content_type():
    """Test to ensure that we do not have two values for `content-type`.

    For more details see https://github.com/getsentry/responses/issues/644
    """

    @responses.activate
    def run():
        responses.get(
            "https://example.org/",
            json={},
            headers={"Content-Type": "application/json"},
        )
        responses.start()

        resp = requests.get("https://example.org/")

        assert resp.headers["Content-Type"] == "application/json"

    run()
    assert_reset()


def test_content_length_error(monkeypatch):
    """
    Currently 'requests' does not enforce content length validation,
    (validation that body length matches header). However, this could
    be expected in next major version, see
    https://github.com/psf/requests/pull/3563

    Now user can manually patch URL3 lib to achieve the same

    See discussion in
    https://github.com/getsentry/responses/issues/394
    """

    @responses.activate
    def run():
        responses.add(
            responses.GET,
            "http://example.com/api/123",
            json={"message": "this body is too large"},
            adding_headers={"content-length": "2"},
        )
        with pytest.raises(ChunkedEncodingError) as exc:
            requests.get("http://example.com/api/123")

        assert "IncompleteRead" in str(exc.value)

    # Type errors here and on 1250 are ignored because the stubs for requests
    # are off https://github.com/python/typeshed/blob/f8501d33c737482a829c6db557a0be26895c5941
    #   /stubs/requests/requests/packages/__init__.pyi#L1
    original_init = getattr(urllib3.HTTPResponse, "__init__")

    def patched_init(self, *args, **kwargs):
        kwargs["enforce_content_length"] = True
        original_init(self, *args, **kwargs)

    monkeypatch.setattr(urllib3.HTTPResponse, "__init__", patched_init)

    run()
    assert_reset()


def test_stream_with_none_chunk_size():
    """
    See discussion in
    https://github.com/getsentry/responses/issues/438
    """

    @responses.activate
    def run():
        responses.add(
            responses.GET,
            "https://example.com",
            status=200,
            content_type="application/octet-stream",
            body=b"This is test",
            auto_calculate_content_length=True,
        )
        res = requests.get("https://example.com", stream=True)
        for chunk in res.iter_content(chunk_size=None):
            assert chunk == b"This is test"

    run()
    assert_reset()


def test_legacy_adding_headers():
    @responses.activate
    def run():
        responses.add(
            responses.GET,
            "http://example.com",
            body="",
            adding_headers={"X-Test": "foo"},
        )
        resp = requests.get("http://example.com")
        assert resp.headers["X-Test"] == "foo"

    run()
    assert_reset()


def test_legacy_adding_headers_with_content_type():
    @responses.activate
    def run():
        with pytest.raises(RuntimeError) as excinfo:
            responses.add(
                responses.GET,
                "http://example.com",
                body="test",
                content_type="text/html",
                adding_headers={"Content-Type": "text/html; charset=utf-8"},
            )
        assert (
            "You cannot define both `content_type` and `headers[Content-Type]`"
            in str(excinfo.value)
        )

    run()
    assert_reset()


def test_auto_calculate_content_length_string_body():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(
            responses.GET, url, body="test", auto_calculate_content_length=True
        )
        resp = requests.get(url)
        assert_response(resp, "test")
        assert resp.headers["Content-Length"] == "4"

    run()
    assert_reset()


def test_auto_calculate_content_length_bytes_body():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(
            responses.GET, url, body=b"test bytes", auto_calculate_content_length=True
        )
        resp = requests.get(url)
        assert_response(resp, "test bytes")
        assert resp.headers["Content-Length"] == "10"

    run()
    assert_reset()


def test_auto_calculate_content_length_json_body():
    @responses.activate
    def run():
        content_type = "application/json"

        url = "http://example.com/"
        responses.add(
            responses.GET,
            url,
            json={"message": "success"},
            auto_calculate_content_length=True,
        )
        resp = requests.get(url)
        assert_response(resp, '{"message": "success"}', content_type)
        assert resp.headers["Content-Length"] == "22"

        url = "http://example.com/1/"
        responses.add(responses.GET, url, json=[], auto_calculate_content_length=True)
        resp = requests.get(url)
        assert_response(resp, "[]", content_type)
        assert resp.headers["Content-Length"] == "2"

    run()
    assert_reset()


def test_auto_calculate_content_length_unicode_body():
    @responses.activate
    def run():
        url = "http://example.com/test"
        responses.add(
            responses.GET, url, body="михољско лето", auto_calculate_content_length=True
        )
        resp = requests.get(url)
        assert_response(resp, "михољско лето", content_type="text/plain; charset=utf-8")
        assert resp.headers["Content-Length"] == "25"

    run()
    assert_reset()


def test_auto_calculate_content_length_doesnt_work_for_buffered_reader_body():
    @responses.activate
    def run():
        url = "http://example.com/test"
        responses.add(
            responses.GET,
            url,
            body=BufferedReader(BytesIO(b"testing")),  # type: ignore
            auto_calculate_content_length=True,
        )
        resp = requests.get(url)
        assert_response(resp, "testing")
        assert "Content-Length" not in resp.headers

    run()
    assert_reset()


def test_auto_calculate_content_length_doesnt_override_existing_value():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(
            responses.GET,
            url,
            body="test",
            headers={"Content-Length": "2"},
            auto_calculate_content_length=True,
        )

        if urllib3.__version__ < "2":
            resp = requests.get(url)
            assert_response(resp, "test")
            assert resp.headers["Content-Length"] == "2"
        else:
            with pytest.raises(ChunkedEncodingError) as excinfo:
                requests.get(url)
            assert "IncompleteRead(4 bytes read, -2 more expected)" in str(
                excinfo.value
            )

    run()
    assert_reset()


def test_multiple_responses():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com", body="test")
        responses.add(responses.GET, "http://example.com", body="rest")
        responses.add(responses.GET, "http://example.com", body="fest")
        responses.add(responses.GET, "http://example.com", body="best")

        resp = requests.get("http://example.com")
        assert_response(resp, "test")

        resp = requests.get("http://example.com")
        assert_response(resp, "rest")

        resp = requests.get("http://example.com")
        assert_response(resp, "fest")

        resp = requests.get("http://example.com")
        assert_response(resp, "best")

        # After all responses are used, last response should be repeated
        resp = requests.get("http://example.com")
        assert_response(resp, "best")

    run()
    assert_reset()


def test_multiple_responses_intermixed():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com", body="test")
        resp = requests.get("http://example.com")
        assert_response(resp, "test")

        responses.add(responses.GET, "http://example.com", body="rest")
        resp = requests.get("http://example.com")
        assert_response(resp, "rest")

        responses.add(responses.GET, "http://example.com", body="best")
        resp = requests.get("http://example.com")
        assert_response(resp, "best")

        # After all responses are used, last response should be repeated
        resp = requests.get("http://example.com")
        assert_response(resp, "best")

    run()
    assert_reset()


def test_multiple_urls():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com/one", body="one")
        responses.add(responses.GET, "http://example.com/two", body="two")

        resp = requests.get("http://example.com/two")
        assert_response(resp, "two")
        resp = requests.get("http://example.com/one")
        assert_response(resp, "one")

    run()
    assert_reset()


def test_multiple_methods():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://example.com/one", body="gotcha")
        responses.add(responses.POST, "http://example.com/one", body="posted")

        resp = requests.get("http://example.com/one")
        assert_response(resp, "gotcha")
        resp = requests.post("http://example.com/one")
        assert_response(resp, "posted")

    run()
    assert_reset()


class TestPassthru:
    def test_passthrough_flag(self, httpserver):
        httpserver.expect_request("/").respond_with_data(
            "OK", content_type="text/plain"
        )
        url = httpserver.url_for("/")

        response = Response(responses.GET, url, body="MOCK")

        @responses.activate
        def run_passthrough():
            responses.add(response)
            resp = requests.get(url)
            assert_response(resp, "OK")

        @responses.activate
        def run_mocked():
            responses.add(response)
            resp = requests.get(url)
            assert_response(resp, "MOCK")

        run_mocked()
        assert_reset()

        response.passthrough = True
        run_passthrough()
        assert_reset()

    def test_passthrough_kwarg(self, httpserver):
        httpserver.expect_request("/").respond_with_data(
            "OK", content_type="text/plain"
        )
        url = httpserver.url_for("/")

        def configure_response(passthrough):
            responses.get(url, body="MOCK", passthrough=passthrough)

        @responses.activate
        def run_passthrough():
            configure_response(passthrough=True)
            resp = requests.get(url)
            assert_response(resp, "OK")

        @responses.activate
        def run_mocked():
            configure_response(passthrough=False)
            resp = requests.get(url)
            assert_response(resp, "MOCK")

        run_mocked()
        assert_reset()

        run_passthrough()
        assert_reset()

    def test_passthrough_response(self, httpserver):
        httpserver.expect_request("/").respond_with_data(
            "OK", content_type="text/plain"
        )
        url = httpserver.url_for("/")

        @responses.activate
        def run():
            responses.add(PassthroughResponse(responses.GET, url))
            responses.add(responses.GET, f"{url}/one", body="one")
            responses.add(responses.GET, "http://example.com/two", body="two")

            resp = requests.get("http://example.com/two")
            assert_response(resp, "two")
            resp = requests.get(f"{url}/one")
            assert_response(resp, "one")
            resp = requests.get(url)
            assert_response(resp, "OK")

            assert len(responses.calls) == 3
            responses.assert_call_count(url, 1)

        run()
        assert_reset()

    def test_passthrough_response_stream(self, httpserver):
        httpserver.expect_request("/").respond_with_data(
            "OK", content_type="text/plain"
        )

        @responses.activate
        def run():
            url = httpserver.url_for("/")
            responses.add(PassthroughResponse(responses.GET, url))
            content_1 = requests.get(url).content
            with requests.get(url, stream=True) as resp:
                content_2 = resp.raw.read()
            assert content_1 == content_2

        run()
        assert_reset()

    def test_passthru_prefixes(self, httpserver):
        httpserver.expect_request("/").respond_with_data(
            "OK", content_type="text/plain"
        )
        url = httpserver.url_for("/")

        @responses.activate
        def run_constructor_argument():
            with responses.RequestsMock(passthru_prefixes=(url,)):
                resp = requests.get(url)
                assert_response(resp, "OK")

        @responses.activate
        def run_property_setter():
            with responses.RequestsMock() as m:
                m.passthru_prefixes = tuple([url])
                resp = requests.get(url)
                assert_response(resp, "OK")

        run_constructor_argument()
        assert_reset()
        run_property_setter()
        assert_reset()

    def test_passthru(self, httpserver):
        httpserver.expect_request("/").respond_with_data(
            "OK", content_type="text/plain"
        )
        url = httpserver.url_for("/")

        @responses.activate
        def run():
            responses.add_passthru(url)
            responses.add(responses.GET, f"{url}/one", body="one")
            responses.add(responses.GET, "http://example.com/two", body="two")

            resp = requests.get("http://example.com/two")
            assert_response(resp, "two")
            resp = requests.get(f"{url}/one")
            assert_response(resp, "one")
            resp = requests.get(url)
            assert_response(resp, "OK")

        run()
        assert_reset()

    def test_passthru_regex(self, httpserver):
        httpserver.expect_request(re.compile("^/\\w+")).respond_with_data(
            "OK", content_type="text/plain"
        )
        url = httpserver.url_for("/")

        @responses.activate
        def run():
            responses.add_passthru(re.compile(f"{url}/\\w+"))
            responses.add(responses.GET, f"{url}/one", body="one")
            responses.add(responses.GET, "http://example.com/two", body="two")

            resp = requests.get("http://example.com/two")
            assert_response(resp, "two")
            resp = requests.get(f"{url}/one")
            assert_response(resp, "one")
            resp = requests.get(f"{url}/two")
            assert_response(resp, "OK")
            resp = requests.get(f"{url}/three")
            assert_response(resp, "OK")

        run()
        assert_reset()

    def test_passthru_does_not_persist_across_tests(self, httpserver):
        """
        passthru should be erased on exit from context manager
        see:
        https://github.com/getsentry/responses/issues/322
        """
        httpserver.expect_request("/").respond_with_data(
            "mocked server", status=969, content_type="text/plain"
        )

        @responses.activate
        def with_a_passthru():
            assert not responses.mock.passthru_prefixes
            responses.add_passthru(re.compile(".*"))

            url = httpserver.url_for("/")
            response = requests.get(url)
            assert response.status_code == 969
            assert response.text == "mocked server"

        @responses.activate
        def without_a_passthru():
            assert not responses.mock.passthru_prefixes
            with pytest.raises(requests.exceptions.ConnectionError):
                requests.get("https://example.com")

        with_a_passthru()
        without_a_passthru()

    def test_passthru_unicode(self):
        @responses.activate
        def run():
            with responses.RequestsMock() as m:
                url = "http://موقع.وزارة-الاتصالات.مصر/"
                clean_url = "http://xn--4gbrim.xn----ymcbaaajlc6dj7bxne2c.xn--wgbh1c/"
                m.add_passthru(url)
                assert m.passthru_prefixes[0] == clean_url

        run()
        assert_reset()

    def test_real_send_argument(self):
        def run():
            # the following mock will serve to catch the real send request from another mock and
            # will "donate" `unbound_on_send` method
            mock_to_catch_real_send = responses.RequestsMock(
                assert_all_requests_are_fired=True
            )
            mock_to_catch_real_send.post(
                "http://send-this-request-through.com", status=500
            )

            with responses.RequestsMock(
                assert_all_requests_are_fired=True,
                real_adapter_send=mock_to_catch_real_send.unbound_on_send(),
            ) as r_mock:
                r_mock.add_passthru("http://send-this-request-through.com")

                r_mock.add(responses.POST, "https://example.org", status=200)

                response = requests.post("https://example.org")
                assert response.status_code == 200

                response = requests.post("http://send-this-request-through.com")
                assert response.status_code == 500

        run()
        assert_reset()


def test_method_named_param():
    @responses.activate
    def run():
        responses.add(method=responses.GET, url="http://example.com", body="OK")
        resp = requests.get("http://example.com")
        assert_response(resp, "OK")

    run()
    assert_reset()


def test_custom_target(monkeypatch):
    requests_mock = responses.RequestsMock(target="something.else")
    std_mock_mock = responses.std_mock.MagicMock()
    patch_mock = std_mock_mock.patch
    monkeypatch.setattr(responses, "std_mock", std_mock_mock)
    requests_mock.start()
    assert len(patch_mock.call_args_list) == 1
    assert patch_mock.call_args[1]["target"] == "something.else"


@pytest.mark.parametrize(
    "url",
    (
        "http://example.com",
        "http://example.com/some/path",
        "http://example.com/other/path/",
    ),
)
def test_request_param(url):  # type: ignore[misc]
    @responses.activate
    def run():
        params = {"hello": "world", "example": "params"}
        responses.add(
            method=responses.GET,
            url=f"{url}?hello=world",
            body="test",
            match_querystring=False,
        )
        resp = requests.get(url, params=params)
        assert_response(resp, "test")
        assert_params(resp, params)

        resp = requests.get(url)
        assert_response(resp, "test")
        assert_params(resp, {})

    run()
    assert_reset()


def test_request_param_with_multiple_values_for_the_same_key():
    @responses.activate
    def run():
        url = "http://example.com"
        params = {"key1": ["one", "two"], "key2": "three"}
        responses.add(
            method=responses.GET,
            url=url,
            body="test",
        )
        resp = requests.get(url, params=params)
        assert_response(resp, "test")
        assert_params(resp, params)

    run()
    assert_reset()


@pytest.mark.parametrize(
    "url", ("http://example.com", "http://example.com?hello=world")
)
def test_assert_call_count(url):  # type: ignore[misc]
    @responses.activate
    def run():
        responses.add(responses.GET, url)
        responses.add(responses.GET, "http://example1.com")

        assert responses.assert_call_count(url, 0) is True

        with pytest.raises(AssertionError) as excinfo:
            responses.assert_call_count(url, 2)
        assert "Expected URL '{}' to be called 2 times. Called 0 times.".format(
            url
        ) in str(excinfo.value)

        requests.get(url)
        assert responses.assert_call_count(url, 1) is True

        requests.get("http://example1.com")
        assert responses.assert_call_count(url, 1) is True

        requests.get(url)
        with pytest.raises(AssertionError) as excinfo:
            responses.assert_call_count(url, 3)
        assert "Expected URL '{}' to be called 3 times. Called 2 times.".format(
            url
        ) in str(excinfo.value)

    run()
    assert_reset()


def test_call_count_with_matcher():
    @responses.activate
    def run():
        rsp = responses.add(
            responses.GET,
            "http://www.example.com",
            match=(matchers.query_param_matcher({}),),
        )
        rsp2 = responses.add(
            responses.GET,
            "http://www.example.com",
            match=(matchers.query_param_matcher({"hello": "world"}),),
            status=777,
        )
        requests.get("http://www.example.com")
        resp1 = requests.get("http://www.example.com")
        requests.get("http://www.example.com?hello=world")
        resp2 = requests.get("http://www.example.com?hello=world")

        assert resp1.status_code == 200
        assert resp2.status_code == 777

        assert rsp.call_count == 2
        assert rsp2.call_count == 2

    run()
    assert_reset()


def test_call_count_without_matcher():
    @responses.activate
    def run():
        rsp = responses.add(responses.GET, "http://www.example.com")
        requests.get("http://www.example.com")
        requests.get("http://www.example.com")
        requests.get("http://www.example.com?hello=world")
        requests.get("http://www.example.com?hello=world")

        assert rsp.call_count == 4

    run()
    assert_reset()


def test_response_calls_indexing_and_slicing():
    @responses.activate
    def run():
        responses.add(responses.GET, "http://www.example.com")
        responses.add(responses.GET, "http://www.example.com/1")
        responses.add(responses.GET, "http://www.example.com/2")

        requests.get("http://www.example.com")
        requests.get("http://www.example.com/1")
        requests.get("http://www.example.com/2")
        requests.get("http://www.example.com/1")
        requests.get("http://www.example.com")

        # Use of a type hints here ensures mypy knows the difference between index and slice.
        individual_call: Call = responses.calls[0]
        call_slice: List[Call] = responses.calls[1:-1]

        assert individual_call.request.url == "http://www.example.com/"

        assert call_slice == [
            responses.calls[1],
            responses.calls[2],
            responses.calls[3],
        ]
        assert [c.request.url for c in call_slice] == [
            "http://www.example.com/1",
            "http://www.example.com/2",
            "http://www.example.com/1",
        ]

    run()
    assert_reset()


def test_response_calls_and_registry_calls_are_equal():
    @responses.activate
    def run():
        rsp1 = responses.add(responses.GET, "http://www.example.com")
        rsp2 = responses.add(responses.GET, "http://www.example.com/1")
        rsp3 = responses.add(
            responses.GET, "http://www.example.com/2"
        )  # won't be requested

        requests.get("http://www.example.com")
        requests.get("http://www.example.com/1")
        requests.get("http://www.example.com")

        assert len(responses.calls) == len(rsp1.calls) + len(rsp2.calls) + len(
            rsp3.calls
        )
        assert rsp1.call_count == 2
        assert len(rsp1.calls) == 2
        assert rsp1.calls[0] is responses.calls[0]
        assert rsp1.calls[1] is responses.calls[2]
        assert rsp2.call_count == 1
        assert len(rsp2.calls) == 1
        assert rsp2.calls[0] is responses.calls[1]
        assert rsp3.call_count == 0
        assert len(rsp3.calls) == 0

    run()
    assert_reset()


def test_fail_request_error():
    """
    Validate that exception is raised if request URL/Method/kwargs don't match
    :return:
    """

    def run():
        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            rsps.add("POST", "http://example1.com")
            rsps.add("GET", "http://example.com")
            rsps.add_passthru("http://other.example.com")

            with pytest.raises(ConnectionError) as excinfo:
                requests.post("http://example.com", data={"id": "bad"})

            msg = str(excinfo.value)
            assert "- POST http://example1.com/ URL does not match" in msg
            assert "- GET http://example.com/ Method does not match" in msg
            assert "Passthru prefixes:\n- http://other.example.com" in msg

    run()
    assert_reset()


@pytest.mark.parametrize(
    "response_params, expected_representation",
    [
        (
            {"method": responses.GET, "url": "http://example.com/"},
            (
                "<Response(url='http://example.com/' status=200 "
                "content_type='text/plain' headers='null')>"
            ),
        ),
        (
            {
                "method": responses.POST,
                "url": "http://another-domain.com/",
                "content_type": "application/json",
                "status": 404,
            },
            (
                "<Response(url='http://another-domain.com/' status=404 "
                "content_type='application/json' headers='null')>"
            ),
        ),
        (
            {
                "method": responses.PUT,
                "url": "http://abcd.com/",
                "content_type": "text/html",
                "status": 500,
                "headers": {"X-Test": "foo"},
                "body": {"it_wont_be": "considered"},
            },
            (
                "<Response(url='http://abcd.com/' status=500 "
                "content_type='text/html' headers='{\"X-Test\": \"foo\"}')>"
            ),
        ),
    ],
)
def test_response_representations(response_params, expected_representation):  # type: ignore[misc]
    response = Response(**response_params)

    assert str(response) == expected_representation
    assert repr(response) == expected_representation


def test_mocked_responses_list_registered():
    @responses.activate
    def run():
        first_response = Response(
            responses.GET,
            "http://example.com/",
            body="",
            headers={"X-Test": "foo"},
            status=404,
        )
        second_response = Response(
            responses.GET, "http://example.com/", body="", headers={"X-Test": "foo"}
        )
        third_response = Response(
            responses.POST,
            "http://anotherdomain.com/",
        )
        responses.add(first_response)
        responses.add(second_response)
        responses.add(third_response)

        mocks_list = responses.registered()

        assert mocks_list == responses.mock.registered()
        assert mocks_list == [first_response, second_response, third_response]

    run()
    assert_reset()


@pytest.mark.parametrize(
    "url,other_url",
    [
        ("http://service-A/foo?q=fizz", "http://service-a/foo?q=fizz"),
        ("http://service-a/foo", "http://service-A/foo"),
        ("http://someHost-AwAy/", "http://somehost-away/"),
        ("http://fizzbuzz/foo", "http://fizzbuzz/foo"),
    ],
)
def test_rfc_compliance(url, other_url):  # type: ignore[misc]
    @responses.activate
    def run():
        responses.add(method=responses.GET, url=url)
        resp = requests.request("GET", other_url)
        assert_response(resp, "")

    run()
    assert_reset()


def test_requests_between_add():
    @responses.activate
    def run():
        responses.add(responses.GET, "https://example.com/", json={"response": "old"})
        assert requests.get("https://example.com/").content == b'{"response": "old"}'
        assert requests.get("https://example.com/").content == b'{"response": "old"}'
        assert requests.get("https://example.com/").content == b'{"response": "old"}'

        responses.add(responses.GET, "https://example.com/", json={"response": "new"})

        assert requests.get("https://example.com/").content == b'{"response": "new"}'
        assert requests.get("https://example.com/").content == b'{"response": "new"}'
        assert requests.get("https://example.com/").content == b'{"response": "new"}'

    run()
    assert_reset()


def test_responses_reuse():
    @responses.activate
    def run():
        url = "https://someapi.com/"
        fail_response = responses.Response(
            method="GET", url=url, body="fail", status=500
        )
        responses.add(responses.GET, url, "success", status=200)
        responses.add(fail_response)
        responses.add(fail_response)
        responses.add(fail_response)
        responses.add(responses.GET, url, "success", status=200)
        responses.add(responses.GET, url, "", status=302)

        response = requests.get(url)
        assert response.content == b"success"

        for _ in range(3):
            response = requests.get(url)
            assert response.content == b"fail"

    run()
    assert_reset()


@pytest.mark.asyncio
async def test_async_calls():  # type: ignore[misc]
    @responses.activate
    async def run():
        responses.add(
            responses.GET,
            "http://twitter.com/api/1/foobar",
            json={"error": "not found"},
            status=404,
        )

        resp = requests.get("http://twitter.com/api/1/foobar")
        assert resp.json() == {"error": "not found"}
        assert responses.calls[0].request.url == "http://twitter.com/api/1/foobar"

    await run()
    assert_reset()


class TestStrictWrapper:
    def test_strict_wrapper(self):
        """Test that assert_all_requests_are_fired could be applied to the decorator."""

        @responses.activate(assert_all_requests_are_fired=True)
        def run_strict():
            responses.add(responses.GET, "https://someapi1.com/", "success")
            responses.add(responses.GET, "https://notcalled1.com/", "success")
            requests.get("https://someapi1.com/")
            assert responses.mock.assert_all_requests_are_fired

        @responses.activate(assert_all_requests_are_fired=False)
        def run_not_strict():
            responses.add(responses.GET, "https://someapi2.com/", "success")
            responses.add(responses.GET, "https://notcalled2.com/", "success")
            requests.get("https://someapi2.com/")
            assert not responses.mock.assert_all_requests_are_fired

        @responses.activate
        def run_classic():
            responses.add(responses.GET, "https://someapi3.com/", "success")
            responses.add(responses.GET, "https://notcalled3.com/", "success")
            requests.get("https://someapi3.com/")
            assert not responses.mock.assert_all_requests_are_fired

        # keep the order of function calls to ensure that decorator doesn't leak to another function
        with pytest.raises(AssertionError) as exc_info:
            run_strict()

        # check that one URL is in uncalled assertion
        assert "https://notcalled1.com/" in str(exc_info.value)

        run_classic()
        run_not_strict()

    @pytest.mark.parametrize("assert_fired", (True, False, None))
    def test_nested_decorators(self, assert_fired):  # type: ignore[misc]
        """Validate if assert_all_requests_are_fired is applied from the correct function.

        assert_all_requests_are_fired must be applied from the function
        where call to 'requests' is done.
        Test matrix of True/False/None values applied to validate different scenarios.
        """

        @responses.activate(assert_all_requests_are_fired=assert_fired)
        def wrapped():
            responses.add(responses.GET, "https://notcalled1.com/", "success")
            responses.add(responses.GET, "http://example.com/1", body="Hello 1")
            assert b"Hello 1" == requests.get("http://example.com/1").content

        @responses.activate(assert_all_requests_are_fired=not assert_fired)
        def call_another_wrapped_function():
            responses.add(responses.GET, "https://notcalled2.com/", "success")
            wrapped()

        if assert_fired:
            with pytest.raises(AssertionError) as exc_info:
                call_another_wrapped_function()

            assert "https://notcalled1.com/" in str(exc_info.value)
            assert "https://notcalled2.com/" in str(exc_info.value)
        else:
            call_another_wrapped_function()


class TestMultipleWrappers:
    """Test to validate that multiple decorators could be applied.

    Ensures that we can call one function that is wrapped with
    ``responses.activate`` decorator from within another wrapped function.

    Validates that mock patch is not leaked to other tests.
    For more detail refer to https://github.com/getsentry/responses/issues/481
    """

    @responses.activate
    def test_wrapped(self):
        responses.add(responses.GET, "http://example.com/1", body="Hello 1")
        assert b"Hello 1" == requests.get("http://example.com/1").content

    @responses.activate
    def test_call_another_wrapped_function(self):
        self.test_wrapped()

    def test_mock_not_leaked(self, httpserver):
        """
        Validate that ``responses.activate`` does not leak to unpatched test.

        Parameters
        ----------
        httpserver : ContentServer
            Mock real HTTP server

        """
        httpserver.expect_request("/").respond_with_data(
            "OK", content_type="text/plain", status=969
        )
        url = httpserver.url_for("/")

        response = requests.get(url)
        assert response.status_code == 969


class TestShortcuts:
    def test_delete(self):
        @responses.activate
        def run():
            responses.delete("http://example.com/1", status=888)
            resp = requests.delete("http://example.com/1")
            assert resp.status_code == 888

        run()
        assert_reset()

    def test_get(self):
        @responses.activate
        def run():
            responses.get("http://example.com/1", status=888)
            resp = requests.get("http://example.com/1")
            assert resp.status_code == 888

        run()
        assert_reset()

    def test_head(self):
        @responses.activate
        def run():
            responses.head("http://example.com/1", status=888)
            resp = requests.head("http://example.com/1")
            assert resp.status_code == 888

        run()
        assert_reset()

    def test_head_with_content_length(self):
        @responses.activate
        def run():
            headers = {"content-length": "1000"}
            responses.head("http://example.com/1", status=200, headers=headers)
            resp = requests.head("http://example.com/1")
            assert resp.status_code == 200
            assert resp.headers["Content-Length"] == "1000"

        run()
        assert_reset()

    def test_options(self):
        @responses.activate
        def run():
            responses.options("http://example.com/1", status=888)
            resp = requests.options("http://example.com/1")
            assert resp.status_code == 888

        run()
        assert_reset()

    def test_patch(self):
        @responses.activate
        def run():
            responses.patch("http://example.com/1", status=888)
            resp = requests.patch("http://example.com/1")
            assert resp.status_code == 888

        run()
        assert_reset()

    def test_post(self):
        @responses.activate
        def run():
            responses.post("http://example.com/1", status=888)
            resp = requests.post("http://example.com/1")
            assert resp.status_code == 888

        run()
        assert_reset()

    def test_put(self):
        @responses.activate
        def run():
            responses.put("http://example.com/1", status=888)
            resp = requests.put("http://example.com/1")
            assert resp.status_code == 888

        run()
        assert_reset()


class TestUnitTestPatchSetup:
    """Validates that ``RequestsMock`` could be used as ``mock.patch``.

    This class is present as example in README.rst

    """

    def setup_method(self):
        self.r_mock = responses.RequestsMock(assert_all_requests_are_fired=True)
        self.r_mock.start()
        self.r_mock.get("https://example.com", status=505)
        self.r_mock.put("https://example.com", status=506)

    def teardown_method(self):
        self.r_mock.stop()
        self.r_mock.reset()

        assert_reset()

    def test_function(self):
        resp = requests.get("https://example.com")
        assert resp.status_code == 505

        resp = requests.put("https://example.com")
        assert resp.status_code == 506


class TestUnitTestPatchSetupRaises:
    """Validate that teardown raises if not all requests were executed.

    Similar to ``TestUnitTestPatchSetup``.

    """

    def setup_method(self):
        self.r_mock = responses.RequestsMock()
        self.r_mock.start()
        self.r_mock.get("https://example.com", status=505)
        self.r_mock.put("https://example.com", status=506)

    def teardown_method(self):
        with pytest.raises(AssertionError) as exc:
            self.r_mock.stop()
        self.r_mock.reset()

        assert "[('PUT', 'https://example.com/')]" in str(exc.value)

        assert_reset()

    def test_function(self):
        resp = requests.get("https://example.com")
        assert resp.status_code == 505


def test_reset_in_the_middle():
    @responses.activate
    def run():
        with responses.RequestsMock() as rsps2:
            rsps2.reset()
        responses.add(responses.GET, "https://example.invalid", status=200)
        resp = requests.request("GET", "https://example.invalid")
        assert resp.status_code == 200

    run()
    assert_reset()


def test_redirect():
    @responses.activate
    def run():
        # create multiple Response objects where first two contain redirect headers
        rsp1 = responses.Response(
            responses.GET,
            "http://example.com/1",
            status=301,
            headers={"Location": "http://example.com/2"},
        )
        rsp2 = responses.Response(
            responses.GET,
            "http://example.com/2",
            status=301,
            headers={"Location": "http://example.com/3"},
        )
        rsp3 = responses.Response(responses.GET, "http://example.com/3", status=200)

        # register above generated Responses in `response` module
        responses.add(rsp1)
        responses.add(rsp2)
        responses.add(rsp3)

        # do the first request in order to generate genuine `requests` response
        # this object will contain genuine attributes of the response, like `history`
        rsp = requests.get("http://example.com/1")
        responses.calls.reset()

        # customize exception with `response` attribute
        my_error = requests.ConnectionError("custom error")
        my_error.response = rsp

        # update body of the 3rd response with Exception, this will be raised during execution
        rsp3.body = my_error

        with pytest.raises(requests.ConnectionError) as exc_info:
            requests.get("http://example.com/1")

        assert exc_info.value.args[0] == "custom error"
        assert rsp1.url in exc_info.value.response.history[0].url
        assert rsp2.url in exc_info.value.response.history[1].url

    run()
    assert_reset()


class TestMaxRetry:
    def set_session(self, total=4, raise_on_status=True):
        session = requests.Session()

        adapter = requests.adapters.HTTPAdapter(
            max_retries=Retry(
                total=total,
                backoff_factor=0.1,
                status_forcelist=[500],
                allowed_methods=["GET", "POST", "PATCH"],
                raise_on_status=raise_on_status,
            )
        )
        session.mount("https://", adapter)
        return session

    def test_max_retries(self):
        """This example is present in README.rst"""

        @responses.activate(registry=registries.OrderedRegistry)
        def run():
            url = "https://example.com"
            rsp1 = responses.get(url, body="Error", status=500)
            rsp2 = responses.get(url, body="Error", status=500)
            rsp3 = responses.get(url, body="Error", status=500)
            rsp4 = responses.get(url, body="OK", status=200)

            session = self.set_session()

            resp = session.get(url)

            assert resp.status_code == 200
            assert rsp1.call_count == 1
            assert rsp2.call_count == 1
            assert rsp3.call_count == 1
            assert rsp4.call_count == 1

        run()
        assert_reset()

    @pytest.mark.parametrize("raise_on_status", (True, False))
    def test_max_retries_exceed(self, raise_on_status):  # type: ignore[misc]
        @responses.activate(registry=registries.OrderedRegistry)
        def run():
            url = "https://example.com"
            rsp1 = responses.get(url, body="Error", status=500)
            rsp2 = responses.get(url, body="Error", status=500)
            rsp3 = responses.get(url, body="Error", status=500)

            session = self.set_session(total=2, raise_on_status=raise_on_status)

            if raise_on_status:
                with pytest.raises(RetryError):
                    session.get(url)
            else:
                resp = session.get(url)
                assert resp.status_code == 500

            assert rsp1.call_count == 1
            assert rsp2.call_count == 1
            assert rsp3.call_count == 1

        run()
        assert_reset()

    def test_max_retries_exceed_msg(self):
        @responses.activate(registry=registries.OrderedRegistry)
        def run():
            url = "https://example.com"
            responses.get(url, body="Error", status=500)
            responses.get(url, body="Error", status=500)

            session = self.set_session(total=1)

            with pytest.raises(RetryError) as err:
                session.get(url)

            assert "too many 500 error responses" in str(err.value)

        run()
        assert_reset()

    def test_adapter_retry_untouched(self):
        """Validate that every new request uses brand-new Retry object"""

        @responses.activate(registry=registries.OrderedRegistry)
        def run():
            url = "https://example.com"
            error_rsp = responses.get(url, body="Error", status=500)
            responses.add(error_rsp)
            responses.add(error_rsp)
            ok_rsp = responses.get(url, body="OK", status=200)

            responses.add(error_rsp)
            responses.add(error_rsp)
            responses.add(error_rsp)
            responses.add(ok_rsp)

            session = self.set_session()
            resp = session.get(url)
            assert resp.status_code == 200

            resp = session.get(url)
            assert resp.status_code == 200

            assert len(responses.calls) == 8

        run()
        assert_reset()


def test_request_object_attached_to_exception():
    """Validate that we attach `request` object to custom exception supplied as body"""

    @responses.activate
    def run():
        url = "https://httpbin.org/delay/2"
        responses.get(url, body=requests.ReadTimeout())

        try:
            requests.get(url, timeout=1)
        except requests.ReadTimeout as exc:
            assert type(exc.request) == requests.models.PreparedRequest

    run()
    assert_reset()
