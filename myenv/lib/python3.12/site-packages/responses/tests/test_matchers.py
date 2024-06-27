import gzip
import re
from typing import Any
from typing import List
from unittest.mock import Mock

import pytest
import requests
from requests.exceptions import ConnectionError

import responses
from responses import matchers
from responses.tests.test_responses import assert_reset
from responses.tests.test_responses import assert_response


def test_body_match_get():
    @responses.activate
    def run():
        url = "http://example.com"
        responses.add(
            responses.GET,
            url,
            body=b"test",
            match=[matchers.body_matcher("123456")],
        )
        resp = requests.get("http://example.com", data="123456")
        assert_response(resp, "test")

    run()
    assert_reset()


def test_body_match_post():
    @responses.activate
    def run():
        url = "http://example.com"
        responses.add(
            responses.POST,
            url,
            body=b"test",
            match=[matchers.body_matcher("123456")],
        )
        resp = requests.post("http://example.com", data="123456")
        assert_response(resp, "test")

    run()
    assert_reset()


def test_query_string_matcher():
    @responses.activate
    def run():
        url = "http://example.com?test=1&foo=bar"
        responses.add(
            responses.GET,
            url,
            body=b"test",
            match=[matchers.query_string_matcher("test=1&foo=bar")],
        )
        resp = requests.get("http://example.com?test=1&foo=bar")
        assert_response(resp, "test")
        resp = requests.get("http://example.com?foo=bar&test=1")
        assert_response(resp, "test")
        resp = requests.get("http://example.com/?foo=bar&test=1")
        assert_response(resp, "test")

    run()
    assert_reset()


def test_request_matches_post_params():
    @responses.activate
    def run(deprecated):
        if deprecated:
            json_params_matcher = getattr(responses, "json_params_matcher")
            urlencoded_params_matcher = getattr(responses, "urlencoded_params_matcher")
        else:
            json_params_matcher = matchers.json_params_matcher
            urlencoded_params_matcher = matchers.urlencoded_params_matcher

        responses.add(
            method=responses.POST,
            url="http://example.com/",
            body="one",
            match=[json_params_matcher({"page": {"name": "first", "type": "json"}})],
        )
        responses.add(
            method=responses.POST,
            url="http://example.com/",
            body="two",
            match=[urlencoded_params_matcher({"page": "second", "type": "urlencoded"})],
        )

        resp = requests.request(
            "POST",
            "http://example.com/",
            headers={"Content-Type": "x-www-form-urlencoded"},
            data={"page": "second", "type": "urlencoded"},
        )
        assert_response(resp, "two")

        resp = requests.request(
            "POST",
            "http://example.com/",
            headers={"Content-Type": "application/json"},
            json={"page": {"name": "first", "type": "json"}},
        )
        assert_response(resp, "one")

    with pytest.deprecated_call():
        run(deprecated=True)
        assert_reset()

    run(deprecated=False)
    assert_reset()


def test_json_params_matcher_not_strict():
    @responses.activate(assert_all_requests_are_fired=True)
    def run():
        responses.add(
            method=responses.POST,
            url="http://example.com/",
            body="one",
            match=[
                matchers.json_params_matcher(
                    {"page": {"type": "json"}},
                    strict_match=False,
                )
            ],
        )

        resp = requests.request(
            "POST",
            "http://example.com/",
            headers={"Content-Type": "application/json"},
            json={
                "page": {"type": "json", "another": "nested"},
                "not_strict": "must pass",
            },
        )
        assert_response(resp, "one")

    run()
    assert_reset()


def test_json_params_matcher_not_strict_diff_values():
    @responses.activate
    def run():
        responses.add(
            method=responses.POST,
            url="http://example.com/",
            body="one",
            match=[
                matchers.json_params_matcher(
                    {"page": {"type": "json", "diff": "value"}}, strict_match=False
                )
            ],
        )

        with pytest.raises(ConnectionError) as exc:
            requests.request(
                "POST",
                "http://example.com/",
                headers={"Content-Type": "application/json"},
                json={"page": {"type": "json"}, "not_strict": "must pass"},
            )
        assert (
            "- POST http://example.com/ request.body doesn't match: "
            "{'page': {'type': 'json'}} doesn't match {'page': {'type': 'json', 'diff': 'value'}}"
        ) in str(exc.value)

    run()
    assert_reset()


def test_failed_matchers_dont_modify_inputs_order_in_error_message():
    json_a = {"array": ["C", "B", "A"]}
    json_b = '{"array" : ["B", "A", "C"]}'
    mock_request = Mock(body=json_b)
    result = matchers.json_params_matcher(json_a)(mock_request)
    assert result == (
        False,
        (
            "request.body doesn't match: {'array': ['B', 'A', 'C']} "
            "doesn't match {'array': ['C', 'B', 'A']}"
        ),
    )


def test_json_params_matcher_json_list():
    json_a = [{"a": "b"}]
    json_b = '[{"a": "b", "c": "d"}]'
    mock_request = Mock(body=json_b)
    result = matchers.json_params_matcher(json_a)(mock_request)
    assert result == (
        False,
        "request.body doesn't match: [{'a': 'b', 'c': 'd'}] doesn't match [{'a': 'b'}]",
    )


def test_json_params_matcher_json_list_empty():
    json_a: "List[Any]" = []
    json_b = "[]"
    mock_request = Mock(body=json_b)
    result = matchers.json_params_matcher(json_a)(mock_request)
    assert result == (True, "")


def test_json_params_matcher_body_is_gzipped():
    json_a = {"foo": 42, "bar": None}
    json_b = gzip.compress(b'{"foo": 42, "bar": null}')
    mock_request = Mock(body=json_b)
    result = matchers.json_params_matcher(json_a)(mock_request)
    assert result == (True, "")


def test_urlencoded_params_matcher_blank():
    @responses.activate
    def run():
        responses.add(
            method=responses.POST,
            url="http://example.com/",
            body="three",
            match=[
                matchers.urlencoded_params_matcher(
                    {"page": "", "type": "urlencoded"}, allow_blank=True
                )
            ],
        )

        resp = requests.request(
            "POST",
            "http://example.com/",
            headers={"Content-Type": "x-www-form-urlencoded"},
            data={"page": "", "type": "urlencoded"},
        )
        assert_response(resp, "three")

    run()
    assert_reset()


def test_query_params_numbers():
    @responses.activate
    def run():
        expected_query_params = {"float": 5.0, "int": 2}
        responses.add(
            responses.GET,
            "https://example.com/",
            match=[
                matchers.query_param_matcher(expected_query_params),
            ],
        )
        requests.get("https://example.com", params=expected_query_params)

    run()
    assert_reset()


def test_query_param_matcher_loose():
    @responses.activate
    def run():
        expected_query_params = {"only_one_param": "test"}
        responses.add(
            responses.GET,
            "https://example.com/",
            match=[
                matchers.query_param_matcher(expected_query_params, strict_match=False),
            ],
        )
        requests.get(
            "https://example.com", params={"only_one_param": "test", "second": "param"}
        )

    run()
    assert_reset()


def test_query_param_matcher_loose_fail():
    @responses.activate
    def run():
        expected_query_params = {"does_not_exist": "test"}
        responses.add(
            responses.GET,
            "https://example.com/",
            match=[
                matchers.query_param_matcher(expected_query_params, strict_match=False),
            ],
        )
        with pytest.raises(ConnectionError) as exc:
            requests.get(
                "https://example.com",
                params={"only_one_param": "test", "second": "param"},
            )

        assert (
            "- GET https://example.com/ Parameters do not match. {} doesn't"
            " match {'does_not_exist': 'test'}\n"
            "You can use `strict_match=True` to do a strict parameters check."
        ) in str(exc.value)

    run()
    assert_reset()


def test_request_matches_empty_body():
    def run():
        with responses.RequestsMock(assert_all_requests_are_fired=True) as rsps:
            # test that both json and urlencoded body are empty in matcher and in request
            rsps.add(
                method=responses.POST,
                url="http://example.com/",
                body="one",
                match=[matchers.json_params_matcher(None)],
            )

            rsps.add(
                method=responses.POST,
                url="http://example.com/",
                body="two",
                match=[matchers.urlencoded_params_matcher(None)],
            )

            resp = requests.request("POST", "http://example.com/")
            assert_response(resp, "one")

            resp = requests.request(
                "POST",
                "http://example.com/",
                headers={"Content-Type": "x-www-form-urlencoded"},
            )
            assert_response(resp, "two")

        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            # test exception raise if matcher body is None but request data is not None
            rsps.add(
                method=responses.POST,
                url="http://example.com/",
                body="one",
                match=[matchers.json_params_matcher(None)],
            )

            with pytest.raises(ConnectionError) as excinfo:
                requests.request(
                    "POST",
                    "http://example.com/",
                    json={"my": "data"},
                    headers={"Content-Type": "application/json"},
                )

            msg = str(excinfo.value)
            assert "request.body doesn't match: {'my': 'data'} doesn't match {}" in msg

        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            rsps.add(
                method=responses.POST,
                url="http://example.com/",
                body="two",
                match=[matchers.urlencoded_params_matcher(None)],
            )
            with pytest.raises(ConnectionError) as excinfo:
                requests.request(
                    "POST",
                    "http://example.com/",
                    headers={"Content-Type": "x-www-form-urlencoded"},
                    data={"page": "second", "type": "urlencoded"},
                )
            msg = str(excinfo.value)
            assert (
                "request.body doesn't match: {'page': 'second', "
                "'type': 'urlencoded'} doesn't match {}"
            ) in msg

    run()
    assert_reset()


def test_request_matches_params():
    @responses.activate
    def run():
        url = "http://example.com/test"
        params = {"hello": "world", "I am": "a big test"}
        responses.add(
            method=responses.GET,
            url=url,
            body="test",
            match=[matchers.query_param_matcher(params)],
            match_querystring=False,
        )

        # exchange parameter places for the test
        params = {
            "I am": "a big test",
            "hello": "world",
        }
        resp = requests.get(url, params=params)

        constructed_url = r"http://example.com/test?I+am=a+big+test&hello=world"
        assert resp.url == constructed_url
        assert resp.request.url == constructed_url

        resp_params = getattr(resp.request, "params")
        assert resp_params == params

    run()
    assert_reset()


def test_fail_matchers_error():
    """
    Validate that Exception is raised if request does not match responses.matchers
        validate matchers.urlencoded_params_matcher
        validate matchers.json_params_matcher
        validate matchers.query_param_matcher
        validate matchers.request_kwargs_matcher
    :return: None
    """

    def run():
        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            rsps.add(
                "POST",
                "http://example.com",
                match=[matchers.urlencoded_params_matcher({"foo": "bar"})],
            )
            rsps.add(
                "POST",
                "http://example.com",
                match=[matchers.json_params_matcher({"fail": "json"})],
            )

            with pytest.raises(ConnectionError) as excinfo:
                requests.post("http://example.com", data={"id": "bad"})

            msg = str(excinfo.value)
            assert (
                "request.body doesn't match: {'id': 'bad'} doesn't match {'foo': 'bar'}"
                in msg
            )

            assert (
                "request.body doesn't match: JSONDecodeError: Cannot parse request.body"
                in msg
            )

        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            rsps.add(
                "GET",
                "http://111.com",
                match=[matchers.query_param_matcher({"my": "params"})],
            )

            rsps.add(
                method=responses.GET,
                url="http://111.com/",
                body="two",
                match=[matchers.json_params_matcher({"page": "one"})],
            )

            with pytest.raises(ConnectionError) as excinfo:
                requests.get(
                    "http://111.com", params={"id": "bad"}, json={"page": "two"}
                )

            msg = str(excinfo.value)
            assert (
                "Parameters do not match. {'id': 'bad'} doesn't match {'my': 'params'}"
                in msg
            )
            assert (
                "request.body doesn't match: {'page': 'two'} doesn't match {'page': 'one'}"
                in msg
            )

        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            req_kwargs = {
                "stream": True,
                "verify": False,
            }
            rsps.add(
                "GET",
                "http://111.com",
                match=[matchers.request_kwargs_matcher(req_kwargs)],
            )

            with pytest.raises(ConnectionError) as excinfo:
                requests.get("http://111.com", stream=True)

            msg = str(excinfo.value)
            assert (
                "Arguments don't match: "
                "{'stream': True, 'verify': True} doesn't match {'stream': True, 'verify': False}"
            ) in msg

    run()
    assert_reset()


@pytest.mark.parametrize(
    "req_file,match_file",
    [
        (b"Old World!", "Old World!"),
        ("Old World!", b"Old World!"),
        (b"Old World!", b"Old World!"),
        ("Old World!", "Old World!"),
        (b"\xacHello World!", b"\xacHello World!"),
    ],
)
def test_multipart_matcher(req_file, match_file):  # type: ignore[misc]
    @responses.activate
    def run():
        req_data = {"some": "other", "data": "fields"}
        responses.add(
            responses.POST,
            url="http://httpbin.org/post",
            match=[
                matchers.multipart_matcher(
                    files={"file_name": match_file}, data=req_data
                )
            ],
        )
        resp = requests.post(
            "http://httpbin.org/post", data=req_data, files={"file_name": req_file}
        )
        assert resp.status_code == 200

        with pytest.raises(TypeError):
            responses.add(
                responses.POST,
                url="http://httpbin.org/post",
                match=[matchers.multipart_matcher(files={})],
            )

    run()
    assert_reset()


def test_multipart_matcher_fail():
    """
    Validate that Exception is raised if request does not match responses.matchers
        validate matchers.multipart_matcher
    :return: None
    """

    def run():
        # different file contents
        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            req_data = {"some": "other", "data": "fields"}
            req_files = {"file_name": b"Old World!"}
            rsps.add(
                responses.POST,
                url="http://httpbin.org/post",
                match=[matchers.multipart_matcher(req_files, data=req_data)],
            )

            with pytest.raises(ConnectionError) as excinfo:
                requests.post(
                    "http://httpbin.org/post",
                    data=req_data,
                    files={"file_name": b"New World!"},
                )

            msg = str(excinfo.value)
            assert "multipart/form-data doesn't match. Request body differs." in msg

            assert (
                r'\r\nContent-Disposition: form-data; name="file_name"; '
                r'filename="file_name"\r\n\r\nOld World!\r\n'
            ) in msg
            assert (
                r'\r\nContent-Disposition: form-data; name="file_name"; '
                r'filename="file_name"\r\n\r\nNew World!\r\n'
            ) in msg

        # x-www-form-urlencoded request
        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            req_data = {"some": "other", "data": "fields"}
            req_files = {"file_name": b"Old World!"}
            rsps.add(
                responses.POST,
                url="http://httpbin.org/post",
                match=[matchers.multipart_matcher(req_files, data=req_data)],
            )

            with pytest.raises(ConnectionError) as excinfo:
                requests.post("http://httpbin.org/post", data=req_data)

            msg = str(excinfo.value)
            assert (
                "multipart/form-data doesn't match. Request headers['Content-Type'] is different."
                in msg
            )
            assert (
                "application/x-www-form-urlencoded isn't equal to multipart/form-data; boundary="
                in msg
            )

        # empty body request
        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            req_files = {"file_name": b"Old World!"}
            rsps.add(
                responses.POST,
                url="http://httpbin.org/post",
                match=[matchers.multipart_matcher(req_files)],
            )

            with pytest.raises(ConnectionError) as excinfo:
                requests.post("http://httpbin.org/post")

            msg = str(excinfo.value)
            assert "Request is missing the 'Content-Type' header" in msg

    run()
    assert_reset()


def test_query_string_matcher_raises():
    """
    Validate that Exception is raised if request does not match responses.matchers
        validate matchers.query_string_matcher
            :return: None
    """

    def run():
        with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
            rsps.add(
                "GET",
                "http://111.com",
                match=[matchers.query_string_matcher("didi=pro")],
            )

            with pytest.raises(ConnectionError) as excinfo:
                requests.get("http://111.com", params={"test": "1", "didi": "pro"})

            msg = str(excinfo.value)
            assert (
                "Query string doesn't match. {'didi': 'pro', 'test': '1'} "
                "doesn't match {'didi': 'pro'}"
            ) in msg

    run()
    assert_reset()


def test_request_matches_headers():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(
            method=responses.GET,
            url=url,
            json={"success": True},
            match=[matchers.header_matcher({"Accept": "application/json"})],
        )

        responses.add(
            method=responses.GET,
            url=url,
            body="success",
            match=[matchers.header_matcher({"Accept": "text/plain"})],
        )

        # the actual request can contain extra headers (requests always adds some itself anyway)
        resp = requests.get(
            url, headers={"Accept": "application/json", "Accept-Charset": "utf-8"}
        )
        assert_response(resp, body='{"success": true}', content_type="application/json")

        resp = requests.get(url, headers={"Accept": "text/plain"})
        assert_response(resp, body="success", content_type="text/plain")

    run()
    assert_reset()


def test_request_header_value_mismatch_raises():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(
            method=responses.GET,
            url=url,
            json={"success": True},
            match=[matchers.header_matcher({"Accept": "application/json"})],
        )

        with pytest.raises(ConnectionError) as excinfo:
            requests.get(url, headers={"Accept": "application/xml"})

        msg = str(excinfo.value)
        assert (
            "Headers do not match: {'Accept': 'application/xml'} doesn't match "
            "{'Accept': 'application/json'}"
        ) in msg

    run()
    assert_reset()


def test_request_headers_missing_raises():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(
            method=responses.GET,
            url=url,
            json={"success": True},
            match=[matchers.header_matcher({"x-custom-header": "foo"})],
        )

        with pytest.raises(ConnectionError) as excinfo:
            requests.get(url, headers={})

        msg = str(excinfo.value)
        assert (
            "Headers do not match: {} doesn't match {'x-custom-header': 'foo'}"
        ) in msg

    run()
    assert_reset()


def test_request_matches_headers_strict_match():
    @responses.activate
    def run():
        url = "http://example.com/"
        responses.add(
            method=responses.GET,
            url=url,
            body="success",
            match=[
                matchers.header_matcher({"Accept": "text/plain"}, strict_match=True)
            ],
        )

        # requests will add some extra headers of its own, so we have to use prepared requests
        session = requests.Session()

        # make sure we send *just* the header we're expectin
        prepped = session.prepare_request(
            requests.Request(
                method="GET",
                url=url,
            )
        )
        prepped.headers.clear()
        prepped.headers["Accept"] = "text/plain"

        resp = session.send(prepped)
        assert_response(resp, body="success", content_type="text/plain")

        # include the "Accept-Charset" header, which will fail to match
        prepped = session.prepare_request(
            requests.Request(
                method="GET",
                url=url,
            )
        )
        prepped.headers.clear()
        prepped.headers["Accept"] = "text/plain"
        prepped.headers["Accept-Charset"] = "utf-8"

        with pytest.raises(ConnectionError) as excinfo:
            session.send(prepped)

        msg = str(excinfo.value)
        assert (
            "Headers do not match: {'Accept': 'text/plain', 'Accept-Charset': 'utf-8'} "
            "doesn't match {'Accept': 'text/plain'}"
        ) in msg

    run()
    assert_reset()


def test_fragment_identifier_matcher():
    @responses.activate
    def run():
        responses.add(
            responses.GET,
            "http://example.com",
            match=[matchers.fragment_identifier_matcher("test=1&foo=bar")],
            body=b"test",
        )

        resp = requests.get("http://example.com#test=1&foo=bar")
        assert_response(resp, "test")

    run()
    assert_reset()


def test_fragment_identifier_matcher_error():
    @responses.activate
    def run():
        responses.add(
            responses.GET,
            "http://example.com/",
            match=[matchers.fragment_identifier_matcher("test=1")],
        )
        responses.add(
            responses.GET,
            "http://example.com/",
            match=[matchers.fragment_identifier_matcher(None)],
        )

        with pytest.raises(ConnectionError) as excinfo:
            requests.get("http://example.com/#test=2")

        msg = str(excinfo.value)
        assert (
            "URL fragment identifier is different: test=1 doesn't match test=2"
        ) in msg
        assert (
            "URL fragment identifier is different: None doesn't match test=2"
        ) in msg

    run()
    assert_reset()


def test_fragment_identifier_matcher_and_match_querystring():
    @responses.activate
    def run():
        url = "http://example.com?ab=xy&zed=qwe#test=1&foo=bar"
        responses.add(
            responses.GET,
            url,
            match_querystring=True,
            match=[matchers.fragment_identifier_matcher("test=1&foo=bar")],
            body=b"test",
        )

        # two requests to check reversed order of fragment identifier
        resp = requests.get("http://example.com?ab=xy&zed=qwe#test=1&foo=bar")
        assert_response(resp, "test")
        resp = requests.get("http://example.com?zed=qwe&ab=xy#foo=bar&test=1")
        assert_response(resp, "test")

    run()
    assert_reset()


class TestHeaderWithRegex:
    @property
    def url(self):  # type: ignore[misc]
        return "http://example.com/"

    def _register(self):
        responses.add(
            method=responses.GET,
            url=self.url,
            body="success",
            match=[
                matchers.header_matcher(
                    {
                        "Accept": "text/plain",
                        "Message-Signature": re.compile(r'signature="\S+",created=\d+'),
                    },
                    strict_match=True,
                )
            ],
        )

    def test_request_matches_headers_regex(self):
        @responses.activate
        def run():
            # this one can not use common _register method because different headers
            responses.add(
                method=responses.GET,
                url=self.url,
                json={"success": True},
                match=[
                    matchers.header_matcher(
                        {
                            "Message-Signature": re.compile(
                                r'signature="\S+",created=\d+'
                            ),
                            "Authorization": "Bearer API_TOKEN",
                        },
                        strict_match=False,
                    )
                ],
            )
            # the actual request can contain extra headers (requests always adds some itself anyway)
            resp = requests.get(
                self.url,
                headers={
                    "Message-Signature": 'signature="abc",created=1243',
                    "Authorization": "Bearer API_TOKEN",
                },
            )
            assert_response(
                resp, body='{"success": true}', content_type="application/json"
            )

        run()
        assert_reset()

    def test_request_matches_headers_regex_strict_match_regex_failed(self):
        @responses.activate
        def run():
            self._register()
            session = requests.Session()
            # requests will add some extra headers of its own, so we have to use prepared requests
            prepped = session.prepare_request(
                requests.Request(
                    method="GET",
                    url=self.url,
                )
            )
            prepped.headers.clear()
            prepped.headers["Accept"] = "text/plain"
            prepped.headers["Message-Signature"] = 'signature="123",created=abc'
            with pytest.raises(ConnectionError) as excinfo:
                session.send(prepped)
            msg = str(excinfo.value)
            assert (
                "Headers do not match: {'Accept': 'text/plain', 'Message-Signature': "
                """'signature="123",created=abc'} """
                "doesn't match {'Accept': 'text/plain', 'Message-Signature': "
                "re.compile('signature=\"\\\\S+\",created=\\\\d+')}"
            ) in msg

        run()
        assert_reset()

    def test_request_matches_headers_regex_strict_match_mismatched_field(self):
        @responses.activate
        def run():
            self._register()
            # requests will add some extra headers of its own, so we have to use prepared requests
            session = requests.Session()
            prepped = session.prepare_request(
                requests.Request(
                    method="GET",
                    url=self.url,
                )
            )
            prepped.headers.clear()
            prepped.headers["Accept"] = "text/plain"
            prepped.headers["Accept-Charset"] = "utf-8"
            # "Accept-Charset" header will fail to match to "Message-Signature"
            with pytest.raises(ConnectionError) as excinfo:
                session.send(prepped)
            msg = str(excinfo.value)
            assert (
                "Headers do not match: {'Accept': 'text/plain', 'Accept-Charset': 'utf-8'} "
                "doesn't match {'Accept': 'text/plain', 'Message-Signature': "
                "re.compile('signature=\"\\\\S+\",created=\\\\d+')}"
            ) in msg

        run()
        assert_reset()

    def test_request_matches_headers_regex_strict_match_mismatched_number(self):
        @responses.activate
        def run():
            self._register()
            # requests will add some extra headers of its own, so we have to use prepared requests
            session = requests.Session()
            # include the "Accept-Charset" header, which will fail to match
            prepped = session.prepare_request(
                requests.Request(
                    method="GET",
                    url=self.url,
                )
            )
            prepped.headers.clear()
            prepped.headers["Accept"] = "text/plain"
            prepped.headers["Accept-Charset"] = "utf-8"
            prepped.headers["Message-Signature"] = 'signature="abc",created=1243'
            with pytest.raises(ConnectionError) as excinfo:
                session.send(prepped)
            msg = str(excinfo.value)
            assert (
                "Headers do not match: {'Accept': 'text/plain', 'Accept-Charset': 'utf-8', "
                """'Message-Signature': 'signature="abc",created=1243'} """
                "doesn't match {'Accept': 'text/plain', 'Message-Signature': "
                "re.compile('signature=\"\\\\S+\",created=\\\\d+')}"
            ) in msg

        run()
        assert_reset()

    def test_request_matches_headers_regex_strict_match_positive(self):
        @responses.activate
        def run():
            self._register()
            # requests will add some extra headers of its own, so we have to use prepared requests
            session = requests.Session()
            prepped = session.prepare_request(
                requests.Request(
                    method="GET",
                    url=self.url,
                )
            )
            prepped.headers.clear()
            prepped.headers["Accept"] = "text/plain"
            prepped.headers["Message-Signature"] = 'signature="abc",created=1243'
            resp = session.send(prepped)
            assert_response(resp, body="success", content_type="text/plain")

        run()
        assert_reset()
